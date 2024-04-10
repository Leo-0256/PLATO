import sys
import copy
import math
import ephem
from astropy.time import Time
import networkx as nx
import numpy as np
from skyfield.api import EarthSatellite, load
import datetime


class TLEGraph:
    def __init__(self, tle_filepath):
        self.__tle_filepath = tle_filepath
        self.__num_orbits = 0
        self.__num_sat_per_orbit = 0
        self.__num_satellites = 0
        self.__distance_matrix = []
        self.__parent_matrix = []
        self.__epoch = Time("2000-01-01 00:00:00", scale="tdb")
        self.__generate_distance_matrix_from_tle_filepath()

    def update_graph_by_new_tle_file(self, tle_filepath):
        self.__tle_filepath = tle_filepath

    def get_all_paths_by_floyd(self, time_str, out_filepath):
        self.__update_distance_matrix(time_str)
        self.__floyd()
        f_out = open(out_filepath, "w+")
        for i in range(self.__num_satellites):
            for j in range(self.__num_satellites):
                print("Path({}-->{}): ".format(i, j), end="", file=f_out)
                self.__print_path(i, j, j, outfile=f_out)
                print(
                    " cost: {}".format(self.__distance_matrix[i][j]),
                    file=f_out,
                    flush=True,
                )

    def __generate_distance_matrix_from_tle_filepath(self):
        with open(self.__tle_filepath, "r") as f_tle:
            firstline = f_tle.readline().strip().split()
            self.__num_orbits = int(firstline[0])
            self.__num_sat_per_orbit = int(firstline[1])
            self.__num_satellites = self.__num_orbits * self.__num_sat_per_orbit
            self.__distance_matrix = [
                [math.inf if i != j else 0 for j in range(self.__num_satellites)]
                for i in range(self.__num_satellites)
            ]
            self.__parent_matrix = [
                [i] * self.__num_satellites for i in range(self.__num_satellites)
            ]

    def __update_distance_matrix(self, time_str):
        """
        Update self.__satellite_list and self.distance_matrix by new tle file, the constellation parameter must be consistent, and we
        do not do the check.
        """
        with open(self.__tle_filepath, "r") as f_tle:
            self.__satellite_list = []
            _ = f_tle.readline()
            for row in f_tle:
                row = row.strip()
                line1 = f_tle.readline().strip()
                line2 = f_tle.readline().strip()
                cur_satellite = ephem.readtle(row, line1, line2)
                self.__satellite_list.append(cur_satellite)

        epoch = self.__epoch
        for i in range(self.__num_orbits):
            for j in range(self.__num_sat_per_orbit):
                cur_satellite_index = i * self.__num_sat_per_orbit + j
                next_sat_in_same_orbit_index = (
                    i * self.__num_sat_per_orbit + (j + 1) % self.__num_sat_per_orbit
                )
                next_sat_in_adjacent_orbit_index = (
                    i + 1
                ) % self.__num_orbits * self.__num_sat_per_orbit + j
                cur_satellite = self.__satellite_list[cur_satellite_index]
                next_sat_in_same_orbit = self.__satellite_list[
                    next_sat_in_same_orbit_index
                ]
                next_sat_in_adjacent_orbit = self.__satellite_list[
                    next_sat_in_adjacent_orbit_index
                ]
                distance_between_sat_in_same_orbit = (
                    self.__distance_m_between_satellites(
                        cur_satellite, next_sat_in_same_orbit, str(epoch), str(time_str)
                    )
                )
                distance_between_sat_in_adjacent_orbit = (
                    self.__distance_m_between_satellites(
                        cur_satellite,
                        next_sat_in_adjacent_orbit,
                        str(epoch),
                        str(time_str),
                    )
                )
                self.__distance_matrix[cur_satellite_index][
                    next_sat_in_same_orbit_index
                ] = self.__distance_matrix[next_sat_in_same_orbit_index][
                    cur_satellite_index
                ] = distance_between_sat_in_same_orbit
                self.__distance_matrix[cur_satellite_index][
                    next_sat_in_adjacent_orbit_index
                ] = self.__distance_matrix[next_sat_in_adjacent_orbit_index][
                    cur_satellite_index
                ] = distance_between_sat_in_adjacent_orbit

    def __distance_m_between_satellites(self, sat1, sat2, epoch_str, date_str):
        """
        Computes the straight distance between two satellites in meters.

        :param sat1:       The first satellite
        :param sat2:       The other satellite
        :param epoch_str:  Epoch time of the observer (string)
        :param date_str:   The time instant when the distance should be measured (string)

        :return: The distance between the satellites in meters
        """

        # Create an observer somewhere on the planet
        observer = ephem.Observer()
        observer.epoch = epoch_str
        observer.date = date_str
        observer.lat = 0
        observer.lon = 0
        observer.elevation = 0

        # Calculate the relative location of the satellites to this observer
        sat1.compute(observer)
        sat2.compute(observer)

        # Calculate the angle observed by the observer to the satellites (this is done because the .compute() calls earlier)
        angle_radians = float(repr(ephem.separation(sat1, sat2)))

        return int(
            math.sqrt(
                sat1.range**2
                + sat2.range**2
                - (2 * sat1.range * sat2.range * math.cos(angle_radians))
            )
        )

    def __floyd(self):
        n = self.__num_satellites
        graph = self.__distance_matrix
        parents = self.__parent_matrix
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if graph[i][k] + graph[k][j] < graph[i][j]:
                        graph[i][j] = graph[i][k] + graph[k][j]
                        parents[i][j] = parents[k][j]

    def __print_path(self, i, j, e, outfile=sys.stdout):
        if i != j:
            self.__print_path(i, self.__parent_matrix[i][j], j, outfile)
        print(j, end="", file=outfile)
        if j != e:
            print("-->", end="", file=outfile)


def compute_path_by_networkx(tles_filepath, time_str, out_filepath):
    epoch = Time("2000-01-01 00:00:00", scale="tdb")
    graph = nx.Graph()
    num_orbits = 0
    num_sat_per_orbit = 0
    num_satellites = 0
    satellite_list = []

    with open(tles_filepath, "r") as f_tle:
        firstline = f_tle.readline().strip().split()
        num_orbits, num_sat_per_orbit = int(firstline[0]), int(firstline[1])
        num_satellites = num_orbits * num_sat_per_orbit
        for row in f_tle:
            row = row.strip()
            line1 = f_tle.readline().strip()
            line2 = f_tle.readline().strip()
            cur_satellite = ephem.readtle(row, line1, line2)
            satellite_list.append(cur_satellite)
    graph.add_nodes_from([i for i in range(num_satellites)])

    for i in range(num_orbits):
        for j in range(num_sat_per_orbit):
            cur_satellite_index = i * num_sat_per_orbit + j
            next_sat_in_same_orbit_index = (
                i * num_sat_per_orbit + (j + 1) % num_sat_per_orbit
            )
            next_sat_in_adjacent_orbit_index = (
                i + 1
            ) % num_orbits * num_sat_per_orbit + j
            cur_satellite = satellite_list[cur_satellite_index]
            next_sat_in_same_orbit = satellite_list[next_sat_in_same_orbit_index]
            next_sat_in_adjacent_orbit = satellite_list[
                next_sat_in_adjacent_orbit_index
            ]
            distance_between_sat_in_same_orbit = distance_m_between_satellites(
                cur_satellite, next_sat_in_same_orbit, str(epoch), str(time_str)
            )
            distance_between_sat_in_adjacent_orbit = distance_m_between_satellites(
                cur_satellite,
                next_sat_in_adjacent_orbit,
                str(epoch),
                str(time_str),
            )
            graph.add_edge(
                cur_satellite_index,
                next_sat_in_same_orbit_index,
                distance=distance_between_sat_in_same_orbit,
            )
            graph.add_edge(
                cur_satellite_index,
                next_sat_in_adjacent_orbit_index,
                distance=distance_between_sat_in_adjacent_orbit,
            )

    f_out = open(out_filepath, "w")
    for i in range(num_satellites):
        for j in range(num_satellites):
            path = nx.shortest_path(graph, i, j, weight="distance")
            length = nx.shortest_path_length(graph, i, j, weight="distance")
            print("Path({}-->{}): ".format(i, j), end="", file=f_out)
            for k in range(len(path)):
                print(path[k], end="", file=f_out)
                if k != len(path) - 1:
                    print("-->", end="", file=f_out)
            print(
                " cost: {}".format(length),
                file=f_out,
                flush=True,
            )


def distance_m_between_satellites(sat1, sat2, epoch_str, date_str):
    """
    Computes the straight distance between two satellites in meters.

    :param sat1:       The first satellite
    :param sat2:       The other satellite
    :param epoch_str:  Epoch time of the observer (string)
    :param date_str:   The time instant when the distance should be measured (string)

    :return: The distance between the satellites in meters
    """

    # Create an observer somewhere on the planet
    observer = ephem.Observer()
    observer.epoch = epoch_str
    observer.date = date_str
    observer.lat = 0
    observer.lon = 0
    observer.elevation = 0

    # Calculate the relative location of the satellites to this observer
    sat1.compute(observer)
    sat2.compute(observer)

    # Calculate the angle observed by the observer to the satellites (this is done because the .compute() calls earlier)
    angle_radians = float(repr(ephem.separation(sat1, sat2)))

    return int(
        math.sqrt(
            sat1.range**2
            + sat2.range**2
            - (2 * sat1.range * sat2.range * math.cos(angle_radians))
        )
    )


if __name__ == "__main__":

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_2x2_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_2x2_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms") # 1.43 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_2x2_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_3x3_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_3x3_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms") # 2.175 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_3x3_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_4x4_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_4x4_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds /1000 , "ms") # 5.102 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_4x4_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_5x5_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_5x5_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds /1000 , "ms") # 10.54 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_5x5_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_10x10_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_10x10_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms") # 319.09 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_10x10_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_16x16_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_16x16_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms") # 991.197 ms
    # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_16x16_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_25x25_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_25x25_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms")
    # # compute_path_by_networkx(tles_filepath, time_str, "./PATH/demo_25x25_nx_path.txt")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/starlink_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # endtime = datetime.datetime.now() # 205.545 ms
    # print((endtime - starttime).microseconds / 1000, "ms")
    # starttime = datetime.datetime.now()
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/starlink_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/kuiper_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms") # 635.933 ms
    # starttime = datetime.datetime.now()
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/kuiper_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_50x50_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_50x50_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms")

    # starttime = datetime.datetime.now()
    # tles_filepath = "./TLE/demo_100x100_tle.txt"
    # tle_graph = TLEGraph(tles_filepath)
    # endtime = datetime.datetime.now() # 369.105 ms
    # print((endtime - starttime).microseconds / 1000, "ms")
    # starttime = datetime.datetime.now()
    # time_str = "2000-01-01 00:00:00"
    # tle_graph.get_all_paths_by_floyd(time_str, "./PATH/demo_100x100_path.txt")
    # endtime = datetime.datetime.now()
    # print((endtime - starttime).microseconds / 1000, "ms")

    pass
