def generate_tles(
    out_filepath,
    constellation_name,
    num_orbits,
    num_sats_per_orbit,
    phase_factor,
    inclination_degree,
    eccentricity,
    arg_of_perigee_degree,
    mean_motion_rev_per_day,
):
    """
    Generate TLE file from parameters.
    Default epoch is 2000-01-01 00:00:00, which is 00001 in yyddd format.
    """

    with open(out_filepath, "w+") as f_out:
        # First line:
        # <number of orbits> <number of satellites per orbit>
        f_out.write("%d %d\n" % (num_orbits, num_sats_per_orbit))

        # Each of the subsequent (number of orbits * number of satellites per orbit) blocks
        # define a satellite as follows:
        # <constallation_name> <global satellite id>
        # <TLE line 1>
        # <TLE line 2>
        satellite_counter = 0
        phase_diff = 360.0 * phase_factor / (num_orbits * num_sats_per_orbit)
        for orbit in range(0, num_orbits):
            raan_degree = orbit * 360 / num_orbits
            orbit_wise_shift = orbit * phase_diff
            for n_sat in range(0, num_sats_per_orbit):
                mean_anomaly_degree = (
                    orbit_wise_shift + (n_sat * 360 / num_sats_per_orbit)
                ) % 360

                tle_line1 = (
                    "1 %05dU 00000ABC 00001.00000000  .00000000  00000-0  00000+0 0    0"
                    % (satellite_counter + 1)
                )

                tle_line2 = "2 %05d %s %s %s %s %s %s    0" % (
                    satellite_counter + 1,
                    ("%8.4f" % inclination_degree),
                    ("%8.4f" % raan_degree),
                    ("%.7f" % eccentricity)[2:],
                    ("%8.4f" % arg_of_perigee_degree),
                    ("%8.4f" % mean_anomaly_degree),
                    ("%11.8f" % mean_motion_rev_per_day),
                )

                tle_line1 = tle_line1 + str(calculate_tle_line_checksum(tle_line1))
                tle_line2 = tle_line2 + str(calculate_tle_line_checksum(tle_line2))

                f_out.write(
                    "%s %d" % (constellation_name, satellite_counter + 1) + "\n"
                )
                f_out.write(tle_line1 + "\n")
                f_out.write(tle_line2 + "\n")

                satellite_counter += 1


def calculate_tle_line_checksum(tle_line_without_checksum):
    if len(tle_line_without_checksum) != 68:
        raise ValueError("TLE line without checksum must have exactly 68 characters")
    checksum = 0
    for i in range(len(tle_line_without_checksum)):
        if tle_line_without_checksum[i].isnumeric():
            checksum += int(tle_line_without_checksum[i])
        if tle_line_without_checksum[i] == "-":
            checksum += 1
    return checksum % 10


if __name__ == "__main__":
    generate_tles("./TLE/starlink_tle.txt", "Starlink", 72, 22, 1, 53.0, 0.0000001, 0, 15.07)
    generate_tles("./TLE/kuiper_tle.txt", "Kuiper", 34, 34, 0, 51.9, 0.0000001, 0, 14.81)

    generate_tles("./TLE/demo_2x2_tle.txt", "Demo-2x2", 2, 2, 0, 53.0, 0.0000001, 0, 15.07)
    generate_tles("./TLE/demo_3x3_tle.txt", "Demo-3x3", 3, 3, 0, 53.0, 0.0000001, 0, 15.07)
    generate_tles("./TLE/demo_4x4_tle.txt", "Demo-4x4", 4, 4, 0, 53.0, 0.0000001, 0, 15.07)
    generate_tles("./TLE/demo_5x5_tle.txt", "Demo-5x5", 5, 5, 0, 53.0, 0.0000001, 0, 15.07)
    generate_tles(
        "./TLE/demo_10x10_tle.txt", "Demo-10x10", 10, 10, 0, 53.0, 0.0000001, 0, 15.07
    )
    generate_tles(
        "./TLE/demo_16x16_tle.txt", "Demo-16x16", 16, 16, 0, 53.0, 0.0000001, 0, 15.07
    )
    generate_tles(
        "./TLE/demo_25x25_tle.txt", "Demo-25x25", 25, 25, 0, 53.0, 0.0000001, 0, 15.07
    )
    generate_tles(
        "./TLE/demo_50x50_tle.txt", "Demo-50x50", 50, 50, 0, 53.0, 0.0000001, 0, 15.07
    )
    generate_tles(
        "./TLE/demo_100x100_tle.txt", "Demo-100x100", 100, 100, 0, 53.0, 0.0000001, 0, 15.07
    )
    pass