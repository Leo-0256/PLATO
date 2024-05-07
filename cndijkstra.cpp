#include <utility>
#include <vector>
#include <queue>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

#define MAXV 5000
#define INF 100000000

struct edge
{
    int to, cost;
};

typedef pair<int, int> P;
int V;

vector<edge> G[MAXV];

int d[MAXV][MAXV];

void dijkstra(int s)
{
    priority_queue<P, vector<P>, greater<P>> que;
    d[s][s] = 0;
    que.push(P(0, s));
    while (!que.empty())
    {
        P p = que.top();
        que.pop();
        int v = p.second;
        if (d[s][v] < p.first)
            continue;
        for (int i = 0; i < G[v].size(); i++)
        {
            edge e = G[v][i];
            if (d[s][e.to] > d[s][v] + e.cost)
            {
                d[s][e.to] = d[s][v] + e.cost;
                que.push(P(d[s][e.to], e.to));
            }
        }
    }
}

extern "C"
{
    void cndijkstra(int num_orbit, int num_sat_per_orbit, int *distanceMatrix)
    {
        int num_sat = num_orbit * num_sat_per_orbit;
        for (int i = 0; i < num_orbit; i++)
        {
            for (int j = 0; j < num_sat_per_orbit; j++)
            {
                int index = i * num_sat_per_orbit + j;
                int index_next_in_same_orbit = i * num_sat_per_orbit + (j + 1) % num_sat_per_orbit;
                int index_next_in_next_orbit = ((i + 1) % num_orbit) * num_sat_per_orbit + j;

                int distance_between_same_orbit = *(distanceMatrix + index * num_sat + index_next_in_same_orbit);
                int distance_between_next_orbit = *(distanceMatrix + index * num_sat + index_next_in_next_orbit);

                G[index].push_back({index_next_in_same_orbit, distance_between_same_orbit});
                G[index].push_back({index_next_in_next_orbit, distance_between_next_orbit});
                G[index_next_in_same_orbit].push_back({index, distance_between_same_orbit});
                G[index_next_in_next_orbit].push_back({index, distance_between_next_orbit});
            }
        }

        // fstream f;
        // f.open("ccc.txt", ios::out);
        // for (int i = 0; i < num_sat; i++)
        // {
        //     for (int j = 0; j < G[i].size(); j++)
        //     {
        //         edge e = G[i][j];
        //         f << e.to << " " << e.cost << "; ";
        //     }
        //     f << endl;
        // }

        for (int i = 0; i < MAXV; i++)
        {
            for (int j = 0; j < MAXV; j++)
            {
                d[i][j] = INF;
            }
        }

        // fstream f;
        // f.open("ccc.txt", ios::out);
        // for (int i = 0; i < num_sat; i++)
        // {
        //     for (int j = 0; j < num_sat; j++)
        //     {
        //         f << d[i][j] << " ";
        //     }
        //     f << endl;
        // }
        // f.close();

        for (int i = 0; i < num_sat; i++)
        {
            dijkstra(i);
        }

        // fstream f;
        // f.open("ccc.txt", ios::out);
        // for (int i = 0; i < num_sat; i++)
        // {
        //     for (int j = 0; j < num_sat; j++)
        //     {
        //         f << d[i][j] << " ";
        //     }
        //     f << endl;
        // }
        // f.close();

        for (int i = 0; i < num_sat; i++)
        {
            for (int j = 0; j < num_sat; j++)
            {
                *(distanceMatrix + i * num_sat + j) = d[i][j];
            }
        }
    }
}

// g++ -shared -Wl,-soname,cndijkstra -o cndijkstra.so -fPIC cndijkstra.cpp