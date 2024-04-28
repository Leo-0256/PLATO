#include <iostream>

using namespace std;

extern "C"
{
    void floyd(int n, int *distanceMatrix)
    {
        for (int k = 0; k < n; k++)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (*(distanceMatrix + i * n + k) + *(distanceMatrix + k * n + j) < *(distanceMatrix + i * n + j))
                    {
                        *(distanceMatrix + i * n + j) = *(distanceMatrix + i * n + k) + *(distanceMatrix + k * n + j);
                    }
                }
            }
        }
    }
}
