#include <stdio.h>
#include "mpi.h"

double sliceArea(double h1, double h2, double d) {
    return (h1 + h2) * d / 2;
}

double f(double x) {
    return x * x;
}

int main(int argc, char **argv)
{
    // integral data
    int d = 9999;
    int a = 1, b = 3;
    double dx = double(b - a) / d;

    // MPI
    MPI_Init(&argc, &argv);
    MPI_Status status;
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Processor task limits
    int chunkSize = d / size;
    int start = chunkSize * rank;
    int end = chunkSize * (rank + 1);

    // chunkSize might be rounded down when if d % size != 1
    if (rank == size - 1)
    {
        end = d;
    }

    // Riman integral
    double area = 0;
    for (int i = start; i < end; i++)
    {
        double lower = f(a + dx * i);
        double higher = f(a + dx * (i + 1));
        area += sliceArea(lower, higher, dx);
    }

    // Fork Join Model, all results should be sent to 0
    if (rank != 0)
    {
        MPI_Send(&area, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
        printf("The result from process %d is %f\n", rank, area);
    }
    else
    {
        for (int i = 0; i < size - 1; i++)
        {
            double message;
            MPI_Recv(&message, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
            area += message;
        }
        printf("The TOTAL integral value is %f", area);
    }
    MPI_Finalize();
    return 0;
}