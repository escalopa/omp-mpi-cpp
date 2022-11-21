#include <stdio.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    int numtask, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtask);

    printf("Hello %d %d \n", rank, numtask);

    MPI_Finalize();
}
