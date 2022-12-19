#include <iostream>
#include "mpi.h"

//// -------------------------------------- MPI INIT  --------------------------------------

//int main(int argc, char **argv)
//{
//    int num_task, rank;
//
//    MPI_Init(&argc, &argv);
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &numtask);
//
//    printf("Hello %d %d \n", rank, num_task);
//
//    MPI_Finalize();
//}

//// -------------------------------------- MPI Send Receive --------------------------------------

//double sliceArea(double h1, double h2, double d) {
//    return (h1 + h2) * d / 2;
//}
//
//double f(double x) {
//    return x * x;
//}
//
//int main(int argc, char **argv)
//{
//    // integral data
//    int d = 9999;
//    int a = 1, b = 3;
//    double dx = double(b - a) / d;
//
//    // MPI
//    MPI_Init(&argc, &argv);
//    MPI_Status status;
//    int rank, size;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    // Processor task limits
//    int chunkSize = d / size;
//    int start = chunkSize * rank;
//    int end = chunkSize * (rank + 1);
//
//    // chunkSize might be rounded down when if d % size != 1
//    if (rank == size - 1)
//    {
//        end = d;
//    }
//
//    // Riman integral
//    double area = 0;
//    for (int i = start; i < end; i++)
//    {
//        double lower = f(a + dx * i);
//        double higher = f(a + dx * (i + 1));
//        area += sliceArea(lower, higher, dx);
//    }
//
//    // Fork Join Model, all results should be sent to 0
//    if (rank != 0)
//    {
//        MPI_Send(&area, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
//        printf("The result from process %d is %f\n", rank, area);
//    }
//    else
//    {
//        for (int i = 0; i < size - 1; i++)
//        {
//            double message;
//            MPI_Recv(&message, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
//            area += message;
//        }
//        printf("The TOTAL integral value is %f", area);
//    }
//    MPI_Finalize();
//    return 0;
//}

//// -------------------------------------- BCast --------------------------------------

//void initData(int n, double* x) {
//    for(int i = 0; i < n; i++) {
//        x[i] = 1.5;
//    }
//}
//
//int main(int argc, char** argv) {
//    // BCast
//    double x[100], totalSum, procSum = 0.0;
//    int procRank, procNum, n = 100;
//    MPI_Status status;
//
//    // init
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
//    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
//
//    // prepare data
//    if (procRank == 0) initData(n, x);
//    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    // calculate partial sum from i1 to i2
//    int k = n / procNum;
//    int i1 = k * procRank;
//    int i2 = k * (procRank+1);
//    for (int i = i1; i < i2; i++) {
//        procSum += x[i];
//    }
//
//    MPI_Reduce(&procSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//    if (procRank == 0) {
//        printf("\nTotal Sum = %10.2f", totalSum);
//    }
//    MPI_Finalize();
//}

//// -------------------------------------- 4.2 -------------------------------------- //

//int main(int argc, char *argv[]) {
//
//    int rank, size;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    int sendarray[100];
//    // init send array
//    for (int i = 0; i < 100; i++) {
//        sendarray[i] = 1000 * rank + i;
//    }
//    int* rbuf = (int *)malloc(size*100*sizeof(int));
//    MPI_Gather(sendarray, 100, MPI_INT, rbuf, 100, MPI_INT, 0, MPI_COMM_WORLD);
//    if(rank == 0){
//        printf("Result in buffer: ");
//        for(int i = 0; i < size*100; i++){
//            printf("%d ", rbuf[i]);
//        }
//        printf("\n");
//    }
//    MPI_Finalize();
//}

//// -------------------------------------- каскадное суммирование -------------------------------------- //

