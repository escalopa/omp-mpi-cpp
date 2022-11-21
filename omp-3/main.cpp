#include <iostream>
#include <omp.h>
#include <ctime>
#include <vector>
#include <cstdio>
#include <valarray>

double Sequential(int Size, int** MatrixA, int** MatrixB, int** MatrixC)
{
    double start = omp_get_wtime();
    for (int i = 0; i < Size; i++)
    {
        for (int j = 0; j < Size; j++)
        {
            for (int k = 0; k < Size; k++)
            {
                MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
    return omp_get_wtime() - start;
}

double Parallel(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
{
    double start = omp_get_wtime();
    int i, j, k;
#pragma omp parallel for private (j, k)
    for (i = 0; i < Size; i++)
    {
        for (j = 0; j < Size; j++)
        {
            for (k = 0; k < Size; k++)
            {
                MatrixC[i * Size + j] += MatrixA[i * Size + k] *
                                         MatrixB[k * Size + j];
            }
        }
    }
    return omp_get_wtime() - start;
}

double TapeSeparation(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
{
    double start = omp_get_wtime();
    int i, j, k;
    int NestedThreadsNum = 2;
    omp_set_nested(true);
    omp_set_num_threads(NestedThreadsNum);
#pragma omp parallel for private (j, k)
    for (i = 0; i < Size; i++)
#pragma omp parallel for private (k)
            for (j = 0; j < Size; j++)
                for (k = 0; k < Size; k++)
                    MatrixC[i * Size + j] += MatrixA[i * Size + k] *
                                             MatrixB[k * Size + j];

    return omp_get_wtime() - start;
}

double Block(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
{
    double start = omp_get_wtime();
    int ThreadNum = 4;
    int GridSize = int(sqrt((double)ThreadNum));
    int BlockSize = Size / GridSize;
    omp_set_num_threads(ThreadNum);
#pragma omp parallel
    {
        int ThreadID = omp_get_thread_num();
        int RowIndex = ThreadID / GridSize;
        int ColIndex = ThreadID % GridSize;
        for (int iter = 0; iter < GridSize; iter++) {
            for (int i = RowIndex * BlockSize;
                 i < (RowIndex + 1) * BlockSize; i++)
                for (int j = ColIndex * BlockSize;
                     j < (ColIndex + 1) * BlockSize; j++)
                    for (int k = iter * BlockSize;
                         k < (iter + 1) * BlockSize; k++)
                        MatrixB[i * Size + j] +=
                                MatrixA[i * Size + k] * MatrixC[k * Size + j];
        }
    }
    return omp_get_wtime() - start;
}

double BlockCache(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
{
    double start = omp_get_wtime();
    int BlockSize = 250;
    int GridSize = int(Size / double(BlockSize));
#pragma omp parallel for
    for (int n = 0; n < GridSize; n++)
        for (int m = 0; m < GridSize; m++)
            for (int iter = 0; iter < GridSize; iter++)
                for (int i = n * BlockSize;i < (n + 1) * BlockSize; i++)
                    for (int j = m * BlockSize;j < (m + 1) * BlockSize; j++)
                        for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++)
                            MatrixC[i * Size + j] +=
                                    MatrixA[i * Size + k] * MatrixB[k * Size + j];\
  return omp_get_wtime() - start;
}


int main()
{
    const int Size = 10;
    int* MatrixA;
    int* MatrixB;
    int* MatrixC;

    MatrixA = new int[Size * Size];
    MatrixB = new int[Size * Size];
    MatrixC = new int[Size * Size];
    for (int i = 0; i < Size; i++)
        for (int j = 0; j < Size; j++)
        {
            MatrixC[i * Size + j] = 0;
            MatrixA[i * Size + j] = rand() % 100;
            MatrixB[i * Size + j] = rand() % 100;
        }

    std::cout << "Parallel: " << Parallel(Size, MatrixA, MatrixB, MatrixC) << "\n";
    std::cout << "TapeSeparation: " << TapeSeparation(Size, MatrixA, MatrixB, MatrixC) << "\n";
    std::cout << "Block: " << Block(Size, MatrixA, MatrixB, MatrixC) << "\n";
    std::cout << "BlockCache: " << BlockCache(Size, MatrixA, MatrixB, MatrixC) << "\n";
}