#include <iostream>
#include "omp.h"

using namespace std;

//// --------------------------------- MAIN  ---------------------------------

int main(){
    cout << "Hello World" << "\n";
    cout << omp_get_num_procs();
}

//// --------------------------------- OMP Init ---------------------------------

//int main()
//{
//    cout << omp_get_num_procs() << endl;
//    int a[100], b[100], c[100], d[100];
//
//#pragma omp parallel for
//    for (int i = 0; i < 100; ++i)
//    {
//        a[i] = i;
//        b[i] = i;
//        d[i] = 2 * i;
//        cout << omp_get_thread_num();
//        cout << omp_get_num_threads();
//    }
//
//#pragma omp parallel for
//    for (int i = 0; i < 100; ++i)
//    {
//        c[i] = a[i] + b[i] + d[i];
//        cout << omp_get_thread_num();
//        cout << omp_get_num_threads();
//    }
//    cout << a[50] << endl;
//}

//// TODO: ADD TASK NAME

// int main()
// {
//     string hw = "hellow world \n";
// #pragma omp parallel
//     cout << hw;
//     return 0;
// }

// int main()
// {
//     int x = 0;
// #pragma omp parallel shared(x) num_threads(5000)
//     x += 1;
//     cout << "x= " << x << endl;
//     return 0;
// }

// int main(int argc, char *argv[])
// {
//     printf("posledovatelniy1\n");
// #pragma omp parallel
//     {
//         printf("parallel\n");
//     }
//     printf("posledovatelniy2\n");
//     return 0;
// }

// int main()
// {
//     int n = 1;
//     printf("n in sequence of fields (started): %d\n", n);
// #pragma omp parallel private(n) num_threads(4)
//     {
//         printf("Value n in thread(on entrance): %d\n", n);
//         int k = 1;
//         printf("Value of k: %d\n", k);
//         n = omp_get_thread_num();
//         printf("Value n in thread(on exit): %d\n", n);
//     }

//     printf("n in sequence of fields (finished): %d\n", n);
//     return 0;
// }

//// --------------------------------- IJK Matrix Multiplication ---------------------------------

//#include <chrono>
//const int n = 3;
//const int trails = 1e5;
//int a[n][n];
//int b[n][n];
//
//void setZero(int *array) {
//    for (int i = 0; i < n * n; i++)
//        array[i] = 0;
//}
//
//void setRandom(int *array) {
//    for (int i = 0; i < n * n; i++)
//        array[i] = (rand() % 1000000);
//}
//
//int calculateElapsed(void (*func)(int[n][n])) {
//    int c[n][n];
//    setZero((int *) c);
//    auto start = std::chrono::system_clock::now();
//
//    func(c);
//
//    auto end = std::chrono::system_clock::now();
//    return (int) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//}
//
//void ijk(int c[n][n]) {
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            for (int k = 0; k < n; k++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void jik(int c[n][n]) {
//    for (int j = 0; j < n; j++) {
//        for (int i = 0; i < n; i++) {
//            for (int k = 0; k < n; k++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void kji(int c[n][n]) {
//    for (int k = 0; k < n; k++) {
//        for (int j = 0; j < n; j++) {
//            for (int i = 0; i < n; i++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void ikj(int c[n][n]) {
//    for (int i = 0; i < n; i++) {
//        for (int k = 0; k < n; k++) {
//            for (int j = 0; j < n; j++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void jki(int c[n][n]) {
//    for (int j = 0; j < n; j++) {
//        for (int k = 0; k < n; k++) {
//            for (int i = 0; i < n; i++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void kij(int c[n][n]) {
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            for (int k = 0; k < n; k++) {
//                c[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
//void elapsedTime(const string &algorithm_name, void (*func)(int[n][n])) {
//    int time_elapsed = 0;
//    for (size_t i = 0; i < trails; i++)
//        time_elapsed += calculateElapsed(func);
//    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
//}

//
//int main() {
//
//    // set random values
//    setRandom((int *) a);
//    setRandom((int *) b);
//
//    // call functions
//    elapsedTime("ijk", ijk);
//    elapsedTime("ikj", ikj);
//    elapsedTime("jik", jik);
//    elapsedTime("jki", jki);
//    elapsedTime("kji", kji);
//    elapsedTime("kij", kij);
//}

//// --------------------------------- Vector Multiplication ---------------------------------

//#include <chrono>
//#include <cmath>
//
// // CONSTANTS
//const int n = 1000;
//const int trails = 1;
//
// // 2D Arrays
//int a[n][n];
//int b[n][n];
//
// // 1D Arrays
//int c[n];
//int d[n];
//
// // -------------------- ARRAY FILLER FUNCTIONS --------------------
//
//void set_zero_by_parallel(int *array, int dimension) {
//#pragma omp for shared(array)
//    {
//        for (int i = 0; i < pow(n, dimension); i++)
//            array[i] = 0;
//    }
//}
//
//void set_random_by_parallel(int *array, int dimension) {
//#pragma omp for shared(array)
//    {
//        for (int i = 0; i < pow(n, dimension); i++)
//            array[i] = (rand() % 1000000);
//    }
//}
//
//
//void set_zero_by_sequence(int *array, int dimension) {
//    for (int i = 0; i < pow(n, dimension); i++)
//        array[i] = 0;
//}
//
//void set_random_by_sequence(int *array, int dimension) {
//    for (int i = 0; i < pow(n, dimension); i++)
//        array[i] = (rand() % 1000000);
//}
//
// // -------------------- PARALLEL CALCULATIONS --------------------
//
//void multiply_two_vectors_by_parallel() {
//    int r[n];
//#pragma omp for
//    {
//        for (int i = 0; i < n; ++i) {
//            r[i] = c[i] * d[i];
//        }
//    }
//}
//
//void multiply_two_matrix_by_parallel() {
//    int r[n][n];
//    set_zero_by_parallel(*r, 2);
//#pragma omp for
//    {
//        for (int i = 0; i < n; i++) {
//            for (int k = 0; k < n; k++) {
//                for (int j = 0; j < n; j++) {
//                    r[i][j] += a[i][k] * b[k][j];
//                }
//            }
//        }
//    }
//}
//
//
// // -------------------- SEQUENCE CALCULATIONS --------------------
//
//void multiply_two_vectors_by_sequence() {
//    int r[n];
//    for (int i = 0; i < n; ++i) {
//        r[i] = c[i] * d[i];
//    }
//}
//
//void multiply_two_matrix_by_sequence() {
//    int r[n][n];
//    set_zero_by_sequence(*r, 2);
//    for (int i = 0; i < n; i++) {
//        for (int k = 0; k < n; k++) {
//            for (int j = 0; j < n; j++) {
//                r[i][j] += a[i][k] * b[k][j];
//            }
//        }
//    }
//}
//
// // --------------------  HELPERS  --------------------
//
//int calculate_time_elapsed(void (*func)()) {
//    auto start = std::chrono::system_clock::now();
//
//    func();
//
//    auto end = std::chrono::system_clock::now();
//    return (int) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//}
//
//void elapsed_time_vector_multiplication(const string &algorithm_name, void (*func)()) {
//    int time_elapsed = 0;
//    for (size_t i = 0; i < trails; i++)
//        time_elapsed += calculate_time_elapsed(func);
//    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
//}
//
//void elapsed_time_matrix_multiplication(const string &algorithm_name, void (*func)()) {
//    int time_elapsed = 0;
//    for (size_t i = 0; i < trails; i++)
//        time_elapsed += calculate_time_elapsed(func);
//    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
//}
// // --------------------  MAIN  --------------------
//
//int main() {
//
//    set_random_by_parallel(*a, 2);
//    set_random_by_parallel(*b, 2);
//    set_random_by_parallel(c, 1);
//    set_random_by_parallel(d, 1);
//
//    elapsed_time_vector_multiplication("VM By SEQ : ", multiply_two_vectors_by_sequence);
//    elapsed_time_vector_multiplication("VM By PAR : ", multiply_two_vectors_by_parallel);
//    elapsed_time_vector_multiplication("MM By SEQ : ", multiply_two_matrix_by_sequence);
//    elapsed_time_vector_multiplication("MM By PAR : ", multiply_two_matrix_by_parallel);
//
//    return 0;
//}

//// --------------------------------- OMP Types ---------------------------------

//#include <iostream>
//#include <omp.h>
//#include <ctime>
//#include <vector>
//#include <cstdio>
//#include <valarray>
//
//double Sequential(int Size, int** MatrixA, int** MatrixB, int** MatrixC)
//{
//    double start = omp_get_wtime();
//    for (int i = 0; i < Size; i++)
//    {
//        for (int j = 0; j < Size; j++)
//        {
//            for (int k = 0; k < Size; k++)
//            {
//                MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
//            }
//        }
//    }
//    return omp_get_wtime() - start;
//}
//
//double Parallel(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
//{
//    double start = omp_get_wtime();
//    int i, j, k;
//#pragma omp parallel for private (j, k)
//    for (i = 0; i < Size; i++)
//    {
//        for (j = 0; j < Size; j++)
//        {
//            for (k = 0; k < Size; k++)
//            {
//                MatrixC[i * Size + j] += MatrixA[i * Size + k] *
//                                         MatrixB[k * Size + j];
//            }
//        }
//    }
//    return omp_get_wtime() - start;
//}
//
//double TapeSeparation(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
//{
//    double start = omp_get_wtime();
//    int i, j, k;
//    int NestedThreadsNum = 2;
//    omp_set_nested(true);
//    omp_set_num_threads(NestedThreadsNum);
//#pragma omp parallel for private (j, k)
//    for (i = 0; i < Size; i++)
//#pragma omp parallel for private (k)
//            for (j = 0; j < Size; j++)
//                for (k = 0; k < Size; k++)
//                    MatrixC[i * Size + j] += MatrixA[i * Size + k] *
//                                             MatrixB[k * Size + j];
//
//    return omp_get_wtime() - start;
//}
//
//double Block(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
//{
//    double start = omp_get_wtime();
//    int ThreadNum = 4;
//    int GridSize = int(sqrt((double)ThreadNum));
//    int BlockSize = Size / GridSize;
//    omp_set_num_threads(ThreadNum);
//#pragma omp parallel
//    {
//        int ThreadID = omp_get_thread_num();
//        int RowIndex = ThreadID / GridSize;
//        int ColIndex = ThreadID % GridSize;
//        for (int iter = 0; iter < GridSize; iter++) {
//            for (int i = RowIndex * BlockSize;
//                 i < (RowIndex + 1) * BlockSize; i++)
//                for (int j = ColIndex * BlockSize;
//                     j < (ColIndex + 1) * BlockSize; j++)
//                    for (int k = iter * BlockSize;
//                         k < (iter + 1) * BlockSize; k++)
//                        MatrixB[i * Size + j] +=
//                                MatrixA[i * Size + k] * MatrixC[k * Size + j];
//        }
//    }
//    return omp_get_wtime() - start;
//}
//
//double BlockCache(int Size, int* MatrixA, int* MatrixB, int* MatrixC)
//{
//    double start = omp_get_wtime();
//    int BlockSize = 250;
//    int GridSize = int(Size / double(BlockSize));
//#pragma omp parallel for
//    for (int n = 0; n < GridSize; n++)
//        for (int m = 0; m < GridSize; m++)
//            for (int iter = 0; iter < GridSize; iter++)
//                for (int i = n * BlockSize;i < (n + 1) * BlockSize; i++)
//                    for (int j = m * BlockSize;j < (m + 1) * BlockSize; j++)
//                        for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++)
//                            MatrixC[i * Size + j] +=
//                                    MatrixA[i * Size + k] * MatrixB[k * Size + j];\
//  return omp_get_wtime() - start;
//}
//
//
//int main()
//{
//    const int Size = 10;
//    int* MatrixA;
//    int* MatrixB;
//    int* MatrixC;
//
//    MatrixA = new int[Size * Size];
//    MatrixB = new int[Size * Size];
//    MatrixC = new int[Size * Size];
//    for (int i = 0; i < Size; i++)
//        for (int j = 0; j < Size; j++)
//        {
//            MatrixC[i * Size + j] = 0;
//            MatrixA[i * Size + j] = rand() % 100;
//            MatrixB[i * Size + j] = rand() % 100;
//        }
//
//    std::cout << "Parallel: " << Parallel(Size, MatrixA, MatrixB, MatrixC) << "\n";
//    std::cout << "TapeSeparation: " << TapeSeparation(Size, MatrixA, MatrixB, MatrixC) << "\n";
//    std::cout << "Block: " << Block(Size, MatrixA, MatrixB, MatrixC) << "\n";
//    std::cout << "BlockCache: " << BlockCache(Size, MatrixA, MatrixB, MatrixC) << "\n";
//}


//int main()
//{
//    printf("Значение OMP_DYNAMIC: %d\n", omp_get_dynamic());
//    omp_set_dynamic(1);
//    printf("Значение OMP_DYNAMIC: %d\n", omp_get_dynamic());
//#pragma omp parallel num_threads(128)
//    {
//#pragma omp master
//        {
//            printf("Параллельная область, %d нитей\n",
//                   omp_get_num_threads());
//        }
//    }
//}

//// --------------------------------- OMP Parallel ---------------------------------

//#include <cstdlib>
//
//void fill(int *array, int size)
//{
//    for (int i = 0; i < size; ++i)
//    {
//        array[i] = rand();
//    }
//}
//
//void print(int *array, int size)
//{
//    for (int i = 0; i < size; ++i)
//    {
//        std::cout << array[i];
//    }
//}
//
//int main(int argc, char *argv[])
//{
//    int N = 100;
//    int a[N];
//    int b[N];
//    int c[N];
//    int d[N];
//
//    fill(b, N);
//    fill(c, N);
//
//#pragma omp parallel
//    {
//        for (int i = 0; i < N; i++)
//            a[i] = b[i] + c[i];
//#pragma omp barrier
//        for (int i = 0; i < N; i++)
//            d[i] = a[i] + b[i];
//    }
//
//    print(a, N);
//    print(d, N);
//}

//#include <stdio.h>
//
//int main(int argc, char *argv[])
//{
//    int count = 0;
//#pragma omp parallel reduction (+: count)`
//    {
//        count++;
//        printf("Текущее значение count: %d\n", count);
//    }
//    printf("Число нитей: %d\n", count);
//}

//#include <stdio.h>
//#include <omp.h>
//int main(int argc, char *argv[])
//{
//#pragma omp parallel
//    {
//        printf("Сообщение 1\n");
//        printf("Сообщение 2\n");
//#pragma omp barrier
//        printf("Сообщение 3\n");
//    }
//}

//#include <stdio.h>
//#include <omp.h>
//int main(int argc, char *argv[])
//{
//    int n;
//#pragma omp parallel
//    {
//#pragma omp critical
//        {
//            n=omp_get_thread_num();
//            printf("Нить %d\n", n);
//        }
//    }
//}

//#include <stdio.h>
//#include <omp.h>
//int main(int argc, char *argv[])
//{
//    int count = 0;
//#pragma omp parallel
//    {
//#pragma omp atomic
//        count++;
//    }
//    printf("Число нитей: %d\n", count);
//}

//// --------------------------------- OMP Course Task ---------------------------------

// int main()
// {
//    int n=0;
//    omp_set_nested(1);
// #pragma omp parallel num_threads(2) reduction(+:n)
//     {
//	n = 0;
//#pragma omp parallel num_threads(2) reduction(+:n)
//	{
//	    n = 1;
//	}
//    }
//    printf("%d ", n);
//    n = 0;
//    omp_set_nested(0);
// #pragma omp parallel num_threads(2) reduction(+:n)
//    {
//	n = 0;
// #pragma omp parallel num_threads(2) reduction(+:n)
//	{
//	    n = 1;
//	}
//    }
//    printf("%d ", n);
// }

//// --------------------------------- OMP Alpha Section ---------------------------------

//int main() {
//    // инициализация
//    int n = 100, rmin = 0, rmax = 1;
//    double h = double(rmax - rmin) / n;
//    double r[n+1], y[n+1], alpha[n+1], beta[n+1], a[n+1], b[n+1], c[n+1], f[n+1], F[n+1], D[n+1];
//
//    // определим r, D, F
//#pragma omp parallel for
//    for (int i = 0; i <= n; i++) {
//        r[i] = h * i + rmin;
//        D[i] = r[i]; // потому что D(r) = r
//        F[i] = 4 * r[i]; // потому что F(r) = 4 * r
//    }
//    // определим a, b, c, f
//#pragma omp parallel sections
//{
//#pragma omp section
//        {
//            for (int i = 0; i <= n; i++) {
//                a[i] = D[i];
//                f[i] = -1 * F[i] * h * h;
//            }
//        }
//#pragma omp section
//        {
//            for (int i = 0; i <= n - 1; i++) {
//                b[i] = D[i + 1];
//                c[i] = D[i] + D[i + 1];
//            }
//        }
//}
// //    alpha[1] = 2 * r[1]; // 0
// //    beta[1] = -1 * (alpha[1] * r[1] * r[1]); // 0
//    alpha[1] = 0;
//    beta[1] = 0;
//    for(int i = 1; i <= n-1; i++) {
//        alpha[i+1] = b[i] / (c[i] - a[i] * alpha[i]);
//        beta[i+1] = (a[i]*beta[i] + f[i]) / (c[i] - a[i]*alpha[i]);
//    }
//    y[n] = 1;
//    for(int i = n-1; i >= 0; i--) {
//        y[i] = alpha[i+1] * y[i+1] + beta[i+1];
//    }
//
//    for (int i = 0; i <= n; i++) {
//        printf("\n y[%f] = %f",r[i]/h, y[i]);
//    }
//
//    return 0;
//}
