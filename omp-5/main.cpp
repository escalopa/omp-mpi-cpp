#include <stdio.h>
#include <omp.h>
#include <cstdlib>
#include <iostream>

void fill(int *array, int size)
{
    for (int i = 0; i < size; ++i)
    {
        array[i] = rand();
    }
}

void print(int *array, int size)
{
    for (int i = 0; i < size; ++i)
    {
        std::cout << array[i];
    }
}

int main(int argc, char *argv[])
{
    int N = 100;
    int a[N];
    int b[N];
    int c[N];
    int d[N];

    fill(b, N);
    fill(c, N);

#pragma omp parallel
    {
        for (int i = 0; i < N; i++)
            a[i] = b[i] + c[i];
#pragma omp barrier
        for (int i = 0; i < N; i++)
            d[i] = a[i] + b[i];
    }

    print(a, N);
    print(d, N);
} 

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
