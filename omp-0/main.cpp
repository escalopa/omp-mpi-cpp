#include "omp.h"
using namespace std;

int main()
{
    cout << omp_get_num_procs() << endl;
    int a[100], b[100], c[100], d[100];

#pragma omp parallel for
    for (int i = 0; i < 100; ++i)
    {
        a[i] = i;
        b[i] = i;
        d[i] = 2 * i;
        cout << omp_get_thread_num();
        cout << omp_get_num_threads();
    }

#pragma omp parallel for
    for (int i = 0; i < 100; ++i)
    {
        c[i] = a[i] + b[i] + d[i];
        cout << omp_get_thread_num();
        cout << omp_get_num_threads();
    }
    cout << a[50] << endl;
}

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
