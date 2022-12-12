#include <cstdio>
int main() {
    // инициализация
    int n = 100, rmin = 0, rmax = 1;
    double h = double(rmax - rmin) / n;
    double r[n+1], y[n+1], alpha[n+1], beta[n+1], a[n+1], b[n+1], c[n+1], f[n+1], F[n+1], D[n+1];

    // определим r, D, F
#pragma omp parallel for
    for (int i = 0; i <= n; i++) {
        r[i] = h * i + rmin;
        D[i] = r[i]; // потому что D(r) = r
        F[i] = 4 * r[i]; // потому что F(r) = 4 * r
    }
    // определим a, b, c, f
#pragma omp parallel sections
{
#pragma omp section
        {
            for (int i = 0; i <= n; i++) {
                a[i] = D[i];
                f[i] = -1 * F[i] * h * h;
            }
        }
#pragma omp section
        {
            for (int i = 0; i <= n - 1; i++) {
                b[i] = D[i + 1];
                c[i] = D[i] + D[i + 1];
            }
        }
}
// TODO: граничный условие alpha, beta не уверени но все вариант работают
//    alpha[1] = 2 * r[1]; // 0
//    beta[1] = -1 * (alpha[1] * r[1] * r[1]); // 0
    alpha[1] = 0;
    beta[1] = 0;
    for(int i = 1; i <= n-1; i++) {
        alpha[i+1] = b[i] / (c[i] - a[i] * alpha[i]);
        beta[i+1] = (a[i]*beta[i] + f[i]) / (c[i] - a[i]*alpha[i]);
    }
    y[n] = 1;
    for(int i = n-1; i >= 0; i--) {
        y[i] = alpha[i+1] * y[i+1] + beta[i+1];
    }

    for (int i = 0; i <= n; i++) {
        printf("\n y[%f] = %f",r[i]/h, y[i]);
    }

    return 0;
}
