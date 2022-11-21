#include <iostream>
#include <chrono>

using namespace std;
const int n = 3;
const int trails = 1e5;
int a[n][n];
int b[n][n];

void setZero(int *array) {
    for (int i = 0; i < n * n; i++)
        array[i] = 0;
}

void setRandom(int *array) {
    for (int i = 0; i < n * n; i++)
        array[i] = (rand() % 1000000);
}

int calculateElapsed(void (*func)(int[n][n])) {
    int c[n][n];
    setZero((int *) c);
    auto start = std::chrono::system_clock::now();

    func(c);

    auto end = std::chrono::system_clock::now();
    return (int) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}

void ijk(int c[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void jik(int c[n][n]) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void kji(int c[n][n]) {
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void ikj(int c[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void jki(int c[n][n]) {
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void kij(int c[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void elapsedTime(const string &algorithm_name, void (*func)(int[n][n])) {
    int time_elapsed = 0;
    for (size_t i = 0; i < trails; i++)
        time_elapsed += calculateElapsed(func);
    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
}

template<class T>
void quickSort(T *a, long n) {

}

int main() {

    // set random values
    setRandom((int *) a);
    setRandom((int *) b);

    // call functions
    elapsedTime("ijk", ijk);
    elapsedTime("ikj", ikj);
    elapsedTime("jik", jik);
    elapsedTime("jki", jki);
    elapsedTime("kji", kji);
    elapsedTime("kij", kij);
}
