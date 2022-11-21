#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

// CONSTANTS
const int n = 1000;
const int trails = 1;

// 2D Arrays
int a[n][n];
int b[n][n];

// 1D Arrays
int c[n];
int d[n];

// -------------------- ARRAY FILLER FUNCTIONS --------------------

void set_zero_by_parallel(int *array, int dimension) {
#pragma omp for shared(array)
    {
        for (int i = 0; i < pow(n, dimension); i++)
            array[i] = 0;
    }
}

void set_random_by_parallel(int *array, int dimension) {
#pragma omp for shared(array)
    {
        for (int i = 0; i < pow(n, dimension); i++)
            array[i] = (rand() % 1000000);
    }
}


void set_zero_by_sequence(int *array, int dimension) {
    for (int i = 0; i < pow(n, dimension); i++)
        array[i] = 0;
}

void set_random_by_sequence(int *array, int dimension) {
    for (int i = 0; i < pow(n, dimension); i++)
        array[i] = (rand() % 1000000);
}

// -------------------- PARALLEL CALCULATIONS --------------------

void multiply_two_vectors_by_parallel() {
    int r[n];
#pragma omp for
    {
        for (int i = 0; i < n; ++i) {
            r[i] = c[i] * d[i];
        }
    }
}

void multiply_two_matrix_by_parallel() {
    int r[n][n];
    set_zero_by_parallel(*r, 2);
#pragma omp for
    {
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    r[i][j] += a[i][k] * b[k][j];
                }
            }
        }
    }
}


// -------------------- SEQUENCE CALCULATIONS --------------------

void multiply_two_vectors_by_sequence() {
    int r[n];
    for (int i = 0; i < n; ++i) {
        r[i] = c[i] * d[i];
    }
}

void multiply_two_matrix_by_sequence() {
    int r[n][n];
    set_zero_by_sequence(*r, 2);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

// --------------------  HELPERS  --------------------

int calculate_time_elapsed(void (*func)()) {
    auto start = std::chrono::system_clock::now();

    func();

    auto end = std::chrono::system_clock::now();
    return (int) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}

void elapsed_time_vector_multiplication(const string &algorithm_name, void (*func)()) {
    int time_elapsed = 0;
    for (size_t i = 0; i < trails; i++)
        time_elapsed += calculate_time_elapsed(func);
    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
}

void elapsed_time_matrix_multiplication(const string &algorithm_name, void (*func)()) {
    int time_elapsed = 0;
    for (size_t i = 0; i < trails; i++)
        time_elapsed += calculate_time_elapsed(func);
    cout << algorithm_name << " " << (double) time_elapsed / (double) trails << " nano" << endl;
}
// --------------------  MAIN  --------------------

int main() {

    set_random_by_parallel(*a, 2);
    set_random_by_parallel(*b, 2);
    set_random_by_parallel(c, 1);
    set_random_by_parallel(d, 1);

    elapsed_time_vector_multiplication("VM By SEQ : ", multiply_two_vectors_by_sequence);
    elapsed_time_vector_multiplication("VM By PAR : ", multiply_two_vectors_by_parallel);
    elapsed_time_vector_multiplication("MM By SEQ : ", multiply_two_matrix_by_sequence);
    elapsed_time_vector_multiplication("MM By PAR : ", multiply_two_matrix_by_parallel);

    return 0;
}
