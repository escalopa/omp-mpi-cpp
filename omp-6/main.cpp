 #include <stdio.h>
 #include <omp.h>

 int main()
 {
    int n=0;
    omp_set_nested(1);
 #pragma omp parallel num_threads(2) reduction(+:n)
     {
	n = 0;
#pragma omp parallel num_threads(2) reduction(+:n)
	{
	    n = 1;
	}
    }
    printf("%d ", n);
    n = 0;
    omp_set_nested(0);
 #pragma omp parallel num_threads(2) reduction(+:n)
    {
	n = 0;
 #pragma omp parallel num_threads(2) reduction(+:n)
	{
	    n = 1;
	}
    }
    printf("%d ", n);
 }