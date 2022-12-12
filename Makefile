omp: 
	mpicc -o omp.o omp.c -fopenmp && ./omp.o

mpi: 
	mpicc -o mpi.o mpi.c && mpirun -np 4 ./mpi.o

.PHONY omp mpi