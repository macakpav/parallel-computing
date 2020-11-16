omp1: omp_bug1.c
	gcc -o omp1.exe omp_bug1.c -fopenmp
	@ ./omp1.exe
omp2: omp_bug2.c
	gcc -o omp2.exe omp_bug2.c -fopenmp
	@ ./omp2.exe
omp3: omp_bug3.c
	gcc -o omp3.exe omp_bug3.c -fopenmp
	@ ./omp3.exe
omp4: omp_bug4.c
	gcc -o omp4.exe omp_bug4.c -fopenmp
	@ ./omp4.sh
omp5: omp_bug5.c
	gcc -o omp5.exe omp_bug5.c -fopenmp
	@ ./omp5.exe
omp6: omp_bug6.c
	gcc -o omp6.exe omp_bug6.c -fopenmp
	@ ./omp6.exe
mpi1: mpi_bug1.c
	mpicc -o mpi1.exe mpi_bug1.c 
	@ mpiexec -np 2 ./mpi1.exe
mpi2: mpi_bug2.c
	mpicc -o mpi2.exe mpi_bug2.c 
	@ mpiexec -np 2  ./mpi2.exe
mpi3: mpi_bug3.c
	mpicc -o mpi3.exe mpi_bug3.c 
	@ mpiexec -np 4  ./mpi3.exe
mpi4: mpi_bug4.c
	mpicc -o mpi4.exe mpi_bug4.c 
	@ mpiexec -np 8  ./mpi4.exe
mpi5: mpi_bug5.c
	mpicc -o mpi5.exe mpi_bug5.c 
	@ mpiexec -np 2  ./mpi5.exe
mpi6: mpi_bug6.c
	mpicc -o mpi6.exe mpi_bug6.c 
	@ mpiexec -np 4  ./mpi6.exe
mpi7: mpi_bug7.c
	mpicc -o mpi7.exe mpi_bug7.c 
	@ mpiexec -np 4  ./mpi7.exe