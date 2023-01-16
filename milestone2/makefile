MyAdvectionSerial: my_advection_program.h my_advection_program.c
	g++ -fopenmp -O3 -flto -march=native -mtune=native my_advection_program.c -o my_advection_program_serial

MyAdvection: my_advection_program.h my_advection_program.c
	g++ -fopenmp -DisParallel -O3 -flto -march=native -mtune=native my_advection_program.c -o my_advection_program_parallel