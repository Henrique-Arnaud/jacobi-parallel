EXEFILE      = jacobi-solver
EXEFILEPARALLEL      = jacobi-solver_parallel
CPUCC     = $(CC)
MPICC			= mpicc
CPPFLAGS  = -std=c11 -m64 -Ofast -Wall
DEFS      =
INCLUDES  =
LIBDIR    =
LIBS     =  -lm
LINK     =  $(LIBDIR) $(LIBS)

CPU_COMPILE = $(CPUCC) $(DEFS) $(INCLUDES) $(CPPFLAGS)
MPI_COMPILE = $(MPICC) $(DEFS) $(INCLUDES) $(CPPFLAGS)


all: jacobi-solver
	$(CPU_COMPILE)	jacobi-solver.o $(LINK) -o $(EXEFILE) $(PAPILIBS)

jacobi-solver:
	$(CPU_COMPILE) -c jacobi-solver.c

runParallel: jacobi-solver_parallel input
	mpirun -np 4 jacobi-solver_parallel < input

buildParallel: jacobi-solver_parallel
	$(MPI_COMPILE)	jacobi-solver_parallel.o $(LINK) -o $(EXEFILEPARALLEL) $(PAPILIBS)

jacobi-solver_parallel: jacobi-solver_parallel.c
	$(MPI_COMPILE) -c jacobi-solver_parallel.c

clean:
	rm *.o jacobi-solver jacobi-solver_parallel
