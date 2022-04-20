
FILE = ./main.c ./fun.c ./matrix_exponential.c ./c8lib.c ./r8lib.c
LDLIBS = libsuperlu_mt_OPENMP.a libcxsparse.a libblas_OPENMP.a \
			cblas_LINUX.a blas_LINUX.a libgsl.a libgslcblas.a -lm

all: compile run

compile:
	gcc -fopenmp -ansi -g $(FILE) $(LDLIBS) -o main

run:
	./main

clean:
	rm -rf *.o
