#include <stdio.h>
#include "./cs.h"

void add_v_bycol(int length, int col_num, double *v, double *Vm);
void get_v(int length, int col_num, double *v, double *Vm);

void set_vector_zero(int length, double *v);

void dense_to_csc(double *dense_matrix_byrow, int length_row, int length_col, double *value, int *rowIdx, int *colPtr, int *num_nonzero);

void create_Imatrix(int length, double *matrix, double value);

void cs_di_create(int i_nzmax, int i_m, int i_n, int* i_p, int* i_i, double* i_x, int i_nz, cs_di* cs_di_p);

void csc_to_arrayByCol(int nnz, int dim_col, int dim_row, double *value, int* rowInd, int* colPtr, double* array);

void read_file_double(char *addr, double *array, int *size);
void read_file_int(char *addr, int *array, int *size);

void create_Imatrix_csc(int dim, double scale_factor, double *val, int *rowIdx, int* colPtr);