#include "./include/fun.h"
#include <stdio.h>
#include <string.h>


void add_v_bycol(int length, int col_num, double *v, double *Vm){
	int i = 0;
	for(i=0; i<length; i++){
		Vm[col_num*length + i] = v[i];
	}
}

void get_v(int length, int col_num, double *v, double *Vm){
	int i = 0;
	for(i=0; i<length; i++){
		v[i] = Vm[col_num*length + i];
	}	

}

void set_vector_zero(int length, double *v){
	int i=0;
	for(i=0; i<length; i++){
		v[i] = 0.0;
	}
}

void dense_to_csc(double *dense_matrix, int length_row, int length_col, double *value, int *rowIdx, int *colPtr, int *num_nonzero){
	int i;
	int j;
	int nnz = 0;
	colPtr[0] = 0;
	for(j=0; j<length_col; j++){
		for(i=0; i<length_row; i++){
			if(dense_matrix[j*length_row+i]!=0){
				value[nnz] = dense_matrix[j*length_row+i];
				rowIdx[nnz] = i;
				nnz++;
			}
		}
		colPtr[j+1] = nnz;
	}
	num_nonzero[0] = nnz;
}

/* create an identical matrix with a scale factor of value */
void create_Imatrix(int length, double *matrix, double value){
	int num = 0;
	int i, j;
	for(i=0; i<length; i++){
		for(j=0; j<length; j++){
			if(j==num){
				matrix[i*length+j] = value;
			}else{
				matrix[i*length+j] = 0;
			}
		}
		num++;
	}
}

/* create an idential matrix with a scale factor, in csc format */
void create_Imatrix_csc(int dim, double scale_factor, double *val, int *rowIdx, int* colPtr){
	int i;
	for(i=0; i<dim; i++){
		val[i] = scale_factor;
		rowIdx[i] = i;
		colPtr[i] = i;
	}
	colPtr[dim] = dim;
}


void cs_di_create(int i_nzmax, int i_m, int i_n, int* i_p, int* i_i, double* i_x, int i_nz, cs_di* cs_di_p) {
	cs_di_p->nzmax = i_nzmax;
	cs_di_p->m = i_m;
	cs_di_p->n = i_n;
	cs_di_p->p = i_p;
	cs_di_p->i = i_i;
	cs_di_p->x = i_x;
	cs_di_p->nz = i_nz;  /* # of entries in triplet matrix, -1 for compressed-col */
}

void csc_to_arrayByCol(int nnz, int dim_col, int dim_row, double *value, int* rowInd, int* colPtr, double* array){
	int i;
	int j;
	for(i=0; i<dim_col;i++){
		for(j=0; j<(colPtr[i+1]-colPtr[i]); j++){
			array[i*dim_row+rowInd[colPtr[i]+j]] = value[colPtr[i]+j];
		}
	}

}


void read_file_double(char *addr, double *array, int *size){
	int i = 0;
	char str [100];  /*the length of each row should smaller than row_length*/
	FILE * fp = fopen (addr, "r");
	if(fp==NULL){
		printf("no such file");
	}else{
		while(!feof(fp)){
			fgets(str, 100, fp);
			array[i] = atof(str);
			i++;
		}
	}
	*size = i-1;
	fclose(fp);
	fp = NULL;
}

void read_file_int(char *addr, int *array, int *size){
	int i = 0;
	char str [100];  /*the length of each row should smaller than row_length*/
	FILE * fp = fopen (addr, "r");
	if(fp==NULL){
		printf("no such file");
	}else{
		while(!feof(fp)){
			fgets(str, 100, fp);
			array[i] = atoi(str);
			i++;
		}
	}
	*size = i-1;
	fclose(fp);
	fp = NULL;
}




