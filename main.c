#include <stdio.h>
#include <stdlib.h>
#include "./include/cblas.h"
#include "./include/slu_mt_ddefs.h"
#include "./include/fun.h"
#include "./include/cs.h"
#include "./include/matrix_exponential.h"

typedef enum
{
    True=1, False=0
}bool;

int main(){
	/*=====================parameter set========================*/
	double t_step = 1.0;
	double gamma = t_step/10;
	double eps = 1e-12;
	int krylov_dim = 50;
	bool converged = False;
	double Cmin = 1e-16;
	int nprocs = 1;   /*number of processor*/

	int i=0;
	int j = 0;

/*=====================input value set======================*/
	/*
	double G_val[13] ={1,-2,-4,-1,5,8,4,2,-3,6,7,4,-5}; 
	int G_rowIdx[13] = {0,1,3,0,1,4,2,3,0,2,3,2,4}; 
	int G_colPtr[6] = {0,3,6,8,11,13}; 

	double C_val[11] ={ 1,1,1,2,8,3,6,6,2,8,1 }; 
	int C_rowIdx[11] = { 0,1,0,1,4,2,3,2,3,1,4 }; 
	int C_colPtr[6] = { 0,2,5,7,9,11 }; 

	double u_t[5] ={0,0,0,0,1}; 
	double u_th[5] ={1,0,0,0,0}; 

	double x_t[5] = {1,0,0,0,0};   
*/
	/* G */
	int data_size = 60000;
	double *Gval;
	int *GrowIdx, *GcolPtr, *Gval_size, *GrowIdx_size, *GcolPtr_size;
	if (!( Gval= doubleMalloc(170000))) SUPERLU_ABORT("Malloc fails for Gval[].");
	if (!( GrowIdx= intMalloc(data_size))) SUPERLU_ABORT("Malloc fails for GrowIdx[].");
	if (!( GcolPtr= intMalloc(data_size))) SUPERLU_ABORT("Malloc fails for GcolPtr[].");
	if (!( Gval_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for Gval_size[].");
	if (!( GrowIdx_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for GrowIdx_size[].");
	if (!( GcolPtr_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for GcolPtr_size[].");
	
	read_file_double("/home/jerry/ei/csc_data/example_ibmpg1t/Gval.txt", Gval, Gval_size);
	read_file_int("/home/jerry/ei/csc_data/example_ibmpg1t/GrowIdx.txt", GrowIdx, GrowIdx_size);
	read_file_int("/home/jerry/ei/csc_data/example_ibmpg1t/GcolPtr.txt", GcolPtr, GcolPtr_size);
	
	/* C */
	double *Cval;
	int *CrowIdx, *CcolPtr, *Cval_size, *CrowIdx_size, *CcolPtr_size;
	if (!( Cval= doubleMalloc(data_size))) SUPERLU_ABORT("Malloc fails for Cval[].");
	if (!( CrowIdx= intMalloc(data_size))) SUPERLU_ABORT("Malloc fails for CrowIdx[].");
	if (!( CcolPtr= intMalloc(data_size))) SUPERLU_ABORT("Malloc fails for CcolPtr[].");
	if (!( Cval_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for Cval_size[].");
	if (!( CrowIdx_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for CrowIdx_size[].");
	if (!( CcolPtr_size= intMalloc(1))) SUPERLU_ABORT("Malloc fails for CcolPtr_size[].");

	read_file_double("/home/jerry/ei/csc_data/example_ibmpg1t/Cval.txt", Cval, Cval_size);
	read_file_int("/home/jerry/ei/csc_data/example_ibmpg1t/CrowIdx.txt", CrowIdx, CrowIdx_size);
	read_file_int("/home/jerry/ei/csc_data/example_ibmpg1t/CcolPtr.txt", CcolPtr, CcolPtr_size);

	int matrix_dim = *CcolPtr_size - 1;   /*matrix dimension of C, G*/

	double *Cmin_val;
	int *Cmin_rowIdx, *Cmin_colPtr;
	if (!(Cmin_val= doubleMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for Cmin_val[].");
	if (!(Cmin_rowIdx = intMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for Cmin_rowIdx[].");
	if (!(Cmin_colPtr = intMalloc(matrix_dim+1))) SUPERLU_ABORT("Malloc fails for Cmin_colPtr[].");

	cs_di Cmin_cs_di, Craw_cs_di;
	cs_di *Cmod_cs_di;  /* C_modify = C + Cmin matrix*/

	create_Imatrix_csc(matrix_dim, Cmin, Cmin_val, Cmin_rowIdx, Cmin_colPtr);
	cs_di_create(matrix_dim, matrix_dim, matrix_dim, Cmin_colPtr, Cmin_rowIdx, Cmin_val, -1, &Cmin_cs_di);
	cs_di_create(*Cval_size, matrix_dim, matrix_dim, CcolPtr, CrowIdx, Cval, -1, &Craw_cs_di);
	Cmod_cs_di = cs_di_add(&Cmin_cs_di, &Craw_cs_di, 1, 1); 
	
	/* u(t), u(t+h) */
	double *u_t, *u_th;
	int *u_t_size, *u_th_size;
	if (!( u_t= doubleMalloc(data_size))) SUPERLU_ABORT("Malloc fails for u_t[].");
	if (!( u_th= doubleMalloc(data_size))) SUPERLU_ABORT("Malloc fails for u_th[].");
	if (!(u_t_size = intMalloc(1))) SUPERLU_ABORT("Malloc fails for u_t_size[].");
	if (!(u_th_size = intMalloc(1))) SUPERLU_ABORT("Malloc fails for u_th_size[].");

	read_file_double("/home/jerry/ei/csc_data/example_ibmpg1t/u_t.txt", u_t, u_t_size);
	read_file_double("/home/jerry/ei/csc_data/example_ibmpg1t/u_th.txt", u_th, u_th_size);

	/* x(t) */
	double * x_t;	
	int *x_t_size;
	if (!(x_t = doubleMalloc(data_size))) SUPERLU_ABORT("Malloc fails for x_t[].");
	if (!(x_t_size = intMalloc(1))) SUPERLU_ABORT("Malloc fails for x_t_size[].");

	read_file_double("/home/jerry/ei/csc_data/example_ibmpg1t/x_t.txt", x_t, x_t_size);

	SuperMatrix C, G, L, U;  
	int info, permc_spec;
	int* C_perm_r;
	int* C_perm_c;

	dCreate_CompCol_Matrix(&C, matrix_dim, matrix_dim, Cmod_cs_di->nzmax, Cmod_cs_di->x, Cmod_cs_di->i, Cmod_cs_di->p, SLU_NC, SLU_D, SLU_GE);
	dCreate_CompCol_Matrix(&G, matrix_dim, matrix_dim, *Gval_size, Gval, GrowIdx, GcolPtr, SLU_NC, SLU_D, SLU_GE);  
	
	if (!(C_perm_r = intMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for C_perm_r[].");
	if (!(C_perm_c = intMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for C_perm_c[].");

	permc_spec = 1;
	get_perm_c(permc_spec, &C, C_perm_c);

/*=====================Arnoldi=======================*/		
	/*
 *	z1 = (C/gamma) * v1
 *	z2 = (I2-gammga J2)^{-1} * v2
 *	wz2 = W*z2
 *	W = [(u(t+h) - u(t))/h   u(t)]
 *	w =(I-gamma* \hat{A})^{-1} * \hat{v}   new vector
 *	w[0:matrix_dim-1] = (C/gamma+G)^{-1} * (C/gamma *v1 + W(I2-gamma*J2)^{-1}*v2)
 *	w[matrix_dim:matrix_dim+1] = z2
 * */
	double *v_buffer, *w, *z1, *z2, *wz2, *v1, *v2;	
	double norm_v1, w_norm, hij;
	int num_nonzero = 0;
	int dim_Hm = 0;
	char trans[1] = {'N'};
	SuperMatrix Z, CaddG, Hm;

	if (!(v_buffer = doubleMalloc(matrix_dim+2))) SUPERLU_ABORT("Malloc fails for v_buffer[].");
	if (!(w = doubleMalloc(matrix_dim+2))) SUPERLU_ABORT("Malloc fails for w[].");
	if (!(z1 = doubleMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for z1[].");
	if (!(z2 = doubleMalloc(2))) SUPERLU_ABORT("Malloc fails for z2[].");
	if (!(wz2 = doubleMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for wz2[].");
	if (!(v1 = doubleMalloc(matrix_dim))) SUPERLU_ABORT("Malloc fails for v1[].");

	double *Hm_val;
	int *Hm_rowIdx, *Hm_colPtr;
	if (!(Hm_val = doubleMalloc(krylov_dim*krylov_dim))) SUPERLU_ABORT("Malloc fails for Hm_val[].");
	if (!(Hm_rowIdx = intMalloc(krylov_dim*krylov_dim))) SUPERLU_ABORT("Malloc fails for Hm_rowIdx[].");
	if (!(Hm_colPtr = intMalloc(krylov_dim+1))) SUPERLU_ABORT("Malloc fails for Hm_colPtr[].");
	Hm_colPtr[0] = 0;

	double *Vm;
	if (!(Vm = doubleMalloc((matrix_dim+2)*krylov_dim))) SUPERLU_ABORT("Malloc fails for Vm[].");

	cblas_dcopy(matrix_dim, x_t, 1, v_buffer, 1);
	v_buffer[matrix_dim] = 0.0;
	v_buffer[matrix_dim +1] = 1.0;
	
	/* get norm and scale*/
	norm_v1 = cblas_dnrm2(matrix_dim+2, v_buffer, 1); 
	cblas_dscal(matrix_dim+2, 1/norm_v1, v_buffer, 1); 
	add_v_bycol(matrix_dim+2, 0, v_buffer, Vm);	

	/*from now on, C = C/gamma */
	cblas_dscal(*Cval_size, 1/gamma, ((NCformat*)C.Store)->nzval, 1);

	cs_di C_cs_di, G_cs_di;
	cs_di_create(((NCformat*)C.Store)->nnz, matrix_dim, matrix_dim, ((NCformat*)C.Store)->colptr, \
				((NCformat*)C.Store)->rowind, ((NCformat*)C.Store)->nzval, -1, &C_cs_di);
	cs_di_create(((NCformat*)G.Store)->nnz, matrix_dim, matrix_dim, ((NCformat*)G.Store)->colptr, \
				((NCformat*)G.Store)->rowind, ((NCformat*)G.Store)->nzval, -1, &G_cs_di);
	
	cs_di *CaddG_cs_di;  /* C/gamma+G */
	CaddG_cs_di = cs_di_add(&C_cs_di, &G_cs_di, 1, 1); 

	dCreate_CompCol_Matrix(&CaddG, matrix_dim, matrix_dim, CaddG_cs_di->nzmax, CaddG_cs_di->x, \
							CaddG_cs_di->i, CaddG_cs_di->p, SLU_NC, SLU_D, SLU_GE);

	/* Arnoldi main loop */
	for (j=0;j<krylov_dim; j++){  
		get_v(matrix_dim+2, j, v_buffer, Vm);	
		z2[0] = v_buffer[matrix_dim] + gamma * v_buffer[matrix_dim+1];
		z2[1] = v_buffer[matrix_dim+1]; 		
		set_vector_zero(matrix_dim, wz2);
		cblas_daxpy(matrix_dim, -(v_buffer[matrix_dim]+gamma*v_buffer[matrix_dim+1])/t_step , u_th, 1, wz2, 1);
		cblas_daxpy(matrix_dim, (v_buffer[matrix_dim]+gamma*v_buffer[matrix_dim+1])/t_step-v_buffer[matrix_dim+1] , u_t, 1, wz2, 1); /*wz2*/
		cblas_dcopy(matrix_dim, v_buffer, 1, v1, 1); 
		sp_dgemm(trans, C.nrow, 1, C.ncol, 1, &C, v1, matrix_dim, 0, z1, matrix_dim); /* z1=C/gamma * v1  */
	
		cblas_daxpy(matrix_dim, 1, z1, 1, wz2, 1); /* wz2:= z1+wz2*/
		dCreate_Dense_Matrix(&Z, matrix_dim, 1, wz2, matrix_dim, SLU_DN, SLU_D, SLU_GE);

		/*Z = (C/gamma+G)^{-1} * (C/gamma *v1 + W(I2-gamma*J2)^{-1}*v2)*/
		pdgssv(nprocs, &CaddG, C_perm_c, C_perm_r, &L, &U, &Z, &info); 

		dCopy_Dense_Matrix(matrix_dim, 1, ((DNformat*)Z.Store)->nzval, 0, w, 0);
		w[matrix_dim] = z2[0];
		w[matrix_dim+1] = z2[1];  

		for (i=0; i<=j ; i++) {
			get_v(matrix_dim+2, i, v_buffer, Vm);
			hij = cblas_ddot(matrix_dim+2, w, 1, v_buffer, 1);
			Hm_val[Hm_colPtr[j]+i] = hij; 
			Hm_rowIdx[Hm_colPtr[j]+i] = i;
			num_nonzero += 1;
			cblas_daxpy(matrix_dim+2, -hij , v_buffer, 1, w, 1);
		}

		w_norm = cblas_dnrm2(matrix_dim+2, w, 1); 
		if (fabs(w_norm) > eps){   
			Hm_val[num_nonzero] = w_norm; 
			Hm_rowIdx[num_nonzero] = j+1;
			num_nonzero += 1;

			cblas_dscal(matrix_dim+2, 1/w_norm, w, 1);
			add_v_bycol(matrix_dim+2, j+1, w, Vm);	
		}else if(fabs(w_norm)<= eps){     /*converged condition*/
			converged = True;
			Hm_colPtr[j+1] = num_nonzero;
			
			break;
		}
		dim_Hm = j+1;
		Hm_colPtr[j+1] = num_nonzero;
	}

	/*===========The result of  krylov subspace approximation=========*/
	dCreate_CompCol_Matrix(&Hm, dim_Hm, dim_Hm, num_nonzero, Hm_val, Hm_rowIdx, Hm_colPtr, SLU_NC, SLU_D, SLU_GE);

	double *matrix_Im;
	if (!(matrix_Im = doubleMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for matrix_Im[].");
	double *Im_val;
	int *Im_rowIdx, *Im_colPtr, *Im_nnz;
	if (!(Im_val= doubleMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for Im_val[].");
	if (!(Im_rowIdx = intMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for Im_rowIdx[].");
	if (!(Im_colPtr = intMalloc(dim_Hm+1))) SUPERLU_ABORT("Malloc fails for Im_colPtr[].");
	if (!(Im_nnz = intMalloc(1))) SUPERLU_ABORT("Malloc fails for Im_nnz[].");

	SuperMatrix Hm_inv;
	int* Hm_perm_r;
	int* Hm_perm_c;
	if (!(Hm_perm_r = intMalloc(dim_Hm))) SUPERLU_ABORT("Malloc fails for Hm_perm_r[].");
	if (!(Hm_perm_c = intMalloc(dim_Hm))) SUPERLU_ABORT("Malloc fails for Hm_perm_c[].");

	SuperMatrix Hm_L, Hm_U;
	double *Hm_inv_val;
	int *Hm_inv_rowIdx, *Hm_inv_colPtr, *Hm_inv_nnz;
	if (!(Hm_inv_val= doubleMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for Hm_inv_val[].");
	if (!(Hm_inv_rowIdx = intMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for Hm_inv_rowIdx[].");
	if (!(Hm_inv_colPtr = intMalloc(dim_Hm+1))) SUPERLU_ABORT("Malloc fails for Hm_inv_colPtr[].");
	if (!(Hm_inv_nnz = intMalloc(1))) SUPERLU_ABORT("Malloc fails for Hm_inv_nnz[].");

	double *I_HmInv;  /* (I-Hm^{-1}) * h/gamma */
	if (!(I_HmInv = doubleMalloc(dim_Hm*dim_Hm))) SUPERLU_ABORT("Malloc fails for I_HmInv[].");

	double *exp_I_HmInv_e1; /* exp((I-Hm^{-1}) * h/gamma)) *e1 */
	if (!(exp_I_HmInv_e1 = doubleMalloc(dim_Hm))) SUPERLU_ABORT("Malloc fails for exp_I_HmInv_e1[].");

	double *x_th; /*x_th*/  
	if (!(x_th = doubleMalloc(matrix_dim+2))) SUPERLU_ABORT("Malloc fails for y[].");

	cs_di Hm_inv_cs_di, Im_cs_di;
	cs_di *I_HmInv_cs_di; /* (I-Hm^{-1}) * h/gamma */
	
	create_Imatrix(dim_Hm, matrix_Im, 1);
	dense_to_csc(matrix_Im, dim_Hm, dim_Hm, Im_val, Im_rowIdx, Im_colPtr, Im_nnz);
	dCreate_Dense_Matrix(&Hm_inv, dim_Hm, dim_Hm, matrix_Im, dim_Hm, SLU_DN, SLU_D, SLU_GE);
	get_perm_c(permc_spec, &Hm, Hm_perm_c);
	
	/* Hm_inv =>>>>> Hm^{-1} */
	pdgssv(nprocs, &Hm, Hm_perm_c, Hm_perm_r, &Hm_L, &Hm_U, &Hm_inv, &info); 

	DNformat *Hmstore;
	Hmstore = (DNformat *)Hm_inv.Store;
	dense_to_csc((double*)Hmstore->nzval, dim_Hm, dim_Hm, Hm_inv_val, \
				Hm_inv_rowIdx, Hm_inv_colPtr, Hm_inv_nnz);

	cs_di_create(*Hm_inv_nnz, dim_Hm, dim_Hm, Hm_inv_colPtr, Hm_inv_rowIdx, Hm_inv_val, -1, &Hm_inv_cs_di);
	cs_di_create(dim_Hm, dim_Hm, dim_Hm, Im_colPtr, Im_rowIdx, Im_val, -1, &Im_cs_di);
	I_HmInv_cs_di = cs_di_add(&Im_cs_di, &Hm_inv_cs_di, t_step/gamma, -t_step/gamma);
	csc_to_arrayByCol(I_HmInv_cs_di->nzmax, I_HmInv_cs_di->m, I_HmInv_cs_di->n, I_HmInv_cs_di->x, \
				I_HmInv_cs_di->i, I_HmInv_cs_di->p, I_HmInv);
	/* exp((I-Hm^{-1}) * h/gamma)) */
	double *exp_I_HmInv = r8mat_expm1(dim_Hm, I_HmInv);
	memcpy(exp_I_HmInv_e1, exp_I_HmInv, dim_Hm*sizeof(double));
	/* x_th = \beta * Vm * exp((I-Hm^{-1}) * h/gamma)) *e1 */
	cblas_dgemv(CblasColMajor, CblasNoTrans, matrix_dim+2, dim_Hm, norm_v1, \
				 Vm, matrix_dim+2, exp_I_HmInv_e1, 1, 0, x_th, 1);

	printf("x_th: \n");
	for(i=0; i<matrix_dim+2; i++){
		printf("%f  \n",x_th[i]);
	}

	SUPERLU_FREE(Gval);
	SUPERLU_FREE(GrowIdx);
	SUPERLU_FREE(GcolPtr);
	SUPERLU_FREE(Cval);
	SUPERLU_FREE(CrowIdx);
	SUPERLU_FREE(CcolPtr);
	SUPERLU_FREE(C_perm_r);
	SUPERLU_FREE(C_perm_c);
	SUPERLU_FREE(x_t);
	SUPERLU_FREE(v_buffer);
	SUPERLU_FREE(Hm_val);
	SUPERLU_FREE(Hm_rowIdx);
	SUPERLU_FREE(Hm_colPtr);
	SUPERLU_FREE(w);
	SUPERLU_FREE(z1);
	SUPERLU_FREE(z2);
	SUPERLU_FREE(wz2);
	SUPERLU_FREE(v1);
	SUPERLU_FREE(Vm);



	


	printf("**********end**********\n");	

	return 0;
}
