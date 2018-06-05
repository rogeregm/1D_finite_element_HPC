/*HPC-Coursework
*	Roger Gonzalez
*	26/03/17
*/

// Common Functions
	void   VEC_print           (int n, double* Vec);
	void   MAT_print           (int m, int n, double* A);
	double MAT_print_col       (int m, int n, double* A);
	void   RowtoCol            (double* A, int n, int m, double* B);
	void   MAT_SCALAR_MULT     (double* A, int size, double alpha);
	void   MAT_VEC_SUM_COL_BAND(double alpha, double* A, int sizeA1, int sizeA2, double beta, double* B, int sizeB, int start);
	void   VEC_SUM             (double* A,double alpha, double* B, double beta, int size);
	void   VEC_VEC_MULT        (double* A, double* B, int size);
	void   INIT_ZERO           (int size, double* A);

	double PARA_PRINT          (int m, int n, double* A, int flag);
	void   row_ones            (double* A, int size, int row);

	void   SolveWithScalapack  (int n, double* H, double* f);

// Element Space functions
	void MAT_stiffness_e       (double* K_e, double A, double E, double I, double l);
	void VEC_force_e           (double* F_e, double l, double qx, double qy);
	void MAT_Mass_e            (double* MAT, double A, double rho, double l);
	void VEC_DYN_force_e       (double* F_e, double l, double qx, double qy, int Nt, double step);

//Global Matrices
	void VEC_force_g           (double* F_g, double* F_e, int N);
	void BAND_stiffness_g      (double* BK_g, double* BK_e, int N, int row1);
	void BAND_stiffness_e      (double* BAND_K_e, double A, double E, double I, double l);
	void MAT_Mass_g            (double* MAT_e, double* MAT_g, int size);
	void VEC_DYN_force_g       (double* F_g, double* F_e, int N, int Nt, double step);
	void PARA_VEC_DYN_force_g  (double* F_g, double* F_e, int N, int Nt, double step);

	void BAND_stiffness_g_PARA (double* BK_g, double* BK_e, int N, int row1);


// BLAS and LAPACK

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dnrm2)(const int& N, const double *X, const int& incX);
    void F77NAME(daxpy)(const int& N, const double& alpha, const double *X,
                             const int& incX, double *Y, const int& incY);
    void F77NAME(dgbmv)(const char& trans, const int& m, const int& n,
                        const int& kl, const int& ku,
                        const double& alpha, const double* a, const int& lda,
                        const double* x, const int& incx, const double& beta,
                    double* y, const int& incy);

    void F77NAME(dgbsv)(const int& n, 
    					const int& kl, 
    					const int& ku, 
                    	const int& nrhs, 
                    	const double * A, 
                    	const int& ldab, 
                    	int * ipiv, 
                    	double * B, 
                    	const int& ldb, 
                    	int* info);
    void F77NAME(pdgbsv)(const int& n, const int& kl, const int& ku, 
                    const int& nrhs, const double * A, const int& ja,
                    const int* desca, int * ipiv, double * B, const int& ib,
                    const int* descb, double* work, const int& lwork, 
                    int* info);

	void Cblacs_get(int, int, int*);
    void Cblacs_pinfo(int*, int*);
    void Cblacs_gridinit(int*, char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_exit(int);
    void Cblacs_gridexit(int);

}
