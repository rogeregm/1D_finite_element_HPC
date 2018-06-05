/*HPC-Coursework
*	Roger Gonzalez
*	26/03/17
*/
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<vector>

#include "cblas.h"
#include "mpi.h"

#include "functions.h"

using namespace std;

/***********************************************************/

// To create the element stiffness matrix
void MAT_stiffness_e(double* K_e, double A, double E, double I, double l){
	
	double  k0 = (A*E)/l;
	double  k1 = (E*I)/l;
	double  k2 = (E*I)/(l*l);
	double  k3 = (E*I)/(l*l*l);
	
	K_e[0]=k0;
	K_e[3]=-k0;
	K_e[7]=12*k3;
	K_e[8]=6*k2;
	K_e[10]=-12*k3;
	K_e[11]=6*k2;
	K_e[13]=6*k2;
	K_e[14]=4*k1;
	K_e[16]=-6*k2;
	K_e[17]=2*k1;
	K_e[18]=-k0;
	K_e[21]=k0;
	K_e[25]=-12*k3;
	K_e[26]=-6*k2;
	K_e[28]=12*k3;
	K_e[29]=-6*k2;
	K_e[31]=6*k2;
	K_e[32]=2*k1;
	K_e[34]=-6*k2;
	K_e[35]=4*k1;

} 

/***********************************************************/
// The make the element stiffness matrix in banded form
void BAND_stiffness_e(double* BAND_K_e, double A, double E, double I, double l){
	
	double  k0 = (A*E)/l;
	double  k1 = (E*I)/l;
	double  k2 = (E*I)/(l*l);
	double  k3 = (E*I)/(l*l*l);


	       /**/
	BAND_K_e[5]=6*k2;
	BAND_K_e[9]=-k0;
	BAND_K_e[10]=-12*k3;
	BAND_K_e[11]=2*k1;
	BAND_K_e[16]=-6*k2;
	BAND_K_e[20]=6*k2;
	BAND_K_e[23]=-6*k2;
		/*Diagonals*/
	BAND_K_e[24]=k0;
	BAND_K_e[25]=12*k3;
	BAND_K_e[26]=4*k1;
	BAND_K_e[27]=k0;
	BAND_K_e[28]=12*k3;
	BAND_K_e[29]=4*k1;
		/**/
	BAND_K_e[31]=6*k2;
	BAND_K_e[34]=-6*k2;
	BAND_K_e[38]=-6*k2;
	BAND_K_e[42]=-k0;
	BAND_K_e[43]=-12*k3;
	BAND_K_e[44]=2*k1;
	BAND_K_e[49]=6*k2;

}

/***********************************************************/
// To create the global stiffness matrix in banded form
void BAND_stiffness_g(double* BK_g, double* BK_e, int N, int row1)
{
	for(int i = 0; i < 3*(N-1) ; i+=3) // Global matrix elements 
  	{
  		if(i!=0)
  		{
  			for(int j = 0 ; j < 3 ; j++)
  			{	
  				BK_g[ (row1+0)*3*(N-1) + i + j] = BK_e[      j+3];
  				BK_g[ (row1+1)*3*(N-1) + i + j] = BK_e[  6 + j+3];
  				BK_g[ (row1+2)*3*(N-1) + i + j] = BK_e[2*6 + j+3];
  				BK_g[ (row1+3)*3*(N-1) + i + j] = BK_e[3*6 + j+3]+BK_e[3*6 + j  ];
  				BK_g[ (row1+4)*3*(N-1) + i + j] = BK_e[4*6 + j  ]+BK_e[4*6 + j  ];
  				BK_g[ (row1+5)*3*(N-1) + i + j] = BK_e[5*6 + j  ]+BK_e[5*6 + j+3];
  				BK_g[ (row1+6)*3*(N-1) + i + j] = BK_e[6*6 + j  ];
  				BK_g[ (row1+7)*3*(N-1) + i + j] = BK_e[7*6 + j  ];
  				BK_g[ (row1+8)*3*(N-1) + i + j] = BK_e[8*6 + j  ];
  			}
  		}
  		else
  		{
  			for(int j = 0 ; j < 3 ; j++)
  			{	
  				BK_g[ (row1+0)*3*(N-1) + i + j] = BK_e[      j  ];
  				BK_g[ (row1+1)*3*(N-1) + i + j] = BK_e[  6 + j  ];
  				BK_g[ (row1+2)*3*(N-1) + i + j] = BK_e[2*6 + j  ];
  				BK_g[ (row1+3)*3*(N-1) + i + j] = BK_e[3*6 + j  ]+BK_e[3*6 + j+3];
  				BK_g[ (row1+4)*3*(N-1) + i + j] = BK_e[4*6 + j+3]+BK_e[4*6 + j+3];
  				BK_g[ (row1+5)*3*(N-1) + i + j] = BK_e[5*6 + j+3]+BK_e[5*6 + j  ];
  				BK_g[ (row1+6)*3*(N-1) + i + j] = BK_e[6*6 + j  ];
  				BK_g[ (row1+7)*3*(N-1) + i + j] = BK_e[7*6 + j  ];
  				BK_g[ (row1+8)*3*(N-1) + i + j] = BK_e[8*6 + j  ];
  			}
  		}
  		if(i==3*(N-2))
  		{
  			for(int j = 0 ; j < 3 ; j++)
  			{	
  				BK_g[ (row1+0)*3*(N-1) + i + j] = BK_e[      j+3];
  				BK_g[ (row1+1)*3*(N-1) + i + j] = BK_e[  6 + j+3];
  				BK_g[ (row1+2)*3*(N-1) + i + j] = BK_e[2*6 + j+3];
  				BK_g[ (row1+3)*3*(N-1) + i + j] = BK_e[3*6 + j+3]+BK_e[3*6 + j  ];
  				BK_g[ (row1+4)*3*(N-1) + i + j] = BK_e[4*6 + j  ]+BK_e[4*6 + j  ];
  				BK_g[ (row1+5)*3*(N-1) + i + j] = BK_e[5*6 + j  ]+BK_e[5*6 + j+3];
  				BK_g[ (row1+6)*3*(N-1) + i + j] = 0;
  				BK_g[ (row1+7)*3*(N-1) + i + j] = 0;
  				BK_g[ (row1+8)*3*(N-1) + i + j] = 0;
  			}
  		}
	}
}

/***********************************************************/
// To create the elemental force vector
void VEC_force_e(double* F_e, double l, double qx, double qy){

	 	// Concentrated force 
	double fx=l*qx/2.0;
	double fy=l*qy/2.0;
	
	F_e[0]=fx;
	F_e[1]=fy;
	F_e[2]=(fy*l*l)/6.0;
	F_e[3]=fx;
	F_e[4]=fy;
	F_e[5]=(-fy*l*l)/6.0;
} 

/***********************************************************/
// To create the global Fore vector
void VEC_force_g(double* F_g, double* F_e, int N)
{

		double Fy=-1000.0;	// Concentrated force

  		for (int i = 0; i < 3*(N-1) ; i+=3)
  		{
  			if(i%2==0)
  			{
  				for(int j = 0; j < 3 ; j++)
  				{
  					F_g[i+j]=F_e[j]+F_e[j+3];

  					if( i+j==0.5*(3*(N-1)-1) ) F_g[i+j]=F_g[i+j]+Fy;	
  				}
  			}
  			else
  			{
  				for(int j = 0; j < 3 ; j++)
  				{
  					F_g[i+j]=F_e[3+j]+F_e[j];

  					if( i+j==0.5*(3*(N-1)-1) ) F_g[i+j]=F_g[i+j]+Fy;	
  				}
  			}			
  		}
}

/*******************************************************/
// To create the Mass elemental matrix ( only diagonal stored)
void MAT_Mass_e(double* MAT, double A, double rho, double l)
{
	double mass=A*rho*l;
	double alpha=1.0/24.0;

	MAT[0]=mass*(1.0/2.0);
	MAT[1]=mass*(1.0/2.0);
	MAT[2]=mass*alpha*l*l;
	MAT[3]=mass*(1.0/2.0);
	MAT[4]=mass*(1.0/2.0);
	MAT[5]=mass*alpha*l*l;
}

// To create global Mass matrix (only diagonal stored)
void MAT_Mass_g(double* MAT_e, double* MAT_g, int size)
{
	for(int i=0 ; i < size ; i+=3)
	{	
		for (int j=0 ; j<3 ; j++)
		{
			MAT_g[i+j]=2*MAT_e[j];
		}
	}
}

// To create the elemental force vector fo the dynamic case
void VEC_DYN_force_e(double* F_e, double l, double qx, double qy, int Nt, double step){

	 	// Concentrated force 
		   qx=l*qx;
		   qy=l*qy;

	double fx=step*qx/2.0;
	double fy=step*qy/2.0;
	
	F_e[0]=fx;
	F_e[1]=fy;
	F_e[2]=(fy*l*l)/6.0;
	F_e[3]=fx;
	F_e[4]=fy;
	F_e[5]=(-fy*l*l)/6.0;

} 

/***********************************************************/
// To create the Global force vector fo the dynamic case
void VEC_DYN_force_g(double* F_g, double* F_e, int N, int Nt, double step)
{

		double Fy=-1000.0;	// Concentrated force
		Fy=step*Fy;

  		for (int i = 0; i < 3*(N-1) ; i+=3)
  		{
  			if(i%2==0)
  			{
  				for(int j = 0; j < 3 ; j++)
  				{
  					F_g[i+j]=F_e[j]+F_e[j+3];

  					if( i+j==0.5*(3*(N-1)-1) ) F_g[i+j] = F_g[i+j]+Fy;
  				}
  			}
  			else
  			{
  				for(int j = 0; j < 3 ; j++)
  				{
  					F_g[i+j]=F_e[3+j]+F_e[j];

  					if( i+j==0.5*(3*(N-1)-1) ) F_g[i+j] = F_g[i+j]+Fy;	
  				}
  			}			
  		}
}

// To create the Global force vector for the case in Parallel
void PARA_VEC_DYN_force_g(double* F_g, double* F_e, int N, int Nt, double step)
{
		double Fy=-1000.0;	// Concentrated force
		Fy=step*Fy;

  		for (int i = 0; i < 3*(N-1) ; i+=3)
  		{
  				for(int j = 0; j < 3 ; j++)
  				{
  					F_g[i+j]=F_e[j]+F_e[j+3];
  				}		
  		}
}

// To create the stiffness matrix for a system processed in parallel
void BAND_stiffness_g_PARA(double* BK_g, double* BK_e, int N, int row1)
{
	for(int i = 0; i < 3*(N-1) ; i+=3) // Global matrix elements 
  	{
  		for(int j = 0 ; j < 3 ; j++)
  		{	
  			BK_g[ (row1+0)*3*(N-1) + i + j] = BK_e[      j+3];
  			BK_g[ (row1+1)*3*(N-1) + i + j] = BK_e[  6 + j+3];
  			BK_g[ (row1+2)*3*(N-1) + i + j] = BK_e[2*6 + j+3];
  			BK_g[ (row1+3)*3*(N-1) + i + j] = BK_e[3*6 + j+3]+BK_e[3*6 + j  ];
  			BK_g[ (row1+4)*3*(N-1) + i + j] = BK_e[4*6 + j  ]+BK_e[4*6 + j  ];
  			BK_g[ (row1+5)*3*(N-1) + i + j] = BK_e[5*6 + j  ]+BK_e[5*6 + j+3];
  			BK_g[ (row1+6)*3*(N-1) + i + j] = BK_e[6*6 + j  ];
  			BK_g[ (row1+7)*3*(N-1) + i + j] = BK_e[7*6 + j  ];
  			BK_g[ (row1+8)*3*(N-1) + i + j] = BK_e[8*6 + j  ];
  		}
  	}
 }

// Function to Solve with Scalapack
// - for a system of N elements
void SolveWithScalapack(int n, double* H, double* f)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nb   = n/size;        // Block size
    const int kl   = 4;             // Number of lower diagonals
    const int ku   = 4;             // Number of upper diagonals
    const int lda  = 1 + 2*kl + 2*ku; // Leading dimension (num of rows)
    const int nrhs = 1;
    const int ja   = 1;             // Offset index in global array (col)
    const int ib   = 1;             // Offset index in global array (row)
    const int lwork= (nb+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku)
                        + max(nrhs*(nb+2*kl+4*ku), 1);
    int info = 0;

    // Set up the CBLACS context
    char order = 'C';
    int nrow = 1;
    int ncol = size;
    int ctx;
    int mype;
    int npe;
    int myrow;
    int mycol;
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get( 0, 0, &ctx );
    Cblacs_gridinit( &ctx, &order, 1, npe );
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);

    // Create descriptors for matrix and RHS vector storage
    int desca[7];
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = ctx;         // Context
    desca[2] = n;           // Problem size
    desca[3] = nb;          // Block size
    desca[4] = 0;           // Process row/column
    desca[5] = lda;         // Local size
    desca[6] = 0;           // Reserved

    int descb[7];
    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = ctx;         // Context
    descb[2] = n;           // Problem size
    descb[3] = nb;          // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = nb;          // Local size
    descb[6] = 0;           // Reserved

    int*  ipiv = new int[n];
    double* wk = new double[lwork];

    // Solve the system H * u' = u (i.e. RHS vector replaced by solution)
    F77NAME(pdgbsv) (n, kl, ku, nrhs, H, ja, desca, ipiv, f, ib, descb,
                        wk, lwork, &info);

    if (info) {
        cout << "Information from computation (PDGBSV): " << info << endl;
    }
    // Free CBLACS context
    Cblacs_gridexit( ctx );

	delete[] wk;
	delete[] ipiv;
}