/*HPC-Coursework
*   Roger Gonzalez
*   26/03/17
*/
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<vector>

#include "cblas.h"
#include "mpi.h"

using namespace std;

double VEC_print(int n, double* Vec){
	for(int i=0; i<n ; i++){
        // cout.precision(10);		
        cout.width(12);
		cout << Vec[i] << " ";
	}
	cout <<  endl;
}

double MAT_print(int m, int n, double* A)
{
	for (int i = 0; i < m ; i++){
		for (int j = 0; j < n ; j++){
			cout.width(12);
			cout << A[i*n+j] << " " ;
		}
		cout << endl;
	}
} 

// To print a column major matrix
double MAT_print_col(int m, int n, double* A)
{
    for (int i = 0; i < m ; i++){
        for (int j = 0; j < n ; j++){
            cout.width(12);
            cout << A[j*m+i] << " " ;
        }
        cout << endl;
    }
} 

// To change  matrix from row major to column major
void RowtoCol(double* A, int n, int m, double* B)
{
    for(int i = 0 ; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            B[i*m+j] = A[j*n+i];
        }
    }
}

// To multiply a Matrix by a scalar
void MAT_SCALAR_MULT(double* A, int size, double alpha)
{
    for(int i = 0 ; i < size ; i++)
    {
        A[i] = alpha*A[i];
    }

}

// To sum a vector to a row in a matrix where the matrix is indexed in column major
void MAT_VEC_SUM_COL_BAND(double alpha, double* A, int sizeA1, int sizeA2, double beta, double* B, int sizeB, int row)
{   
        int j = 0;
    for(int i = row; i < sizeA1*sizeA2 ; i += sizeA1)
    {
        A[i] = alpha*A[i] + beta*B[j];
        j++;
    }
}

// Sum two vectors
void VEC_SUM(double* A,double alpha, double* B, double beta, int size)
{
    for (int i = 0 ; i < size ; i++)
    {
        A[i] = alpha*A[i] + beta*B[i];
    }
}

// To multiply two vectors
void VEC_VEC_MULT(double* A, double* B, int size)
{
    for (int i = 0 ; i < size ; i++)
    {
        A[i] = A[i] * B[i];
    }
}

// To Initialise an array to zeros( safety net)
void INIT_ZERO(int size, double* A)
{
    for(int i = 0 ; i < size ; i++)
    {
        A[i] = 0;
    }
    
}

// To print from only a single processor when running MPI
// - can pring vectors, and matrices in row major or column major indexing
double PARA_PRINT(int m, int n, double* A, int flag)
{
    int p;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(rank == 1 || rank==p-1)
    {
        cout << "Printing done by process: " << rank << endl;

        if(n == 0 && flag != 1)VEC_print(m, A);

        if(n != 0 && flag != 1)MAT_print(m, n, A);

        if(flag==1)MAT_print_col(m, n, A);
    }
}

// To fill the row in a matrix with ones
void row_ones(double* A, int size, int row)
{
    for(int i = 0 ; i < size ; i++)
    {
        A[row*size + i] = 1;
    }
}