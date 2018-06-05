/*HPC-Coursework
*	Roger Gonzalez
*	26/03/17
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include "cblas.h"
#include "mpi.h"
#include "functions.h"

using namespace std;


int main( int argc, char *argv[] )
{
  // INITIALISING MPI
  int retval;
  retval = MPI_Init(&argc, &argv);
  if (retval!=MPI_SUCCESS)
  {
    cout << "An error occurres in initialising MPI" << endl;
  }
  else
  {
    cout << "MPI succesful" <<endl; 
  }

// Declaring input variables
	const double  L=atof(argv[1]); 
	const int     N=atoi(argv[2]); //if(N%2!=0){cout<< "N needs to be even, Default value use : 24 "<<endl; N=24;}
	const double  A=atof(argv[3]);
	const double  I=atof(argv[4]);
	const double  E=atof(argv[5]);
	const double  l=L/N;	//segment lengths
	const double  T      =atof(argv[6]);
	const double  Nt     =atof(argv[7]);
	const int     Not    =atoi(argv[7]);
	const double  rho    =atof(argv[8]);
	const int     choice1=atoi(argv[9]);
	const int     choice2=atoi(argv[10]);

	const double del_T  =(T/Nt);
	const double timefac=1/(del_T*del_T);

	const double  qx=0;
	const double  qy=-1000.0;// N/m
  const double  Fy=-1000.0;

	cout << "----------------------------Inputs---------------------" << endl;
	cout << "lenght (m)          = " << L << endl;
	cout << "Elements            = " << N << endl;
	cout << "Area   (m^2)        = " << A << endl;
	cout << "Inertia(m^4)        = " << I << endl;
	cout << "Young's Modulus (Pa)= " << E << endl;
	cout << "Time                = " << T << endl;
	cout << "Number of Timesteps = " << Nt << endl;
	cout << "Timestep            = " << del_T <<endl;
  cout << "Timefactor          = " << timefac <<endl;
	cout << "Density (kg/m^3)    = " << rho << endl;
	cout << "Mass                = " << A*rho*l <<endl;
	cout << endl;

	const int g_size = 3*(N-1);
	
/************************************************************/
  // Stiffness elemental matrix	
	double* K_e = new double[6*6]();
	MAT_stiffness_e( K_e, A, E, I, l);
  // Banded Stiffness elemental Matrix
	double* BAND_K_e = new double[9*6](); 
	BAND_stiffness_e(BAND_K_e, A, E, I, l);

/**********************************************************/
// Force elemental vector
	double * F_e = new double [6]();
	VEC_force_e(F_e, l, qx, qy);

/************************* STATIC SOLUTION *************************************/
 
if(choice1==0)
{
  // Global matrix in banded form
  double* BAND_K_g = new double [13*(g_size)]();
  BAND_stiffness_g(BAND_K_g, BAND_K_e, N, 4);
  double* Col_K_g = new double [13*(g_size)]();
  RowtoCol(BAND_K_g, g_size, 13, Col_K_g);
  /**************************************************************************/
  // Global force vector
  double* F_g = new double [g_size]();
  VEC_force_g(F_g, F_e, N);

  /***************************************************************************/
  // STATIC Solving linear system Ax=B
 
  int     info = 0;
  int* 	  ipiv = new int [g_size];
  F77NAME(dgbsv) ( g_size, 4, 4, 1, Col_K_g, 13, ipiv, F_g, g_size, &info) ;

  /*************************************************************************************************/
  // SAVING NUMERICAL SOLUTION
  double* Numerical = new double [N+1]();
  int counter =1;
  for(int i=1 ; i<g_size ; i+=3)
  {
  	  	Numerical[counter]=F_g[i];
  	  	counter++;
  }
/**********************************************************************/
  // Analytical Solution
  		double* x= new double [N+1];
  		double An_distributed;//= new double [N];
  		double An_concentrated;//= new double [N];
  		double* Analytical= new double [N];
      for(int i = 0 ; i < N+1 ; i++) x[i] = i*l;

      int xcounter = 0;
  		for (int i = 0 ; i < 0.5*N ; i++)
  		{
  			An_distributed      = (qy*x[i]*x[i]) * (L-x[i]) * (L-x[i]) / (24*E*I);
  			An_concentrated     = (Fy*x[i]*x[i]) * (3*L - 4*x[i]) / (48*E*I);
  			Analytical[xcounter]= An_distributed + An_concentrated; 
        xcounter++ ;   
  		}
      for (int i=0.5*N ; i>-1 ; i--)
      {
        An_distributed      =(qy*x[i]*x[i]) * (L-x[i]) * (L-x[i]) / (24*E*I);
        An_concentrated     =(Fy*x[i]*x[i]) * (3*L-4*x[i]) / (48*E*I);
        Analytical[xcounter]= An_distributed + An_concentrated;
        xcounter++;     
      }

  //Writing Analytical and Numerical Solutions to file
	ofstream file1("Static_displacements.txt");
    if (!file1.good()) 
    {
        cout << "Error: unable to open output file: Static_displacements.txt" << endl;
    }
    else 
    {
    	for (int i = 0 ; i < N+1 ; ++i) 
    	{
    		file1.width(10);
    		file1 << x[i] << "," << Analytical[i] << "," << Numerical[i] << "," << endl;
    	}
    	file1.close();
        cout << "Written file: Static_displacements.txt" << endl;
    }
} // END OF STATIC

  /***************************************************************************/
 /******************************* DYNAMIC ***********************************/
/***************************************************************************/
if(choice1==1)
{
  // Building Mass Elemental Matrix
 	double* M_e= new double[6]();
 	MAT_Mass_e( M_e, A,  rho, l);
  // Building Mass Global Matrix
 	double* M_g= new double[g_size]();
 	MAT_Mass_g(M_e, M_g, g_size);

  if(choice2==0)
 	{
 		double* BAND_K_g = new double [9*(g_size)]();
 		double* Col_K_g  = new double [9*(g_size)]();
    // Building Global Stiffness matrix
 		BAND_stiffness_g( BAND_K_g, BAND_K_e, N, 0);
 		RowtoCol(BAND_K_g, g_size, 9, Col_K_g);
 /*****************************************************************************************************/
 		MAT_SCALAR_MULT( M_g, g_size, timefac);
 /************************************ EXPLICIT *************************************************/
 	  MAT_VEC_SUM_COL_BAND(1, Col_K_g, 9, g_size, -2, M_g, g_size, 4);
 /*************************************************************************************************/
 		// Dynamic System
 		double* F_g    = new double[g_size]();
 		double* u_n    = new double[g_size]();
 		double* u_n_1  = new double[g_size]();
 		double* temp1  = new double[g_size]();
 		double* temp2  = new double[g_size]();
 		double* storage= new double[(N+1)*Not]();
 		
 		INIT_ZERO(g_size,u_n);
 		INIT_ZERO(g_size,u_n_1);
    int       info = 0;

  // double* osc = new double [100*10000];
  // for(int out=1 ; out<100 ; out++)
  // {
    double load_time_factor = 1;
    double load_factor = 1/load_time_factor;
 		for (int count = 0 ; count*del_T < T ; count++)
 		{	
 		 // Building Force vector with variable load factor
      if(count*del_T < load_time_factor*T)
      {
        INIT_ZERO(g_size,F_g);
 			  VEC_DYN_force_e(F_e,   l, qx, qy, Nt, del_T*count*load_factor);
 			  VEC_DYN_force_g(F_g, F_e,  N, Nt, del_T*count*load_factor);
      }
      else
      {
        INIT_ZERO(g_size,F_g);
        VEC_force_e(F_e,   l, qx, qy);
        VEC_force_g(F_g, F_e,  N);
      }
      // ReBuilding Mass Matrix
      INIT_ZERO(g_size,M_g);
      MAT_Mass_g(M_e, M_g, g_size);
      /************************************BUILIDNG MASS/timesept MATRIX*****************************/
      MAT_SCALAR_MULT( M_g, g_size, timefac);
      /************************************BUILDING STIFFNESS MATRIX********************************/
      INIT_ZERO(g_size*9,BAND_K_g);
      INIT_ZERO(g_size*9, Col_K_g);
      BAND_stiffness_g( BAND_K_g, BAND_K_e, N, 0);
      RowtoCol(BAND_K_g, g_size, 9, Col_K_g);
      /*****************************************      K-M      **************************************/
      MAT_VEC_SUM_COL_BAND(1, Col_K_g, 9, g_size, (-2), M_g, g_size, 4); //rewrites Col_K_g, but keeps M_g
      /**************************************STARTING CALCULATIONS**********************************/
      INIT_ZERO(g_size,temp2);
      //BUILDING RHS
      VEC_VEC_MULT(M_g, u_n_1, g_size);//rewrites M_g keeps un-1
 			F77NAME(dgbmv)('N', g_size, g_size, 4, 4, 1, Col_K_g, 9, u_n, 1, 0, temp2, 1); //writes to temp2 and keeps u_n
 			VEC_SUM(temp2, -1,   M_g, (-1), g_size);// writes to temp2
 			VEC_SUM(  F_g, 1, temp2, 1, g_size);//writes to F_g : now RHS of system, and keeps temp2
      /****************************RE-BUILDING MASS MATRIX***********************************/
      INIT_ZERO(g_size,M_g);
      MAT_Mass_g(M_e, M_g, g_size);
      /************************************BUILIDNG MASS/timefactor MATRIX*****************************/
      MAT_SCALAR_MULT( M_g, g_size, timefac);
      /********************************************************************************************/
 			// Solve System Ax=B : A is DM1 , B = Fn-DM2un-DM1un-1 , x=un+1
  		  	int* 	  ipiv = new int [g_size];
        	F77NAME(dgbsv) ( g_size, 0, 0, 1, M_g, 1, ipiv, F_g, g_size, &info);
      //Updating
      INIT_ZERO(g_size,u_n_1); 
 			VEC_SUM(u_n_1, 0, u_n , 1, g_size);
      INIT_ZERO(g_size,u_n); 
 			VEC_SUM(  u_n, 0, F_g , 1, g_size);
     	/***************************Storing Variables*****************************************/   
        	int counter = 1;
        	for(int i=1 ; i<g_size ; i+=3)
        	{
        		storage[count*(N+1)+counter]=F_g[i];
        		counter++;
        	}
      // osc[10000*out + count] = u_n[3*(N/2) -2];
 		}
  //     INIT_ZERO(g_size,u_n);
  //     INIT_ZERO(g_size,u_n_1);
  // }
  //   ofstream filetemp("Task2_Osc.txt");
  //   for(int i= 0 ; i<100 ; i++){
  //     for(int j= 0 ; j<10000 ; j++){
  //       filetemp.width(15);
  //       filetemp << osc[i*10000 + j] << ",";
  //     }
  //     filetemp<< endl;
  //   }
  //   filetemp.close();

    // Writing Displacement of every node at every timestep into text file
 		ofstream file2("Dynamic_Explicit_displacements.txt");
    	if (!file2.good())
    	{
     	   cout << "Error: unable to open output file: Dynamic_Explicit_displacements.txt" << endl;
    	}
    	else 
    	{	
    		for (int i = 0; i < N+1 ; i++)
    		{
				  for (int j = 0; j < Not ; j++)
				  {
					 file2.width(15);
					 file2 << storage[j*(N+1)+i] << ",";
				  }
				  file2 << endl;
			  }
			  cout << endl;
    		file2.close();
        cout << "Written file: Dynamic_Explicit_displacements.txt" << endl;
    	}
    }// END of EXPLICIT

  /******************************************************************************/
 /***************************** IMPLICIT ***************************************/
/******************************************************************************/
  	if(choice2 == 1)
  	{
     // Constants for Newmarks method
 		 double beta  = 1.0/4.0;
 		 double gamma = 1.0/2.0;
 		 double beta_timefac = timefac/beta;
 		
 		 double* BAND_K_g = new double [13*(g_size)]();
 		 double* Col_K_g  = new double [13*(g_size)]();
     // Building Effective stiffness matrix
     MAT_SCALAR_MULT( M_g, g_size, beta_timefac);
 		 BAND_stiffness_g( BAND_K_g, BAND_K_e, N, 4);
 		 RowtoCol(BAND_K_g, g_size, 13, Col_K_g);
 		 MAT_VEC_SUM_COL_BAND(1, Col_K_g, 13, g_size, 1, M_g, g_size, 8);
 		/***************************************************************************************/
 		// Declaring variables for solution
 		 double* F_g      = new double[g_size]();

 		 double* u_n      = new double[g_size]();
 		 double* accel    = new double[g_size]();
 		 double* vel      = new double[g_size](); 
     
     // Declaring copies
     double*   u_n_copy1 = new double[g_size]();
 		 double* accel_copy1 = new double[g_size]();
 		 double*   vel_copy1 = new double[g_size]();
     double*   u_n_copy2 = new double[g_size]();
 		 double* accel_copy2 = new double[g_size]();
 		 double*   vel_copy2 = new double[g_size]();
     double*   u_n_copy3 = new double[g_size]();
     double* accel_copy3 = new double[g_size]();
     double*   vel_copy3 = new double[g_size]();
     double*    M_g_copy = new double[g_size]();
     double*   u_n_p1_copy1 = new double[g_size]();
     double* accel_p1_copy1 = new double[g_size]();

     // Declaring storage matrix
 		 double* storage_u     =new double[Not*(N+1)]();

     // Variables necessary for DGBSV
 		 int     info = 0;
 		 int* 	 ipiv = new int [g_size];
    // double* osc = new double [100*10000];
    // for(int out=1 ; out<100 ; out++)
    // {
     // loading time variables for oscillation test
     double load_time_factor=1;
     double load_factor=1/load_time_factor;
 		 for (int count = 0 ; count*del_T < T ; count++)
 		 {
        /********************COPIES AT TIME N************************************/
        INIT_ZERO( g_size,   u_n_copy1);
        INIT_ZERO( g_size, accel_copy1);
        INIT_ZERO( g_size,   vel_copy1);
        INIT_ZERO( g_size,   u_n_copy2);
        INIT_ZERO( g_size, accel_copy2);
        INIT_ZERO( g_size,   vel_copy2);
        INIT_ZERO( g_size,   u_n_copy3);
        INIT_ZERO( g_size, accel_copy3);
        INIT_ZERO( g_size,   vel_copy3);
        VEC_SUM(  u_n_copy1, 0,   u_n, 1, g_size); 
        VEC_SUM(accel_copy1, 0, accel, 1, g_size); 
        VEC_SUM(  vel_copy1, 0,   vel, 1, g_size);
        VEC_SUM(  u_n_copy2, 0,   u_n, 1, g_size); 
        VEC_SUM(accel_copy2, 0, accel, 1, g_size); 
        VEC_SUM(  vel_copy2, 0,   vel, 1, g_size);
        VEC_SUM(  u_n_copy3, 0,   u_n, 1, g_size); 
        VEC_SUM(accel_copy3, 0, accel, 1, g_size); 
        VEC_SUM(  vel_copy3, 0,   vel, 1, g_size); 
        /*****************************Global Mass Matrix**************************/
        INIT_ZERO( g_size,   M_g);
        MAT_Mass_g(M_e, M_g, g_size);
        VEC_SUM(  M_g_copy, 0,   M_g, 1, g_size);
        /*************************Effective Stiffness MATRIX***************************/
        INIT_ZERO( g_size,   BAND_K_g);
        INIT_ZERO( g_size,    Col_K_g);
        MAT_SCALAR_MULT( M_g_copy, g_size, beta_timefac);
        BAND_stiffness_g( BAND_K_g, BAND_K_e, N, 4);
        RowtoCol(BAND_K_g, g_size, 13, Col_K_g);
        MAT_VEC_SUM_COL_BAND(1, Col_K_g, 13, g_size, 1, M_g_copy, g_size, 8);
        /***************************** Force *********************************/
        if(count*del_T < load_time_factor*T)
        {
          INIT_ZERO(g_size,F_g);
          VEC_DYN_force_e(F_e,   l, qx, qy, Nt, del_T*(count+1)*load_factor);
          VEC_DYN_force_g(F_g, F_e,  N, Nt, del_T*(count+1)*load_factor);
        }
        else
        {
          INIT_ZERO(g_size,F_g);
          VEC_force_e(F_e,   l, qx, qy);
          VEC_force_g(F_g, F_e,  N);
        }
 			  /***************************-Displacement-******************************/
        // Building RHS
 			  MAT_SCALAR_MULT( accel_copy1, g_size,  ((1/(2*beta))-1) );     
        MAT_SCALAR_MULT(   vel_copy1, g_size,  (1/(del_T*beta)) );     
        MAT_SCALAR_MULT(   u_n_copy1, g_size,  beta_timefac  );     
 			  VEC_SUM(vel_copy1, 1, accel_copy1, 1, g_size );
 			  VEC_SUM(u_n_copy1, 1, vel_copy1, 1, g_size );
 			  VEC_VEC_MULT(M_g, u_n_copy1, g_size);// writes to mass matrix
 			  VEC_SUM(F_g, 1, M_g, 1, g_size );
        // Solving System
 			  F77NAME(dgbsv)( g_size, 4, 4, 1, Col_K_g, 13, ipiv, F_g, g_size, &info);
        // Updating un+1 with solution
        VEC_SUM(u_n, 0, F_g, 1, g_size); 
        INIT_ZERO( g_size, u_n_p1_copy1);
        VEC_SUM( u_n_p1_copy1, 0,   u_n, 1, g_size); 

        // Storing Displacement
          int counter=1;
          for(int i=1 ; i<g_size ; i+=3)
          {
            storage_u[count*(N+1)+counter]=u_n[i];
            counter++;
          }
 			/******************************-Acceleration-*************************/
 			  //Building RHS
 			  MAT_SCALAR_MULT(accel_copy2, g_size, ((1/(2*beta))-1));	  
 			  MAT_SCALAR_MULT(  vel_copy2, g_size, 1/(beta*del_T));	    
 			  VEC_SUM        (u_n_p1_copy1, 1, u_n_copy2, -1, g_size);  
 			  MAT_SCALAR_MULT(u_n_p1_copy1, g_size, beta_timefac);
 			  VEC_SUM(  vel_copy2, -1, accel_copy2, -1, g_size);
 			  VEC_SUM(u_n_p1_copy1, 1, vel_copy2, 1, g_size);

 			  VEC_SUM(accel, 0 , u_n_p1_copy1, 1, g_size);//copying new accel into accel
        // copy of new acceleration for use in finding vel
        INIT_ZERO( g_size, accel_p1_copy1);
        VEC_SUM(  accel_p1_copy1, 0,   accel, 1, g_size);

 			  /***************************-Velocity-********************************/
 			  // Building RHS
        MAT_SCALAR_MULT(accel_p1_copy1, g_size, del_T*gamma);
 			  MAT_SCALAR_MULT(accel_copy3, g_size, del_T*(1-gamma));
 			  VEC_SUM(accel_copy3, 1, accel_p1_copy1, 1, g_size);
 			  VEC_SUM(vel_copy3, 1, accel_copy3, 1, g_size);

        VEC_SUM(vel, 0 , vel_copy3, 1, g_size);// velocity being updated
        // osc[10000*out + count] = u_n[3*(N/2) -2];
 		  }
    //   INIT_ZERO(g_size,u_n);
    //   INIT_ZERO(g_size,accel);
    //   INIT_ZERO(g_size,vel);
    // }
    // ofstream filetemp("Task3_Osc.txt");
    // for(int i= 0 ; i<100 ; i++){
    //   for(int j= 0 ; j<10000 ; j++){
    //     filetemp.width(15);
    //     filetemp << osc[i*10000 + j] << ",";
    //   }
    //   filetemp<< endl;
    // }
    // filetemp.close();
     
     // Writing to file, displacement of each node at each timestep
 		 ofstream file3("Dynamic_Implicit_displacements.txt");
    		if (!file3.good())
    		{
     	   		cout << "Error: unable to open output file: Dynamic_Implicit_displacements.txt" << endl;
    		}
    		else 
    		{	
    			for (int i = 0; i < N+1 ; i++)
    			{
					 for (int j = 0; j < Not ; j++)
            {
              file3.width(15);
              file3 << storage_u[j*(N+1)+i] << ",";
            }
           file3 << endl;
				  }
				  cout << endl;
    			file3.close();
        	cout << "Written file: Dynamic_Implicit_displacements.txt" << endl;
        }
  	}// END OF IMPLICIT
}

  /***********************************************************************************************/
 /********************************************PARALLEL*******************************************/
/***********************************************************************************************/
if(choice1==2)
{
  // Characterising running conditions
  int p;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  const int pN = N/p + 2 ;  // number of local gridpoint for solvable matrix
        int pg_size = 3*pN; // size of vectors in each processor
  if(p<2){pg_size=3*(N-1);}


  MPI_Barrier(MPI_COMM_WORLD);
  if(choice2==0)
  {
    /****************************EXPLICIT*************************************************/
    //Building Mass MAtrix
    double* pM_e = new double[6];
    double* pM_g = new double[pg_size];
    double* pM_g_inv = new double[pg_size];

    MAT_Mass_e(pM_e, A, rho, l);
    MAT_Mass_g( pM_e, pM_g, pg_size);

    // Building Stiffness Matrix for each process
    double* pK_g     = new double[9*pg_size];
    double* pCol_K_g = new double[9*pg_size];
    BAND_stiffness_g(pK_g, BAND_K_e, pN+1, 0);
    
    for(int i=0; i<pg_size ; i++)
    {
      pM_g[i]*=timefac;
      pM_g_inv[i]=1/pM_g[i];
      pK_g[4*pg_size + i]+=(-2*pM_g[i]);
    }
    RowtoCol(pK_g, pg_size, 9, pCol_K_g);

    double* u_n   = new double [pg_size];
    double* u_n_1 = new double [pg_size];
    double* inbox = new double [pg_size];
    double* post  = new double [      6];
    double* temp2 = new double [pg_size];
    double* F_g   = new double [pg_size];
    double* result= new double [ g_size];

    INIT_ZERO(pg_size,u_n);
    INIT_ZERO(pg_size,u_n_1);
    INIT_ZERO(pg_size,inbox);
    INIT_ZERO(pg_size,temp2);
    INIT_ZERO(pg_size,result);

    double load_time_factor=1.0;
    double load_factor=1.0/load_time_factor;
    int info;

    for (int count = 0 ; count*del_T < T ; count++)
    {
      // Building Force Vector
      INIT_ZERO(pg_size,F_g);
      VEC_DYN_force_e(F_e,   l, qx, qy, Nt, del_T*(count)*load_factor);
      if(p!=2){VEC_DYN_force_g(F_g, F_e, N,  Nt, del_T*(count)*load_factor);} // for single Process
      else
      {
        PARA_VEC_DYN_force_g(F_g, F_e, pN+1, Nt, del_T*(count)*load_factor);  // multiple processes
        if(rank==0)F_g[pg_size - 8] = F_g[pg_size - 8]+Fy*del_T*(count)*load_factor;
        if(rank==1)F_g[          7] = F_g[          7]+Fy*del_T*(count)*load_factor;
      }

      //ReBuilding Mass Matrix
      INIT_ZERO(pg_size,pM_g);
      MAT_Mass_g( pM_e, pM_g, pg_size);
      MAT_SCALAR_MULT(pM_g, pg_size, timefac);

      // Building RHS
      VEC_VEC_MULT(pM_g, u_n_1, pg_size);//rewrites M_g keeps un-1
      F77NAME(dgbmv)('N', pg_size, pg_size, 4, 4, 1, pCol_K_g, 9, u_n, 1, 0, temp2, 1); //writes to temp2 and keeps un
      VEC_SUM(temp2,  -1,   pM_g,   -1, pg_size);// writes to temp2
      VEC_SUM(  F_g,  1,  temp2,  1, pg_size);// writes to F_g : now RHS of system, and keeps temp2

      //SOLVING AND UPDATING
      INIT_ZERO(pg_size,u_n_1); 
      VEC_SUM(u_n_1, 0, u_n , 1, pg_size);// updates un-1
      VEC_VEC_MULT(F_g, pM_g_inv, pg_size); //Solving System
      INIT_ZERO(pg_size,u_n);
      VEC_SUM(u_n, 0, F_g , 1, pg_size); //updates u_n

      // Sharing results dependent on number of processes used
      if(p==2)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==0)
        {
          for (int i = 0 ; i < 6 ; i++)
          {
            post[i]= u_n[pg_size-15+i];
          }
          MPI_Send(post, 6, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
          MPI_Recv(inbox, 6, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          for (int i = 0 ; i < 6 ; i++)
          {
            u_n[ pg_size-6 + i ] = inbox[i];
          }
        }  
        if(rank==1)
        {
          for (int i = 0 ; i < 6 ; i++)
          {
            post[i] = u_n[ 9 + i];
          }
          MPI_Recv(inbox, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(post, 6, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

          for (int i = 0 ; i < 6 ; i++)
          {
            u_n[i]= inbox[i];
          }
        }
      }
    } // End of time loop

    // Passing results dependent on number of processes used and writes to file
    if(p==2)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==1){MPI_Send(u_n, pg_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);}

      if(rank==0)
      {
        MPI_Recv(inbox, pg_size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for(int i = 1; i < g_size ; i++)
        {
          if(i<=pg_size){result[i]=u_n[i];}

          if(i>pg_size ){result[i]=inbox[i-(pg_size-15)];}
        }
        ofstream file4("Parallel_Explicit_displacements.txt");
        if (!file4.good()) 
        {
          cout << "Error: unable to open output file: Parallel_Explicit_displacements.txt" << endl;
        }
        else 
        {
          file4 << 0 << "," << endl;
          for (int i = 1 ; i < g_size+1 ; i+=3) 
          {
            file4.width(10);
            file4 << result[i]<< "," << endl;
          }
          file4 << 0 << "," << endl;
          file4.close();
          cout << "Written file: Parallel_Explicit_displacements.txt" << endl;
        }
      }
    }
    else
    {
      ofstream file4("Parallel_Explicit_displacements_SINGLE.txt");
        if (!file4.good()) 
        {
          cout << "Error: unable to open output file: Parallel_Explicit_displacements_SINGLE.txt" << endl;
        }
        else 
        {
          file4 << 0 << "," << endl;
          for (int i = 1 ; i < g_size+1 ; i+=3) 
          {
            file4.width(10);
            file4 << u_n[i]<< "," << endl;
          }
          file4 << 0 << "," << endl;
          file4.close();
          cout << "Written file: Parallel_Explicit_displacements_SINGLE.txt" << endl;
        }
    }
  }

  /***********************************************************************************************************************/

  if(choice2==1)
  {
    //characterising problem size for scalapack
    const int full_size = 3*N;
    const int ku = 4;
    const int kl = 4;
          int blk_size = 3*(N/p);
    const int odd_size = 3*(N/p -1);
    const int B_width  = 1+(2*ku)+(2*kl);

    const double beta  = 1.0/4.0;
    const double gamma = 1.0/2.0;
    const double beta_timefac = timefac/beta; 
    
    // Declaring mass elemental
    double* M_e= new double[6]();
    MAT_Mass_e( M_e, A,  rho, l);

    double* F_g = new double [blk_size]();

    double* u_n   = new double[blk_size]();
    double* vel   = new double[blk_size]();
    double* accel = new double[blk_size]();

    // Declaring copies
    double*   u_n_copy1 = new double[blk_size]();
    double* accel_copy1 = new double[blk_size]();
    double*   vel_copy1 = new double[blk_size]();
    double*   u_n_copy2 = new double[blk_size]();
    double* accel_copy2 = new double[blk_size]();
    double*   vel_copy2 = new double[blk_size]();
    double*   u_n_copy3 = new double[blk_size]();
    double* accel_copy3 = new double[blk_size]();
    double*   vel_copy3 = new double[blk_size]();
    
    double*    M_g_copy = new double[blk_size]();

    double*   u_n_p1_copy1 = new double[blk_size]();
    double* accel_p1_copy1 = new double[blk_size]();


    // Variables for testing oscillations
    double load_time_factor = 1;
    double load_factor = 1/load_time_factor;
    int    info;
    
    MPI_Barrier(MPI_COMM_WORLD); // Waiting for all processes to enter loop
    for (int count = 0 ; count*del_T < T ; count++)
    {
      /********************COPIES AT TIME N************************************/
        INIT_ZERO( blk_size,   u_n_copy1);
        INIT_ZERO( blk_size, accel_copy1);
        INIT_ZERO( blk_size,   vel_copy1);
        INIT_ZERO( blk_size,   u_n_copy2);
        INIT_ZERO( blk_size, accel_copy2);
        INIT_ZERO( blk_size,   vel_copy2);
        INIT_ZERO( blk_size,   u_n_copy3);
        INIT_ZERO( blk_size, accel_copy3);
        INIT_ZERO( blk_size,   vel_copy3);
        VEC_SUM(  u_n_copy1, 0,   u_n, 1, blk_size); 
        VEC_SUM(accel_copy1, 0, accel, 1, blk_size); 
        VEC_SUM(  vel_copy1, 0,   vel, 1, blk_size);
        VEC_SUM(  u_n_copy2, 0,   u_n, 1, blk_size); 
        VEC_SUM(accel_copy2, 0, accel, 1, blk_size); 
        VEC_SUM(  vel_copy2, 0,   vel, 1, blk_size);
        VEC_SUM(  u_n_copy3, 0,   u_n, 1, blk_size); 
        VEC_SUM(accel_copy3, 0, accel, 1, blk_size); 
        VEC_SUM(  vel_copy3, 0,   vel, 1, blk_size);

      /*****************************-BUILDING SYSTEM -************************************/ 
      // BUILDING Force
      VEC_DYN_force_e(F_e, l, qx, qy, Nt, del_T*(count+1)*load_factor);
      if(p < 2)
      {
        VEC_DYN_force_g(F_g, F_e, N,  Nt, del_T*(count+1)*load_factor);
      }
      else
      {
        if(rank < p-1)
        {
          PARA_VEC_DYN_force_g(F_g, F_e, (N/p+1), Nt, del_T*(count+1)*load_factor);
          if(rank == p/2-1){ F_g[blk_size - 2] = F_g[blk_size - 2] + Fy*del_T*(count+1)*load_factor;}
        }
        if(rank == p-1)
        {
         PARA_VEC_DYN_force_g(F_g, F_e, (N/p), Nt, del_T*(count+1)*load_factor);
        }
      }

      /************************* BUILDING MASS *******************************/
      double* M_g= new double [blk_size]();
      if(rank<p-1 || p==1)
      {
        MAT_Mass_g(M_e, M_g, blk_size);
      }
      if(rank==p-1 && p>1)
      {
        MAT_Mass_g(M_e, M_g, odd_size);
      }
      MAT_SCALAR_MULT( M_g, blk_size, beta_timefac);

      /***************** BUILDING EFFECTIVE STIFFNESS MATRIX *****************/
      double*     K_g = new double[B_width * blk_size]();
      double* Col_K_g = new double[B_width * blk_size]();
      if(rank < p-1)
      {
        if(p==1){BAND_stiffness_g_PARA( K_g, BAND_K_e,  (N/p), 8);}
        else    {BAND_stiffness_g_PARA( K_g, BAND_K_e,  (N/p+1), 8);}
        for(int i=0; i<blk_size ; i++)
            {K_g[12*blk_size + i] += M_g[i];}
        RowtoCol(K_g, blk_size, B_width, Col_K_g);
      }
      if(rank==p-1 && p>1)
      {
        double* temp1 = new double[B_width*odd_size]();
        double* temp2 = new double[B_width*odd_size]();
        BAND_stiffness_g_PARA( temp1, BAND_K_e,  N/p, 8);

        for(int i=0; i<odd_size ; i++)
            {temp1[12*odd_size + i] += M_g[i];}
        RowtoCol(temp1, odd_size, B_width, temp2);
        row_ones(K_g, blk_size, 12);
        RowtoCol(K_g, blk_size, B_width, Col_K_g);
        for(int i=0 ; i < B_width*odd_size ; i++)
            { Col_K_g[i] = temp2[i]; }
      }
      /*************************** REBUILDING MASS ***************************/
      INIT_ZERO(blk_size, M_g);
      if(rank<p-1)
      {
        MAT_Mass_g(M_e, M_g, blk_size);
      }
      if(rank==p-1 && p>1)
      {
        MAT_Mass_g(M_e, M_g, odd_size);
      }
      /******************************-Displacement-*************************/
      // Building RHS
      MAT_SCALAR_MULT( accel_copy1, blk_size,  ((1/(2*beta))-1) );     
      MAT_SCALAR_MULT(   vel_copy1, blk_size,  (1/(del_T*beta)) );
      MAT_SCALAR_MULT(   u_n_copy1, blk_size,  beta_timefac  );     
      VEC_SUM(vel_copy1, 1, accel_copy1, 1, blk_size );
      VEC_SUM(u_n_copy1, 1,   vel_copy1, 1, blk_size );
      VEC_VEC_MULT(M_g, u_n_copy1, blk_size);// writes to mass matrix
      VEC_SUM(F_g, 1, M_g, 1, blk_size );

      MPI_Barrier(MPI_COMM_WORLD);
      SolveWithScalapack(full_size, Col_K_g, F_g);
      VEC_SUM(u_n, 0, F_g, 1, blk_size); // updating un+1 with solution

      INIT_ZERO( blk_size, u_n_p1_copy1);
      VEC_SUM( u_n_p1_copy1, 0,   u_n, 1, blk_size); 

      /******************************-Acceleration-*************************/
      //Building RHS
      MAT_SCALAR_MULT(accel_copy2, blk_size, ((1/(2*beta))-1));   
      MAT_SCALAR_MULT(  vel_copy2, blk_size, 1/(beta*del_T));     
      VEC_SUM        (u_n_p1_copy1, 1, u_n_copy2, -1, blk_size);  
      MAT_SCALAR_MULT(u_n_p1_copy1, blk_size, beta_timefac);
      VEC_SUM(  vel_copy2, -1, accel_copy2, -1, blk_size);
      VEC_SUM(u_n_p1_copy1, 1, vel_copy2, 1, blk_size);
      VEC_SUM(accel, 0 , u_n_p1_copy1, 1, blk_size);//copying new accel into accel
      
      INIT_ZERO( blk_size, accel_p1_copy1);
      VEC_SUM(  accel_p1_copy1, 0,   accel, 1, blk_size);// copy of new acceleration for use in finding vel

      /***************************-Velocity-********************************/
      // Building RHS
      MAT_SCALAR_MULT(accel_p1_copy1, blk_size, del_T*gamma);
      MAT_SCALAR_MULT(accel_copy3, blk_size, del_T*(1-gamma));
      VEC_SUM(accel_copy3, 1, accel_p1_copy1, 1, blk_size);
      VEC_SUM(vel_copy3, 1, accel_copy3, 1, blk_size);
      VEC_SUM(vel, 0 , vel_copy3, 1, blk_size);// velocity being updated     
    }//End of time loop

    // PARA_PRINT(blk_size,0,u_n,0);

    // Gathering results from Scalapack into process 0
    if(rank > 0 && p>1)
      {
        MPI_Send(u_n, blk_size, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
      }

    if(rank == 0 && p>1)
    {
      double* inbox1  = new double [blk_size]();
      double* inbox2  = new double [blk_size]();
      double* inbox3  = new double [blk_size]();
              MPI_Recv(inbox1, blk_size, MPI_DOUBLE, 1, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      if(p>2){MPI_Recv(inbox2, blk_size, MPI_DOUBLE, 2, 2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}
      if(p>2){MPI_Recv(inbox3, blk_size, MPI_DOUBLE, 3, 3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}
   
    
      double* storage = new double [g_size]();
        for(int i=1 ; i<blk_size ; i+=3)
        {
                storage[             i] = u_n   [i];
        if(p>1){storage[  blk_size + i] = inbox1[i];}
        if(p>2){storage[2*blk_size + i] = inbox2[i];}
        if(p>2){storage[3*blk_size + i] = inbox3[i];}
        }
    
      ofstream file5("PARA_Implicit_displacements.txt");
      if (!file5.good())
      {
        cout << "Error: unable to open output file: PARA_Implicit_displacements.txt" << endl;
      }
      else 
      { 
        file5.width(15);
        file5 << 0 << "," << endl;
        for (int i = 1; i < g_size ; i+=3)
        {
          file5.width(15);
          file5 << storage[i] << "," << endl;
        }
        file5.width(15);
        file5 << 0 << "," << endl;
        cout << endl;
        file5.close();
        cout << "Written file: PARA_Implicit_displacements.txt by RANK :  " << rank << endl;
      }
    }
    else if(p==1)
    {
      ofstream file5("PARA_Implicit_displacements.txt");
      if (!file5.good())
      {
        cout << "Error: unable to open output file: PARA_Implicit_displacements.txt" << endl;
      }
      else 
      { 
        file5.width(15);
        file5 << 0 << "," << endl;
        for (int i = 1; i < g_size ; i+=3)
        {
          file5.width(15);
          file5 << u_n[i] << "," << endl;
        }
        file5.width(15);
        file5 << 0 << "," << endl;
        cout << endl;
        file5.close();
        cout << "Written file: PARA_Implicit_displacements.txt by RANK :  " << rank << endl;
      }
    }
  }
}

MPI_Finalize();
return 0;
}