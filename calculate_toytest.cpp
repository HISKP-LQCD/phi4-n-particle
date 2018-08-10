#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/unistd.h>

#include <iostream>

#include "types.h"
#include "action.h"
#include "parameters.h"
#include "scalar.h"
#include "metropolis.h"
#include "correlators.h"



//###########################################################################################//
// (I.)                                                                                      //
//             Function for the calculation of the parameters LAMBDA and KAPPA               //
//                      (needed for the calculation of the action S):                        //
//                                                                                           //
//###########################################################################################//

double KAPPA, LAMBDA;

void calculate_parameters() {

  LAMBDA = (4*lambda_c - (8+m2_0)*(-8 -m2_0 + sqrt(8*lambda_c + (8+m2_0)*(8+m2_0)))  )/(8*lambda_c);
  KAPPA = (-8 - m2_0 + sqrt(8*lambda_c + (8.0+m2_0)*(8.0+m2_0)) )/(4*lambda_c);

  if(KAPPA<0 || KAPPA>1) {
    
    LAMBDA = (4*lambda_c + (8+m2_0)*(8 +m2_0 + sqrt(8*lambda_c + (8+m2_0)*(8+m2_0)))  )/(8*lambda_c);
    KAPPA = (-8 - m2_0 - sqrt(8*lambda_c + (8.0+m2_0)*(8.0+m2_0)) )/(4*lambda_c);
    
  }
}



//###########################################################################################//
// (II.)                                                                                     //
//             MAIN FUNCTION FOR THE CALCULATION OF SCALAR FIELD CONFIGURATIONS:             //
//                                                                                           //
//      The field configurations "scalar_X_Y_Z_T_(n_conf).txt" are stored in the folder      //
//   located at "path_read" (selectable in "parameters.h") which are read in in the course   //
//         of the calculation of the correlation functions (see "calculate_corr.cpp")        //
//                                                                                           //
//###########################################################################################//

int main(int argc, char *argv[]) {

  double time0=omp_get_wtime( ), time0_b, time1_b;

  scalar_field phi, phi2;                            // scalar_field defined in "types.h"
  double action, action_next, acceptance, phase;
  int i,j=0;
  int n_tot;
  complex corr[T], corr_im[T];
  complex caux, caux2;
  complex corr_mean[T], vev_mean;
  complex corr_mean_im[T];
  char *endptr;    
  int nthreads, tid;

  
  //=========================================================================================//
  // (II.A)                                                                                  //
  // Print the number of saved phi field configurations:                                     //
  //                                                                                         //
  //=========================================================================================//

  printf("\n=====================================================\n");
  printf("\nNumber of saved phi field configurations: %d\n", n_save); // new
 

  //=========================================================================================//
  // (II.B)                                                                                  //
  // Create a folder "test" located at "path_read", where the scalar field configurations    //
  // are printed (II.I,J), from which the correlation functions can be calculated. The       //
  // corresponding path "path_read" can be chosen in "parameters.h".                         //
  //                                                                                         //
  //=========================================================================================//

  // stat checks if the directory exists:
  
  struct stat st = {0};

  // mkdir creates the new directory and the option 0700 allows the owner to read, write and
  // execute the files:
  
  if (stat(path_read, &st)==-1) {
    mkdir(path_read, 0700);
  }


  //=========================================================================================//
  // (II.C)                                                                                  //
  // Parallel Computing via Open MP:                                                         //
  //                                                                                         //
  //=========================================================================================//
  
#pragma omp parallel private(nthreads, tid)
  {
    // Obtain thread number
    tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);

    // Only master thread does this
    if (tid == 0) {
      
	nthreads = omp_get_num_threads();
	printf("Number of threads = %d\n", nthreads);

	printf("T = %d\n", T);
  
	if(nthreads > T/2 || T%nthreads!=0) {
	  printf("Number of threads must be smaller then T/2 and T multiple of nthreads  \n");
	  exit(0);
	}	
     }
  }   // All threads join master thread and disband

  nthreads = omp_get_num_threads();


  //=========================================================================================//
  // (II.D)                                                                                  //
  // Create a file "action.out" ("w") and open it for writing:                               //
  //                                                                                         //
  //=========================================================================================//
  
  FILE * faction;
  faction = fopen("action.out", "w");
  fprintf(faction, "step of %d\n", nprint_field);

  
  //=========================================================================================//
  // (II.E)                                                                                  //
  // Calculate the parameters Lambda and Kappa:                                              //
  //                                                                                         //
  //=========================================================================================//
  
  calculate_parameters();
  
  printf("START \n");
  printf("Lambda_c = %f, mass^2_0 = %f \n", lambda_c, m2_0);
  printf("Lambda = %f, Kappa = %f \n", LAMBDA, KAPPA);

  
  srand (clock());


  //=========================================================================================//
  // (II.F)                                                                                  //
  // Initialize the 2 fields needed in the course of the metropolis algorithm (The           //
  // initialization process is controlled by the function initialize_field() defined in      //
  // "scalar.cpp"):                                                                          //
  //                                                                                         //
  //=========================================================================================//
  
  initialize_field(&phi);
  initialize_field(&phi2);


  //=========================================================================================//
  // (II.G)                                                                                  //
  // If start_random is set to 0 (this is done in "parameters.h") the file start_conf (see   //
  // "parameters.h") is read in as start configuration called phi (The operating function    //
  // fread_field is defined in "scalar.cpp"):                                                //
  //                                                                                         //
  //=========================================================================================//
  
  if(start_random == 0) {
    
    printf("start config for phi is %s \n", start_conf);
    fread_field(start_conf, &phi);
    
    printf("start action = %f\n", eval_action_nogauge(phi));
    printf("Lambda is %f and Kappa is %f \n", LAMBDA,KAPPA);    

  }


  //=========================================================================================//
  // (II.H)                                                                                  //
  // If start_random is set to 1 (this is done in "parameters.h") a hot start is performed,  //
  // where a disordered, random configuration is utilized for the calculation of the field   //
  // configurations (II.I):                                                                  //
  //                                                                                         //
  //=========================================================================================//
  
  
  if(start_random == 1) {
    
    printf("start action = %f\n", eval_action_nogauge(phi));
    printf("\n=====================================================\n");
    
    acceptance = metropolis(&phi, n_term_field, faction);
    printf("acceptance: %f \n", acceptance);
  
  }


  clock_t endconf,startconf;
  double time_spentconf;


  //=========================================================================================//
  // (II.I)                                                                                  //
  // Create n_save configuration files "scalar_X_Y_Z_T_(n_conf).txt" which are read in in    //
  // the ROUTINE FOR THE CALCULATION OF THE CORRELATION FUNCTIONS (see "calculate_corr.cpp")://
  //                                                                                         //
  //=========================================================================================//
  
  for(i=0;i<n_save;i++) { // n_save is the number of configuration files saved at "path_read"

    time0_b = omp_get_wtime( );
    
    acceptance = metropolis(&phi, n_term_save, faction);
    printf("out acceptance: %f \n", acceptance);

    
    //=======================================================================================//
    // (II.J)                                                                                //
    // Print the ith field configuration (real and imaginary part of phi at all possible     //
    // field points (t,x,y,z) of the lattice) into a file opened by fprint_field()           //
    // (see "scalar.cpp"):                                                                   //
    //                                                                                       //
    //=======================================================================================//
    
    fprint_field(phi, (long long) (i+1)*n_term_save + (long long) start_random_conf*(1-start_random));
    
    
    time1_b = omp_get_wtime( );
    printf("Duration %f seconds \n", time1_b-time0_b);
  }
  
  double time1=omp_get_wtime( );
  printf("Duration %f seconds \n", time1-time0);
}
