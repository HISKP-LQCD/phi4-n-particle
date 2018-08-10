#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/unistd.h>

#include "complex.h"
#include "types.h"
#include "action.h"
#include "parameters.h"
#include "scalar.h"
#include "metropolis.h"
#include "correlators.h"



double KAPPA,LAMBDA;

static FILE* save_fopen(const char filename[]) {
  
  FILE* f = fopen(filename, "w");
  if(f == NULL){
    printf("Failed to open %s\n", filename);
    exit(1);
  }
  
return f;
}



//###########################################################################################//
// (I.)                                                                                      //
//              Function for the calculation of the parameters Lambda and Kappa              //
//(needed for the calculation of the action S; here: metadata for the correlation functions)://
//                                                                                           //
//###########################################################################################//

void calculate_parameters() {
  
  LAMBDA = (4*lambda_c - (8+m2_0)*(-8 -m2_0 + sqrt(8*lambda_c + (8+m2_0)*(8+m2_0))))/(8*lambda_c);
  KAPPA = (-8 - m2_0 + sqrt(8*lambda_c + (8.0+m2_0)*(8.0+m2_0)))/(4*lambda_c);
  
  if(KAPPA<0 || KAPPA>1) {
    
    LAMBDA = (4*lambda_c + (8+m2_0)*(8 +m2_0 + sqrt(8*lambda_c + (8+m2_0)*(8+m2_0))))/(8*lambda_c);
    KAPPA = (-8 - m2_0 - sqrt(8*lambda_c + (8.0+m2_0)*(8.0+m2_0)))/(4*lambda_c);    
  }
}



//###########################################################################################//
// (II.)                                                                                     //
//                 ROUTINE FOR THE CALCULATION OF THE CORRELATION FUNCTIONS:                 //
//                                                                                           //
//###########################################################################################//

int main() {
  
  double time0=omp_get_wtime( );
  
  char corr_filename_n[40] = "";
  char filename_sc[70] = "";
  
  //int n = n_fields; // This is only needed if one omits the for loop over n
  long n;
  int i,j;
  complex corr_n[T];
  
  scalar_field phi;
  
  clock_t begin = clock();
  
  
  //=========================================================================================//
  // (II.A)                                                                                  //
  // Create a folder "analysis" (if it does not exist) located at "path_corr", where the     //
  // output files (correlation functions) are stored in. The corresponding path "path_corr"  //
  // can be chosen in  "parameters.h".                                                       //
  //                                                                                         //
  //=========================================================================================//
  
  // stat checks if the directory already exists:
  struct stat st = {0};
  
  // mkdir creates the new folder analysis (path chosen in parameters.h) and the option 0700
  // allows the owner to read, write and execute the files:
  if (stat(path_corr, &st)==-1) {
    mkdir(path_corr, 0700);
  }
  
  
  //=========================================================================================//
  // (II.B)                                                                                  //
  // Parallel Computing via Open MP:                                                         //
  //                                                                                         //
  //=========================================================================================//
  
  int nthreads, tid;
  // Fork a team of threads giving them their own copies of variables
#pragma omp parallel private(nthreads, tid)
  {
    // Obtain thread number
    tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);
    
    // Only master thread does this
    if (tid == 0) {
      
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n", nthreads);
    }
  } // All threads join master thread and disband
 
  
  //=========================================================================================//
  // (II.C)                                                                                  //
  // Initialize field (see "scalar.cpp") and calculate KAPPA, LAMBDA:                        //
  //                                                                                         //
  //=========================================================================================//
  
  initialize_field(&phi);
  calculate_parameters();

  
  //=========================================================================================//
  // (II.D)                                                                                  //
  // Create (fopen(...,"w")) and open the metadata file for writing in the folder "analysis" //
  // located at path_corr (see "parameters.h"). The file ending ".tsv" stands for            //
  // "tab-separated values". Then, print the metadata file including Kappa, Lambda as well   //
  // as T,X,Y,Z,n_analyse:                                                                   //
  //                                                                                         //
  //=========================================================================================//

  char path_corr_aux[80];
  sprintf(path_corr_aux, path_corr);
  
  char meta_file[30]  = "metadata_conf.tsv";
  char path_meta[120];
  strcpy(path_meta, path_corr_aux);
  strcat(path_meta, meta_file);

  FILE *mconf  = save_fopen(path_meta);

  fprintf(mconf,"%f %f %f \n", LAMBDA,KAPPA);
  fprintf(mconf,"%d %d %d %d %d", X,Y,Z,T,n_analyse);
  
  
  //=========================================================================================//
  // (II.E)                                                                                  //
  // Here, one creates (fopen(...,"w")) and opens n_fields-times a correlator file for the   //
  // n particle correlation function (1<=n<=n_fields) as well as the metadata file for       //
  // writing in the folder "analysis" located at path_corr (see "parameters.h").             //
  // The file ending ".tsv" stands for "tab-separated values":                               //
  //                                                                                         //
  //=========================================================================================//

  for(n=0;n<n_fields;n++) {
    
    char path_corr_aux[80];
    
    sprintf(path_corr_aux, path_corr);
    
    char path_corr_n[120];
    sprintf(corr_filename_n, "correlators_%d_phi_phi4p.tsv", (n+1));
    strcpy(path_corr_n, path_corr_aux);
    strcat(path_corr_n, corr_filename_n);
    
    
    FILE *fconf_n = save_fopen(path_corr_n);
    
    
    //=======================================================================================//
    // (II.F)                                                                                //
    // Print parameters KAPPA, LAMBDA as well as T,X,Y,Z,n_analyse into the header of        //
    // correlation function output files:                                                    //
    //                                                                                       //
    //=======================================================================================//
    
    fprintf(fconf_n,"# Number of particles n_fields=%d \n", (n+1));
    fprintf(fconf_n,"# LAMBDA=%f KAPPA=%f %f \n", LAMBDA,KAPPA);
    fprintf(fconf_n,"# X=%d Y=%d Z=%d T=%d n_analyse=%d \n", X,Y,Z,T,n_analyse);
    fprintf(fconf_n,"Point Re Im \n"); 
    
    
    //=======================================================================================//
    // (II.G)                                                                                //
    // n_analyse-times txt-files containing the start configurations for the fields phi      //
    // called "scalar_X_Y_Z_T_(i+1)*n_term_save" are read in, which are located at           //
    // "path_read" (the number n_analyse and the path, where the start configurations are    //
    // stored ("path_read") can be chosen in "parameters.h").                                // 
    // By means of those fields phi, the correlation functions are build:                    //
    //                                                                                       //
    //=======================================================================================//
    
    for(i=0;i<n_analyse;i++) {
      if(i%n_restrict==0) {
	
	sprintf(filename_sc,"%sscalar_%d_%d_%d_%d_%li.txt", path_read, X,Y,Z,T, (long long) (i+1)*n_term_save);
	fread_field(filename_sc,&phi);      
	
#pragma omp parallel for private(j)
	
	// Since the field configuration arrays phi are provided by "fread_field", the n particle
	// correlation functions "correlator_n()" that are defined in "correlators.cpp" can be
	// calculated:
	
	for(j=0;j<T/2+1;j++) {

	  // Calculate the n particle correlation function:
	  if(correlator == 0) {
	    corr_n[j] = correlator_n(phi,(n+1),j,0,0,0);
	  }

	  // Calculate the derivative of the n particle correlation function:
	  if(correlator == 1) {
	    corr_n[j] = sub_complex(correlator_n(phi,(n+1),j,0,0,0), correlator_n(phi,(n+1),(j+1),0,0,0));
	  }
	}
	
	for(j=0;j<T/2+1;j++) {
	  fprintf(fconf_n,"%d %e %e \n", j, corr_n[j].re, corr_n[j].im);
	}
	
	for(j=T/2+1;j<T;j++) {
	  fprintf(fconf_n,"%d %e %e \n", j, corr_n[T-j].re, corr_n[T-j].im);
	}
	
	if((i+1)%100==0) {
	  printf("Analysed conf number %d \n",i+1);    
	}
      }
      else{}
    }
    
    
    //=======================================================================================//
    // (II.H)                                                                                //
    // Close the files opened in (II.E):                                                     //
    //                                                                                       //
    //=======================================================================================//
    
    fclose(fconf_n);

    
  } // end of the for loop over n
  
  double time1=omp_get_wtime( );
  printf("Duration %f seconds \n",time1-time0);
}
