#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#include <random>
#include <iostream>

#include "action.h"
#include "parameters.h"
#include "types.h"
#include "generator_singleton.h"



//*******************************************************************************************//
//                                                                                           //
// NOTE:                                                                                     //
//                                                                                           //
// Lattice points:                                                                           //
// A (complex) scalar field (e.g. p_aux in (II.)) is only defined at discrete lattice points //
// "lattice_point(t,x,y,z)" (see (I.)).                                                      //
//                                                                                           //
// Field points:                                                                             //
// The complex "scalar_field" value "aux[lattice_point(t,x,y,z)].re" is the real part and    //
// "aux[lattice_point(t,x,y,z)].im" the imag. part of the field p_aux at the lattice point   //
// "lattice_point(t,x,y,z)". For the typedef of the complex "scalar_field" see "types.h".    //
//                                                                                           //
//*******************************************************************************************//



//###########################################################################################//
// (I.)                                                                                      //
//     The integer "lattice_point" contains the complete information of the lattice point    //
//    4-vector with components t,x,y,z and is utilized to indicate discrete scalar fields    //
//                                  (e.g. p_aux in (II.)):                                   //
//                                                                                           //
//###########################################################################################//

int lattice_point(int t, int x, int y, int z) {
  
  int lattice_point;
  int t_aux, x_aux, y_aux, z_aux;

  t_aux = (T+t)%T;
  x_aux = (X+x)%X;
  y_aux = (Y+y)%Y;
  z_aux = (Z+z)%Z;

  lattice_point = x_aux * Y*Z*T + y_aux * Z*T + z_aux * T + t_aux;
  
  return lattice_point;
}



//###########################################################################################//
// (II.)                                                                                     //
//           Random initialization in "metropolis.cpp", "calculate_toytest.cpp" and          //
//     "calculate_corr.cpp" of the lattice field phi (and phi2) consisting of two arrays     //
//            containing the real and imaginary components for all lattice points            //
//                          "lattice_point(t,x,y,z)", respectively:                          //
//                                                                                           //
//###########################################################################################//

int initialize_field(scalar_field *p_aux) {   // typedef of scalar_field in types.h
  
  //=========================================================================================//
  //                                                                                         //
  // By means of malloc, the memory for the array size volume = X*Y*Z*T is reserved. This    //
  // means that it can outlive the scope where it has been created and therefore also can    //
  // be used in other scopes of this program.                                                //
  //                                                                                         //
  //=========================================================================================//
  
  *p_aux = (scalar_field)  malloc(volume * sizeof(complex));
  
                                                               
  scalar_field aux = *p_aux; // typedef of scalar_field in "types.h"
  double random;
  int t,x,y,z;

  
  // "ZeroOne_distribution()" generates real numbers uniformly distr. between 0 and 1:
  auto ZeroOne_distribution = std::uniform_real_distribution<double>(0.0,1.0);
  
  
  for(t=0;t<T;t++) {
    for(x=0;x<X;x++) {
      for(y=0;y<Y;y++) {
	for(z=0;z<Z;z++) {

	  // get() is the thread number dependend mersenne twister defined in the class
	  // "GeneratorSingleton" in "generator_singleton.h":
	  random = 2 * PI * ZeroOne_distribution(GeneratorSingleton::get());
	  
	  aux[lattice_point(t,x,y,z)].re =  cos(random) ; // array containing real components
	  aux[lattice_point(t,x,y,z)].im =  sin(random) ; // array containing imag. components

	}
      }
    }
  }
}



//###########################################################################################//
// (III.)                                                                                    //
//        This function copies one single complex field point by taking a scalar field       //
//        point p_old and replacing it with p_new (This is applied in the MC algorithm       // 
//       "metropolis.cpp" (I.) and (II.) when an updated field point phi2 is accepted):      //
//                                                                                           //
//###########################################################################################//

int copy_field_point(scalar_field *p_old, scalar_field *p_new, int t, int x, int y, int z) {

  scalar_field aux   = *p_old; // typedef of scalar_field in "types.h"
  scalar_field aux2  = *p_new;	    

  aux[lattice_point(t,x,y,z)].re = aux2[lattice_point(t,x,y,z)].re;
  aux[lattice_point(t,x,y,z)].im = aux2[lattice_point(t,x,y,z)].im;
}



//###########################################################################################//
// (IV.)                                                                                     //
//    This function copies an entire field since the for loops over t,x,y,z are included.    //
//       It replaces all scalar field points in space-time p_old and with p_new. This        //
//     function is applied in the MC algorithm in the course of the initialization of the    //
//                         entire field (see "metropolis.cpp" (III.))                        //
//                                                                                           //
//###########################################################################################//

int copy_field(scalar_field *p_old, scalar_field *p_new) {

  scalar_field aux  = *p_old; // typedef of scalar_field in "types.h"
  scalar_field aux2 = *p_new;
  
  int t,x,y,z;
  double temp, temp2;
  
  //  #pragma omp parallel for private(t,x,y,z) shared(aux,aux2,temp,temp2)
  for(t=0;t<T;t++) {
    for(x=0;x<X;x++) {
      for(y=0;y<Y;y++) {
	for(z=0;z<Z;z++) {

	  temp  = aux2[lattice_point(t,x,y,z)].re;
	  temp2 = aux2[lattice_point(t,x,y,z)].im;
	  aux[lattice_point(t,x,y,z)].re = temp;
	  aux[lattice_point(t,x,y,z)].im = temp2;
	}
      }
    }
  }
}



//###########################################################################################//
// (V.)                                                                                      //
//         This function updates one single field point by calculating a small change        //
//               of the field point *p_aux (It is applied in "metropolis.cpp"):              //
//      "ZeroOne_distribution" generates real numbers uniformly distr. between 0 and 1.      //
//                                                                                           //
//###########################################################################################//

void update_field_point(scalar_field *p_aux, int t, int x, int y, int z) {

  int i;
  double phase, module;
  complex random, caux;
  scalar_field aux = *p_aux; // typedef of scalar_field in "types.h"

  
  // "ZeroOne_distribution()" generates real numbers uniformly distr. between 0 and 1:
  auto ZeroOne_distribution = std::uniform_real_distribution<double>(0.0,1.0);

  // get() is the thread number dependend mersenne twister defined in the class
  // "GeneratorSingleton" in "generator_singleton.h":
  random.re = aux[lattice_point(t,x,y,z)].re - deltarho + 2*deltarho * ZeroOne_distribution(GeneratorSingleton::get());  
  random.im = aux[lattice_point(t,x,y,z)].im - deltarho + 2*deltarho * ZeroOne_distribution(GeneratorSingleton::get());
  
  aux[lattice_point(t,x,y,z)].re = random.re;
  aux[lattice_point(t,x,y,z)].im = random.im;

}



//###########################################################################################//
// (VI.)                                                                                     //
//      In "calculate_toytest.cpp" (II.I,J) this function builds the field configuration     //
//    files "scalar_X_Y_Z_T_(n_conf).txt" which are read in in "calculate_corr.cpp" (II.).   //
//      It prints all possible combinations of the lattice components t,x,y,z as well as     //
//     the real and imaginary part of the scalar field "phiaux" at these points into the     //
//                                opened file at "path_read":                                //  
//                                                                                           //
//###########################################################################################//

int fprint_field(scalar_field phiaux, long long n_conf) {

  char filename[40]="";
  int x,y,z,t;
  double real, imaginary;
  FILE * fs;

  //=========================================================================================//
  //                                                                                         //
  // Create a file called "filename" ("w") located at "path_read" and open it for writing:   //
  //                                                                                         //
  //=========================================================================================//
  
  sprintf(filename,"%s/scalar_%d_%d_%d_%d_%li.txt", path_read, X,Y,Z,T,n_conf);          
  printf("%s\n", filename);
  fs = fopen(filename, "w");

  for(t=0;t<T;t++) {
    for(x=0;x<X;x++) {
      for(y=0;y<Y;y++) {
	for(z=0;z<Z;z++) {

	  real = phiaux[lattice_point(t,x,y,z)].re;
	  imaginary = phiaux[lattice_point(t,x,y,z)].im;

	    //===============================================================================//
	    //                                                                               //
	    // Print line by line x,y,z,t,real,imaginary into the opened file:               //
	    //                                                                               //
	    //===============================================================================//
	    
	    fprintf(fs,"%d %d %d %d %f %f\n",x,y,z,t,real,imaginary);

	}
      }
    }
  }
  fclose(fs);
  return 0;
}


//###########################################################################################//
// (VII.)                                                                                    //
//  This functions is utilized in "calculate_toytest.cpp" to read in the start configuration //
//       (start_conf) stored in a txt-file whose path can be chosen in "parameters.h".       //
//     Furthermore it is used in "calculate_corr.cpp" to read in the configuration files     //
//            created with fprint_field() or alternatively provided configuration.           //
//                                                                                           //
//###########################################################################################//

int fread_field(const char *filename, scalar_field *p_phi) {

  FILE* file = fopen(filename, "r");   // Opens the file "filename" which should be read in
  char line[100];
  //char line [BUFSIZ];
  int x,y,z,t;
  double real, imaginary;
  scalar_field phi = *p_phi;

  while(fgets(line, sizeof(line), file)) { // The function fgets() reads in a line of the
                                           // textfile "file" and stores this string in an
                                           // array of chars pointed to by "line". It will
                                           // stop when (n-1) characters are read in.
    
    sscanf(line,"%d %d %d %d %lf %lf",&x,&y,&z,&t,&real,&imaginary);
    phi[lattice_point(t,x,y,z)].re = real;
    phi[lattice_point(t,x,y,z)].im = imaginary;
    
  }
  fclose(file);
  return 0;
}
