#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "complex.h"
#include "types.h"
#include "action.h"
#include "parameters.h"
#include "scalar.h"



/* The declarations of complex n particle correlation function correlator_n(...)
   is embedded into "operators.h". */

//####################################################################################//
// (I.)                                                                               //
//     Calculate the n particle correlation function "C_{nphi}(t) = correlator_n".    //
//   (The integer n_aux runs from 1 up to n_fields (see "calculate_corr.cpp"(II.E))   //
//         and is the number of particles in the finite volume, respectively):        //
//                                                                                    //
//####################################################################################//


complex correlator_n(scalar_field phi, long n_aux, int dt, int px, int py, int pz) {
  
  int t1,x,y,z;
  
  complex phi_t_sink, phi_t_source;
  complex exp_ipx;
  complex aux_source, aux_sink;
  complex operator_sink, operator_source;
  complex corr;

  corr.re = 0.;
  corr.im = 0.;

  for(t1=0;t1<T;t1++) {
    
    phi_t_source.re = 0.;
    phi_t_source.im = 0.;
    phi_t_sink.re = 0.;
    phi_t_sink.im = 0.;
    
    for(x=0;x<X;x++) {
      for(y=0;y<Y;y++) {
	for(z=0;z<Z;z++) {
	  
	  exp_ipx.re = cos(2*PI*((double) px*x/X + (double) py*y/Y + (double) pz*z/Z));
	  exp_ipx.im = sin(2*PI*((double) px*x/X + (double) py*y/Y + (double) pz*z/Z));

	  
	  //==========================================================================//
	  // (I.A)                                                                    //
	  // Calculate a single summand of the sink field phi(k_1) without            //
	  // the prefactor 1/T, i.e. "phi(\vec{x}_1,t1) exp(...)":                    //
	  //                                                                          //
	  //==========================================================================//
	  
	  aux_sink   = prod_complex(exp_ipx, phi[lattice_point(t1,x,y,z)]);
	  phi_t_sink.re +=   aux_sink.re/(X*Y*Z);
	  phi_t_sink.im +=   aux_sink.im/(X*Y*Z);	  

	  
	  //==========================================================================//
	  // (I.B)                                                                    //
	  // Calculate a single summand of the source field phi(k_2) without          //
	  // the prefactor 1/T, i.e. "phi(\vec{x}_2,t1+dt) exp(...)":                 //
	  //                                                                          //
	  //==========================================================================//
	  
	  aux_source  = prod_complex(exp_ipx, phi[lattice_point((t1+dt)%T,x,y,z)]);
	  phi_t_source.re +=  aux_source.re/(X*Y*Z);
	  phi_t_source.im +=  aux_source.im/(X*Y*Z);
	  
	}
      }
    }

    //================================================================================//
    // (I.C)                                                                          //
    // Calculate the n particle sink and source operator:                             //
    //                                                                                //
    //================================================================================//
    
    operator_sink   = pow(phi_t_sink, n_aux);
    operator_source = pow(phi_t_source, n_aux);
    
    
    //================================================================================//
    // (I.D)                                                                          //
    // Calculate the n particle correlation function (real and imaginary part):       //
    //                                                                                //
    //================================================================================//
    
    corr.re += (prod_complex(operator_sink, conjugate(operator_source)).re/T);
    corr.im += (prod_complex(operator_sink, conjugate(operator_source)).im/T);

    }
  
  return corr;
}
