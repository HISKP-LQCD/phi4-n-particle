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
#include "scalar.h"
#include "types.h"
#include "generator_singleton.h"



// If the function "metropolis()" (see (II.)) is applied, then the metropolis algorithm is
// calculated combined for odd and even t (depending on the argument "core" in the called
// function "metropolis_core()").

//###########################################################################################//
// (I.)                                                                                      //
//                        METROPOLIS ALGORITHM FOR ALL ODD OR EVEN t:                        //
//                                                                                           //
//   At the beginning of each MC step, the coordinates t,x,y,z of the initial lattice point  //
//    lattice_point(t,x,y,z) are chosen randomly (if core==0: only odd t, if core==1: only   //
//  even t). By means of these space-time coordinates, a small change of the field point is  //
//    calculated via update_field_point(), see "scalar.cpp". The resulting field point is    //
//    called phi2. From the initial phi and updated phi2, the change in action is computed   //
//  ("action.cpp") and if exp(-deltaS)>random with 0<=random<=1, the new field point phi2 is //
//   adopted via copy_field_point() (see "scalar.cpp"). If this condition is not satisfied,  //
//                              the initial field point is kept.                             //
//                                                                                           //
//###########################################################################################//

int metropolis_core(scalar_field *p_phi, scalar_field *p_phi2, int core) {
  
  scalar_field phi  = *p_phi;  // Initial scalar field point phi
  scalar_field phi2 = *p_phi2; // Updated field point phi2 after small change
  
  int update = 0;
  int x,y,z,t;
  double deltaS, random;

  // Parallel Computing via Open MP:
#pragma omp parallel private(x, y, z, t, deltaS, random)
  {
    int i=0;
    int nthreads = omp_get_num_threads();  
    int pid = omp_get_thread_num();
    
    
    auto x_distribution       = std::uniform_int_distribution<int>(0,X-1);
    auto y_distribution       = std::uniform_int_distribution<int>(0,Y-1);
    auto z_distribution       = std::uniform_int_distribution<int>(0,Z-1);
    auto t_distribution       = std::uniform_int_distribution<int>(0,(T/nthreads)-1);
    auto ZeroOne_distribution = std::uniform_real_distribution<double>(0.0,1.0);
    
    
    //=======================================================================================//
    // (I.A)                                                                                 //
    // This is the for loop over all n_metropolis MC steps (n_metropolis is the number of    //
    // local updates in a single Metropolis step):                                           //
    //                                                                                       //
    //=======================================================================================//
    
    for(i=0;i<n_metropolis;i++) {
      
      
      //=====================================================================================//
      // (I.B)                                                                               //
      // Spatial coordinates of the ith lattice point are created randomly for the           //
      // ith update:                                                                         //
      //                                                                                     //
      //=====================================================================================//

      // get() is the thread number dependend mersenne twister defined in the class
      // "GeneratorSingleton" in "generator_singleton.h":
      x = x_distribution(GeneratorSingleton::get());
      y = y_distribution(GeneratorSingleton::get());
      z = z_distribution(GeneratorSingleton::get());
      
      
      //=====================================================================================//
      // (I.C)                                                                               //
      // Odd temporal coordinate of the lattice point is created randomly for the ith update://
      //                                                                                     //
      //=====================================================================================//

      if(core == 0) {
	do{
	  t = t_distribution(GeneratorSingleton::get()) + pid*T/nthreads;
	}
	while(t%2==0);
      }

      //=====================================================================================//
      // (I.D)                                                                               //
      //Even temporal coordinate of the lattice point is created randomly for the ith update://
      //                                                                                     //
      //=====================================================================================//
      
      if(core == 1) {
	do{
	  t = t_distribution(GeneratorSingleton::get()) + pid*T/nthreads;
	}
	while(t%2==1);
      }
      
      //=====================================================================================//
      // (I.E)                                                                               //
      // The function update_field_point() (see scalar.cpp) effects a small change of the    //
      // initial field point depending on the coordinates chosen randomly above. The new     //
      // scalar field point is p_aux = phi2:                                                 //
      //=====================================================================================//

      update_field_point(&phi2,t,x,y,z);
      
      //=====================================================================================//
      // (I.F)                                                                               //
      // From the initial field point phi_old=phi and the updated field point phi_new=phi2,  //
      // the change in action \Delta S_E is calculated via delta_action_nogauge()            //
      // (see action.cpp):                                                                   //
      //                                                                                     //
      //=====================================================================================//
      
      deltaS = delta_action_nogauge(phi,phi2,x,y,z,t);

      //=====================================================================================//
      // (I.G)                                                                               //
      //   Generate a probability randomly with 0 <= P <= 1 (real numbers uniformly distr.   //
      //      between 0 and 1) and accept the updated field point phi2 with probability      //
      //            exp(-\Delta S_E) > P (which includes the case \Delta S_E <=0):           //
      //                                                                                     //
      //=====================================================================================//

      random = ZeroOne_distribution(GeneratorSingleton::get());
      
      if(exp(-deltaS) > random) {
	
	copy_field_point(&phi,&phi2,t,x,y,z);  // The function copy_field_point() (see
	                                  // scalar.cpp) replaces the field point "phi" with
	                                  // "phi2" which means that the updated field point
	                                  // is accepted.
	
	if(pid==0 && i==0) {               // Only if we have a new scalar field point
	  update = 1;                      // (pid==0 && i==0 means new master thread),
	}                                  // "update" is set to 1.
      }

      //=====================================================================================//
      // (I.H)                                                                               //
      // If the condition exp(-deltaS) > random is not satisfied, then the old field point   //
      // "phi" is kept ("reject"):                                                           //
      //                                                                                     //
      //=====================================================================================//
      
      copy_field_point(&phi2,&phi,t,x,y,z);
    }
  }
  return update;
}



//###########################################################################################//
// (II.)                                                                                     //
//                               COMPLETE METROPOLIS ALGORITHM:                              //
//                                                                                           //
//    This function unifies the even and the odd metropolis core provided by the function    //
//                             "metropolis_core()" given in (I.).                            //
//                                                                                           //
//###########################################################################################//

double metropolis(scalar_field *p_phi, int n_field, FILE *faction) {

  int i=0;
  int core_odd, core_even;
  int n_acc = 0;
  double acc=0;
  scalar_field phi = *p_phi;
  scalar_field phi2;
  double action;
  

  // In "parameters.h": Use any non-zero integer as a seed
  
 
  initialize_field(&phi2);
  copy_field(&phi2,&phi);   

  //=========================================================================================//
  //                                                                                         //
  // Combine the metropolis algorithm for odd t and even t:                                  //
  //                                                                                         //
  //=========================================================================================//
  
  while(n_acc<n_field) {

    core_odd  = metropolis_core(p_phi,&phi2,0);
    core_even = metropolis_core(p_phi,&phi2,1);

    i += 2;
    n_acc+= (core_odd + core_even);

  }
  
  phi = *p_phi;
  action = eval_action_nogauge(phi);
  printf("=====================================================\n");
  printf("Accepted step = %d, action = %f \n",n_acc, action);
  fprintf(faction, "%e\n",action);
  acc = (double) n_acc/i;
    
  printf("acceptance = %f \n", (double) n_acc/i);
  
  return acc;
}
