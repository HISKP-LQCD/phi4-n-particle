#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "parameters.h"
#include "scalar.h"
#include "types.h"



//###########################################################################################//
// (I.)                                                                                      //
//                                                                                           //
//###########################################################################################//

double eval_action_nogauge(scalar_field phi) {   // scalar_field defined in "types.h" as
                                                 // complex;
                                                 // phi defined in "calculate_toytest.cpp".
  double action = 0.;
  int x,y,z,t,gamma;
  int ipt;
  float add=0;
  
  for(t=0;t<T;t++) {
    for(x=0;x<X;x++) {
      for(y=0;y<Y;y++) {
	for(z=0;z<Z;z++) {

	  // LAMBDA*(phi^2 - 1)^2 + phi^2
	  add=0;
          add = prod_complex( conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,y,z)]).re;
	  action += LAMBDA*( add -1 )*( add -1 ) + add; 

	  // -KAPPA  phi_x^* U_x,mu phi_x+mu +cc)
	  add=0;
	  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point((t+1)%T,x,y,z)]).re;
	  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,(x+1)%X,y,z)]).re;
	  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,(y+1)%Y,z)]).re;
	  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,y,(z+1)%Z)]).re;
	  
	  action += -KAPPA*2*add;
	  
	}
      }
    }
  }
  return action;
};



//###########################################################################################//
// (II.)                                                                                     //
//     Function for the  calculation of the change in action \Delta S, which is needed in    //
//                                       "metropolis.cpp":                                   //
//                                                                                           //
//###########################################################################################//

double delta_action_nogauge(scalar_field phi, scalar_field phi_new, int x, int y, int z, int t) {
  
  double action = 0.;
  double action_new = 0.;
  int ipt;
  float add=0;

  //=========================================================================================//
  //                                                                                         //
  // Calculation of the Euclidean action S_E:                                                //
  //                                                                                         //
  //=========================================================================================//
  
  // LAMBDA*(phi^2 - 1)^2 + phi^2
  add=0;
  add = prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,y,z)]).re;
  action += LAMBDA*( add -1 )*( add -1 ) + add; 
  
  // -KAPPA  phi_x phi_x 
  add=0;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point((t+1)%T,x,y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,(x+1)%X,y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,(y+1)%Y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,y,(z+1)%Z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point((T+t-1)%T,x,y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,(X+x-1)%X,y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,(Y+y-1)%Y,z)]).re;
  add += prod_complex(conjugate(phi[lattice_point(t,x,y,z)]),phi[lattice_point(t,x,y,(Z+z-1)%Z)]).re;
  
  action += -1*KAPPA*add*2;

  //=========================================================================================//
  //                                                                                         //
  // Calculation of the new action S_E,new:                                                  //
  //                                                                                         //
  //=========================================================================================//
  
  // LAMBDA*(phi^2 - 1)^2 + phi^2    
  add=0;
  add = prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,x,y,z)]).re;
  action_new += LAMBDA*( add -1 )*( add -1 ) + add; 
  
  // -KAPPA  phi_x phi_x 
  add=0;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point((t+1)%T,x,y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,(x+1)%X,y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,x,(y+1)%Y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,x,y,(z+1)%Z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point((T+t-1)%T,x,y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,(X+x-1)%X,y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,x,(Y+y-1)%Y,z)]).re;
  add += prod_complex(conjugate(phi_new[lattice_point(t,x,y,z)]),phi_new[lattice_point(t,x,y,(Z+z-1)%Z)]).re;
  
  action_new += -1*KAPPA*add*2;

  return (action_new - action);
};
