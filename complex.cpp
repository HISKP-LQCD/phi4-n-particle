#include <assert.h>

#include <cmath>

#include "complex.h"



//###########################################################################################//
// (I.)                                                                                      //
//                    Function for the calculation of the complex product:                   //
//                                                                                           //
//###########################################################################################//

complex prod_complex(complex z1, complex z2) {
  return static_cast<std::complex<double>>(z1) *
    static_cast<std::complex<double>>(z2);
};



//###########################################################################################//
// (II.)                                                                                     //
//                   Function for the calculation of the complex conjugate:                  //
//                                                                                           //
//###########################################################################################//

complex conjugate(complex z1) {
  complex result;
  result.re = z1.re;
  result.im = -1*z1.im;
  return result;
};



//###########################################################################################//
// (III.)                                                                                    //
//                                     Complex power function:                               //
//                                                                                           //
//###########################################################################################//

complex pow(complex base, complex power) {
  return std::pow(static_cast<std::complex<double>>(base), static_cast<std::complex<double>>(power));
}

complex pow(complex base, long power) {
  return std::pow(static_cast<std::complex<double>>(base), power);
}

std::ostream &operator<<(std::ostream &os, complex const &c) {
  os << static_cast<std::complex<double> const>(c);
  return os;
}

/*
complex pow(complex base, long power) {
  assert(power >= 0);
  
  complex result = base;
  
  for (int number = 1; number < power; ++number) {
    result = prod_complex(base, result);
  }
  return result;
};
*/

//###########################################################################################//
// (IV.)                                                                                     //
//              Function for the substraction of two complex numbers z1 and z2:              //
//                                                                                           //
//###########################################################################################//

complex sub_complex(complex z1, complex z2){
  complex result;
  result.re = z1.re - z2.re;
  result.im = z1.im - z2.im;
  return result;
};
