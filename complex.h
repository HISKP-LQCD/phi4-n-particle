#pragma once

#include <complex>
#include <iosfwd>

// Here, the data type "complex" is created:
class complex {
  public:
  double re;
  double im;

  operator std::complex<double>() const {
    return std::complex<double>(re, im);
  }

  complex() {}
  complex(double r, double i) : re(r), im(i) {}

  complex(std::complex<double> c): re(c.real()), im(c.imag()) {}

};

std::ostream &operator<<(std::ostream &os, complex const &c);


complex prod_complex(complex z1, complex z2);
complex conjugate(complex z1);

complex pow(complex base, complex power);
complex pow(complex base, long power);

complex sub_complex(complex z1, complex z2);