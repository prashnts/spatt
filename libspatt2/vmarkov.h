/**
    Copyright 2006 Mark Hoebeke, Vincent Miele & Gregory Nuel.

    This file is part of SPatt 2

    SPatt 2 is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SPatt 2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SPatt 2; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**/
#ifndef VMARKOV_H
#define VMARKOV_H 1

#include <cstdio>
#include <vector>
#include <iterator>
#include <math.h>
#include "markov.h"

#define MU_TOL 1e-10

namespace spatt {

class vmarkov {

 public:
  
  unsigned short _alphabet_size;
  unsigned short _m; // markov order
  std::vector<double> _pi;
  std::vector<double> _mu;
  
  unsigned long _alpha;
  unsigned long _K;
  unsigned long _Len;

  double *_M;
  double **_Sigma;

  vmarkov(markov &M);

  ~vmarkov();

  // compute y = x * Pi
  void xPi(std::vector<double> &x,std::vector<double> &y);
  void xPi(double *x,double *y);
  // compute y = Pi * x
  void Pix(std::vector<double> &x,std::vector<double> &y);
  void Pix(double *x,double *y);

  // compute _mu and _alpha
  void compute_mu();

  void print_mu();

  void compute_Sigma(unsigned long ell,unsigned long start=0);

  void print_Sigma();

  void dump_Sigma(const char *file);
  
};
};
#endif
