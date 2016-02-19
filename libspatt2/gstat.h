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
#ifndef GSTAT_H
#define GSTAT_H 1

#include <cstdio>
#include <vector>
#include <iterator>
#include <math.h>
#include <limits>
#include "sequence.h"
#include "stat.h"

#define MU_TOL 1e-10
#define ITER_MAX 1000

namespace spatt {

class gstat : public stat {

 public:

  std::vector<double> _mu;
  unsigned long _alpha;

  double _esp;
  double _var;
  double _zscore;

  gstat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool presence,bool verbose);
  gstat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool presence,bool verbose);
  
  virtual void compute(bool interrupt=true,int offset=0,bool verbose=false);

  virtual void compute(std::vector<double> &mu0,bool interrupt=true,int offset=0,bool verbose=false);

  virtual void print(const char *label,std::string &format_float);

  void compute_mu(bool verbose=false);

  inline double mean() {
    return _esp;
  };

  inline double var() {
    return _var;
  };

  inline double zscore() {
    return _zscore;
  };

  inline double sd() {
    return sqrt(_var);
  };

    
};

};
#endif
