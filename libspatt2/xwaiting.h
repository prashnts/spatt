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
#ifndef XWAITING_H
#define XWAITING_H 1

#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include "sequence.h"
#include "dfa.h"
#include "pmc.h"

#define OVER 0
#define UNDER 1

#define TOL 1e-10

namespace spatt {

class xwaiting {

 public:
  
  int _rep;
  dfa *_Dfa;
  pmc *_Y;
  sequence *_Seq;
  unsigned short _m;

  std::vector<std::vector<double> > _u;
  double _lambda;
  unsigned long _t0;

  std::vector<unsigned long> _pos;
  std::vector<unsigned long> _obs;
  std::vector<bool> _final;
  std::vector<double> _pvalue;


  // rep==OVER or rep==UNDER
  xwaiting(dfa &D,pmc &Y,sequence &S,int rep,bool verbose=false,bool debug=false);
  
  // fills _u and compute _lambda and _t0
  void precompute(bool debug=false);

  void compute(bool debug=false);

  // returm sum_{i<=t} P(tau_from=i)
  double cumsum_tau(unsigned long from,unsigned long t);

  // return sum_{i>t} P(tau_from=i) 
  double tail_tau(unsigned long from,unsigned long t);

  void print();



};
};
#endif
