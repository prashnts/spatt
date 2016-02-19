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
#ifndef STAT_H
#define STAT_H 1

#include <cstdio>
#include <vector>
#include <iterator>
#include "sequence.h"
#include "dfa.h"
#include "pmc.h"

#define OVER 0
#define UNDER 1

namespace spatt {

typedef struct {
  unsigned long start;
  unsigned long length;
} seq;

class stat {

 public:

  int _rep;
  dfa *_Dfa;
  pmc *_Y;
  unsigned long sequence_length;
  sequence *_Seq;
  unsigned short _m;
  unsigned long _nobs;
  unsigned long _npresence;
  std::vector<seq> _seq;
  double _pvalue;

  // rep==OVER or rep==UNDER
  stat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool verbose);
  stat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool verbose);

  virtual void compute(bool interrupt=true,bool verbose=false);

  virtual void compute(std::vector<double> &mu0,bool interrupt=true,bool verbose=false);

  virtual void print(const char *label);

  inline double pvalue(){
    return _pvalue;
  };

  inline unsigned long nobs(){
    return _nobs;
  };

};
};
#endif
