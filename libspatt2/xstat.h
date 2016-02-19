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
#ifndef XSTAT_H
#define XSTAT_H 1

#include <cstdio>
#include <vector>
#include <list>
#include <iterator>
#include "stat.h"

namespace spatt {

class xstat : public stat {

 public:

  bool _presence;

  unsigned long _c; // = _nobs if _rep=OVER, _nobs+1 if _rep=UNDER;

  std::vector<std::vector<double> > _x; // computation vector

  // rep==OVER or rep==UNDER
  xstat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool presence,bool verbose);
  xstat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool presence,bool verbose);

  virtual void compute(bool interrupt=true,int offset=0,bool verbose=false);

  virtual void compute(std::vector<double> &mu0,bool interrupt=true,int offset=0,bool verbose=false);

  virtual void print(const char *label,std::string &format_float);

  void print_distribution(std::string &format_float);

};
};

#endif
