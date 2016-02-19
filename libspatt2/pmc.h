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
#ifndef PMC_H
#define PMC_H 1

#include <map>
#include <vector>
#include <iterator>
#include "dfa.h"

namespace spatt {

typedef struct {
  unsigned long from;
  unsigned long to;
  unsigned long index;
  double proba;
} term;

class pmc {

 public:

  /* statistics */
  unsigned long _nstates;
  unsigned short _alphabet_size;
  unsigned short _m; // order >=0
  
  std::vector<unsigned long> _recode;

  /* transitions */
  std::vector<term> _Pt; // regular transitions
  std::vector<term> _Qt; // counting transitions
  std::vector<double> _SQ; // sum of _Qt (by columns)
  
  /* start */
  std::vector<unsigned long> _start; // starting positions
  std::vector<unsigned long> _final; // final states


  /* functions */
  
  /* constructor from a dfa and a Markov parameter */
  pmc(dfa &A,std::vector<double> &markov_param,bool verbose,bool debug);

  /* destructor */
  ~pmc();

  void print();

  void sci_export(const char *file);

  void indexed_export(const char *file);

  // set a vector to 0
  inline void zero(std::vector<double> &x){
  for (std::vector<double>::iterator it=x.begin(); it!=x.end(); it++)
     *it=0.0;
  };

  // add sparse matrix s vector x to res
  inline void add_Mx(std::vector<term> &M,std::vector<double> &x,std::vector<double> &res){
    for (std::vector<term>::iterator it=M.begin(); it!=M.end(); it++)
      res[it->from]+=x[it->to]*it->proba;
  };

  // add vector x sparse matrix s to res
  inline void add_xM(std::vector<double> &x,std::vector<term> &M,std::vector<double> &res){
    for (std::vector<term>::iterator it=M.begin(); it!=M.end(); it++)
      res[it->to]+=x[it->from]*it->proba;
  };

  // add scalar vector x sparse matrix s to res
  inline void add_sxM(double s,std::vector<double> &x,std::vector<term> &M,std::vector<double> &res){
    for (std::vector<term>::iterator it=M.begin(); it!=M.end(); it++)
      res[it->to]+=s*x[it->from]*it->proba;
  };

  inline double scalprod(std::vector<double> &x,std::vector<double> &y){
    double res=0.0;
    for (unsigned long i=0; i<_nstates; i++)
      res+=x[i]*y[i];
    return res;
  };


  // example: compute y=x*(P+exp(t)*Q)
  // zero(y); add_xM(x,_Pt,y); add_sxM(exp(t),x,_Qt,y);

  // return the pmc state corresponding to p 
  // -1 if p have been removed in the pmc
  inline unsigned long recode(unsigned long p) {
    return _recode[p]-1;
  };


};

};

#endif 
