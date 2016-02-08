/*
    This file is part of SPatt.

    SPatt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SPatt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SPatt; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Copyright 2004, 2005 Grégory Nuel, Mark Hoebeke.
*/

/*********************************************************************/   
/*  								     */
/*  Class stat: class for statistics computations                    */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef TRANSITION_H
#define TRANSITION_H

#include "spattparameters.h"
#include "alphabet.h"
#include "sequence.h"
#include "count.h"
#include "markov.h"
#include "input.h"
#include <string>
#include <map>

#define MINUS_INF -1e300

using namespace std;

namespace spatt {

  struct __transition {
    long number;
    double proba;
    struct __transition *to;
  };
  typedef struct __transition trans;
  
  class transition {
    
  public:
    transition(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    ~transition();
    void build(word *w);
    void build(word *w1,word *w2);
    void build(pattern *patt);
    // inline stuff
    inline const long get_order() const { return _order;}
    // large deviations stuff
    void multiply_Q_by(double t); // this function multiply all counting proba by t
    void transition_by_vector(double *x,double *res); // compte res = transition *  x
    void fill_fortran(double *res); // fill a fortran matrix with the transition
    // exact stuff
    void P_times(double *u,double *res); // res = P * u
    void P_times_plus_Q_times(double *u,double *v,double *res); // res = P * u + Q * v
    void sum_Q(double *res); // res = sum of Q by row (sum of row1, sum of row2, ...)
    void log_P_times(double *u,double *res); // res = log (P * exp(u))
    void log_P_times_plus_Q_times(double *u,double *v,double *res); // res = log( P * exp(u) + Q * exp(v) )

    inline double mylog(double x) {
      if (x<=0)
	return MINUS_INF;
      else
	return log(x);
    };

    inline double myexp(double x) {
      if (x<=MINUS_INF)
	return 0.0;
      else
	return exp(x);
    };

  private:
    const spattparameters *_params;
    const alphabet *_alpha;
    const sequence *_seq;
    const count *_occ;
    const markov *_M;
    const input *_In;
    int _h,_k;
    map <string,trans *> _common_states,_states;
    trans *_common_trans_table,*_trans_table;
    long _allocated_size;
    long _n_common_states,_n_states;
    long _order,_trans_size;    
    long *_positions;
    double *_probas;
    long _n_counting_trans;
    long _counting_pos_size;
    long *_counting_pos;
    long _counted;
     
    void init();
    void add(string s);
    void tag(); 
    void connect();
    string longuest_suffix(string s);
    void print_trans();
    void print_fixed_trans();
    void fix_trans();
    void counting(string s);
    void print_counting();   
  };
};
#endif
