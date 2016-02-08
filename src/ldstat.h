/* $Id: ldstat.h 968 2006-11-20 09:57:54Z gnuel $ */
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

#ifndef LDSTAT_H
#define LDSTAT_H

#include "ldspattparameters.h"
#include "sstat.h"
#include "transition.h"
#include "fortran.h"
//#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

#define GGOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define CHEBYSHEV_ORDER 20
#define CHEBYSHEV_HALF_SPAN 1e-4
#define LDSTAT_PI 3.14159
#define DELTA 1

static double mmaxarg1,mmaxarg2;
#define DDMAX(a,b) (mmaxarg1=(a),mmaxarg2=(b),(mmaxarg1) > (mmaxarg2) ? (mmaxarg1) : (mmaxarg2))

extern "C" {
  double ldstat_sparse_F(double,void*);
  double ldstat_full_F(double,void*);
  double ldstat_full_rate(double,void*);
  double ldstat_sparse_rate(double,void*);
}

using namespace std;

namespace spatt {
  
  class ldstat: public sstat {
    
  public:
    ldstat(spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    virtual ~ldstat();
    /* compute the stat (here nobs) for simple word, word and ic or pattern */
    virtual void compute(word *w);
    virtual void compute(word *w1,word *w2);
    virtual void compute(pattern *patt);
    ///* if called with no previous call to compute, displays the output format */
    //virtual void print_format(FILE *fout);
    //virtual void print_regular(FILE *fout);
    //virtual void print_normalized(FILE *fout);
    double sparse_F(double t,void *not_used); // L(t)-a*t in the sparse case
    double full_F(double t,void *not_used); // L(t)-a*t in the full case
    double full_rate(double a,void *not_used); // I(a) in the full case
    double sparse_rate(double a,void *not_used); // I(a) in the sparse case
    
  private:
    transition *_T;
    long _order;
    bool _precise;
    gsl_cheb_series *_rate,*_drate,*_ddrate;
    double _V,_eta;
    int _h,_k;
    long _n;
    dnaupd_par *_sparse_par;
    double _a,_ori_a;
    double _In;
    dgeev_par *_full_par;
    double _AA[MAX_FULL_DIM*MAX_FULL_DIM];
    

    void comp();
    void mnbrack(double *ax,double *bx,double *cx, double *fa,double *fb, double *fc,
		 double (*func)(double,void*));
    double rate(double (*func)(double,void*));



  protected:
    /* added here new members */
    
  };
};
#endif
