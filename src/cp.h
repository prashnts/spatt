/* $Id: cp.h 583 2005-12-09 15:54:10Z gnuel $ */
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

#ifndef CP_H
#define CP_H

#include <cstdio>
#include <cstring>
//#include <gsl/gsl_sf_hyperg.h>

#include "spattglobals.h"
#include "word.h"
#include "alphabet.h"
#include "sequence.h"
#include "count.h"
#include "markov.h"
#include "cdf.h"
#include "timer.h"

using namespace std;

namespace spatt {

  const int LOGFACT_BLOCK_SIZE = 1024;
  const int CP_MP_NDIGIT       =   15;
  const int CP_MP_NOUTPUT      =    4;

  const double CP_RELATIVE_PRECISION = 1e-10;
  const double LOGEPSILON            = -29.93361;
  const double LOGZERO = -690.0;

/* return log(n factorial) */
/* use a static table */
  double logfact(int n);

  /* return the cpstat of a pattern */
  /* mu = expectation */
  /* A = overlap parameter */
  /* sstat = simple statistic of the pattern */
  /* N = number of occurrences */
  /* linear complexity new implementation */
  double newcpstat(double expectation,double A,double sstat,long N);
  
  /* return A for a word w */
  double cpA(const char *w,const markov *M);
  
  double pgeopois(double x,double lambda,double theta,int lowertail,int logscale);
  void inline nextL(double *Lcurrent, double *Lprec, double *Lprecprec, long n,double z, double theta);
  void inline updateSum(double *A, double *S,double logadd);

  // compute P(N=n) for N compound poisson with geometric tail
  double complete_dcpoi(int n,double lambda,double theta, int alpha,double *p);
  double barbour_dcpoi(int n,double lambda,double theta, int alpha,double *p);
  double fast_dcpoi(int n,double lambda,double theta, int alpha,double *p);

  // compute P(N<=n) for N compound poisson with geometric tail
  double complete_pcpoi(int n,double lambda,double theta, int alpha,double *p);
  double barbour_pcpoi(int n,double lambda,double theta, int alpha,double *p);
  double fast_pcpoi(int n,double lambda,double theta, int alpha,double *p);

  // compute P(N>=n) for N compound poisson with geometric tail
  double complete_qcpoi(int n,double lambda,double theta, int alpha,double *p);
  double barbour_qcpoi(int n,double lambda,double theta, int alpha,double *p);
  double fast_qcpoi(int n,double lambda,double theta, int alpha,double *p);


  // compute log P(N=n) for N compound poisson with geometric tail
  double complete_log_dcpoi(int n,double lambda,double theta, int alpha,double *p);
  double fast_log_dcpoi(int n,double lambda,double theta, int alpha,double *p);
  
}
#endif
