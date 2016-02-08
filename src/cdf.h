/* $Id: cdf.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Provides support for Poisson and Gaussian cdf with special	     */
/*  care for very small p-value (log of p-value are given here)	     */
/*  								     */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef CDF_H
#define CDF_H

#define EPS 3.0e-7
#define ITMAX 10000
#define FPMIN 1.0e-30
#define MAXIT 100
#define MAGNITUDE_LIMIT 37.5
#define PI 3.14159

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_cdf.h>

  /* returns incomplete gamma function P(a,x) */
  double gammp(double a, double x);

  /* returns incomplete gamma function Q(a,x) = 1 - P(a,x) */
  double gammq(double a, double x);

  /* returns incomplete gamma function P(a,x) through series */
  /* output: gamser contain the evaluation and gln ln(Gamma(a)) */
  void gser(double *gamser, double a, double x, double *gln);

  /* returns incomplete gamma function Q(a,x) through continued fraction */
  /* output: gamcf contain the evaluation and gln ln(Gamma(a)) */
  void gcf(double *gammcf, double a, double x, double *gln);

  /* returns the value of ln(Gamma(xx)) for xx>0 */
  double gammln(double xx);

  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X<=k) with X Poisson of mean x */
  /* returns also log of that value */
  double pcdfpoi(double k,double x,double *m,double *e);

  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X>=k) with X Poisson of mean x */
  /* returns also log of that value */
  double qcdfpoi(double k,double x,double *m,double *e);

  /* return the value n! using a static array */
  double factrl(int n);

  /* returns the incomplete beta function Ix(a,b) */
  double betai(double a, double b, double x);

  /* evaluates continued fraction for incomplete beta */
  /* function by modified Lentz's method */
  double betacf(double a, double b, double x);

  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X<=k) with X Binomial (n,p) */
  /* returns also log of that value */
  /* uses Poisson approximation when necessary */
  double pcdfbin(double k,long n,double p,double *m,double *e);

  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X>=k) with X Binomial (n,p) */
  /* returns also log of that value */
  /* uses Poisson approximation when necessary */
  double qcdfbin(double k,long n,double p,double *m,double *e);

  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X<=x) with X Gaussian N(0,1) */
  /* returns also log of that value */
  double pcdfnor(double x,double *m,double *e);
  
  /* returns mantiss (in m) and exposant (in e) */
  /* of the P(X>=x) with X Gaussian N(0,1) */
  /* returns also log of that value */
  double qcdfnor(double x,double *m,double *e);

#ifdef __cplusplus
}
#endif

#endif
