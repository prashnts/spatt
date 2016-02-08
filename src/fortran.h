/* $Id: fortran.h 553 2005-11-24 13:11:19Z gnuel $ */
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
/*  Class fortran parameter: help to call fortran routines           */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef FORTRAN_H
#define FORTRAN_H

#include <cstdio>
#include <cstdlib>

using namespace std;

/* complex16 type */
typedef struct {
  double re;
  double im;
} complex16;


extern "C" {
  void dgeev_(char&,char&,int&,double*,int&,double*,double*,double*,int&,double*,int&,double*,int&,int&);
  void zgeev_(char&,char&,int&,complex16*,int&,complex16*,complex16*,int&,complex16*,int&,complex16*,int&,double*,int&);
  void dnaupd_(int&,char&,int&,char*,int&,double&,double*,int&,double*,int&,int*,int*,double*,double*,int&,int&);
  void dneupd_(int&,char&,int*,double*,double*,double*,int&,double&,double&,double*,char&,int&,char*,int&,double&,double*,int&,double*,int&,int*,int*,double*,double*,int&,int&);
}

namespace spatt {

  const int   MAX_FULL_DIM       =   20;
  const float DNAUPD_TOL         = 1e-6;
  const int DNAUPD_MAXITER       =  200;
  const int DNAUPD_KRYLOVDIM     =   20; /* must be smaller than MAX_FULL_DIM*/
  const int DNEUPD_GET_VECTORS   =    1;
  const int DNEUPD_MAX_NEV       =    2;     


class dgeev_par {

 private:
  double *_data;
  long _size;

 public:

  char _jobvl;
  char _jobvr;
  int _n;
  int _lda;
  double *_wr; 
  double *_wi; 
  double *_vl; 
  int _ldvl;
  double *_vr;
  int _ldvr;
  int _lwork;
  int _info;
  double *_work;

  /* allocate and connect memory */
  dgeev_par(int n,int compute_left,int compute_right);
  /* do nothing, connect, or allocate more */
  void reset(int n,int compute_left,int compute_right,int workl);
  void alloc();
  void connect();
  /* print param */
  void print();
  /* free all memory */
  ~dgeev_par();

};

  /* remark: zgeev should be replaced by dgeev here and in markov.cc */
class zgeev_par {

 private:
  complex16 *_data;
  long _size;

 public:

  char _jobvl;
  char _jobvr;
  int _n;
  //complex16 *_a;
  int _lda;
  complex16 *_w; 
  complex16 *_vl; 
  int _ldvl;
  complex16 *_vr;
  int _ldvr;
  int _lwork;
  double *_rwork; 
  int _info;
  complex16 *_work;

  /* allocate and connect memory */
  zgeev_par(int n,int compute_left,int compute_right);
  /* do nothing, connect, or allocate more */
  void reset(int n,int compute_left,int compute_right,int workl);
  void alloc();
  void connect();
  /* print param */
  void print();
  /* free all memory */
  ~zgeev_par();

};

class dnaupd_par {

 private:

  double *_data;
  long _size;

 public:

  int _ido;
  char _bmat;
  int _n;
  char _which[3];
  int _nev;
  double _tol;
  double *_resid;
  int _ncv;
  double *_v;
  int _ldv;
  int _iparam[11];
  int _ipntr[14];
  double *_workd;
  double *_workl;
  int _lworkl;
  int _info;

  /* allocate and connect memory */
  dnaupd_par(int n,int ncv);
  /* do nothing, connect, or allocate more */
  void reset(int n,int ncv);
  void alloc();
  void connect();
  /* print param */
  void print();
  /* free all memory */
  ~dnaupd_par();
  
};

class dneupd_par {

 private:
  int _mode;
  double *_data;
  long _size;

 public:
  
  int _rvec;
  char _howmny;
  int _select[DNAUPD_KRYLOVDIM];
  double _dr[DNEUPD_MAX_NEV+1];
  double _di[DNEUPD_MAX_NEV+1];
  double *_z;
  int _ldz;
  double _sigmar;
  double _sigmai;
  double _workev[3*DNAUPD_KRYLOVDIM];
  int _n;
  int _nev;

  /* allocate and connect memory */
  dneupd_par(int n,int nev,int mode);
 /* do nothing, connect, or allocate more */
  void reset(int n,int nev,int mode);
  void alloc();
  void connect();
 /* print */
  void print();
  ~dneupd_par();
  
};

};
#endif
