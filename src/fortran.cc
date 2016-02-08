/* $Id: fortran.cc 553 2005-11-24 13:11:19Z gnuel $ */
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

#include "fortran.h"

namespace spatt {

  zgeev_par::zgeev_par(int n,int compute_left,int compute_right) {

    //printf("sizeof(double)=%i\n",sizeof(double));
    //printf("sizeof(complex16)=%i\n",sizeof(complex16));

    /* set default values */
    _lwork=-1;
    _ldvr=1;
    _ldvl=1;

    /* set values */
    _n=n;
    _lda=_n;
    if (compute_left==0) {
      _jobvl='N';
      _ldvl=1;
    } else {
      _jobvl='V';
      _ldvl=_n;
    }
    if (compute_right==0) {
      _jobvr='N';
      _ldvr=1;
    } else {
      _jobvr='V';
      _ldvr=_n;
    }
    _size=(_n+_ldvl+_ldvr+1)*_n+1;
    _size*=sizeof(complex16);
    /* alloc */
    alloc();
    /* connect */
    connect();
 
  };

  dgeev_par::dgeev_par(int n,int compute_left,int compute_right) {

    /* set default values */
    _lwork=-1;
    _ldvr=1;
    _ldvl=1;

    /* set values */
    _n=n;
    _lda=_n;
    if (compute_left==0) {
      _jobvl='N';
      _ldvl=1;
    } else {
      _jobvl='V';
      _ldvl=_n;
    }
    if (compute_right==0) {
      _jobvr='N';
      _ldvr=1;
    } else {
      _jobvr='V';
      _ldvr=_n;
    }
    _size=(2+_ldvl+_ldvr)*_n+1;
    _size*=sizeof(double);
    /* alloc */
    alloc();
    /* connect */
    connect();
 
  };

  void zgeev_par::alloc() {

    _data=NULL;
    _data=(complex16 *)malloc(_size);
    //printf("zgeev_par::alloc _data=%p\n",_data);
    if (_data==NULL) {
      fprintf(stderr,"not enough memory to allocate zgeev_par\n");
      exit(EXIT_FAILURE);
    }

  };

  void dgeev_par::alloc() {

    _data=NULL;
    _data=(double *)malloc(_size);
    if (_data==NULL) {
      fprintf(stderr,"not enough memory to allocate zgeev_par\n");
      exit(EXIT_FAILURE);
    }

  };

  zgeev_par::~zgeev_par() {
  
    //printf("zgeev_par::~zgeev_par _data=%p\n",_data);
    if (_data!=NULL)
      free(_data);

  };

  dgeev_par::~dgeev_par() {
  
    if (_data!=NULL)
      free(_data);

  };

  void zgeev_par::print() {

    printf("jobvl='%c'\tjobvr='%c'\tn=%i\n",_jobvl,_jobvr,_n);

  };

  void dgeev_par::print() {

    printf("jobvl='%c'\tjobvr='%c'\tn=%i\tlwork=%i\n",_jobvl,_jobvr,_n,_lwork);

  };

  void zgeev_par::connect(){

    complex16 *pos;
    pos=_data;
    //_a=pos;
    //pos+=_n*_n;
    _w=pos;
    pos+=_n*_n;
    _vl=pos;
    pos+=_ldvl*_n;
    _vr=pos;
    pos+=_ldvr*_n;
    _rwork=(double *)pos;
    pos+=_n;
    _work=pos;    

  };

  void dgeev_par::connect(){
    
    double *pos;
    pos=_data;
    _wr=pos;
    pos+=_n;
    _wi=pos;
    pos+=_n;
    _vl=pos;
    pos+=_ldvl*_n;
    _vr=pos;
    pos+=_ldvr*_n;
    _work=pos;    

  };

  void zgeev_par::reset(int n,int compute_left,int compute_right,int lwork) {
  
    /* update values */
    _n=n;
    _lda=_n;
    if (compute_left==0) {
      _jobvl='N';
      _ldvl=1;
    } else {
      _jobvl='V';
      _ldvl=_n;
    }
    if (compute_right==0) {
      _jobvr='N';
      _ldvr=1;
    } else {
      _jobvr='V';
      _ldvr=_n;
    }
    _lwork=lwork;
    /* check needed memory*/
    long new_size;
    new_size=(_n+_ldvl+_ldvr+1)*_n;
    if (_lwork>0) {
      new_size+=_lwork;
    } else {
      new_size++;
    }
    new_size*=sizeof(complex16);
    if (new_size<=_size) {
      /* no alloc needed */
    } else {
      _size=new_size;
      free(_data);
      alloc();
    }
    connect();
  };
  
  void dgeev_par::reset(int n,int compute_left,int compute_right,int lwork) {
  
    /* update values */
    _n=n;
    _lda=_n;
    if (compute_left==0) {
      _jobvl='N';
      _ldvl=1;
    } else {
      _jobvl='V';
      _ldvl=_n;
    }
    if (compute_right==0) {
      _jobvr='N';
      _ldvr=1;
    } else {
      _jobvr='V';
      _ldvr=_n;
    }
    _lwork=lwork;
    /* check needed memory*/
    long new_size;
    new_size=(2+_ldvl+_ldvr)*_n;
    if (_lwork>0) {
      new_size+=_lwork;
    } else {
      new_size++;
    }
    new_size*=sizeof(double);
    if (new_size<=_size) {
      /* no alloc needed */
    } else {
      _size=new_size;
      free(_data);
      alloc();
    }
    connect();
  };

  dneupd_par::dneupd_par(int n,int nev,int mode){

    /* set default values */
    _rvec=1;
    _howmny='A';
  
    /* set values */
    _n=n;
    _nev=nev;
    _ldz=_n;
    _mode=mode;
  
    _data=NULL;
    if (_mode==DNEUPD_GET_VECTORS) {
      _size=(long)sizeof(double)*(_n*(_nev+1));
      //printf("n=%i\tnev=%i\n",n,nev);
      //printf("%i = %i x %i\n",_size,sizeof(double),_n*(_nev+1));
    } else {
      _size=0;
    }

    alloc();
    connect();
  
  };

  void dneupd_par::alloc(){
    _data=NULL;
    if (_size>0) {
      _data=(double *)malloc(_size);
      if (_data==NULL) {
	fprintf(stderr,"Not enough memory for dneupd_par allocation\n");
	exit(EXIT_FAILURE);
      }
    }
  };

  void dneupd_par::connect(){
    _z=_data;
  };

  void dneupd_par::print(){
    printf("size=%li\n",_size);
    printf("data=%p\tz=%p\n",_data,_z);
    printf("n=%i\tnev=%i\trvec=%i\thowmny='%c'\tldz=%i\n",_n,_nev,_rvec,_howmny,_ldz);
  };

  void dneupd_par::reset(int n,int nev,int mode){

    if (_mode!=DNEUPD_GET_VECTORS && mode==DNEUPD_GET_VECTORS) {
      /* need to allocate memory */
      _size=(long)sizeof(double)*(_n*(_nev+1));
      alloc();
      connect();
    }
    if (_mode==DNEUPD_GET_VECTORS && mode!=DNEUPD_GET_VECTORS) {
      /* need to free memory */
      _size=0;
      free(_data);
      _data=NULL;
      connect();
    }
    _mode=mode;
    if (_n!=n || _nev!=nev) {
      /* update values */
      _n=n;
      _nev=nev;
      _ldz=_n;
      /* realloc if necessary */
      int new_size=n*(nev+1);
      if (new_size>_size) {
	free(_data);
	_size=new_size;
	alloc();
	connect();
      }
    }
  
  };

  dneupd_par::~dneupd_par(){
    if (_size>0)
      free(_data);
  }

  dnaupd_par::dnaupd_par(int n,int ncv){

    /* set default values */
    _ido=0;
    _bmat='I';
    _which[0]='L';
    _which[1]='M';
    _which[2]='\0';
    _nev=1;
    _tol=DNAUPD_TOL;
    _iparam[0]=1;
    _iparam[2]=DNAUPD_MAXITER;
    _iparam[3]=1;
    _iparam[6]=1; 
    _info=0;
  
    /* set values */
    _n=n;
    _ncv=ncv;
    _ldv=_n;
    _lworkl=_n+_ncv*_n+3*_n+3*_ncv*_ncv+6*_ncv;

    /* allocate */
    _size=(long)sizeof(double)*(_lworkl);
    //printf("%i * %i = %li\n",sizeof(double),_lworkl,_size);
    alloc();
    /* connect */
    connect();
  
  };

  void dnaupd_par::alloc() {
    _data=NULL;
    _data=(double *)malloc(_size);
    if (_data==NULL) {
      fprintf(stderr,"not enough memory to allocate dnaupd_par\n");
      exit(EXIT_FAILURE);
    }
  };

  void dnaupd_par::connect() {
    double *pos=_data;
    _resid=pos;
    pos+=_n;
    _v=pos;
    pos+=(_ncv*_n);
    _workd=pos;
    pos+=3*_n;
    _workl=pos;
  };

  void dnaupd_par::reset(int n,int ncv){

    //printf("dnaupd_par::reset(%i,%i)\n",n,ncv);
    /* set default values */
    _ido=0;
    _bmat='I';
    _which[0]='L';
    _which[1]='M';
    _nev=1;
    _tol=DNAUPD_TOL;
    _iparam[0]=1;
    _iparam[2]=DNAUPD_MAXITER;
    _iparam[3]=1;
    _iparam[4]=1;
    _iparam[6]=1;  

    /* set values */
    _ido=0;
  
    if (_n==n && _ncv==ncv) {
      /* nothing to do */
    } else {
      //printf("reallocating dnaupd_par\n");
      //long new_size;
      //new_size=n+ncv*n+3*n+3*ncv*ncv+6*ncv;
      _n=n;
      _ncv=ncv;
      _ldv=_n;
      _lworkl=_n+_ncv*_n+3*_n+3*_ncv*_ncv+6*_ncv;
      long new_size=(long)sizeof(double)*(_lworkl);
      if (new_size>_size) {
	/* realloc */
	free(_data);
	_size=new_size;
	alloc();
      }
      connect();
    }
  };

  void dnaupd_par::print(){
    printf("size=%li\n",_size);
    printf("ido=%i\tbmat='%c'\tn=%i\twhich=\"%s\"\tnev=%i\ttol=%.2e\tncv=%i\tldv=%i\n",_ido,_bmat,_n,_which,_nev,_tol,_ncv,_ldv);
    printf("iparam=");
    for (int i=0; i<11; i++)
      printf("[%i]:%i\t",i+1,_iparam[i]);
    printf("\n");
    printf("lworkl=%i\tinfo=%i\n",_lworkl,_info);
  }

  dnaupd_par::~dnaupd_par(){
    if (_data!=NULL)
      free(_data);
  };

};
