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
#include "vmarkov.h"

using namespace std;

namespace spatt {

vmarkov::vmarkov(markov &model) {
  
  _alphabet_size=model._alphabet_size;
  _m=model._m;
  _pi=model._param;
  _K=(unsigned long)pow((double)_alphabet_size,(double)_m);
  _Len=_K+(unsigned long)pow((double)_alphabet_size,(double)_m+1);

  _M=NULL;
  _Sigma=NULL;
};

vmarkov::~vmarkov() {

  if (_M) 
    delete[] _M;
  if (_Sigma) {
    delete[] _Sigma[0];
    delete[] _Sigma;
  }
};

void vmarkov::xPi(std::vector<double> &x,std::vector<double> &y){
 unsigned long aux;
  // y=0  
  for (unsigned long i=0; i<_Len; i++)
    y[i]=0.0;
  for (unsigned long i=0; i<_K; i++)
    for (unsigned short b=0; b<_alphabet_size; b++) {
      aux=_alphabet_size*i+b;
      for (unsigned short p=0; p<=_alphabet_size; p++)
	y[_K+aux]+=x[i+p*_K]*_pi[aux];
    }
};

void vmarkov::xPi(double *x,double *y){
 unsigned long aux;
  // y=0  
  for (unsigned long i=0; i<_Len; i++)
    y[i]=0.0;
  for (unsigned long i=0; i<_K; i++)
    for (unsigned short b=0; b<_alphabet_size; b++) {
      aux=_alphabet_size*i+b;
      for (unsigned short p=0; p<=_alphabet_size; p++)
	y[_K+aux]+=x[i+p*_K]*_pi[aux];
    }
};

void vmarkov::Pix(std::vector<double> &x,std::vector<double> &y){
  
  unsigned long aux;
  // y=0  
  for (unsigned long i=0; i<_Len; i++)
    y[i]=0.0;
  for (unsigned long i=0; i<_K; i++)
    for (unsigned short b=0; b<_alphabet_size; b++) {
      aux=_alphabet_size*i+b;
      for (unsigned short p=0; p<=_alphabet_size; p++)
	y[i+p*_K]+=x[_K+aux]*_pi[aux];
    }
};

void vmarkov::Pix(double *x,double *y){
  
  unsigned long aux;
  // y=0  
  for (unsigned long i=0; i<_Len; i++)
    y[i]=0.0;
  for (unsigned long i=0; i<_K; i++)
    for (unsigned short b=0; b<_alphabet_size; b++) {
      aux=_alphabet_size*i+b;
      for (unsigned short p=0; p<=_alphabet_size; p++)
	y[i+p*_K]+=x[_K+aux]*_pi[aux];
    }
};

void vmarkov::compute_mu(){
  { // initialize _mu
    _mu.resize(_Len);
    double aux=1.0/(double)_Len;
    for (vector<double>::iterator it=_mu.begin(); it!=_mu.end(); it++)
      *it=aux;
  } // end of initialization
  _alpha=0;
  double test=1.0;
  vector<double> aux(_Len);
  while (test>MU_TOL) {
    //printf("alpha=%i mu[%i]=%e\n",_alpha,_Len-1,_mu[_Len-1]);
    // old _mu is saved in aux
    aux.swap(_mu);
    // _mu is updated
    xPi(aux,_mu);
    _alpha++;
    test=fabs(_mu[0]-aux[0]);
    for (unsigned long i=1; i<_Len; i++) {
      double tmp=fabs(_mu[i]-aux[i]);
      if (tmp>test)
	test=tmp;
    }      
  }
};

void vmarkov::print_mu() {
  printf("alpha = %i\n",_alpha);
  printf("mu = [ ");
  for (unsigned long i=0; i<_Len; i++)
    printf("%e ",_mu[i]);
  printf("]\n");
};

void vmarkov::compute_Sigma(unsigned long ell,unsigned long start) {

  if (ell<(_m+_alpha)) {
    fprintf(stderr,"vmarkov::compute_Sigma: ell too small. Aborting.\n");
    exit(EXIT_FAILURE);
  } 
  
  // allocate memory
  _M = new double[_Len];
  _Sigma = new double*[_Len];
  _Sigma[0] = new double[_Len*_Len];
  {
    double *pos=_Sigma[0];
    for (unsigned long i=0; i<_Len; i++) {
      _Sigma[i]=pos;
      pos+=_Len;
    }
  }
  double *x = new double[_Len];
  double *y = new double[_Len];
  double **A = new double*[_alpha];
  A[0] = new double[_alpha*_Len];
  {
    double *pos=A[0];
    for (unsigned long i=0; i<_alpha; i++) {
      A[i]=pos;
      pos+=_Len;
    }
  }
  double **B = new double*[_alpha];
  B[0] = new double[_alpha*_Len];
  {
    double *pos=B[0];
    for (unsigned long i=0; i<_alpha; i++) {
      B[i]=pos;
      pos+=_Len;
    }
  }
  // end of memory allocation

  // initialization
  for (unsigned long i=0; i<_Len; i++) {
    _M[i]=0.0;
    for (unsigned long a=0; a<_alpha; a++)
      A[a][i]=0.0;
    for (unsigned long j=0; j<_Len; j++)
      _Sigma[i][j]=0.0;      
  }
  A[0][start]=1.0;
  _M[start]=1.0;
  for (unsigned long d=1; d<_alpha; d++) {
    xPi(A[d-1],A[d]);
    for (unsigned long i=0; i<_Len; i++) 
      _M[i]+=A[d][i];
  } // end d loop
  //printf("A = [ ");
  //for (unsigned long d=0; d<_alpha; d++)
  //  printf("%e ",A[d][_Len-1]);
  //printf("]\n");
  for (unsigned long i=0; i<_Len; i++) 
    _M[i]+=(ell-_m-_alpha+1)*_mu[i];
  for (unsigned long i=0; i<_K; i++)
    for (unsigned long p=1; p<=_alphabet_size; p++)
      _M[i]+=_M[i+p*_K];
  // _M is now computed
  for (unsigned long j=0; j<_K; j++) {
    for (unsigned long q=0; q<=_alphabet_size; q++) {
      unsigned long j0=j+q*_K;
      // x=y=0 B[0]=e_j0
      for (unsigned long i=0; i<_Len; i++) {
	x[i]=0.0;
	y[i]=0.0;
	B[0][i]=0.0;
      }
      B[0][j0]=1.0;
      unsigned long dmax=ell-_m-_alpha;
      if (_alpha-1<dmax)
	dmax=_alpha-1;
      for (unsigned long d=1; d<=dmax; d++) {
	Pix(B[d-1],B[d]);
	for (unsigned long i=0; i<_Len; i++) {
	  x[i]+=(ell-_m-_alpha+1.0-d)*B[d][i];
	  y[i]+=B[d][i];
	}
      } // end d loop
      if (ell>=(_m+2*_alpha)) {
	for (unsigned long i=0; i<_Len; i++) {
	  x[i]+=(ell-_m-2.0*_alpha+1.0)*(ell-_m-2.0*_alpha+2.0)*_mu[j0]/2.0;
	  y[i]+=(ell-_m-2.0*_alpha+1.0)*_mu[j0];
	}
      }
      for (unsigned long i=0; i<_Len; i++)
	x[i]*=_mu[i];
      for (long i=_m+_alpha-1; i>=_m; i--) {
	if ((ell-i)<_alpha) {
	  for (unsigned long ii=0; ii<_Len; ii++)
	    y[ii]+=B[ell-i][ii];
	} else {
	  for (unsigned long ii=0; ii<_Len; ii++)
	    y[ii]+=_mu[j0];
	}
	for (unsigned long ii=0; ii<_Len; ii++)
	  x[ii]+=y[ii]*A[i-_m][ii];
      } // end i loop
      for (unsigned long i=0; i<_K; i++) {
	for (unsigned long p=0; p<=_alphabet_size; p++) {
	  unsigned long i0=i+p*_K;
	  _Sigma[i][j]+=x[i0];
	  _Sigma[j][i]+=x[i0];
	  if (p>0) {
	    _Sigma[i0][j]+=x[i0];
	    _Sigma[j][i0]+=x[i0];
	  }
	  if (q>0) {
	    _Sigma[i][j0]+=x[i0];
	    _Sigma[j0][i]+=x[i0];
	  }
	  if (p*q>0) {
	    _Sigma[i0][j0]+=x[i0];
	    _Sigma[j0][i0]+=x[i0];
	  }
	  _Sigma[i0][j0]-=_M[i0]*_M[j0];
	  if (i==j) {
	    if (q==0) {
	      _Sigma[i0][j0]+=_M[i0];
	    } else if ((p==0)||(p==q)) {
	      _Sigma[i0][j0]+=_M[j0];
	    }
	  }
	} // end p loop
      } // end i loop
    } // end q loop
  } // end j loop

  // free local memory
  delete[] x;
  delete[] y;
  delete[] A[0];
  delete[] A;
  delete[] B[0];
  delete[] B;
      


};

void vmarkov::print_Sigma() {
  printf("M = [ ");
  for (unsigned long i=0; i<_Len; i++)
    printf("%.2f ",_M[i]);
  printf("]\n");
  for (unsigned long i=0; i<_Len; i++) {
    printf("Sigma[%i] = [ ",i);
    for (unsigned long j=0; j<_Len; j++)
      printf("%.2f ",_Sigma[i][j]);
    printf("]\n");    
  }
};

void vmarkov::dump_Sigma(const char *file){
  FILE *out=NULL;
  out=fopen(file,"w");
  if (out==NULL) {
    fprintf(stderr,"vmarkov::dump_Sigma(): cannot write on file \"%s\". Aborting.\n",file);
  }
  for (unsigned long i=0; i<_Len; i++)
    fprintf(out,"%e ",_M[i]);
  fprintf(out,"\n");
  for (unsigned long i=0; i<_Len; i++) {
    for (unsigned long j=0; j<_Len; j++)
      fprintf(out,"%e ",_Sigma[i][j]);
    fprintf(out,"\n");    
  }
  fclose(out);
};
};
