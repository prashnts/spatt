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
#include "gstat.h"

using namespace std;

namespace spatt {

gstat::gstat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool presence,bool verbose)
  : stat(D,Y,S,rep,nobs,verbose) {
  
  if (verbose)
    fprintf(stderr,"start of gstat::gstat()\n");

  if (presence) {
    fprintf(stderr,"presence option not implemented in gstat is ignored\n");
  }
    
  _esp=-1.0;
  _var=-1.0;

  compute_mu(verbose);

};

gstat::gstat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool presence,bool verbose)
  : stat(D,Y,sequence_length,rep,nobs,verbose) {

  if (verbose)
    fprintf(stderr,"start of gstat::gstat()\n");

  if (presence) {
    fprintf(stderr,"presence option not implemented in gstat is ignored\n");
  }

  _esp=-1.0;
  _var=-1.0;

  compute_mu(verbose);

};

void gstat::compute(bool interrupt,int offset,bool verbose){
  vector<double> mu0;
  compute(mu0,interrupt,offset,verbose);
}

void gstat::compute(vector<double> &mu0,bool interrupt,int offset,bool verbose){
  
  if (verbose)
    printf(">>> call gstat::compute()\n"); 
    

  // gather information
  unsigned long sum_ell=0;
  map<unsigned long,vector<unsigned long> > starts;
  for (vector<seq>::iterator it=_seq.begin(); it!=_seq.end(); it++) {
    starts[it->start].push_back(it->length);
    sum_ell+=it->length;
  }

  if (verbose) {
    for (map<unsigned long,vector<unsigned long> >::iterator it=starts.begin(); it!=starts.end(); it++) {
      printf("start = %i, ell = [ ",it->first);
      for (vector<unsigned long>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++)
	printf("%i ",*it2);
      printf("]\n");
    }
  }
  vector<unsigned long> final;
  for (vector<unsigned long>::iterator it=_Dfa->_final.begin(); it!=_Dfa->_final.end(); it++)
    final.push_back(_Y->recode(*it));
  unsigned long F=final.size();
  if (verbose) {
    printf("final (F=%i) = [ ",F);
    for (vector<unsigned long>::iterator it=final.begin(); it!=final.end(); it++)
      printf("%i ",*it);
    printf("]\n");
  }
  unsigned long L=_Y->_nstates;

  //printf("ok\n");

  // if ignore interruption treat data as a single sequence
  if (!interrupt) {
    unsigned long first_start=_seq.begin()->start;
    sum_ell-=offset*_seq.size();
    starts.clear();
    starts[first_start].push_back(sum_ell);
  }

  // allocate working space
  vector<vector<double> > A(F);
  vector<vector<double> > B(F);
  for (unsigned long f=0; f<F; f++) {
    A[f].resize(_alpha);
    B[f].resize(_alpha);
  }
  vector<double> x(L);
  vector<double> aux(L);
  vector<double> y(F);

  
  // initialization
  _esp=0;
  _var=0;
  double muF=0.0;
  for (unsigned long f=0; f<F; f++)
    muF+=_mu[final[f]];
  { // fill all B[f] with e_f Pi^d e_F' for 0<=d<alpha
    // x=e_F
    for (unsigned long i=0; i<L; i++)
      x[i]=0.0;
    for (unsigned long f=0; f<F; f++)
      x[final[f]]=1.0;
    if (_alpha>0) {
      for (unsigned long f=0; f<F; f++)
	B[f][0]=x[final[f]];
    }
    for (unsigned long d=1; d<_alpha; d++) {
      // x = Pi * x
      aux.swap(x);
      _Y->zero(x);
      _Y->add_Mx(_Y->_Pt,aux,x);
      _Y->add_Mx(_Y->_Qt,aux,x);
      // fill B
      for (unsigned long f=0; f<F; f++)
	B[f][d]=x[final[f]];
    } // end d loop      
  } // end compute B
  //if (verbose) {
  //  for (unsigned long f=0; f<F; f++) {
  //    printf("B[%i] = [ ",final[f]);
  //    for (unsigned long d=0; d<_alpha; d++)
  //	printf("%.4f ",B[f][d]);
  //    printf("]\n");
  //  }
  //}

  // loop on starting states
  for (map<unsigned long,vector<unsigned long> >::iterator it=starts.begin(); it!=starts.end(); it++) {
    { // compute A[f] = mu0 Pi^d e_f'
      // x=mu0
      for (unsigned long i=0; i<L; i++)
	x[i]=0.0;
      if (mu0.empty()) {
	x[it->first]=1.0;
      } else {
	unsigned long imax=_Y->_start.size();
	for (unsigned long i=0; i<imax; i++)
	  x[_Y->_start[i]]=mu0[i];
      }
      if (_alpha>0) {
	for (unsigned long f=0; f<F; f++)
	  if (final[f]==it->first)
	    A[f][0]=x[final[f]];
      }
      for (unsigned long d=1; d<_alpha; d++) {
	// x = x * Pi 
	aux.swap(x);
	_Y->zero(x);
	_Y->add_xM(aux,_Y->_Pt,x);
	_Y->add_xM(aux,_Y->_Qt,x);
	// fill A
	for (unsigned long f=0; f<F; f++)
	  A[f][d]=x[final[f]];
      } // end d loop      
    } // end compute A
    //printf("compute A: done\n");
    //if (verbose) {
    //  for (unsigned long f=0; f<F; f++) {
    //	printf("A[%i] = [ ",final[f]);
    //	for (unsigned long d=0; d<_alpha; d++)
    //	  printf("%.4f ",A[f][d]);
    //	printf("]\n");
    //  }
    //}
    // loop on sequence lengths
    unsigned long ell=0;
    for (vector<unsigned long>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++) {
      //if (verbose)
      //	printf("ell=%i\n",*it2);
      // do not treat sequences shorter than _m !!
      ell=*it2;
      if (ell>=_m) {
      double lesp=0.0; // local esperance
      double lsum=0.0; // local sum of C(f,F)
      // y=0
      for (unsigned long f=0; f<F; f++)
	y[f]=0.0;
      if (ell<(_m+_alpha)) {
	for (unsigned long d=0; d<=ell-_m; d++)
	  for (unsigned long f=0; f<F; f++)
	    lesp+=A[f][d];
	for (long i=ell-1; i>=_m; i--) {
	  for (unsigned long f=0; f<F; f++) {
	    y[f]+=B[f][ell-i];
	    lsum+=y[f]*A[f][i-_m];
	  }	    
	}	  
      } else {
	for (unsigned long d=0; d<_alpha; d++)
	  for (unsigned long f=0; f<F; f++)
	    lesp+=A[f][d];
	lesp+=(ell-_m-_alpha+1)*muF;
	unsigned long dmax=ell-_m-_alpha;
	if ((_alpha-1)<dmax)
	  dmax=_alpha-1;
	for (unsigned long d=1; d<=dmax; d++) {
	  for (unsigned long f=0; f<F; f++) {
	    y[f]+=B[f][d];
	    lsum+=(ell-_m-_alpha+1-d)*_mu[final[f]]*B[f][d];
	  }
	} // end loop on d
	if (ell>=(_m+2*_alpha)) {
	  for (unsigned long f=0; f<F; f++)
	    y[f]+=(ell-_m-2*_alpha+1)*muF;
	  lsum+=(ell-_m-2*_alpha+1.0)*(ell-_m-2*_alpha+2.0)/2.0*muF*muF;
	} // end if
	for (long i=_m+_alpha-1; i>=_m; i--) {
	  for (unsigned long f=0; f<F; f++) {
	    if ((ell-i)<_alpha)
	      y[f]+=B[f][ell-i];
	    else
	      y[f]+=muF;
	    lsum+=y[f]*A[f][i-_m];
	  }
	} // end i loop
      }
      if (verbose) {
	printf("start=%i, length=%i, esp=%f, var=%f\n",it->first,ell,lesp,2.0*lsum+lesp-lesp*lesp);
      }
      _esp+=lesp;
      _var+=2.0*lsum+lesp-lesp*lesp;
      }
    } // end loop on sequence lengths
  } // end loop on starting states
  if (verbose) {
    printf("mean=%f, sd=%f\n",_esp,sqrt(_var));
  }
  
  _zscore=(_nobs-_esp)/sqrt(_var);
  
  if (_rep==OVER) {
    if (_zscore>0)
      _pvalue=0.5*erfc(M_SQRT1_2*_zscore);
    else
      _pvalue=0.5+0.5*erf(-M_SQRT1_2*_zscore);
  } else {
    if (_zscore<0)
      _pvalue=0.5*erfc(-M_SQRT1_2*_zscore);
    else
      _pvalue=0.5+0.5*erf(M_SQRT1_2*_zscore);
  }
  //printf("_pvalue=%e\tapprox=%e\n",_pvalue,1/fabs(_zscore)*exp(-0.5*_zscore*_zscore)*M_2_SQRTPI*0.5*M_SQRT1_2);

};

void gstat::compute_mu(bool verbose) {

  if (verbose)
    fprintf(stderr,"gstat::compute_mu()\n");
  _alpha=0;  
  // initialize mu
  _mu.resize(_Y->_nstates);
  {
    double aux=1.0/(double)_Y->_nstates;
    for (unsigned long i=0; i<_Y->_nstates; i++)
      _mu[i]=aux;
  }
  double test=1.0;
  vector<double> aux(_Y->_nstates);
  int iter=0;
  while ((test>MU_TOL)&&(iter<ITER_MAX)) {
    iter++;
    if (verbose)
      printf("test=%e\n",test);
     // old _mu is saved in aux
    aux.swap(_mu);
    // compute _mu = aux * Pi
    _Y->zero(_mu);
    _Y->add_xM(aux,_Y->_Pt,_mu);
    _Y->add_xM(aux,_Y->_Qt,_mu);
    // update _alpha
    _alpha++;
    // update test
    test=fabs(_mu[0]-aux[0]);
    for (unsigned long i=1; i<_Y->_nstates; i++) {
      double tmp=fabs(_mu[i]-aux[i]);
      if (tmp>test)
	test=tmp;
    }      
  } // end while loop
  if (iter==ITER_MAX) {
    _alpha==numeric_limits<long>::max();
  }

  if (verbose) {
    printf("alpha=%i\n",_alpha);
    printf("mu = [ ");
    for (unsigned long i=0; i<_Y->_nstates; i++) 
      printf("%.2e ",_mu[i]);
    printf("]\n");
  }
};

void gstat::print(const char *label,std::string &format_float) {
  printf("pattern=%s\tNobs=%i\tmean=",label,_nobs);
  printf(format_float.c_str(),_esp);
  printf("\tsd=");
  printf(format_float.c_str(),sqrt(_var));
  printf("\tz-score=");
  printf(format_float.c_str(),_zscore);
  printf("\t");
  if (_rep==OVER)
    printf("P(N>=Nobs)=");
  else
    printf("P(N<=Nobs)=");
  printf(format_float.c_str(),_pvalue);
  printf("\n");
};

};
