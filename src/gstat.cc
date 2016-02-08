/* $Id: gstat.cc 504 2005-10-19 11:53:32Z mhoebeke $ */
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

#include "gstat.h"

namespace spatt {

  gstat::gstat(const gspattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
	       const markov &M,const input &I) : stat(params,alpha,seq,occ,M,I) //,_toto(42) 
  {
    /* add here specific stuff */
    _sd=0.0;
    _power=NULL;
    _size_power=0;
    _n_power=0;

  };

  void gstat::compute(word *w){
    this->stat::compute(w);
    precomp(w);
    comp();
  };
  
  void gstat::compute(word *w1,word *w2){
    fprintf(stderr,"gstat not yet available for degenerate patterns\n");
    exit(EXIT_FAILURE);
  };
  
  void gstat::compute(pattern *patt){
    if (patt->get_optimized_word_list_size()>1 || _params->get_both_strands()) {
      /* case of a degenerate pattern */
      fprintf(stderr,"gstat not yet available for degenerate patterns\n");
      exit(EXIT_FAILURE);
    } else {
      /* case of a simple word */
      this->stat::compute(patt);
      {
	precomp(patt->get_optimized_word_list()[0]);
	comp();
      }
    }
  };

  void gstat::precomp(word *w) {
    _natt=_M->expect(*w);
    comp_overlap(w);
    comp_main(w);
//    printf("overlap=%f\n",comp_overlap(w));
//    printf("main=%f\n",comp_main(w));
//    printf("_natt=%f\n",_natt);
//    printf("variance=%f\n",_main+_overlap-_natt*_natt);
//    printf("sd=%f\n",sqrt(_main+_overlap-_natt*_natt));
    _sd=sqrt(_main+_overlap-_natt*_natt);
  }

  void gstat::comp(){

    double Z=(_nobs-_natt)/_sd;
    _zscore=Z;

    //printf("E=%f\n",_natt);
    //printf("sd=%f\n",_sd);      
    //printf("_zscore=%f\n",_zscore);
    
    /* compute statistic from Z score */
    {
      double m,e;
      if (Z>0) {
	_stat=-qcdfnor(Z,&m,&e);
      } else {
	_stat=pcdfnor(Z,&m,&e);
      }
    }
    //printf("stat=%f\n",_stat);
    post_comp();
  };

  gstat::~gstat() {
    if (_power!=NULL) {
      free(_power);
    }

  };
  
  void gstat::print_format(FILE *fout) {
    // fixme: better idea ?
    if (dynamic_cast<gspattparameters *>(const_cast<spattparameters *>(_params))->get_output_zscore())
      fprintf(fout,"pattern\tnobs\tnatt\tsd\tzscore\n");
    else
      fprintf(fout,"pattern\tnobs\tnatt\tsd\tstat\n");
  };

  void gstat::print_regular(FILE *fout) {
    // fixme: better idea ?
    if (dynamic_cast<gspattparameters *>(const_cast<spattparameters *>(_params))->get_output_zscore())
      fprintf(fout,"%s\t%li\t%.2f\t%.2f\t%f\n",_label,_nobs,_natt,_sd,_zscore);
    else
      fprintf(fout,"%s\t%li\t%.2f\t%.2f\t%f\n",_label,_nobs,_natt,_sd,_stat);
  };

  void gstat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%.2f\t%e\n",_label,_nobs,_natt,_sd,_stat);
  }

  double gstat::comp_overlap(word *w) {
    return comp_overlap(w,w);
  }

  double gstat::comp_overlap(word *v,word *w) {
    _g=v->get_size();
    _h=w->get_size();
    _ell=_occ->get_n();
    //printf("_ell=%i\n",_ell);
    
    double res=0.0;
    if (_h>_g) {
      res=comp_overlap(w,v);
    } else {
	for (int i=1; i<=(_g-1); i++) {
	  res+=(_ell-_h-i+1)*overlap(v,w,i);
	}
	for (int i=1; i<=(_h-1); i++) {
	  res+=(_ell-_g-i+1)*overlap(w,v,i);
	}
	for (int i=1; i<=(_g-_h+1); i++) {
	  res+=(_ell-_g+1)*fit(v,w,i);
	}
    }
    _overlap=res;
    return _overlap;
  }

  double gstat::comp_main(word *w) {
    return comp_main(w,w);
  }

  double gstat::comp_main(word *v,word *w) {
    _g=v->get_size();
    _h=w->get_size();
    _ell=_occ->get_n();
    
    double res=0.0;
    if (_h>_g) {
      res=comp_main(w,v);
    } else {
      /* computations here */
      if (_M->get_order()<=0) {
	/* simple case */
	//	res=(_ell-_ell*_g-_ell*_h+_g*_g-2*_g+_h*_h-2*_h+_g*_h+1)*_M->mu(v->get_label())*_M->mu(w->get_label());
	//printf("_ell-_g-_h+1=%i\n",_ell-_g-_h+1);
	//printf("_M->mu(v->get_label())=%f\n",_M->mu(v->get_label()));
	res=(double)(_ell-_g-_h+1)*(double)(_ell-_g-_h+2)*_M->mu(v->get_label())*_M->mu(w->get_label());
      } else {
	/* Markov case */
	/* evaluate _rank */
	{
	  double res;
	  /* temp evaluation through MACHINE precision */
	  res=log(GSTAT_MACHINE_RELATIVE_PRECISION)/log(_M->get_secondmag());
//
//	  res=log(GSTAT_TARGET_RELATIVE_PRECISION);
//	  res+=log(1-_M->get_secondmag());
//	  {
//	    double aux;
//	    aux=_ell;
//	    aux*=(_nobs-_natt)*(_nobs-_natt)/_natt/_natt;
//	    {
//	      double p,m;
//	      string s(w->get_label());
//	      p=_M->mu(&s);
//	      string ss(s,0,_M->get_order()-1);
//	      m=_M->mu(&ss);
//	      aux*=p*p/m;
//	    }
//	    res+=log(aux);
//	  }
//	  res/=log(_M->get_secondmag());
//	  res-=1.0;
	  //printf("res=%f\n",res);
	  _rank=(long)res+1+SECURITY_MARGIN;
	} // end rank evaluation
	//printf("_rank=%i\n",_rank);
	/* alloc memory for Pi^1 ... Pi^_rank */
	{
	  _order=_M->get_dim()*_alpha->get_size();
	  _one_size=_order*_order;
	  long new_size=_one_size*_rank;
	  if (new_size>_size_power) {
	    /* alloc new one */
	    double *new_power;
	    new_power=(double *)malloc(sizeof(double)*new_size);
	    if (new_power==NULL) {
	      fprintf(stderr,"not enough memory in gstat\n");
	      exit(EXIT_FAILURE);
	    }
	    /* copy last content */
	    for (long i=0; i<_size_power; i++) {
	      new_power[i]=_power[i];
	    }
	    /* free last power */
	    if (_power!=NULL)
	      free(_power);
	    _power=new_power;
	    _size_power=new_size;
	  }
	}
	/* fill power */
	if (_rank>_n_power) {
	  if (_n_power==0) {
	    /* fill the first matrix with Pi */
	    //printf("todo: fill _power with Pi^1\n");
	    /* fill matrix a in the fortran way */
	    {
	      int _k=_alpha->get_size();
	      long _dim=_M->get_dim();
	      long _n=_dim*_k;
	      double **_model=_M->get_model();
	      double *a=_power;
	      int shift=0;
	      int jrot=0;
	      int irot=0;
	      int pos=0;
	      for (int j=0; j<_n; j++) {
		int mpos=0;
		// column j
		for (int i=0; i<_n; i++) {
		  //printf("A[%i][%i]=",i,j);
		  // line i
		  if (irot==shift) {
		    a[pos]=_model[mpos][jrot];
		    shift+=_dim;
		  } else {
		    a[pos]=0.0;
		  }
		  mpos++;
		  pos++;
		  irot++;
		}
		jrot++;
		if (jrot==_k) {
		  jrot=0;
		  shift++;
		}
	      }
	      // verif
	      //for (int i=0; i<_n; i++) {
	      //  for (int j=0; j<_n; j++) {
	      //    printf("%f ",a[j*_n+i]);
	      //  }
	      //  printf("\n");
	      //}
	    } /* end fill matrix a */
	    _n_power=1;
	  }
	  for (long i=(_n_power+1); i<=(_rank); i++) {
	    /* compute Pi^i */
	    //printf("todo: fill _power with Pi^%i\n",i);
	    {
	      double *one=_power;
	      double *last=_power;
	      double *current=_power;
	      last+=(i-2)*_one_size;
	      current+=(i-1)*_one_size;
	      //{
	      //  double *aux=one;
	      //  printf("Pi^1=\n");
	      //  for (int i=0; i<_order; i++) {
	      //    for (int j=0; j<_order; j++) {
	      //      printf("%f ",aux[j*_order+i]);
	      //    }
	      //    printf("\n");
	      //  }
	      //}
	      //{
	      //  double *aux=last;
	      //  printf("Pi^%i=\n",i-1);
	      //  for (int i=0; i<_order; i++) {
	      //    for (int j=0; j<_order; j++) {
	      //      printf("%f ",aux[j*_order+i]);
	      //    }
	      //    printf("\n");
	      //  }
	      //}
	      // now compute Pi^i 
	      {
		for (int i=0; i<_order; i++) {
		  for (int j=0; j<_order; j++) {
		    // current[i][j]
		    double aux=0.0;
		    for (int k=0; k<_order; k++) {
		      // one[i][k]*last[k][j]
		      aux+=one[k*_order+i]*last[j*_order+k];
		    }
		    current[j*_order+i]=aux;
		  }
		}
	      }
	      //{
	      //  double *aux=current;
	      //  printf("Pi^%i=\n",i);
	      //  for (int i=0; i<_order; i++) {
	      //    for (int j=0; j<_order; j++) {
	      //      printf("%f ",aux[j*_order+i]);
	      //    }
	      //    printf("\n");
	      //  }
	      //}
	    }
	  }
	  _n_power=_rank;
	} // power filled
	/* now compute the main part */
	{
	  double pv=_M->mu(v->get_label());
	  double pw=_M->mu(w->get_label());
	  res=(double)(_ell-_g-_h-_rank+1)*(double)(_ell-_g-_h-_rank)*pv*pw;
	  // compute code for end and start of words
	  // todo
	  long codev1=0,codev2=0; // code of v[1-m] and v[_g-m+1-_g]
	  long codew1=0,codew2=0; // code of v[1-m] and v[_g-m+1-_g]
	  double *aux=_power;
	    for (int i=1; i<=_rank; i++) {
	      // Pi^i[codev2][codew1]+Pi^i[codew2][codew1]
	      res+=(_ell-_g-_h-i)*(aux[codev2+codew1*_order]+aux[codew2+codew1*_order])*pv*pw;
	      aux+=_one_size;
	    }
	}
      } // end Markov case
    } // end if (_h>_g) 
    _main=res;
    return _main;
  }

  double gstat::overlap(word *v,word *w,int d){
    _g=v->get_size();
    _h=w->get_size();
    char *sv=v->get_label();
    char *sw=w->get_label();
    if (d<1 || d>=_g) {
      fprintf(stderr,"third parameter out of range in gstat::overlap(%s,%s,%i)\n",sv,sw,d);
      exit(EXIT_FAILURE);
    }
    //printf("*****************\noverlap(%s,%s,%i)\n",sv,sw,d);
    /* check if v[2-g] == w [1-(g-d)] */
    bool test=true;
    for (int i=1; i<=(_g-d); i++) {
      if (sv[d+i-1]!=sw[i-1]) {
	test=false;
	break;
      }
    }
    if (test) {
      string s1(sw,_g-d,_h-1);
      string s2(sv);
      s1=s2+s1;
      double res;
      int m=_M->get_order();
      if (m<0)
	m=1;
      if (_g>=m) {
	res=_M->mu(&s1);
	return res;
      } else {
	fprintf(stderr,"overlap(word,word,int) not yet implemented when Markov order is greater than a word size\n");
	exit(EXIT_FAILURE);
      }
      
    } else {
      /* return 0.0 */
      //printf("incompatible overlap\n");
      return 0.0;
    }    
  }

  double gstat::fit(word *v,word *w,int d){
    _g=v->get_size();
    _h=w->get_size();
    char *sv=v->get_label();
    char *sw=w->get_label();
    //printf("*****************\nfit(%s,%s,%i)\n",sv,sw,d);
    if (d<1 || d>(_g-_h+1)) {
      fprintf(stderr,"third parameter out of range in gstat::fit(%s,%s,%i)\n",sv,sw,d);
      exit(EXIT_FAILURE);
    }
    if (_h>_g) {
      fprintf(stderr,"'%s' must be longer than '%s' out of range in gstat::fit(%s,%s,%i)\n",sv,sw,sv,sw,d);
      exit(EXIT_FAILURE);
    }
    /* check if w == v [d-(d+h-1)] */
    bool test=true;
    for (int i=1; i<=_h; i++) {
      if (sw[i-1]!=sv[d-1+i-1]) {
	test=false;
	break;
      }
    }
    if (test) {
      //printf("res=%f\n",_M->mu(sw));
      return _M->mu(sw);
    } else {
      /* does not fit in */
      //printf("incompatible overlap\n");
      return 0.0;
    }

  }


};
