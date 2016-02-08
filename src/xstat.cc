/* $Id: xstat.cc 486 2005-10-10 09:29:31Z gnuel $ */
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

#include "xstat.h"

namespace spatt {

  xstat::xstat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
		const markov &M,const input &I) : cpstat(params,alpha,seq,occ,M,I)//,_toto(42) 
  {
    _power=NULL;
    _size_power=0;
    _n_power=0;

    /* add here specific stuff */
    /* check for markov model order */
    if (_M->get_order()>0) {
      //fprintf(stderr,"gstat not yet available for Markov model of order > 0\n");
      //exit(EXIT_FAILURE);
      comp_rank();
      update_power();
    }
    _data=NULL;
    _size_data=0;


    _W=new WordFam();

  };
  
  void xstat::compute(word *w){
    this->cpstat::compute(w);
    _k=w->get_size();

    if (_M->get_order()>0) {
      _W->SetWord(w);
      x_comp();
    } else {
      init();
      comp();
    }

    post_comp();
  };
  
  void xstat::compute(word *w1,word *w2){
    this->cpstat::compute(w1,w2);
    _k=w1->get_size();

    if (_M->get_order()>0) {
      _W->SetWord(w1,w2);
      x_comp();
    } else {
      init();
      comp();
    }

    post_comp();
  };

  void xstat::compute(pattern *patt){
    if (patt->get_optimized_word_list_size()>1 || _params->get_both_strands()) {
      /* case of a degenerate pattern */
      /* not yet working */
      if (_M->get_order()>0) {
      //if (false) {
	this->sstat::compute(patt);
	_k=this->sstat::get_pattern_length();
	_W->SetWord(patt);
	x_comp();
      } else {
	fprintf(stderr,"xstat not yet available for degenerate patterns and order 0 Markov model\n");
	exit(EXIT_FAILURE);
      }      
    } else {
      /* case of a simple word */
      this->cpstat::compute(patt);
      _k=this->sstat::get_pattern_length();

      if (_M->get_order()>0) {
	_W->SetWord(patt->get_optimized_word_list()[0]);
	x_comp();
      } else {
	init();
	comp();
      }

      post_comp();
    }
  };  

  void xstat::x_comp() {
    // verif 
    //_W->Print();
    comp_ell();
    //printf("old _stat: %f\n",_stat);
    _P_W= new PAppearFast(_M,_power,_stat,_W,_ell,_ell0) ;
    {
      int Nocc;
      double Probas;//,GaussStats;
      _P_W->GetStatsFromSeq(&Nocc,&Probas,_nobs);//,&GaussStats);
      //printf("Nocc=%i\tProbas=%e\n",Nocc,Probas);
      //printf("_stat=%f\n",_stat);
      if (_stat<0) {
	_stat=log(-Probas)/log(10.0);
      } else {
	_stat=-log(Probas)/log(10.0);
      }
      //printf("_stat=%f\n",_stat);
    }
    
    delete _P_W;
  }

  void xstat::comp_ell() {
    _ell=_occ->get_n();
    _ell0=_ell;
    if (_stat>0.0)
      _nmax=_nobs;
    else
      _nmax=_nobs+1;
    /* check for tail sum */
    if (_stat>=(log(XSTAT_MACHINE_RELATIVE_PRECISION)/log(10.0)-log(XSTAT_TARGET_RELATIVE_PRECISION)/log(10.0))) {
      _tail_sum=false;
    } else {
      _tail_sum=true;
    }
    
    if (_tail_sum) {
      /* compute _ell0 */
      // fixme: 1000 could be too short for some cases
      _ell0=_ell+1000; // to replace by better later
    }
    //printf("ell=%i\tell0=%i\n",_ell,_ell0);
  }

  void xstat::init() {

    comp_ell();

    /* check memory requirements */
    long new_size_data=2*(_ell0+1)+_k;
    
    /* alloc if necessary */
    if (new_size_data>_size_data) {
      if (_data!=NULL)
	free(_data);
      _size_data=new_size_data;
      _data=(double *)malloc(sizeof(double)*_size_data);
      if (_data==NULL) {
	fprintf(stderr,"not enough memory to allocate xstat\n");
	exit(EXIT_FAILURE);
      }      
    }

    /* connect */
    {
      double *pos=_data;
      _old=pos;
      pos+=_ell0+1;
      _current=pos;
      pos+=_ell0+1;
      _ac=pos;
    }

  }
  
  void xstat::comp() {

    /* get stationary distribution */
    double *mu=_M->get_stationary();

    /* compute auto-correlation */
    {
      // rewrite pattern as array of int
      int *w=(int *)malloc(sizeof(int)*(_k+1));
      if (w==NULL) {
	fprintf(stderr,"not enough memory in xstat\n");
	exit(EXIT_FAILURE);
      }
      w[0]=_k;
      for (int i=1; i<=_k; i++) {
	w[i]=_alpha->char2code(_label[i-1]);
      }
	
      // compute mu(pattern)
      double aux=1.0;
      for (int i=1; i<=_k; i++) {
	aux*=mu[w[i]];
      }
      _ac[0]=aux;
      
      // fill _ac
      aux=mu[w[1]];
      if (w[1]==w[_k])
	_ac[_k-1]=_ac[0]/aux;
      else
	_ac[_k-1]=0.0;
      for (int i=2; i<_k; i++) {
	aux*=mu[w[i]];
	int test=1;
	for (int j=1; j<=i; j++)
	  if (w[j]!=w[_k-i+j]) {
	    test=0;
	    break;
	  }
	if (test==1)
	  _ac[_k-i]=_ac[0]/aux;
	else
	  _ac[_k-i]=0.0;
      }
      free(w);
    } // end compute ac    
    
    /* verify _ac */
    if (_params->get_debug_level()>1) {
      printf("_ac:\n");
      for (int i=0; i<_k; i++) {
	printf("\t_ac[%i]=%f\n",i,_ac[i]);
      }
    }

    /* recurrence */
    {
      //double temp=1.0;
      double *aux;
      double sum;
      /* compute the first line */
      for (long x=0; x<_k; x++) {
	_current[x]=0.0;
      }
      _current[_k]=_ac[0];
      //temp-=_current[_k];
      //printf("x=%i\tp(x)=%e\t1.0-cumsum=\n",_k,_current[_k],temp);
      for (long x=_k+1; x<=_ell0; x++) {
	_current[x]=_current[x-1]-_ac[0]*_current[x-_k];
	sum=0.0;
	for (long z=x-_k+1; z<=x-1; z++)
	  if (_ac[x-z]!=0.0)
	    sum+=(_current[z]-_current[z-1])*_ac[x-z];
	_current[x]-=sum;
	//temp-=_current[x];
	//printf("x=%i\tp(x)=%e\t1.0-cumsum=%e\n",x,_current[x],temp);
	//std::cout<<"_current["<<x<<"]->"<<_current[x]<<endl;
      }
      if (_params->get_debug_level()>1) {
	//std::cout<<"n=1 ->"<<_current[_k]<<endl;
	//std::cout<<"n=1 ->"<<_current[_ell0]<<endl;
	printf("n=1/%i\n",_nmax);
      }
      
      /* compute all other lines */
      for (long n=2; n<=_nmax; n++) {
	aux=_old;
	_old=_current;
	_current=aux;
	for (long x=0; x<_k+n-1; x++)
	  _current[x]=0.0;
	for (long x=_k+n-1; x<=_ell0; x++) {
	  _current[x]=_current[x-1]+_ac[0]*(_old[x-_k]-_current[x-_k]);
	  sum=0.0;
	  for (long z=x-_k+1; z<=x-1; z++)
	    if (_ac[x-z]!=0.0)
	      sum+=(_old[z]+_current[z-1]-_current[z]-_old[z-1])*_ac[x-z];
	  _current[x]+=sum;
	}
	if (_params->get_debug_level()>1) {
	  //std::cout<<"n="<<n<<" ->"<<_current[_k+n-1]<<endl;
	  //std::cout<<"n="<<n<<" ->"<<_current[_ell0]<<endl;
	  printf("n=%i/%i\n",n,_nmax);
	}
      }
    } // end recurrence

    /* time to sum */
    {
      double res,p,q;
      if (_tail_sum) {
	/* _tail_sum true */
	q=0.0;
	for (long x=_ell+1; x<=_ell0; x++)
	  q+=_current[x];
	res=log(q)/log(10.0);
	printf("sum_%i^%i=%e\tlast term=%e\n",_ell+1,_ell0,q,_current[_ell0]);
      } else {
	/* _tail_sum false */
	if (_stat<0.0) {
	  // compute under-representation
	  p=0.0;
	  q=1.0;    
	  for (long x=0; x<=_ell; x++)
	    p+=_current[x];
	  q-=p;
	  //std::cout<<"p="<<p<<endl;
	  //std::cout<<"q="<<q<<endl;
	  res=log(q)/log(10.0);
	} else {
	  // compute over-representation
	  p=0.0;
	  if (_nobs==0) {
	    p=1.0;
	  } else {
	    for (long x=0; x<=_ell; x++)
	      p+=_current[x];
	  }
	  res=-log(p)/log(10.0);
	}      
      } // end if _tail_sum
      _stat=res;
    } //end summing part
  }

  xstat::~xstat() {
    if (_data!=NULL)
      free(_data);
    if (_power!=NULL) {
      free(_power);
    }
    delete _W;
  };

  void xstat::comp_rank() {
    /* evaluate _rank */
    double res;
    /* temp evaluation through MACHINE precision */
    res=log(XSTAT_MACHINE_RELATIVE_PRECISION)/log(_M->get_secondmag());
    //printf("res=%f\n",res);
    _rank=(long)res+1+XSTAT_SECURITY_MARGIN;
  }

  void xstat::update_power(){
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

  }


  void xstat::print_format(FILE *fout) {
    fprintf(fout,"pattern\tnobs\tnatt\tstat\n");
  };

  void xstat::print_regular(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%f\n",_label,_nobs,_natt,_stat);
  };

  void xstat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%e\n",_label,_nobs,_natt,_stat);
  }

};
