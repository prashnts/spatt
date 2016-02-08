/* $Id: ldstat.cc 891 2006-09-27 07:12:56Z gnuel $ */
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

#include "ldstat.h"

double ldstat_sparse_F(double t,void *that) {
  //printf("ldstat_sparse_F(%f,%p)\n",t,that);
  spatt::ldstat *current=(spatt::ldstat *)that;
  return current->sparse_F(t,NULL);
}

double ldstat_full_F(double t,void *that) {
  //printf("ldstat_full_F(%f,%p)\n",t,that);
  spatt::ldstat *current=(spatt::ldstat *)that;
  return current->full_F(t,NULL);
}

double ldstat_full_rate(double a,void *that) {
  spatt::ldstat *current=(spatt::ldstat *)that;
  return current->full_rate(a,NULL);
}

double ldstat_sparse_rate(double a,void *that) {
  spatt::ldstat *current=(spatt::ldstat *)that;
  return current->sparse_rate(a,NULL);
}

namespace spatt {

  ldstat::ldstat(spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
	       const markov &M,const input &I) : sstat(params,alpha,seq,occ,M,I) //,_toto(42) 
  {
    // fixme: is it necessary to remove const for params to do this dynamic cast
    ldspattparameters *sparams=dynamic_cast<ldspattparameters *>(&params);      
    /* set default values */
    _rate=NULL;
    _drate=NULL;
    _ddrate=NULL;
    _precise=sparams->get_precise();
    _sparse_par=NULL;
    _full_par=NULL;
 
    _T=new transition(params,alpha,seq,occ,M,I);

    // allocate chebyshev series if necessary
    if (_precise) {
      _rate=gsl_cheb_alloc (CHEBYSHEV_ORDER);
      _drate=gsl_cheb_alloc (CHEBYSHEV_ORDER);
      _ddrate=gsl_cheb_alloc (CHEBYSHEV_ORDER);
    }
  };


  void ldstat::compute(word *w){
    this->sstat::compute(w);
    _T->build(w);
    comp();
  };
  
  void ldstat::compute(word *w1,word *w2){
    this->sstat::compute(w1,w2);
    _T->build(w1,w2);
    comp();
  };
  
  void ldstat::compute(pattern *patt){
    this->sstat::compute(patt);   
    _T->build(patt);
    comp();
  };

  void ldstat::comp() {
    //printf("calling of ldstat::comp()\n");
    _order=_T->get_order();
    _In=0.0;
    //printf("get_pattern_length()=%i\n",get_pattern_length());
    _n=_occ->get_n()-get_pattern_length()+1;
    _a=(double)_nobs/(double)_n;
    //print_transition();
    //print_fixed_transition();
    //print_counting();
    /* all the work is now done here */
    if (_order<=MAX_FULL_DIM) {
      /* full case */
      //fprintf(stderr,"warning ! full case not yet implemented in ldstat, using sstat instead\n");
      /* alloc or adjust full_par */
      if (_full_par==NULL) {
	_full_par=new dgeev_par(_order,0,1);
	if (_full_par==NULL) {
	  fprintf(stderr,"Not enough memory to alloc full_par in ldstat\n");
	  exit(EXIT_FAILURE);
	}
      } else {
	//printf("in comp(): _order=%i\n",_order);
	_full_par->reset(_order,0,1,0);
      }
      /* rate function evaluation */
      _In=rate(&ldstat_full_F);

      //printf("eta=%e\tI=%e\n",_eta,_In);

      /* precise case */
      //printf("_precise=%i && _nobs=%i\n",_precise,_nobs);
      if (_precise && _nobs!=0) {
	// store original value of _a
	
	_ori_a=_a;	
	_a=0.0;
	gsl_function F;
	F.function = &ldstat_full_F;
	F.params = this;
	gsl_cheb_init(_rate,&F,_eta-CHEBYSHEV_HALF_SPAN,_eta+CHEBYSHEV_HALF_SPAN);
	_a=_ori_a;
	gsl_cheb_calc_deriv(_drate,_rate);
	gsl_cheb_calc_deriv(_ddrate,_drate);
	//printf("verification (full): L'(%e)=%e\ta=%e\n",_eta,gsl_cheb_eval(_drate,_eta),_ori_a);
	_V=gsl_cheb_eval(_ddrate,_eta);
	
	//printf("_V=%f\n",_V);
	if (_V<=0)
	  fprintf(stderr,"ignored precise correction factor V(%e)=%e < 0 for pattern '%s'\n",_a,_V,_label);
	
	  
      } // end precise

      
    } else {
      /* sparse case */
      //fprintf(stderr,"warning ! sparse case not yet implemented in ldstat, using sstat instead\n");
      /* alloc or adjust sparse_par */
      if (_sparse_par==NULL) {
	_sparse_par=new dnaupd_par(_order,MAX_FULL_DIM);
	if (_sparse_par==NULL) {
	  fprintf(stderr,"Not enough memory to alloc sparse_par in ldstat\n");
	  exit(EXIT_FAILURE);
	}
      } else {
	_sparse_par->reset(_order,MAX_FULL_DIM);
      }
      /* rate function evaluation */
      _In=rate(&ldstat_sparse_F);

      if (_precise && _nobs!=0) {
	// store original value of _a
	_ori_a=_a;	
	_a=0.0;
	gsl_function F;
	F.function = &ldstat_sparse_F;
	F.params = this;
	gsl_cheb_init(_rate,&F,_eta-CHEBYSHEV_HALF_SPAN,_eta+CHEBYSHEV_HALF_SPAN);
	_a=_ori_a;
	gsl_cheb_calc_deriv(_drate,_rate);
	gsl_cheb_calc_deriv(_ddrate,_drate);
	//printf("verification (sparse): L'(%e)=%e\ta=%e\n",_eta,gsl_cheb_eval(_drate,_eta),_ori_a);	
	_V=gsl_cheb_eval(_ddrate,_eta);
	
	if (_V<=0)
	  fprintf(stderr,"ignored precise correction factor V(%e)=%e < 0 for pattern '%s'\n",_a,_V,_label);
	
	  
      } // end precise

    }
    /* computing stat from _In */
    //printf("precise: V=%e\n",_V);
    if (_stat>0) {
      _stat=-_n*_In/log(10.0);
      if (_precise && _nobs!=0 && _V>0.0)
	_stat+=0.5*log(2.0*LDSTAT_PI*_n*_V)/log(10.0)+log((1-exp(-fabs(_eta))))/log(10.0);
      {
	if (_stat<0.0)
	  _stat=0.0;
      }
    } else {
      _stat=_n*_In/log(10.0);
      if (_precise && _nobs!=0 && _V>0.0)
	_stat-=0.5*log(2.0*LDSTAT_PI*_n*_V)/log(10.0)+log((1-exp(-fabs(_eta))))/log(10.0);
      {
	if (_stat>0.0)
	  _stat=0.0;
      }
    }
    //printf("stat=%e\n",_stat);
    post_comp();
    ///* check normalize */
    //if (_params->get_normalize()) {
    //  _stat/=(double)_n;
    //}
    
  }

  double ldstat::sparse_F(double t,void *not_used) {
   
    double res;
    _sparse_par->reset(_order,MAX_FULL_DIM);

    /* multiply by exp(t) all counting parameters */
    double aux=exp(t);
    _T->multiply_Q_by(aux);

    /* main loop */
  main:
    dnaupd_(_sparse_par->_ido,_sparse_par->_bmat,_sparse_par->_n,_sparse_par->_which,_sparse_par->_nev,_sparse_par->_tol,_sparse_par->_resid,_sparse_par->_ncv,_sparse_par->_v,_sparse_par->_ldv,_sparse_par->_iparam,_sparse_par->_ipntr,_sparse_par->_workd,_sparse_par->_workl,_sparse_par->_lworkl,_sparse_par->_info);
    //printf("dnaupd_call done\n");
    if (_sparse_par->_ido==-1 || _sparse_par->_ido==1) {
      /* compute product Y=OP*X where */
      /* ipntr[0]-1 is the pointer into workd for X */
      /* ipntr[1]-1 is the pointer into workd for Y */ 
      //printf("matrix x vector product\n");
      double *x=&_sparse_par->_workd[_sparse_par->_ipntr[0]-1];
      double *y=&_sparse_par->_workd[_sparse_par->_ipntr[1]-1];
      _T->transition_by_vector(x,y);
      goto main; 
    } /* end if ido */

    //printf("implicitly restarted algorithm returned:\n");
    //printf("\tinfo=%i\n",_sparse_par->_info);
    //printf("\tnumber of iter=%i\n",_sparse_par->_iparam[2]);
    //printf("\tnumber of converged ritzvalues=%i\n",_sparse_par->_iparam[4]);
    //printf("\tnumber of prod=%i\n",_sparse_par->_iparam[8]);
    //printf("\tnumber of re-orth steps=%i\n",_sparse_par->_iparam[10]);    

    /* get the result */
    /* there is no check here that the first one is the right one */
    /* but this seems to be allways the case */
    {
      //printf("ipntr:");
      //for (int i=0; i<14; i++)
      //  printf("[%i]:%i\t",i+1,_sparse_par->_ipntr[i]);
      //printf("\n");
      double *re,*im;
      re=&_sparse_par->_workl[_sparse_par->_ipntr[5]-1];
      //im=&_sparse_par->_workl[_sparse_par->_ipntr[6]-1];
      //printf("eigenvalues:\n");
      //for (int i=0; i<_sparse_par->_ncv; i++) {
      //  printf("\tre=%f\tim=%f\n",re[i],im[i]);
      //}
      res=re[0];
    }

    /* restore initial values */
    _T->multiply_Q_by(1.0/aux);

    //printf("sparse_F(%f)=log(%f)-%e*%f=%e\n",t,res,_a,t,log(res)-_a*t);
    
    return log(res)-_a*t;
  }

  double ldstat::full_rate(double a,void *not_used){
    double res=0.0;
    _a=a;
    res=-rate(&ldstat_full_F);
    return res;
  }

  double ldstat::sparse_rate(double a,void *not_used){
    double res=0.0;
    _a=a;
    res=-rate(&ldstat_sparse_F);
    return res;
  }

  double ldstat::full_F(double t,void *not_used) {
   
    double res=1.0;

    _full_par->reset(_order,0,1,0);

    //printf("full_F(%f,NULL) [_order=%i]\n",t,_order);
    /* multiply by exp(t) all counting parameters */
    double aux=exp(t);
    _T->multiply_Q_by(aux);

    /* filling matrix _AA */
    _T->fill_fortran(_AA);
    /* get the result */
    /* workspace query */
    _full_par->_lwork=-1;
    //_full_par->print();
    //printf("_lda=%i\n",_full_par->_lda);
    dgeev_(_full_par->_jobvl,_full_par->_jobvr,_full_par->_n,_AA,_full_par->_lda,_full_par->_wr,_full_par->_wi,_full_par->_vl,_full_par->_ldvl,_full_par->_vr,_full_par->_ldvr,_full_par->_work,_full_par->_lwork,_full_par->_info);
    /* realloc */
    _full_par->reset(_order,0,1,(int)_full_par->_work[0]);
    //_full_par->print();
    //printf("_lda=%i\n",_full_par->_lda);
    /* real call */
    dgeev_(_full_par->_jobvl,_full_par->_jobvr,_full_par->_n,_AA,_full_par->_lda,_full_par->_wr,_full_par->_wi,_full_par->_vl,_full_par->_ldvl,_full_par->_vr,_full_par->_ldvr,_full_par->_work,_full_par->_lwork,_full_par->_info);
    /* get res */
    {
      res=-1.0;
      //printf("eigenvalues:\n");
      for (int i=0; i<_order; i++) {
        //printf("\tre=%f\tim=%f\n",_full_par->_wr[i],_full_par->_wi[i]);	
	// check for real eigenvalue
	if (fabs(_full_par->_wi[i])<1e-10) {
	  if (_full_par->_wr[i]>res)
	    res=_full_par->_wr[i];
	}	
      }
      //printf("res=%f\n",res);
      if (res<0) {
	fprintf(stderr,"no real and >0 largest eigenvalue find in full_F\n");
	exit(EXIT_FAILURE);
      }
    }
    /* verification */
    //printf("A:\n");
    //for (int i=0; i<_order; i++) {
    //  printf("\t");
    //  for (int j=0; j<_order; j++)
    //    printf("%.2f\t",_AA[i+_order*j]);
    //  printf("\n");
    //}
    
    /* restore initial values */
    _T->multiply_Q_by(1.0/aux);
    
    //printf("full_F(%f)=log(%f)-%e*%f=%e\n",t,res,_a,t,log(res)-_a*t);
    return log(res)-_a*t;
  }

  void ldstat::mnbrack(double *ax,double *bx,double *cx, double *fa,double *fb, double *fc,
		       double (*func)(double,void*)){
  
    double ulim,u,r,q,fu,dum;
    
    *fa=(*func)(*ax,this);
    *fb=(*func)(*bx,this);
    if (*fb>*fa) {
      SHFT(dum,*ax,*bx,dum);
      SHFT(dum,*fa,*fb,dum);
    }
    *cx=(*bx)+GGOLD*(*bx-*ax);
    *fc=(*func)(*cx,this);
    while (*fb>*fc) {
      r=(*bx-*ax)*(*fb-*fc);
      q=(*bx-*cx)*(*fb-*fa);
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(DDMAX(fabs(q-r),TINY),q-r));
      ulim=(*bx)+GLIMIT*(*cx-*bx);
      if ((*bx-u)*(u-*cx) > 0.0) {
	fu=(*func)(u,this);
	if (fu < *fc) {
	  *ax=(*bx);
	  *bx=u;
	  *fa=(*fb);
	  *fb=fu;
	  return;
	} else if (fu > *fb) {
	  *cx=u;
	  *fc=fu;
	  return;
	}
	u=(*cx)+GGOLD*(*cx-*bx);
	fu=(*func)(u,this);
      } else if ( (*cx-u)*(u-ulim) > 0.0 ) {
	fu=(*func)(u,this);
	if (fu < *fc) {
	  SHFT(*bx,*cx,u,*cx+GGOLD*((*cx-*bx)));
	  SHFT(*fb,*fc,fu,(*func)(u,this));
	}
      } else if ( (u-ulim)*(ulim-*cx) >= 0.0 ) {
	u=ulim;
	fu=(*func)(u,this);
      } else {
	u=(*cx)+GGOLD*(*cx-*bx);
	fu=(*func)(u,this);
      }
      SHFT(*ax,*bx,*cx,u);
      SHFT(*fa,*fb,*fc,fu);
    } /* end while */
  };

  double ldstat::rate(double (*func)(double,void*)){
    //  double ldstat::sparse_rate(){
  
    double ax,bx,cx,fa,fb,fc;
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;


    //F.function = &ldstat_sparse_F;
    F.function = func;
    F.params = this;

    //printf("calling rate\n");
    /* special case first */
    if (_a==0.0) {
      return  (*F.function)(-1000,this);
    } else {
      /* brackting the miniminum */
      ax=-1.0; bx=1.0;
      mnbrack(&ax,&bx,&cx,&fa,&fb,&fc,(*F.function));
      //printf("[%e,%e,%e] with function values = [%e,%e,%e]\n",ax,bx,cx,fa,fb,fc);
      /* brent minimizer */
      T = gsl_min_fminimizer_brent;
      s = gsl_min_fminimizer_alloc (T);
      if (ax<cx)
	gsl_min_fminimizer_set_with_values (s,&F,bx,fb,ax,fa,cx,fc);
      else
	gsl_min_fminimizer_set_with_values (s,&F,bx,fb,cx,fc,ax,fa);
      /* main loop */
      do
	{
	  iter++;
	  status = gsl_min_fminimizer_iterate (s);

	  bx = gsl_min_fminimizer_x_minimum (s);
	  ax = gsl_min_fminimizer_x_lower (s);
	  cx = gsl_min_fminimizer_x_upper (s);
	  fb = gsl_min_fminimizer_f_minimum (s);
	  fa = gsl_min_fminimizer_f_lower (s);
	  fc = gsl_min_fminimizer_f_upper (s);
	//printf("[%e,%e,%e] with function values = [%e,%e,%e]\n",ax,bx,cx,fa,fb,fc);	
	//	  printf("rel errors: %e\t%e\n",(fa-fb)/fb,(fc-fb)/fb);
	  status=gsl_min_test_interval(ax, cx, 0.001, 0.0);

	//  if (status == GSL_SUCCESS)
	//    printf ("Converged:\n");
	//  
	 // printf("%i\tF(%e)=%e\n",iter,bx,fb);
	}
      while (status == GSL_CONTINUE && iter < max_iter);
    }
    
    _eta=bx;
    return fb;
  };


  ldstat::~ldstat() {
    delete _T;
    if (_sparse_par!=NULL)
      delete _sparse_par;
    if (_full_par!=NULL)
      delete _full_par;
    // free chebyshev series if necessary
    if (_rate!=NULL)
      gsl_cheb_free(_rate);
    if (_drate!=NULL)
      gsl_cheb_free(_drate);
    if (_ddrate!=NULL)
      gsl_cheb_free(_ddrate);

  };
  
};
