/* $Id: newxstat.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "newxstat.h"

namespace spatt {

  newxstat::newxstat(spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
		     const markov &M,const input &I) : sstat(params,alpha,seq,occ,M,I)//,_toto(42) 
  {
    _workspace=NULL;
    _workspace_size=0;
    _T = new transition(params,alpha,seq,occ,M,I);
    _n=_seq->get_valid_char();

  };
  

  newxstat::~newxstat() {
    delete _T;
    if (_workspace!=NULL)
      free(_workspace);
  };


  void newxstat::compute(word *w){
    this->sstat::compute(w);
    _T->build(w);
    if (fabs(_stat)>LIMIT_LOG_COMP)
      log_comp();
    else
      comp();
  };
  
  void newxstat::compute(word *w1,word *w2){
    this->sstat::compute(w1,w2);
    _T->build(w1,w2);
    if (fabs(_stat)>LIMIT_LOG_COMP)
      log_comp();
    else
      comp();
  };
  
  void newxstat::compute(pattern *patt){
    this->sstat::compute(patt);   
    _T->build(patt);
    if (fabs(_stat)>LIMIT_LOG_COMP)
      log_comp();
    else
      comp();
  };

  void newxstat::comp() {
    _order=_T->get_order();
    if (_stat>0) {
      _nblock=_nobs;
    } else {
      _nblock=_nobs+1;
    }
    _dim=_order*_nblock;
    /* check memory requirements */
    {
      unsigned long new_workspace_size=2*_dim+_order;
      if (new_workspace_size>_workspace_size) {
	/* realloc space */
	if (_workspace!=NULL) {
	  free(_workspace);	  
	}
	_workspace_size=new_workspace_size;
	_workspace=(double *)malloc(sizeof(double)*_workspace_size);
	if (_workspace==NULL) {
	  fprintf(stderr,"not enough memory for workspace allocation in xstat\n");
	  exit(EXIT_FAILURE);
	}
      }
    } /* end memory requirements */
    /* connect stuff */
    {
      double *dpos=_workspace;
      _current_u=dpos;
      dpos+=_dim;
      _last_u=dpos;
      dpos+=_dim;
      _sum=dpos;
    }
    /* two cases: over and under */
    if (_stat>0) {
      /* over case */
      // initialization
      for (long i=0; i<_dim; i++) {
	_current_u[i]=0.0;
      }
      _T->sum_Q(_current_u);
      // verif
      //for (long i=0; i<_dim; i++) {
      //  printf("u[%i]=%f\n",i,_current_u[i]);
      //}      
      for (long i=0; i<_order; i++) {
	_sum[i]=0.0;
      }
      double *aux=NULL;
      double *pos_in_current,*pos_in_last1,*pos_in_last2;
      // main loop
      //int percent=0,newpercent=0;
      for (long i=2; i<=_n-1; i++) {
	//newpercent=(100*i)/_n;
	//if (newpercent>percent) {
	//  percent=newpercent;
	//  fprintf(stderr,"%i %%\r",(int)(100*(double)i/_n));
	//  fflush(stderr);
	//}
	// switch _last_u and _current_u;
	aux=_last_u;
	_last_u=_current_u;
	_current_u=aux;
	// put fmci transition times _last_u in _current_u
	pos_in_current=_current_u;
	pos_in_last1=_last_u;
	// block 0
	// _current_u_0 = P * _last_u_0
	_T->P_times(pos_in_last1,pos_in_current);		
	// loop on blocks
	for (long j=1; j<_nblock; j++) {
	  // new positions
	  pos_in_current+=_order;
	  pos_in_last2=pos_in_last1;
	  pos_in_last1+=_order;
	  // _current_u_j = P * _last_u_j + Q * _last_u_j-1
	  _T->P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	} // end blocks loop
	// _sum = _sum + last block of current_u (pointed by pos_in_current)
	for (long j=0; j<_order; j++) {
	  _sum[j]+=pos_in_current[j];
	}
      }
      // get the final result
      {
	//printf("stationary size = %i\n",_M->get_stationary_size());
	int imax=_M->get_stationary_size();
	_pvalue=0.0;
	double *_stationary=_M->get_stationary();
	for (int i=0; i<imax; i++) {
	  _pvalue+=_stationary[i]*_sum[i];
	}
	_stat=-log(_pvalue)/log(10.0);
      }
      /* end of over case */
    } else {
      /* under case */
      // initialization
      for (long i=0; i<_dim; i++) {
	_current_u[i]=1.0;
      }
      double *aux=NULL;
      double *pos_in_current,*pos_in_last1,*pos_in_last2;
      //double last=1.0;
      // main loop
      for (long i=2; i<=_n; i++) {
	// switch _last_u and _current_u;
	aux=_last_u;
	_last_u=_current_u;
	_current_u=aux;
	// put fmci transition times _last_u in _current_u
	pos_in_current=_current_u;
	pos_in_last1=_last_u;
	// block 0
	// _current_u_0 = P * _last_u_0
	//_T->P_times(pos_in_last1,pos_in_current);		
	_T->P_times(pos_in_last1,pos_in_current);		
	// loop on blocks
	for (long j=1; j<_nblock; j++) {
	  // new positions
	  pos_in_current+=_order;
	  pos_in_last2=pos_in_last1;
	  pos_in_last1+=_order;
	  // _current_u_j = P * _last_u_j + Q * _last_u_j-1
	  //_T->P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	  _T->P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	} // end blocks loop
	//printf("i=%i\tu[0]=%e\n",i,pos_in_current[0]);
	//{
	//  _pvalue=0.0;
	//  double *_stationary=_M->get_stationary();
	//  for (int ii=0; ii<_alpha->get_size(); ii++) {
	//    _pvalue+=_stationary[ii]*pos_in_current[ii];
	//  }
	//  printf("i=%i\tpvalue=%e\tdiff=%e\n",i,_pvalue,last-_pvalue);
	//  last=_pvalue;
	//}

      } // end main loop
      // get the final result
      {
	//printf("stationary size = %i\n",_M->get_stationary_size());
	int imax=_M->get_stationary_size();
	_pvalue=0.0;
	double *_stationary=_M->get_stationary();
	for (int i=0; i<imax; i++) {
	  _pvalue+=_stationary[i]*pos_in_current[i];
	}
	_stat=log(_pvalue)/log(10.0);
      }
      /* end of under case */
    }

  };

  void newxstat::log_comp() {
    _order=_T->get_order();
    if (_stat>0) {
      _nblock=_nobs;
    } else {
      _nblock=_nobs+1;
    }
    _dim=_order*_nblock;
    /* check memory requirements */
    {
      unsigned long new_workspace_size=2*_dim+_order;
      if (new_workspace_size>_workspace_size) {
	/* realloc space */
	if (_workspace!=NULL) {
	  free(_workspace);	  
	}
	_workspace_size=new_workspace_size;
	_workspace=(double *)malloc(sizeof(double)*_workspace_size);
	if (_workspace==NULL) {
	  fprintf(stderr,"not enough memory for workspace allocation in xstat\n");
	  exit(EXIT_FAILURE);
	}
      }
    } /* end memory requirements */
    /* connect stuff */
    {
      double *dpos=_workspace;
      _current_u=dpos;
      dpos+=_dim;
      _last_u=dpos;
      dpos+=_dim;
      _sum=dpos;
    }
    /* two cases: over and under */
    if (_stat>0) {
      /* over case */
      // initialization
      for (long i=0; i<_dim; i++) {
	_current_u[i]=0.0;
      }
      _T->sum_Q(_current_u);
      // _current_u = log (_current_u )
      for (long i=0; i<_dim; i++) {
	_current_u[i]=log(_current_u[i]);
      }
      for (long i=0; i<_order; i++) {
	_sum[i]=log(0.0);
      }
      double *aux=NULL;
      double *pos_in_current,*pos_in_last1,*pos_in_last2;
      // main loop
      for (long i=2; i<=_n-1; i++) {
	// switch _last_u and _current_u;
	aux=_last_u;
	_last_u=_current_u;
	_current_u=aux;
	// put fmci transition times _last_u in _current_u
	pos_in_current=_current_u;
	pos_in_last1=_last_u;
	// block 0
	// _current_u_0 = P * _last_u_0
	_T->log_P_times(pos_in_last1,pos_in_current);		
	// loop on blocks
	for (long j=1; j<_nblock; j++) {
	  // new positions
	  pos_in_current+=_order;
	  pos_in_last2=pos_in_last1;
	  pos_in_last1+=_order;
	  // _current_u_j = P * _last_u_j + Q * _last_u_j-1
	  _T->log_P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	} // end blocks loop
	// get largest term of _sum and add
	double z=log(0.0);
	for (long j=0; j<_order; j++) {
	  if (_sum[j]>z)
	    z=_sum[j];
	}
	for (long j=0; j<_order; j++) {
	  if (pos_in_current[j]>z)
	    z=pos_in_current[j];
	}
	if (z==log(0.0))
	  z=0;
	//printf("z=%e\n",z);
	// _sum = _sum + last block of current_u (pointed by pos_in_current)	
	for (long j=0; j<_order; j++) {
	  _sum[j]=exp(_sum[j]-z)+exp(pos_in_current[j]-z);
	}
	for (long j=0; j<_order; j++) {
	  _sum[j]=log(_sum[j])+z;
	}
	//printf("i=%i\tlog10(P|1)=%f\n",i,_sum[0]/log(10.0));
      }
      // get the final result
      {
	int imax=_M->get_stationary_size();
	_pvalue=0.0;
	// find largest term
	double z=log(0.0);
	for (int i=0; i<imax; i++) {
	  if (_sum[i]>z)
	    z=_sum[i];
	}
	// compute sum
	double *_stationary=_M->get_stationary();
	for (int i=0; i<imax; i++) {
	  _pvalue+=_stationary[i]*exp(_sum[i]-z);
	}
	_pvalue=log(_pvalue)+z;
	_stat=-_pvalue/log(10.0);
	_pvalue=exp(_pvalue);
      }
      /* end of over case */
    } else {
      /* under case */
      // initialization
      for (long i=0; i<_dim; i++) {
	_current_u[i]=0.0;
      }
      double *aux=NULL;
      double *pos_in_current,*pos_in_last1,*pos_in_last2;
      //double last=1.0;
      // main loop
      for (long i=2; i<=_n; i++) {
	// switch _last_u and _current_u;
	aux=_last_u;
	_last_u=_current_u;
	_current_u=aux;
	// put fmci transition times _last_u in _current_u
	pos_in_current=_current_u;
	pos_in_last1=_last_u;
	// block 0
	// _current_u_0 = P * _last_u_0
	//_T->P_times(pos_in_last1,pos_in_current);		
	_T->log_P_times(pos_in_last1,pos_in_current);		
	// loop on blocks
	for (long j=1; j<_nblock; j++) {
	  // new positions
	  pos_in_current+=_order;
	  pos_in_last2=pos_in_last1;
	  pos_in_last1+=_order;
	  // _current_u_j = P * _last_u_j + Q * _last_u_j-1
	  //_T->P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	  _T->log_P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
	} // end blocks loop
	//printf("i=%i\tu[0]=%e\n",i,pos_in_current[0]);
	//{
	//  _pvalue=0.0;
	//  double *_stationary=_M->get_stationary();
	//  for (int ii=0; ii<_alpha->get_size(); ii++) {
	//    _pvalue+=_stationary[ii]*pos_in_current[ii];
	//  }
	//  printf("i=%i\tpvalue=%e\tdiff=%e\n",i,_pvalue,last-_pvalue);
	//  last=_pvalue;
	//}

      } // end main loop
      // get the final result
      {
	int imax=_M->get_stationary_size();
	_pvalue=0.0;
	// find largest term
	double z=log(0.0);
	for (int i=0; i<imax; i++) {
	  if (pos_in_current[i]>z)
	    z=pos_in_current[i];
	}
	// compute sum
	double *_stationary=_M->get_stationary();
	for (int i=0; i<imax; i++) {
	  _pvalue+=_stationary[i]*exp(pos_in_current[i]-z);
	}
	_pvalue=log(_pvalue)+z;
	_stat=_pvalue/log(10.0);
	_pvalue=exp(_pvalue);
      }
      /* end of under case */
    }

  };

  void newxstat::print_format(FILE *fout) {
    fprintf(fout,"pattern\tnobs\tnatt\tstat\n");
  };

  void newxstat::print_regular(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%f\n",_label,_nobs,_natt,_stat);
  };

  void newxstat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%e\n",_label,_nobs,_natt,_stat);
  }

};
