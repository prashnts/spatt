/* $Id: markov.cc 864 2006-09-13 13:21:10Z gnuel $ */
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
/*  Class count: counting occurrences of words                       */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/


#include "markov.h"

#include <cstdio>
#include <cmath>

using namespace std;

namespace spatt {

/* markov from a markovfile */
  markov::markov(const spattparameters &params,
	       const alphabet &alpha) : 
  _params(&params), _alpha(&alpha), _model(NULL), 
  _augmented_model(NULL), _stationary(NULL) 
  {
  
  if (_params->use_markov_file()==0) {
    fprintf(stderr,"markov constructor from file invoqued without file\n");
    exit(EXIT_FAILURE);
  }


  _model=NULL;
  _stationary=NULL;
  _secondmag=0.0;
  _augmented_model=NULL;

  _filename[0]='\0';

  _order=_params->get_markov_order();

  _k=_alpha->get_size();
  _dim=(long)pow((float)_k,_order-1);
  _n=_dim*_k;
  if (_order<=0)
    _n=_k;

  sprintf(_filename,_params->get_markov_filename());

  if (_order==-1) {
    fprintf(stderr,"Markov order -1 not usable with a Markov file\n");
    exit(EXIT_FAILURE);
  }
  
  /* opens stream */
  _stream=fopen(_filename,"r");
  if (_stream==NULL) {
    fprintf(stderr,"Cannot open markov file \"%s\"\n",_filename);
    exit(EXIT_FAILURE);
  }

  //printf("%s successfully opened in read mode\n",_filename);

  if (_order==0) {
    printf("Markov order 0 from a Markov file not yet implemented\n");
    exit(EXIT_FAILURE);
  } else {
    char buffer[MAX_STRING_LENGTH];
    char *pointfield;
    long line=0;
    long col=0;
    double sum=0;

    // allocating _model
    {
      _model=(double **)malloc(sizeof(double)*_n);
      if (_model==NULL) {
	fprintf(stderr,"Failure in memory alloc in markov::markov\n");
	exit(EXIT_FAILURE);
      }
      _model[0]=(double *)malloc(sizeof(double)*_n*_k);
      if (_model[0]==NULL) {
	fprintf(stderr,"Failure in memory alloc in markov::markov\n");
	exit(EXIT_FAILURE);
      }
      double *pos=_model[0];
      for (long i=0; i<_n; i++) {
	_model[i]=pos;
	pos+=_k;
      }    
    } // end alloc block

    // reading parameters from the file   
    col=0;
    sum=0;
    while (fgets(buffer,MAX_STRING_LENGTH,_stream)!=NULL) {
      if (buffer[0]!=COMMENT_CHAR) {
	// splitting line collecting parameters
	pointfield = (char *)strtok(buffer, " \t\n,:;") ;
	while (pointfield!=NULL)
	  {
	    if (col>=_k) {
	      fprintf(stderr,"Error reading Markov file line %li: too many columns\n",line+1);
	      exit(EXIT_FAILURE);
	    }
	    //printf("line=%i col=%i reading=%s\n",line,col,pointfield);
	    if (line>=_n) {
	      fprintf(stderr,"Error reading Markov file: too many parameters lines\n");
	      exit(EXIT_FAILURE);
	    }
	    _model[line][col]=atof(pointfield);
	    sum+=_model[line][col];
	    col++ ;
	    if (col>=_k) {
	      col=0;
	      // normalizing parameters
	      if (sum!=0.0) {
		for (long i=0; i<_k; i++)
		  _model[line][i]/=sum;
	      }
	      line++;
	      sum=0.0;
	    }
	    pointfield = (char *)strtok(NULL, " \t\n,:;") ;
	  }
      } // end if buffer[0]!=COMMENT_CHAR
    } // end while fgets
    
    // testing number of lines treated
    if (line!=_n || col!=0) {
      fprintf(stderr,"Error reading Markov file: too few parameters\n");
      exit(EXIT_FAILURE);
    }

    // computing stationary distribution
    //compute_stationary();
    compute_stationary();

  }

  /* output model if necessary */
  if (_params->use_model_file()==1) {
    dump_model(_params->get_model_filename());
  }
  
  //printf("markov from a file not yet implemented\n");
  //exit(EXIT_FAILURE);

  //printf("_model=%p\n",_model);
  //printf("_stationary=%p\n",_stationary);
  

  fclose(_stream);
};


/* compute augmented model */
void markov::compute_augmented(int augmented_order){

  int done=0;

  /* checking order value */
  if (_order>augmented_order) {
    fprintf(stderr,"Critical error in markov::compute_augmented\n");
    fprintf(stderr,"augmented order supplied is incompatible with markov order\n");
    exit(EXIT_FAILURE);
  }

  /* checking presence of an anterior augmented model */
  if (_augmented_model!=NULL) {
    if (_augmented_order==augmented_order) {
      /* everything is done */
      done=1;
    } else {
      /* free _augmented_model */
      if (_augmented_model[0]!=NULL)
	free(_augmented_model[0]);
      free(_augmented_model);
    }
  }

  if (done==0) {
    _augmented_order=augmented_order;
    _DIM=(long)pow((float)_k,_augmented_order);
    /* allocating augmented_model */
    _augmented_model=(double **)malloc(sizeof(double)*_DIM);
    if (_augmented_model==NULL) {
      fprintf(stderr,"Failure in memory alloc in markov::operator=\n");
      exit(EXIT_FAILURE);
    }
    _augmented_model[0]=(double *)malloc(sizeof(double)*_DIM*_k);
    if (_augmented_model[0]==NULL) {
      fprintf(stderr,"Failure in memory alloc in markov::operator=\n");
      exit(EXIT_FAILURE);
    }
    double *pos=_augmented_model[0];
    for (long i=0; i<_DIM; i++) {
      _augmented_model[i]=pos;
      pos+=_k;
    }    
    /* filling augmented_model with model */
    /* easy, fill _augmented_model[0] with repeats */
    /* of _model_[0] or _stationary if _order<=0 */
    if (_order<=0) {
      long jmax=_DIM;
      long shift=0;
      for (long j=0; j<jmax; j++) {
	for (long i=0; i<_k; i++)
	  _augmented_model[0][shift+i]=_stationary[i];
	shift+=_k;
      }
    } else {
      long jmax=_DIM/_n;
      long shift=0;
      for (long j=0; j<jmax; j++) {
	for (long i=0; i<_n*_k; i++)
	  _augmented_model[0][shift+i]=_model[0][i];
	shift+=_n*_k;
      }
    }
    
  }
  

};

/* markov from a sequence */
  markov::markov(const spattparameters &params,
		 const alphabet &alpha,
		 const count &occ) :
    _params(&params), _alpha(&alpha), _occ(&occ),
    _model(NULL), _augmented_model(NULL), _stationary(NULL)
  {
    
    _seq=_occ->get_seq();

  _secondmag=0.0;
  _order=_params->get_markov_order();

  _filename[0]='\0';

  _k=_alpha->get_size();
  _dim=(long)pow((float)_k,_order-1);
  _n=_dim*_k;
  if (_order<=0)
    _n=_k;
  
  if (_occ->get_length() < _order +1 ) {
    fprintf(stderr,"Critical error in markov:markov ! length >= order + 1 is needed for Markov estimation.\n");
    exit(EXIT_FAILURE);
  }
  
  if (_order>-2) {

    if (_order<=0) {
      //fprintf(stderr,"Markov model with no memory (order 0 and -1) not yet implemented\n");
      //exit(EXIT_FAILURE);
      /* allocation */
      _stationary=(double *)malloc(sizeof(double)*_k);
      if (_stationary==NULL) {
	fprintf(stderr,"Critical error in markov:markov ! Not enough memory.\n");      
	exit(EXIT_FAILURE);
      }
      /* estimation */
      if (_order==-1) {
	for (int i=0; i<_k; i++)
	  _stationary[i]=1/(double)_k;
      } else {
	double sum=0;
	for (int i=0; i<_k; i++) {
	  _stationary[i]=_occ->get(1,i);
	  sum+=_stationary[i];
	}
	if (sum!=0.0) {
	  for (int i=0; i<_k; i++) 
	    _stationary[i]/=sum;
	}
      }
      /* verification */
      //for (int i=0; i<_k; i++)
      //  printf("_s[%i]=%f ",i,_stationary[i]);
      //printf("\n");
      if (_params->use_stationary_file()==1) {
	dump_stationary(_params->get_stationary_filename());
      }
    }// end if (_order <=0 )
    else {
          
    //printf("markov: _k=%i _order=%i _dim=%i _n=%i\n",_k,_order,_dim,_n);
    /* _model memory allocation */
    {
      double *pos;
      _model=(double **)malloc(sizeof(double *)*_n);
      if (_model==NULL) {
	fprintf(stderr,"Critical error in markov:markov ! Not enough memory.\n");
	exit(EXIT_FAILURE);
      }
      _model[0]=(double *)malloc(sizeof(double)*_n*_k);
      if (_model[0]==NULL) {
	fprintf(stderr,"Critical error in markov:markov ! Not enough memory.\n");
	exit(EXIT_FAILURE);
      }    
      pos=_model[0];
      for (long i=0; i<_n; i++) {
	_model[i]=pos;
	pos+=_k;
      }    
    }
    /* _model estimation */
    {
      long code=0;
      double sum=0;
      for (long i=0; i<_n; i++){
	sum=0;
	for (int j=0; j<_k; j++) {
	  _model[i][j]=_occ->get(_order+1,code);
	  sum+=_model[i][j];
	  code++;
	}
	if (sum>0) {
	  for (int j=0; j<_k; j++)
	    _model[i][j]/=sum;
	} else {
	  /* do nothing */
	}
      }
    }
    /* output _model */
    //for (long i=0; i<_n; i++){
    //  for (int j=0; j<_k; j++) 
    //    printf("%f ",_model[i][j]);
    //  printf("\n");
    //}


    /* compute stationary */
    //compute_stationary();
    compute_stationary();

    /* output _stationary */
    //for (long i=0; i<_n; i++)
    //  printf("%f ",_stationary[i]);
    //printf("\n");
    }

    /* output model if necessary */
    if (_params->use_model_file()==1) {
      dump_model(_params->get_model_filename());
    }
    
  } // end if (_order>-2)
};

void markov::compute_stationary(){
  
  double **V;
  double **H;
  double *f;
  double *aux;

  if (_model!=NULL) {
    
    /* memory allocation */

    /* memory allocation of f */
    f=(double *)malloc(sizeof(double)*_n);
    //printf("compute_stationary f=%p\n",f);
    if (f==NULL) {
      fprintf(stderr,"Critical error in markov::compute_stationary ! Not enough memory.\n");
      exit(EXIT_FAILURE);
    }

    //printf("n=%i limit=%i\n",_n,MAX_FULL_DIM);
    if (_n<=MAX_FULL_DIM) {
      /* full case */
      //printf("full\n");
      zgeev_par zpar(_n,1,0);
      //char jobvl='v'; /* left eigenvector needed to be computed */
      //char jobvr='n'; /* but not right ones */
      //int n=_n;
      complex16 a[MAX_FULL_DIM*MAX_FULL_DIM]; /* dim n x n */
      //complex16 tmpwork[MAX_FULL_DIM*MAX_FULL_DIM];
      //int lda=_n;
      //complex16 w[MAX_FULL_DIM*MAX_FULL_DIM]; /* dim n */
      //complex16 vl[MAX_FULL_DIM*MAX_FULL_DIM]; /* dim n x n */
      //int ldvl=_n;
      //complex16 vr[MAX_FULL_DIM*MAX_FULL_DIM]; /* dim n x n */
      //int ldvr=_n;
      //complex16 tempwork;
      //complex16 *work=NULL;
      //int lwork=-1;
      //double rwork[2*MAX_FULL_DIM]; /* dim 2n */
      //int info;

      /* fill matrix a */
      {
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
	      a[pos].re=_model[mpos][jrot];
	      a[pos].im=0.0;
	      shift+=_dim;
	    } else {
	      a[pos].re=0.0;
	      a[pos].im=0.0;
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
      } /* end fill matrix a */

      // verif
      //for (int i=0; i<_n; i++) {
      //  for (int j=0; j<_n; j++) {
      //    printf("%f ",a[j*_n+i]);
      //  }
      //  printf("\n");
      //}

      /* calls zgeev for a dimension query */
      //zgeev_(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vr,&ldvr,work,&lwork,rwork,&info);
      //zgeev_(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,&tempwork,lwork,rwork,info);
      //zpar.print();
      zgeev_(zpar._jobvl,zpar._jobvr,zpar._n,a,zpar._lda,zpar._w,zpar._vl,zpar._ldvl,zpar._vr,zpar._ldvr,zpar._work,zpar._lwork,zpar._rwork,zpar._info);
      int lwork=(int)zpar._work[0].re;
//      zgeev_(zpar._jobvl,zpar._jobvr,zpar._n,a,zpar._lda,zpar._w,zpar._vl,zpar._ldvl,zpar._vr,zpar._ldvr,tmpwork,zpar._lwork,zpar._rwork,zpar._info);
//      int lwork=(int)tmpwork[0].re;
      //printf("after work space query, info=%i and workspace=%i\n",zpar._info,lwork);
      //work=(complex16 *)malloc(sizeof(complex16)*lwork);
      //if (work==NULL) {
      //  fprintf(stderr,"Not enough memory to allocate workspace (%i octets required)\n",lwork);
      //  exit(EXIT_FAILURE);
      //}
      zpar.reset(_n,1,0,lwork);
      //zpar.print();
      zgeev_(zpar._jobvl,zpar._jobvr,zpar._n,a,zpar._lda,zpar._w,zpar._vl,zpar._ldvl,zpar._vr,zpar._ldvr,zpar._work,zpar._lwork,zpar._rwork,zpar._info);
      //      zgeev_(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info);
      //free(work);
      //printf("second call to zgeev returned  info=%i\n",zpar._info);

      /* looking for the position of eigenvalue 1.0 and for the second magnitude*/
      {
	complex16 *w=zpar._w;
	int firstpos=-1;
	int secondpos=-1;
	double secondmag=0.0;
	double currentmag=0.0;
	
	for (int i=0; i<_n; i++) {
	  currentmag=sqrt(w[i].re*w[i].re+w[i].im*w[i].im);
	  //printf("vp[%i]-> re=%f\tim=%f\tmagnitude=%f\n",i,w[i].re,w[i].im,currentmag);
	  if (fabs(currentmag-1.0)<1e-10) {
	    if (fabs(w[i].im)<1e-10) {
	      firstpos=i;
	    } else {
	      if (currentmag>secondmag) {
		secondmag=currentmag;
		secondpos=i;
	      }
	    }	      
	  } else {
	    if (currentmag>secondmag) {
	      secondmag=currentmag;
	      secondpos=i;
	    }	    
	  }
	}
	
	//printf("firstpos=%i\tsecondpos=%i\tsecondmag=%f\n",firstpos,secondpos,secondmag);

	{
	  complex16 *vl=zpar._vl;
	  complex16 tempwork;
	  for (int i=0; i<_n; i++) {
	    tempwork=vl[i+firstpos*_n];
	    //printf("mu[%i]-> re=%f\tim=%f\n",i,tempwork.re,tempwork.im);
	    if (fabs(tempwork.im)>1e-10) {
	      fprintf(stderr,"Numerical inconsistence while computing the stationary distribution: %e != 0.0\n",fabs(tempwork.im));
	      exit(EXIT_FAILURE);
	    }
	    f[i]=tempwork.re;
	  }
	}
	_secondmag=secondmag;

      } /* end searching block */
    } else {
      /* sparse case */
      //printf("sparse\n");
      /* create a dnaupd_par */
      dnaupd_par dpar(_n,DNAUPD_KRYLOVDIM);
      dpar._nev=2; // we need 2 eigenvalues
      //dpar.print();
      
    main:
      dnaupd_(dpar._ido,dpar._bmat,dpar._n,dpar._which,dpar._nev,dpar._tol,dpar._resid,dpar._ncv,dpar._v,dpar._ldv,dpar._iparam,dpar._ipntr,dpar._workd,dpar._workl,dpar._lworkl,dpar._info);
      //printf("dnaupd_call done\n");
      if (dpar._ido==-1 || dpar._ido==1) {
	/* compute product Y=OP*X where */
	/* ipntr[0]-1 is the pointer into workd for X */
	/* ipntr[1]-1 is the pointer into workd for Y */ 
	//printf("matrix x vector product\n");
	double *x=&dpar._workd[dpar._ipntr[0]-1];
	double *y=&dpar._workd[dpar._ipntr[1]-1];
	int shift=0;
	int jrot=0;
	for (long j=0; j<_n; j++) {
	  /* compute y[j] */
	  //printf("y[%i]=",j);
	  int mpos=0;
	  y[j]=0.0;
	  for (long i=shift; i<_n; i+=_dim) {
	    //printf("x[%i]*_model[%i][%i]+",i,i,jrot);
	    y[j]+=x[i]*_model[i][jrot];
	  }
	  //printf("\n");
	  jrot++;
	  if (jrot==_k) {
	    jrot=0;
	    shift++;
	  }
	}
	goto main; 
      } /* end if ido */

      //printf("implicitly restarted algorithm returned:\n");
      //printf("\tinfo=%i\n",dpar._info);
      //printf("\tnumber of iter=%i\n",dpar._iparam[2]);
      //printf("\tnumber of converged ritzvalues=%i\n",dpar._iparam[4]);
      //printf("\tnumber of prod=%i\n",dpar._iparam[8]);
      //printf("\tnumber of re-orth steps=%i\n",dpar._iparam[10]);    
      {
	//printf("ipntr:");
	//for (int i=0; i<14; i++)
	//  printf("[%i]:%i\t",i+1,dpar._ipntr[i]);
	//printf("\n");
	double *re,*im;
	re=&dpar._workl[dpar._ipntr[5]-1];
	im=&dpar._workl[dpar._ipntr[6]-1];
	//printf("eigenvalues:\n");
	//for (int i=0; i<dpar._ncv; i++) {
	//  printf("\tre=%f\tim=%f\n",re[i],im[i]);
	//}
	
	/* call of dneupd */
	//printf("calling dneupd_par\n");
	dneupd_par ddpar(dpar._n,dpar._nev,DNEUPD_GET_VECTORS);
	//ddpar.print();
	//dpar.print();
	//printf("ipntr:");
	//for (int i=0; i<14; i++)
	//  printf("[%i]:%i\t",i+1,dpar._ipntr[i]);
	//printf("\n");
	dneupd_(ddpar._rvec,ddpar._howmny,ddpar._select,ddpar._dr,ddpar._di,ddpar._z,ddpar._ldz,ddpar._sigmar,ddpar._sigmai,ddpar._workev,
		dpar._bmat,dpar._n,dpar._which,dpar._nev,dpar._tol,dpar._resid,dpar._ncv,dpar._v,dpar._ldv,dpar._iparam,dpar._ipntr,dpar._workd,dpar._workl,dpar._lworkl,dpar._info);      

	/* print results */
	//printf("eigenvalues:\n");
	//for (int i=0; i<ddpar._nev; i++) {
	//  printf("\tre=%f\tim=%f\n",ddpar._dr[i],ddpar._di[i]);
	//}
	_secondmag=sqrt(ddpar._dr[1]*ddpar._dr[1]+ddpar._di[1]*ddpar._di[1]);
	//printf("eigenvector:\n");      
	for (int i=0; i<ddpar._n; i++) {
	  //printf("\t%f\n",ddpar._z[i]);
	  f[i]=ddpar._z[i];
	}
      }      

    } /* end sparse case */

    /* get stationary distrib */
    _stationary=f;

    // verification
    //printf("stationary:\n");
    //for (long i=0; i<_n; i++)
    //  printf("%f\t",_stationary[i]);
    //printf("\n");
    
    /* normalize */
    {
      double sum=0;
      for (long i=0; i<_n; i++)
	sum+=_stationary[i];
      for (long i=0; i<_n; i++)
	_stationary[i]/=sum;      
    }

    /* output stationary if necessary */
    if (_params->use_stationary_file()==1) {
      dump_stationary(_params->get_stationary_filename());
    }

    /* free unused memory */
  } /* end if (_model!=NULL) */
};


double markov::expect(word &w) const{

  if (_order<=0) {
    int c;
    /* very simple */
    //printf("expect(%s)\n",w.get_label());
    w.set_expected(1.0);
    for (int i=0; i<w.get_size(); i++) {
      c=_alpha->char2code(w.get_label()[i]);
      //printf("%c[%i]=%f ",w.get_label()[i],c,_stationary[c]);
      w.set_expected(w.get_expected()*_stationary[c]);
    }
  } else {
    /* more complicated */
    if (w.get_size()<=_order) {
      /* word is shorter than markov order */
      /* lets complete the word with all */
      /* possible end */
      w.set_expected(0);
      int g=_order-w.get_size();
      long dim_g=(long)pow((float)_alpha->get_size(),g);
      long c=w.get_code()*dim_g;
      for (int i=0; i<dim_g; i++) 
	w.set_expected(w.get_expected()+_stationary[c+i]);
    } else {
      /* the word must be parsed */
      //fprintf(stderr,"warning ! general case not yet implemented in markov:expect\n");
      /* compute code of W1 ... Wm where m=_order */
      long c=0;
      long power=1;
      for (int i=_order-1; i>=0; i--) {
	c+=_alpha->char2code(w.get_label()[i])*power;
	power*=_alpha->get_size();
      }
      //printf("code[0]=%li\n",c);
      power/=_alpha->get_size();
      /* add mu(W1 ... Wm) */
      w.set_expected(_stationary[c]);
      long cc;
      for (int i=_order; i<w.get_size(); i++) {
	/* remove first letter */
	cc=c-_alpha->char2code(w.get_label()[i-_order])*power;
	/* shift to the left */
	cc*=_alpha->get_size();
	/* add new letter */
	cc+=_alpha->char2code(w.get_label()[i]);
	//printf("code[%i]=%li\n",i-_order+1,cc);
	/* update _expected */
	w.set_expected(w.get_expected()*_model[c][_alpha->char2code(w.get_label()[i])]);
	/* switch c and cc */
	c=cc;
      }
    }
  }
  //w.set_expected(w.get_expected()*(_occ->get_n()-w.get_size()+1));
  //printf("valid_char=%i nseq=%i and n=%i\n",_seq->_valid_char,_seq->_nseq,_seq->_valid_char-_seq->_nseq*(w.get_size()-1));
  w.set_expected(w.get_expected()*(_seq->_valid_char-_seq->_nseq*(w.get_size()-1)));
  return w.get_expected();
};

double markov::mu(long code,int size) const {
  string *w=_alpha->code2string(code,size);
  double res=mu(_alpha->code2string(code,size));
  delete w;
  return res;
}

double markov::mu(string *s) const {
  return mu(s->c_str());
}

double markov::mu(const char *label) const {
  //printf("expect(%s)\n",label);
  int size=strlen(label);
  
  double expected;
  if (_order<=0) {
    int c;
    /* very simple */
    expected=1.0;
    for (int i=0; i<size; i++) {
      c=_alpha->char2code(label[i]);
      expected*=_stationary[c];
      //printf("mu(%c)=%f\n",label[i],_stationary[c]);
    }
  } else {
    /* more complicated */
    if (size<=_order) {
      /* word is shorter than markov order */
      /* lets complete the word with all */
      /* possible end */
      expected=0;
      int g=_order-size;
      long dim_g=(long)pow((float)_alpha->get_size(),g);
      long code;
      /* code word from the label */
      {
	code=0;
	int power=1;
	int c;
	for (int i=size-1; i>=0; i--) {
	  c=_alpha->char2code(label[i]);
	  if (c<0 || c>=_alpha->get_size()) {
	    fprintf(stderr,"Critical error in markov::expect(%s,%i) ! \"%s\" uses invalid char\n",label,size,label);
	    exit(EXIT_FAILURE);
	  }
	  code+=power*c;
	  power*=_alpha->get_size();
	}
      }
      long c=code*dim_g;
      for (int i=0; i<dim_g; i++) 
	expected+=_stationary[c+i];
    } else {
      /* the word must be parsed */
      /* compute code of W1 ... Wm where m=_order */
      long c=0;
      long power=1;
      for (int i=_order-1; i>=0; i--) {
	c+=_alpha->char2code(label[i])*power;
	power*=_alpha->get_size();
      }
      //printf("code[0]=%li\n",c);
      power/=_alpha->get_size();
      /* add mu(W1 ... Wm) */
      expected=_stationary[c];
      long cc;
      for (int i=_order; i<size; i++) {
	/* remove first letter */
	cc=c-_alpha->char2code(label[i-_order])*power;
	/* shift to the left */
	cc*=_alpha->get_size();
	/* add new letter */
	cc+=_alpha->char2code(label[i]);
	/* update _expected */
	expected*=_model[c][_alpha->char2code(label[i])];
	/* switch c and cc */
	c=cc;
      }
    }
  }
  return expected;
};

double markov::expect(char *label) const{

//  //printf("expect(%s)\n",label);
//  int size=strlen(label);
//
//  double expected;
//  if (_order<=0) {
//    int c;
//    /* very simple */
//    expected=1.0;
//    for (int i=0; i<size; i++) {
//      c=_alpha->char2code(label[i]);
//      expected*=_stationary[c];
//      //printf("mu(%c)=%f\n",label[i],_stationary[c]);
//    }
//    expected*=_occ->_n-size+1;
//  } else {
//    /* more complicated */
//  if (size<=_order) {
//    /* word is shorter than markov order */
//    /* lets complete the word with all */
//    /* possible end */
//    expected=0;
//    int g=_order-size;
//    long dim_g=(long)pow((float)_alpha->get_size(),g);
//    long code;
//    /* code word from the label */
//    {
//      code=0;
//      int power=1;
//      int c;
//      for (int i=size-1; i>=0; i--) {
//	c=_alpha->char2code(label[i]);
//	if (c<0 || c>=_alpha->get_size()) {
//	  fprintf(stderr,"Critical error in markov::expect(%s,%i) ! \"%s\" uses invalid char\n",label,size,label);
//	  exit(EXIT_FAILURE);
//	}
//	code+=power*c;
//	power*=_alpha->get_size();
//      }
//    }
//    long c=code*dim_g;
//    for (int i=0; i<dim_g; i++) 
//      expected+=_stationary[c+i];
//  } else {
//    /* the word must be parsed */
//    /* compute code of W1 ... Wm where m=_order */
//    long c=0;
//    long power=1;
//    for (int i=_order-1; i>=0; i--) {
//      c+=_alpha->char2code(label[i])*power;
//      power*=_alpha->get_size();
//    }
//    //printf("code[0]=%li\n",c);
//    power/=_alpha->get_size();
//    /* add mu(W1 ... Wm) */
//    expected=_stationary[c];
//    long cc;
//    for (int i=_order; i<size; i++) {
//      /* remove first letter */
//      cc=c-_alpha->char2code(label[i-_order])*power;
//      /* shift to the left */
//      cc*=_alpha->get_size();
//      /* add new letter */
//      cc+=_alpha->char2code(label[i]);
//      /* update _expected */
//      expected*=_model[c][_alpha->char2code(label[i])];
//      /* switch c and cc */
//      c=cc;
//    }
//  }
  //return (_occ->get_n()-strlen(label)+1)*mu(label);
  //printf("valid_char=%i nseq=%i and n=%i\n",_seq->_valid_char,_seq->_nseq,_seq->_valid_char-_seq->_nseq*(strlen(label)-1));
  return (_seq->_valid_char-_seq->_nseq*(strlen(label)-1))*mu(label);
  //  }
};

/* compute P(fromto|from) */
double markov::trans(const char *from, const char *to) const{

  double res=1.0;
  int lfrom=strlen(from);
  int lto=strlen(to);
  char *label=NULL;

  label=(char*)malloc(sizeof(char)*(lto+lfrom+1));
  strcpy (label,from);
  strcat (label,to);
  // verif
  //printf("from=%s|to=%s->fromto->%s",from,to,label);

  res=mu(label)/mu(from);
  free(label);
  return res;

//  //printf("trans(from=%s,to=%s)\n",from,to);
//  if (_order<=0) {
//    /* order 0 model */
//    int c;
//    for (int i=0; i<lto; i++) {
//      c=_alpha->char2code(to[i]);    
//      res*=_stationary[c];
//    }    
//  } else {
//    /* higher order models */
//    fprintf(stderr,"markov::trans not yet implemented for order > 0\n");
//    exit(EXIT_FAILURE);
//    /* get code of the last m letters of from */
//    
//  }
//  return res;
}


/* dump _model in a file */
void markov::dump_model(const char* filename){
  FILE *stream;
    
  stream=fopen(filename,"w");
  if (stream==NULL) {
    fprintf(stderr,"Cannot open %s in write mode\n",filename);
    exit(EXIT_FAILURE);
  }

  //printf("dump_model: _order=%i\n",_order);
  if ( (_model==NULL && _order>0) || (_stationary==NULL && _order<=0) ) {
    fprintf(stderr,"Warning ! dump_model called but there is nothing to dump !\n");
  } else {
    fprintf(stream,"# file produced by %s %s\n",PACKAGE,VERSION);
    fprintf(stream,"# Markov model parameters\n");
    fprintf(stream,"# alphabet = {%s}\n",_alpha->get_label().c_str());
    fprintf(stream,"# order = %i\n",_order);
    if (_order>0) {
      fprintf(stream,"# number of lines = %li\n",_n);
      fprintf(stream,"# number of parameters = %li\n",_n*_k);
    } else {
      fprintf(stream,"# number of lines = 1\n");
      fprintf(stream,"# number of parameters = %li\n",_n);
    }
    if (_order>0) {
      for (long i=0; i<_n; i++) {
	for (long j=0; j<_k; j++)
	  fprintf(stream,"%.10f\t",_model[i][j]);
	fprintf(stream,"\n");
      }
    } else {
      for (long j=0; j<_k; j++)
	fprintf(stream,"%.10f\t",_stationary[j]);
      fprintf(stream,"\n");
    }
    fprintf(stream,"# end of file\n");    
  }


  fclose(stream);
}

/* dump stationary in a file */
void markov::dump_stationary(const char* filename){
  FILE *stream;
    
  stream=fopen(filename,"w");
  if (stream==NULL) {
    fprintf(stderr,"Cannot open %s in write mode\n",filename);
    exit(EXIT_FAILURE);
  }

  if (_stationary==NULL) {
    fprintf(stderr,"Warning ! dump_stationary called but there is nothing to dump !\n");
  } else {
    fprintf(stream,"# file produced by %s %s\n",PACKAGE,VERSION);
    fprintf(stream,"# Stationary distribution\n");
    fprintf(stream,"# alphabet = {%s}\n",_alpha->get_label().c_str());
    fprintf(stream,"# order = %i\n",_order);
    fprintf(stream,"# second magnitude = %f\n",_secondmag);    
    fprintf(stream,"# number of parameters = %li\n",_n);
    for (long i=0; i<_n; i++) {
      fprintf(stream,"%.10f\n",_stationary[i]);
    }
    fprintf(stream,"# end of file\n");    
  }


  fclose(stream);
}


/* empty constructor */
  markov::markov() : 
    _params(NULL), _alpha(NULL), _seq(NULL), _occ(NULL),
    _stream(NULL), _model(NULL), _augmented_model(NULL),
    _order(-2), _stationary(NULL), _n(-1), _dim(-1), _k(-1) 
  {
    
  };

/* affectation operator */
markov & markov::operator=(const markov &source){
  //fprintf(stderr,"affectation operator call for object markov\n");
  //exit(EXIT_FAILURE);    
  // simple copy
  if (this != &source) {
    _params=source._params;
    _alpha=source._alpha;
    _seq=source._seq;
    _occ=source._occ;
    _stream=source._stream;
    _order=source._order;
    _n=source._n;
    _dim=source._dim;
    _k=source._k;
    sprintf(_filename,source._filename);
    
    // copy with alloc : _stationary
    if (source._stationary!=NULL) {
      _stationary=(double *)malloc(sizeof(double)*_n);
      if (_stationary==NULL) {
	fprintf(stderr,"Failure in memory alloc in markov::operator=\n");
	exit(EXIT_FAILURE);
      }
      for (long i=0; i<_n; i++)
	_stationary[i]=source._stationary[i];
    }
    
    // copy with alloc : _model
    if (source._model!=NULL) {
      _model=(double **)malloc(sizeof(double)*_n);
      if (_model==NULL) {
	fprintf(stderr,"Failure in memory alloc in markov::operator=\n");
	exit(EXIT_FAILURE);
      }
      _model[0]=(double *)malloc(sizeof(double)*_n*_k);
      if (_model[0]==NULL) {
	fprintf(stderr,"Failure in memory alloc in markov::operator=\n");
	exit(EXIT_FAILURE);
      }
      double *pos=_model[0];
      for (long i=0; i<_n; i++) {
	_model[i]=pos;
	pos+=_k;
      }    
      for (long i=0; i<_n*_k; i++)
	_model[0][i]=source._model[0][i];
    }
    
    // _augmented_model not copied
    _augmented_model=NULL;
  } 
  return *this;

};

/* close stream */
markov::~markov(){
  if (_stationary)
    free(_stationary);
  if (_model) {
    if (_model[0])
      free(_model[0]);
    free(_model);
  }
  if (_augmented_model) {
    if (_augmented_model[0])
      free(_augmented_model[0]);
    free(_augmented_model);
  }
};

};
