#include "amarkov.h"

using namespace std;
using namespace spatt;

amarkov::amarkov(const char *sfile,const char *mfile,bool verbose,bool debug) {

  // read mfile first
  if (verbose)
    printf("reading Markov transition from file \"%s\"\n",mfile);
  // open file
  FILE *input=fopen(mfile,"r");
  if (input==NULL) {
    fprintf(stderr,"can not read \"%s\". Aborting.\n",mfile);
    exit(EXIT_FAILURE);
  }
  char buffer[ABUFFER_SIZE];
  char *current=NULL;
  // get first line
  while ( fgets(buffer,ABUFFER_SIZE,input) && buffer[0]==ACOMMENT_CHAR ) {}
  // treat first line
  current=strtok(buffer,DELIM);
  if (current==NULL) {
    fprintf(stderr,"empty lines forbidden in file \"%s\". Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _k=atoi(current);
  if (_k<=0) {
    fprintf(stderr,"negative alphabet size forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  current=strtok(NULL,DELIM);
  if (current==NULL) {
    fprintf(stderr,"Markov order is missing. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _m=atoi(current);
  if (_m<0) {
    fprintf(stderr,"negative Markov order forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  // end first line treatment
  // compute number of starting states
  _nstart=(int)pow((double)_k,(double)_m);
  // allocate memory for start index
  _start=(int*)malloc(sizeof(int)*_nstart);
  if (_start==NULL) {
    fprintf(stderr,"not enough memory to allocate sindex. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  // read _start from file
  {
    int i=0; 
    while ( fgets(buffer,ABUFFER_SIZE,input) &&  i<_nstart ) {
      if (buffer[0]!=ACOMMENT_CHAR) {
	_start[i]=atoi(buffer);
	i++;
      } // end if buffer[0]
    } // end while
    if (i!=_nstart) {
      fprintf(stderr,"wrong number of starting states in \"%s\". Aborting.\n",mfile);
      exit(EXIT_FAILURE);
    }
  } // end _start reading
  // get the next line
  while ( fgets(buffer,ABUFFER_SIZE,input) && buffer[0]==ACOMMENT_CHAR ) {}
  // read the number of final states
  _nfinal=atoi(buffer);
  if (verbose) {
    printf("k=%i\tm=%i\tnstart=%i\tnfinal=%i\n",_k,_m,_nstart,_nfinal);
  }
  // allocate _final and _h
  _final=(int *)malloc(sizeof(int)*2*_nfinal);
  if (_final==NULL) {
    fprintf(stderr,"not enough memory to allocate _final. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _h=_final+_nfinal;
  { // read _final from file
    int i=0; 
    while ( fgets(buffer,ABUFFER_SIZE,input) &&  i<_nfinal ) {
      if (buffer[0]!=ACOMMENT_CHAR) {
	_final[i]=atoi(buffer);
	i++;
      } // end if buffer[0]
    } // end while
    if (i!=_nfinal) {
      fprintf(stderr,"wrong number of final states in \"%s\". Aborting.\n",mfile);
      exit(EXIT_FAILURE);
    }
  } // end _final reading
  // get the next line
  while ( fgets(buffer,ABUFFER_SIZE,input) && buffer[0]==ACOMMENT_CHAR ) {}
  // treat it
  current=strtok(buffer,DELIM);
  if (current==NULL) {
    fprintf(stderr,"empty lines forbidden in file \"%s\". Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _nstate=atoi(current);
  if (_nstate<0) {
    fprintf(stderr,"negative number of states forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  current=strtok(NULL,DELIM);
  if (current==NULL) {
    fprintf(stderr,"missing parameter. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _nz=atoi(current);
  if (_nz<=0) {
    fprintf(stderr,"negative number of non zero terms forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  current=strtok(NULL,DELIM);
  if (current==NULL) {
    fprintf(stderr,"missing parameter. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _Pnz=atoi(current);
  if (_Pnz<0) {
    fprintf(stderr,"negative number of non zero terms forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  current=strtok(NULL,DELIM);
  if (current==NULL) {
    fprintf(stderr,"missing parameter. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _Qnz=atoi(current);
  if (_Qnz<0) {
    fprintf(stderr,"negative number of non zero terms forbidden. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  // allocate main memory segment
  _data_size=_nz*(2*sizeof(int)+sizeof(double));
  _data=(int*)malloc(_data_size);
  if (_data==NULL) {
    fprintf(stderr,"not enough memory to allocate data. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  // connecting pointers on _data
  {
    { // main part
      int *pos=_data;
      _from=pos;
      pos+=_nz;
      _to=pos;
      pos+=_nz;
      _proba=(double *)pos;      
    } // end main part
    { // P and Q part
      int *ipos=_from;
      _Pfrom=ipos;
      ipos+=_Pnz;
      _Qfrom=ipos;
      ipos=_to;
      _Pto=ipos;
      ipos+=_Pnz;
      _Qto=ipos;
      double *dpos=_proba;
      _Pproba=dpos;
      dpos+=_Pnz;
      _Qproba=dpos;
    } // end P and Q part
  } // end pointer connnection
  { // read data from file
    int i=0; 
    while ( fgets(buffer,ABUFFER_SIZE,input) &&  i<_nz ) {
      if (buffer[0]!=ACOMMENT_CHAR) {
	//printf("%s",buffer);
	// parse line
	current=strtok(buffer,DELIM);
	if (current==NULL) {
	  fprintf(stderr,"empty lines forbidden in file \"%s\". Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	_from[i]=atoi(current);
	if (_from[i]<0) {
	  fprintf(stderr,"wrong index. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	current=strtok(NULL,DELIM);
	if (current==NULL) {
	  fprintf(stderr,"missing parameter. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	_to[i]=atoi(current);
	if (_to[i]<0) {
	  fprintf(stderr,"wrong index. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	current=strtok(NULL,DELIM);
	if (current==NULL) {
	  fprintf(stderr,"missing parameter. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	_proba[i]=atof(current);
	if (_proba[i]<0 || _proba[i]>1) {
	  fprintf(stderr,"wrong proba. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	//printf("Pi(%i,%i)=%.10e\n",_from[i],_to[i],_proba[i]);
	i++;
      } // end if buffer[0]
    } // end while
    if (i!=_nz) {
      fprintf(stderr,"wrong number of transitions in \"%s\". Aborting.\n",mfile);
      exit(EXIT_FAILURE);
    }
  }  // end read data
  // end line treatment
  if (verbose) {
    printf("nstate=%i\tnz=%i\tPnz=%i\tQnz=%i\n",_nstate,_nz,_Pnz,_Qnz);
    if (debug) {
      printf("P:\n");
      for (int i=0; i<_Pnz; i++) {
	printf(" P(%i,%i)=%.10e\n",_Pfrom[i]+1,_Pto[i]+1,_Pproba[i]);
      }
      printf("Q:\n");
      for (int i=0; i<_Qnz; i++) {
	printf(" Q(%i,%i)=%.10e\n",_Qfrom[i]+1,_Qto[i]+1,_Qproba[i]);
      }
    }
  }
  // closing file
  fclose(input);
  // end mfile part

  // allocate _stationary
  _stationary=(double *)malloc(sizeof(double)*_nstate);
  if (_stationary==NULL) {
    fprintf(stderr,"not enough memory to allocate _stationary. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  
  // allocate _dstart
  _dstart=(double *)malloc(sizeof(double)*_nstart);
  if (_dstart==NULL) {
    fprintf(stderr,"not enough memory to allocate dstart. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  // read sfile if provided
  if ( strlen(sfile)>0 ) {
    if (verbose)
      printf("reading starting distribution from file \"%s\"\n",sfile);
    // open file
    input=fopen(sfile,"r");
    if (input==NULL) {
      fprintf(stderr,"can not read file \"%s\". Aborting.\n",sfile);
      exit(EXIT_FAILURE);
    }
    int i=0;
    while (  fgets(buffer,ABUFFER_SIZE,input) &&  i<_nstart ) {
      if (buffer[0]!=ACOMMENT_CHAR) {
	_dstart[i]=atof(buffer);
	if (_dstart[i]<0 || _dstart[i]>1) {
	  fprintf(stderr,"wrong transition in \"%s\". Aborting.\n",sfile);
	  exit(EXIT_FAILURE);
	}
	i++;
      } // end if buffer[0]
    } // end while
    // close file
    fclose(input);
  } else {
    if (verbose)
      printf("using uniform starting distribution\n");
    for (int i=0; i<_nstart; i++)
      _dstart[i]=1.0/(double)_nstart;
  } // end sfile part

  normalize();

  compute_h();

  compute_stationary();

  if (debug) {
    printf("stationary:\n");
    for (int i=0; i<_nstate; i++)
      printf("S(%i)=%e\n",i,_stationary[i]);    
  }

};

void amarkov::normalize() {

  { // _dstart normalization
    double sum=0.0;
    for (int i=0; i<_nstart; i++) 
      sum+=_dstart[i];
    for (int i=0; i<_nstart; i++) 
      _dstart[i]/=sum;    
  } // end _dstart normalization
  { // transition normalization
    // initialization
    for (int i=0; i<_nstate; i++)
      _stationary[i]=0.0;
    // store  Pi * 1 in _stationary
    for (int i=0; i<_nz; i++) {
      //printf("Pi(%i,%i)=%.10e\n",_from[i],_to[i],_proba[i]);      
      _stationary[_from[i]]+=_proba[i];
      //printf("_stationary[%i]=%f\n",_from[i],_stationary[_from[i]]);
    }
    // normalization
    for (int i=0; i<_nz; i++) {
      //printf("from=%i\tto=%i\tproba=%e\n",_Pfrom[i],_Pto[i],_Pproba[i]);
      //printf("_stationary[%i]=%f\n",_from[i],_stationary[_from[i]]);
      _proba[i]/=_stationary[_from[i]];
      //printf("from=%i\tto=%i\tproba=%e\n",_Pfrom[i],_Pto[i],_Pproba[i]);
    }
  } // end transition normalization
};


void amarkov::leftprod(double *x,double t,double *res){
  // initalization
  for (int i=0; i<_nstate; i++)
    res[i]=0.0;
  // loop on P
  for (int i=0; i<_Pnz; i++)
    res[_Pto[i]]+=x[_Pfrom[i]]*_Pproba[i];
  // loop on Q
  for (int i=0; i<_Qnz; i++)
    res[_Qto[i]]+=t*x[_Qfrom[i]]*_Qproba[i];
};

void amarkov::rightprod(double *x,double t,double *res){
  // initalization
  for (int i=0; i<_nstate; i++)
    res[i]=0.0;
  // loop on P
  for (int i=0; i<_Pnz; i++) {
    //printf("from=%i\tto=%i\tproba=%e\n",_Pfrom[i],_Pto[i],_Pproba[i]);
    //printf("res[%i]+=%f*%f\n",_Pfrom[i],x[_Pto[i]],_Pproba[i]);
    res[_Pfrom[i]]+=x[_Pto[i]]*_Pproba[i];
  }
  // loop on Q
  for (int i=0; i<_Qnz; i++)
    res[_Qfrom[i]]+=t*x[_Qto[i]]*_Qproba[i];
};


void amarkov::compute_stationary(bool verbose){

  // fixme: for nstate < DNAUPD_KRYLOVDIM, arnoldi sometimes
  // fails. Should be far better to switch to full linear algebra
  // in those case.

  // fixme: One should add a verification at the end of this function
  // to check if _stationary is right or not

  // simple test
  //double *aux=(double*)malloc(sizeof(double)*_nstate);
  //double *swap;
  //for (int i=0; i<_nstate; i++)
  //  aux[i]=1.0/_nstate;
  //for (int i=0; i<10; i++) {
  //  leftprod(aux,1.0,_stationary);
  //  swap=aux; aux=_stationary; _stationary=swap;
  //  printf("i=%i\taux[1]=%e\n",i,aux[1]);
  //}

  /* create a dnaupd_par */
  // fixme: treat the case _nstate<DNAUPD_KRYLOVDIM
  // with full linear algebra rather than sparse one
  int krylovdim=DNAUPD_KRYLOVDIM;
  if (_nstate<krylovdim)
    krylovdim=_nstate;
  if (verbose)
    printf("krylovdim=%i\n",krylovdim);
  dnaupd_par dpar(_nstate,krylovdim);
  dpar._nev=2; // we need 2 eigenvalues
  //dpar.print();

  /* get an initial random vector */
  srand(time(0));  
  for (int i=0; i<dpar._n; i++) {
    dpar._resid[i]=rand()/(double)RAND_MAX;
  }
  dpar._info=1;
  
 main:
  dnaupd_(dpar._ido,dpar._bmat,dpar._n,dpar._which,dpar._nev,dpar._tol,dpar._resid,dpar._ncv,dpar._v,dpar._ldv,dpar._iparam,dpar._ipntr,dpar._workd,dpar._workl,dpar._lworkl,dpar._info);
  //printf("dnaupd_call done\n");
  printf("dnaupd ido=%i\n",dpar._ido);
  if (dpar._ido==-1 || dpar._ido==1) {
    /* compute product Y=OP*X where */
    /* ipntr[0]-1 is the pointer into workd for X */
    /* ipntr[1]-1 is the pointer into workd for Y */ 
    //printf("matrix x vector product\n");
    double *x=&dpar._workd[dpar._ipntr[0]-1];
    double *y=&dpar._workd[dpar._ipntr[1]-1];
    leftprod(x,1.0,y);    
    goto main; 
  } /* end if ido */
  if (verbose) {
    printf("implicitly restarted algorithm returned:\n");
    printf("\tinfo=%i\n",dpar._info);
    printf("\tnumber of iter=%i\n",dpar._iparam[2]);
    printf("\tnumber of converged ritzvalues=%i\n",dpar._iparam[4]);
    printf("\tnumber of prod=%i\n",dpar._iparam[8]);
    printf("\tnumber of re-orth steps=%i\n",dpar._iparam[10]);    
  }
  { // post treatment
    //printf("ipntr:");
    //for (int i=0; i<14; i++)
    //  printf("[%i]:%i\t",i+1,dpar._ipntr[i]);
    //printf("\n");
    double *re,*im;
    re=&dpar._workl[dpar._ipntr[5]-1];
    im=&dpar._workl[dpar._ipntr[6]-1];
    if (verbose) {
      printf("eigenvalues:\n");
      for (int i=0; i<dpar._iparam[4]; i++) {
    	printf("\tre=%e\tim=%e\n",re[i],im[i]);
      }
    }
    
    /* call of dneupd */
    printf("calling dneupd_par\n");
    dneupd_par ddpar(dpar._n,dpar._nev,DNEUPD_GET_VECTORS);
    //dneupd_par ddpar(dpar._n,dpar._iparam[4],DNEUPD_GET_VECTORS);
    //ddpar.print();
    //dpar.print();
    //printf("ipntr:");
    //for (int i=0; i<14; i++)
    //  printf("[%i]:%i\t",i+1,dpar._ipntr[i]);
    //printf("\n");    
    dneupd_(ddpar._rvec,ddpar._howmny,ddpar._select,ddpar._dr,ddpar._di,ddpar._z,ddpar._ldz,ddpar._sigmar,ddpar._sigmai,ddpar._workev,
	    dpar._bmat,dpar._n,dpar._which,dpar._nev,dpar._tol,dpar._resid,dpar._ncv,dpar._v,dpar._ldv,dpar._iparam,dpar._ipntr,dpar._workd,dpar._workl,dpar._lworkl,dpar._info);      
    
    printf("returned info:%i\n",dpar._info);
    /* print results */
    if (verbose) {
      printf("eigenvalues:\n");
      for (int i=0; i<ddpar._nev; i++) {
        printf("\tre=%f\tim=%f\n",ddpar._dr[i],ddpar._di[i]);
      }
    }
    if ( fabs(ddpar._dr[0]-1.0)>EPSILON || fabs(ddpar._di[0]>EPSILON) ) {
      fprintf(stderr,"eigenvalue 1.0 not found in spectrum. Check Markov parameters or try another run.\n");
      exit(EXIT_FAILURE);
    }
    printf("eigenvector:\n");      
    {
      double sum=0.0;
      for (int i=0; i<ddpar._n; i++)
	sum+=ddpar._z[i];	
      for (int i=0; i<ddpar._n; i++) {
	if (verbose) {
	  printf("%e\n",ddpar._z[i]/sum);
	  //printf("%e (%e)\n",ddpar._z[i]/sum,aux[i]);
	}
	_stationary[i]=ddpar._z[i]/sum;
      }
    } // end stationary block
    if (verbose)
      printf("secondmag=%f\n",_secondmag);

  } // end post treatment
};

double amarkov::nexp(long n) {

  double res=0.0;
  for (int i=0; i<_nfinal; i++) {
    res+=(n-_h[i]+1)*_stationary[_final[i]];
  }
  return res;
};


double amarkov::bstat(long n,long nobs,double nexp) {

  double m,e,r,res;
  //printf("nobs=%i\tnexp=%e\n",nobs,nexp);
  if (nobs>nexp) {
    /* pattern is over-represented */
    r=qcdfbin(nobs,n-_min_h+1,nexp/(n-_min_h+1),&m,&e);
    res=-e-log(m)/log(10.0);
  } else {
    /* pattern is under-represented */
    r=pcdfbin(nobs,n-_min_h+1,nexp/(n-_min_h+1),&m,&e);
    res=e+log(m)/log(10.0);
  }
  return res;
};

double amarkov::cpstat(long n,long nobs,double nexp,bool verbose,bool debug) {

  double res=0.0;

  // allocate A
  double **A=(double **)malloc(sizeof(double*)*_nfinal);
  if (!A) {
    fprintf(stderr,"not enough memory to allocate A in bstat. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  A[0]=(double *)malloc(sizeof(double)*_nfinal*_nfinal);
  if (!A[0]) {
    fprintf(stderr,"not enough memory to allocate A[0] in bstat. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  { // connect
    double *pos=A[0];
    for (int i=0; i<_nfinal; i++) {
      A[i]=pos;
      pos+=_nfinal;
    } // end for
  } // end connect
  // end allocate A

  // allocate aux vectors x (size nstate) and y (size nfinal)
  double *data=(double *)malloc(sizeof(double)*(2*_nstate+_nfinal));
  if (!data) {
    fprintf(stderr,"not enough memory to x in bstat. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  double *aux=data;
  double *x=data+_nstate;
  double *y=data+2*_nstate;
  // end aux allocation

  // main loop on final states
  for (int f=0; f<_nfinal; f++) {
    // init
    for (int i=0; i<_nfinal; i++)
      y[i]=0.0;
    for (int i=0; i<_nstate; i++)
      x[i]=0.0;
    x[_final[f]]=1.0;
    int h=0;
    bool finished=false;
    while (!finished) {
      // x = x * Pi
      leftprod(x,1.0,aux);
      double *swap=x;
      x=aux;
      aux=swap;
      h++;
      finished=true;
      // y = y + x(F) and x(F)=0
      for (int i=0; i<_nfinal; i++) {
	//printf("x[_final[%i]]=%f\n",i,x[_final[i]]);
	if (h<_h[i]) {
	  y[i]+=x[_final[i]];
	  finished=false;
	}
	x[_final[i]]=0.0;
      }
    } // end while loop
    // A(s,:)=y
    for (int i=0; i<_nfinal; i++)
      A[f][i]=y[i];
  } // and loop on final state
			
  if (verbose) {
    printf("A:\n");
    for (int i=0; i<_nfinal; i++) {
      for (int j=0; j<_nfinal; j++)
	printf("%e\t",A[i][j]);
      printf("\n");
    }
  }

//  for (int i=0; i<_nfinal; i++) {
//    printf("%i\n",i);
//    printf("x[%i]=%e\n",i,x[i]);
//    printf("y[%i]=%e\n",i,y[i]);
//    printf("aux[%i]=%e\n",i,aux[i]);
//  }
  
  //double nexp=0.0;
  if (_nfinal>0) {
    double esp=0.0;
    // compound poisson case
    if (verbose)
      printf("compound Poisson\n");
    // fill mass
    double mass[MASS_KMAX];
    double totmass=0.0;
    // put (I-A)*1 in y
    for (int i=0; i<_nfinal; i++) {
      y[i]=1.0;
      for (int j=0; j<_nfinal; j++) {
	y[i]-=A[i][j];
      }
      totmass+=y[i]*_stationary[_final[i]];
    }
    printf("totmass=%e\n",totmass);
    // put (I-A)*y in x 
    for (int i=0; i<_nfinal; i++) {
      x[i]=y[i];
      for (int j=0; j<_nfinal; j++) {
	x[i]-=A[i][j]*y[j];
      }
      if (verbose)
	printf("((I-A)^2 . 1)[%i]=%e\n",i,x[i]);
    }
    // put  (n-h+1)mu in y
    for (int i=0; i<_nfinal; i++) {
      //y[i]=(n-_h[i]+1)*_stationary[_final[i]];
      y[i]=_stationary[_final[i]];
      esp+=y[i];
      //nexp+=y[i];
      if (verbose)
	printf("mu[%i]=%e\n",i,y[i]);
    }
    esp*=n;
    double sum=0.0;
    double *swap;
    int k=0;
    while ( (1.0-sum)>MASS_TOL && k<MASS_KMAX) {
      // mass[k] = y * x
      mass[k]=0;
      for (int i=0; i<_nfinal; i++) {
	mass[k]+=y[i]*x[i];
      }
      mass[k]/=totmass;
      sum+=mass[k];
      //esp+=(k+1)*mass[k];
      if (verbose) {
	printf("mass[%i]=%e (test=%e)\n",k,mass[k],(1.0-sum));
      }
      // compute y = y * A
      for (int i=0; i<_nfinal; i++) {
	aux[i]=0.0;
	for (int j=0; j<_nfinal; j++) {
	  aux[i]+=y[j]*A[j][i];
	  //printf("x[%i] * A[%i][%i] = %e * %e = %e\n",j,j,i,x[j],A[j][i],x[j]*A[j][i]); 
	}
	//printf("aux[%i]=%e\n",i,aux[i]);
      }
      swap=y;
      y=aux;
      aux=swap;
      // iter k
      k++;      
    } // end while
    if (k==MASS_KMAX) {
      fprintf(stderr,"Warning ! lambda_k has not yet converged, %e is still missing.\n",(1.0-sum));
    }
    // correct lambda to fit with expectation
    double lambda=n*totmass;//(nexp/esp);
    if (verbose) {
      printf("k=%i\tesp=%e\tlambda=%e\n",k,esp,lambda);
    }
    if (k==1) {
      // it is a simple binomial case
      return bstat(n,nobs,nexp);
    }
    // reduce k for computations
    if (k>nobs)
      k=nobs+1;
    // allocate working space
    //printf("nobs=%i\n",nobs);
    double **work=(double **)malloc(sizeof(double *)*k);
    if (work==NULL) {
      fprintf(stderr,"not enough memory to compute cpstat. Aborting.\n");
      exit(EXIT_FAILURE);
    }
    work[0]=(double *)malloc(sizeof(double)*k*(nobs+1));
    if (work[0]==NULL) {
      fprintf(stderr,"not enough memory to compute cpstat. Aborting.\n");
      exit(EXIT_FAILURE);
    }
    { // connect
      double *dpos=work[0];
      for (int i=0; i<k; i++) {
	work[i]=dpos;
	dpos+=(nobs+1);
      }
    } // end connect
    // initialize work[0][j]
    double m,e;
    if (nobs < nexp) {
      if (verbose)
	printf("use cdf\n");
      // use P(<=j)
      for (long j=0; j<(nobs+1); j++)
	//work[0][j]=exp( pcdfpoi((double)j,lambda*mass[0],&m,&e) )*exp(lambda*(mass[0]-1.0));
	// new log computations below
	work[0][j]=pcdfpoi((double)j,lambda*mass[0],&m,&e)+(lambda*(mass[0]-1.0));
    } else {
      // fixme: many issues with ccdf. A tail sum must be done.
      // current implementation should not be considered reliable for ccdf
      if (verbose)
	printf("use ccdf\n");
      // use P(>=j)      
      for (long j=0; j<(nobs+1); j++) {
	//printf("test=%f\n",exp( qcdfpoi((double)j,lambda*mass[0],&m,&e) ));
	//work[0][j]=exp( qcdfpoi((double)j,lambda*mass[0],&m,&e) )*exp(lambda*(mass[0]-1.0));      
	// new log computations below
	work[0][j]=qcdfpoi((double)j,lambda*mass[0],&m,&e)+(lambda*(mass[0]-1.0));      
      }
    }
    // end work[1][j] initialization
    // main loop
    for (int i=1; i<k; i++) {
      // compute work[i][j]
      for (long j=0; j<(nobs+1); j++) {
	long lmax=j/(i+1);
	//printf("lmax=%i/%i=%i\n",j,i,lmax);
	//work[i][j]=0.0;
	//for (long l=0; l<=lmax; l++) {
	//  work[i][j]+=pow(lambda*mass[i],(double)l)/factrl(l)*work[i-1][j-l*(i+1)];
	//} // end l loop
	//if (nobs >= nexp) {
	//  work[i][j]+=exp(qcdfpoi((double)(lmax+1),lambda*mass[i],&m,&e))*exp(lambda*mass[i])*work[i-1][0];
	//}
	// log computations below
	work[i][j]=1.0;
	for (long l=1; l<=lmax; l++) {
	  work[i][j]+=pow(lambda*mass[i],(double)l)/factrl(l)*exp(work[i-1][j-l*(i+1)]-work[i-1][j]);
	}
	if (nobs >= nexp) {
	  work[i][j]+=exp(qcdfpoi((double)(lmax+1),lambda*mass[i],&m,&e))*exp(lambda*mass[i]);
	}
	work[i][j]=work[i-1][j]+log(work[i][j]);
      } // end j loop
    } // end main loop
    if (debug) {
      printf("work:\n",k,nobs+1);
      for (int i=0; i<k; i++) {
	for (long j=0; j<(nobs+1); j++) {
	  printf("work[%i][%i]=%e\t",i,j,work[i][j]);
	}
	printf("\n");
      }
    }
    // get the pvalue
    //    res=-log(work[k-1][nobs])/log(10.0);
    // log computation below
        res=-work[k-1][nobs]/log(10.0);
    // free working space
    free(work[0]);
    free(work);
    // end compound poisson case
  } else {
    // geometric poisson case
    if (verbose)
      printf("geometric Poisson\n");
    if (nobs>nexp)
      res=-pgeopois(nobs,(1.0-A[0][0])*nexp,1.0-A[0][0],0,1)/log(10.0);
    else
      res=-pgeopois(nobs,(1.0-A[0][0])*nexp,1.0-A[0][0],1,1)/log(10.0);      
    // geometric poisson
  }

  // free memory
  free(A[0]);
  free(A);
  free(data);

  if (nobs>nexp) {
    return res;
  } else {
    return -res;
  }
};

bool amarkov::is_null(double *x,int size){
  bool res=true;
  for (int i=0; i<size; i++){
    if (x[i]!=0.0) {
      res=false;
      break;
    }
  }
  return res;
};

void amarkov::compute_h(bool verbose){

  double *data,*x,*y,*swap;
  
  _min_h=10000;

  if (verbose) {
    printf("compute_h: ");
  }
  
  // memory allocation
  data=(double *)malloc(sizeof(double)*2*_nstate);
  if (data==NULL) {
    fprintf(stderr,"not enough memory to compute h. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  x=data;
  y=data+_nstate;

  // init x to dirac in 0
  for (int i=0; i<_nstate; i++) 
    x[i]=0.0;
  x[0]=1.0;
  
  // init _h
  for (int i=0; i<_nfinal; i++)
    _h[i]=-1;

  // main loop
  bool finished=false;
  int h=0;
  while (!finished) {
    // x = x * P
    leftprod(x,1.0,y);
    swap=x;
    x=y;
    y=swap;
    h++;
    // get min size
    finished=true;
    for (int i=0; i<_nfinal; i++) {
      if (_h[i]==-1) {
	finished=false;
	if (x[_final[i]]>0) {
	  _h[i]=h;
	  if (h<_min_h)
	    _min_h=h;
	}
      } // end if _h[i]
    } // loop on i

  } // end while
  
  // free memory
  free(data);

  if (verbose) {
    for (int i=0; i<_nfinal; i++) {
      printf("h[%i]=%i\t",i,_h[i]);
    }
    printf("\n");
  }
};

//double amarkov::xstat(long n,long nobs, double nexp,bool verbose,bool debug){
//
//  // fixme : must code here log_P_times and  log_P_times_plus_Q_times
//  // xstat not yet sucessfully compiled
//
//
//  long order=_nstates;
//  long nblock;
//  if (nobs>nexp) {
//    nblock=nobs;
//  } else {
//    nblock=nobs+1;
//  }
//  long dim=order*nblock;
//  // memory allocation
//  double *workspace=(double *)malloc(sizeof(double)*(2*dim+order));
//  if (workspace==NULL) {
//    fprintf(stderr,"not enough memory for workspace allocation in xstat\n");
//    exit(EXIT_FAILURE);
//  }
//  double *current_u=workspace;
//  double *last_u=workspace+dim;
//  double *sum=workspace+2*dim;
//  /* two cases: over and under */
//  if (nobs>nexp) {
//    /* over case */
//    // initialization
//    for (long i=0; i<dim; i++) {
//      current_u[i]=0.0;
//    }
//    // put sum of Q (per line) in current_u
//    for (long i=0; i<order; i++) {
//      current_u[i]=0.0;
//    }
//    for (int i=0; i<_Qnz; i++) {
//      current_u[_Qfrom[i]]+=_Qproba[i];
//    }
//    // current_u = log (current_u )
//    for (long i=0; i<dim; i++) {
//      current_u[i]=log(current_u[i]);
//    }
//    for (long i=0; i<order; i++) {
//      sum[i]=log(0.0);
//    }
//    double *aux=NULL;
//    double *pos_in_current,*pos_in_last1,*pos_in_last2;
//    // main loop
//    for (long i=2; i<=_n-1; i++) {
//      // switch last_u and current_u;
//      aux=last_u;
//      last_u=current_u;
//      current_u=aux;
//      // put fmci transition times last_u in current_u
//      pos_in_current=current_u;
//      pos_in_last1=last_u;
//      // block 0
//      // current_u_0 = P * last_u_0
//      _T->log_P_times(pos_in_last1,pos_in_current);		
//      // loop on blocks
//      for (long j=1; j<nblock; j++) {
//	// new positions
//	pos_in_current+=order;
//	pos_in_last2=pos_in_last1;
//	pos_in_last1+=order;
//	// current_u_j = P * last_u_j + Q * last_u_j-1
//	_T->log_P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
//      } // end blocks loop
//	// get largest term of sum and add
//      double z=log(0.0);
//      for (long j=0; j<order; j++) {
//	if (sum[j]>z)
//	  z=sum[j];
//      }
//      for (long j=0; j<order; j++) {
//	if (pos_in_current[j]>z)
//	  z=pos_in_current[j];
//      }
//      if (z==log(0.0))
//	z=0;
//      //printf("z=%e\n",z);
//      // sum = sum + last block of current_u (pointed by pos_in_current)	
//      for (long j=0; j<order; j++) {
//	sum[j]=exp(sum[j]-z)+exp(pos_in_current[j]-z);
//      }
//      for (long j=0; j<order; j++) {
//	sum[j]=log(sum[j])+z;
//      }
//      //printf("i=%i\tlog10(P|1)=%f\n",i,sum[0]/log(10.0));
//    }
//    // get the final result
//    {
//      int imax=_M->get_stationary_size();
//      _pvalue=0.0;
//      // find largest term
//      double z=log(0.0);
//      for (int i=0; i<imax; i++) {
//	if (sum[i]>z)
//	  z=sum[i];
//      }
//      // compute sum
//      double *_stationary=_M->get_stationary();
//      for (int i=0; i<imax; i++) {
//	_pvalue+=_stationary[i]*exp(sum[i]-z);
//      }
//      _pvalue=log(_pvalue)+z;
//      _stat=-_pvalue/log(10.0);
//      _pvalue=exp(_pvalue);
//    }
//    /* end of over case */
//  } else {
//    /* under case */
//    // initialization
//    for (long i=0; i<dim; i++) {
//      current_u[i]=0.0;
//    }
//    double *aux=NULL;
//    double *pos_in_current,*pos_in_last1,*pos_in_last2;
//    //double last=1.0;
//    // main loop
//    for (long i=2; i<=_n; i++) {
//      // switch last_u and current_u;
//      aux=last_u;
//      last_u=current_u;
//      current_u=aux;
//      // put fmci transition times last_u in current_u
//      pos_in_current=current_u;
//      pos_in_last1=last_u;
//      // block 0
//      // current_u_0 = P * last_u_0
//      //_T->P_times(pos_in_last1,pos_in_current);		
//      _T->log_P_times(pos_in_last1,pos_in_current);		
//      // loop on blocks
//      for (long j=1; j<nblock; j++) {
//	// new positions
//	pos_in_current+=order;
//	pos_in_last2=pos_in_last1;
//	pos_in_last1+=order;
//	// current_u_j = P * last_u_j + Q * last_u_j-1
//	//_T->P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
//	_T->log_P_times_plus_Q_times(pos_in_last1,pos_in_last2,pos_in_current);
//      } // end blocks loop
//	//printf("i=%i\tu[0]=%e\n",i,pos_in_current[0]);
//	//{
//	//  _pvalue=0.0;
//	//  double *_stationary=_M->get_stationary();
//	//  for (int ii=0; ii<_alpha->get_size(); ii++) {
//	//    _pvalue+=_stationary[ii]*pos_in_current[ii];
//	//  }
//	//  printf("i=%i\tpvalue=%e\tdiff=%e\n",i,_pvalue,last-_pvalue);
//	//  last=_pvalue;
//	//}
//      
//    } // end main loop
//      // get the final result
//    {
//      int imax=_M->get_stationary_size();
//      _pvalue=0.0;
//      // find largest term
//      double z=log(0.0);
//      for (int i=0; i<imax; i++) {
//	if (pos_in_current[i]>z)
//	  z=pos_in_current[i];
//      }
//      // compute sum
//      double *_stationary=_M->get_stationary();
//      for (int i=0; i<imax; i++) {
//	_pvalue+=_stationary[i]*exp(pos_in_current[i]-z);
//      }
//      _pvalue=log(_pvalue)+z;
//      _stat=_pvalue/log(10.0);
//      _pvalue=exp(_pvalue);
//    }
//    /* end of under case */
//  }
// 
//};


amarkov::~amarkov() {

  if (_start!=NULL)
    free(_start);
  if (_data!=NULL)
    free(_data);
  if (_stationary!=NULL)
    free(_stationary);
  if (_final!=NULL)
    free(_final);

}; 

