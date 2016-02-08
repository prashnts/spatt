/* $Id: cp.cc 616 2006-01-05 09:05:59Z gnuel $ */
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

#include "cp.h"

#ifndef INFINITY
#define INFINITY 1e300
#endif

using namespace std;

namespace spatt {
double logfact(int n){

  static double *table=NULL;
  static int ntop=-1;
  static int nblock=0;

  /* if n<0 then free memory */
  if (n<0) {
    if (table!=NULL)
      free(table);
    return -1;
  } else {
    /* if ntop==-1 then alloc and compute first block */
    if (ntop==-1) {
      /* alloc */
      table=(double *)malloc(sizeof(double)*(LOGFACT_BLOCK_SIZE+1));
      /* and compute */
      ntop=LOGFACT_BLOCK_SIZE;
      nblock=1;
      table[0]=0.0;
      for (int i=1; i<=ntop; i++)
	table[i]=table[i-1]+log((double)i);
    }
    /* if requested value not yet computed, extend table */
    if (n>ntop) {
      nblock+=(n-ntop)/LOGFACT_BLOCK_SIZE+1;
      realloc(table,sizeof(double)*(nblock*LOGFACT_BLOCK_SIZE+1));
      for (int i=ntop+1; i<=nblock*LOGFACT_BLOCK_SIZE; i++)
	table[i]=table[i-1]+log((double)i);
      ntop=nblock*LOGFACT_BLOCK_SIZE;
    }
    /* return result */
    return table[n];
  }
}

double cpA(const char *w,const markov *M){
  int i,j,h;
  int *is_period=NULL;
  double A=0.0;

  h=strlen(w);
  
  // is_period initialization
  is_period=(int *)malloc(sizeof(int)*(h+1));
  is_period[0]=h-1;
  for (i=1; i<h; i++)
    is_period[i]=1;

  // loop on possible periods
  for (i=1; i<h; i++) {
    // check if i is a period
    // printf("checking %i ...\n",i);
    j=0;
    while (is_period[i]==1 && j<(h-i)) {
      if (w[j]!=w[i+j]) {
	is_period[i]=0;
	is_period[0]--;
      } else {
	j++;
      }
    }
    // if it is, remove multiple
    // and do optionnal action
    if (is_period[i]==1) {
      // remove multiple
      j=2*i;
      while (j<h) {
	is_period[j]=0;
	is_period[0]--;
	j+=i;
      }
      // action
      //printf("%i is a principal period\n",i);
      A+=M->trans(w,&w[h-i]);
    }
  }
  // end loop on possible periods
  
  //printf("total number of principal periods: %i\n",is_period[0]);
  
  return A;
}

double newcpstat(double expectation,double A,double sstat,long N){

  if (sstat>0) {
    return -pgeopois(N,(1.0-A)*expectation,1.0-A,0,1)/log(10.0);
  } else {
    return pgeopois(N,(1.0-A)*expectation,1.0-A,1,1)/log(10.0);
  }
}

/* If N ~ GP(lambda,theta) geometric Poisson distribution of parameter   */
/* lambda for the Poisson part and theta in ]0,1[ for the geometric part */
/* this function compute:                                                */
/* - P(N>=x) if lowertail=0 and logscale=0                               */
/* - log P(N>=x) if lowertail=0 and logscale!=0                          */
/* - P(N<=x) if lowertail!=0 and logscale=0                              */
/* - log P(N<=x) if lowertail!=0 and logscale!=0                         */
double pgeopois(double x,double lambda,double theta,int lowertail,int logscale){

  long n0,n;
  double z,res;
  double S,A,Lcurrent,Lprec,Lprecprec;
  int converged=0;

  /* check parameters */
  if ( theta <=0 || theta>=1 || lambda<=0 || x<0 ) {
    fprintf(stderr,"Check parameter values !\n");
    exit(EXIT_FAILURE);
  }

  /* compute n0 from x */
  n0=(long)floor(x);

  /* compute z */
  z=lambda*theta/(1-theta);

  /* degenerate cases */
  if (n0==0) {
    if (lowertail)
      res=-lambda;
    else
      res=0.0;
    goto final;
  }
  if (n0==1) {
    if (lowertail)
      res=-lambda+log(1.0+(1.0-theta)*z);
    else
      res=log(-expm1(-lambda));
    goto final;
  }      

  /* from now n0>=2 */
  
  /* compute L0 and L1 */
  Lprecprec=-lambda;
  Lprec=-lambda+log( (1.0-theta)*z );  

  if (lowertail) {
    /* initialize cumsum */
    A=Lprecprec;
    S=1.0+(1.0-theta)*z;
    for (n=2; n<=n0; n++) {
      /* log recurrence */
      nextL(&Lcurrent,&Lprec,&Lprecprec,n,z,theta);
      /* update cumsum */
      updateSum(&A,&S,Lcurrent);
    }
    res=A+log(S);
  }

  if (!lowertail) {
    /* compute first n0 terms */
    for (n=2; n<=n0; n++) {
      /* log recurrence */
      nextL(&Lcurrent,&Lprec,&Lprecprec,n,z,theta);
    }    
    /* initialize cumsum */
    A=Lprec;
    S=1.0;
    converged=0;
    n=n0+1;
    while (!converged) {
      /* log recurrence */
      nextL(&Lcurrent,&Lprec,&Lprecprec,n,z,theta);
      /* check convergence */
      if (Lcurrent<LOGEPSILON+A+log(S))
	converged=1;
      /* update cumsum */
      updateSum(&A,&S,Lcurrent);
      /* increment n */
      n++;
    }
    res=A+log(S);
  }

  /* Returning the results */
 final:
  if (logscale)
    return res;
  else
    return exp(res);
};

void nextL(double *Lcurrent, double *Lprec, double *Lprecprec, long n,double z, double theta){
    *Lcurrent=(2.0*n-2+z)/n*(1.0-theta);
    *Lcurrent+=(2.0-n)/n*(1.0-theta)*(1.0-theta)*exp(*Lprecprec-*Lprec);
    *Lcurrent=*Lprec+log(*Lcurrent);
    *Lprecprec=*Lprec;
    *Lprec=*Lcurrent;  
};

void updateSum(double *A, double *S,double logadd){
  static double newS;
  /* compute next S */
  newS=*S+exp(logadd-*A);
  /* check range */
  if (newS>=INFINITY || newS==0) {
    //printf("change scale A=%f -> A=%f\n",*A,*A+log(*S));
    *A=*A+log(*S);
    *S=1+exp(logadd-*A);
  } else {
    *S=newS;
  }
};

  double complete_dcpoi(int n,double lambda,double theta, int alpha,double *p){

    
    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation

    double **D;
    double *P;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(n+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // n+2 (and not n+1) for P
      D[0]=(double *)malloc(sizeof(double )*(n+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=n; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      P=pos;
    } // end memory allocation

    // precompute all P[i]
    for (int i=1; i<=n; i++) {
      if (i<=alpha) {
	P[i]=p[i-1];
      } else {
	P[i]=C*pow(theta,i-alpha);
      }
    }

    { // initialization
      D[1][0]=exp(-lambda);
      double aux=P[1]*lambda;
      for (int j=1; j<=n; j++) {
	D[1][j]=D[1][j-1]*aux/(double)j;
	//printf("D[%i][%i]=%f\n",1,j,log(D[1][j]));
      }
    } // end initialization
    
    
    { // main loop
      //printf("lambda*P[1]=%e\n",lambda*P[1]);
      for (int i=2; i<=n; i++) {
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux=1.0;
	  for (int k=0; k<=kmax; k++) {
	    D[i][j]+=aux*D[i-1][j-i*k];
	    aux*=lambda*P[i]/(k+1.0);
	  } // end k loop
	} // end j loop
      } // end i loop
    } // end main loop

    // get result
    double res=D[n][n];

    //printf("complete_dcpoi:\n");
    //for (int j=1; j<=n-alpha; j++) { 
    //  printf("%e\t",D[alpha][n-alpha-j]);
    //}
    //printf("\n");
    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;

  }

  double barbour_dcpoi(int n,double lambda,double theta, int alpha,double *p) {

    //timer T;

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation

    double **D;
    { // memory allocation
      //malloc(10000000);
      D=(double **)malloc(sizeof(double *)*(n+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      D[0]=(double *)malloc(sizeof(double )*(n+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=0; i<=n; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
    } // end memory allocation

    { // initialization
      for (int j=1; j<=n; j++) {
	D[j][j]=1.0/(double)j;
	// put lambdaj in D[0][j]
	D[0][j]=0.0;
	if (j<=alpha) {
	  D[0][j]=lambda*p[j-1];
	} else {
	  D[0][j]=lambda*C*pow(theta,j-alpha);
	}
      }
    } // end initialization

    //T.start();
    { // main loop
      //int nop=0;
      for (int i=1; i<=n; i++) {	
	//printf("i=%i\tnop=%i\tt=%.2fs for main loop\n",i,nop,T.elapsed_time());     
	//nop=0;
	for (int j=1; j<i; j++) {
	  //printf("C_{%i,%i}:\n",i,j);
	  // compute c_{i,j}
	  D[i][j]=0.0;
	  //double aux1=0.0;
	  //double *aux2=D[0];
	  for (int k=1; k<=i-j; k++) {
	    //printf("%i * lambda_%i / %i * C_{%i,%i}\n",k,k,i,i-k,j);
	    D[i][j]+=k/(double)i*D[0][k]*D[i-k][j];
	    //aux1+=k/(double)i*aux2[k]*D[i-k][j];
	    //nop++;
	  } // end k loop
	  //D[i][j]=aux1;
	} // end loop on j
      } // end loop on i
    } // end main loop
    //printf("\n%.2fs for main loop\n",T.elapsed_time());     

    // get result
    double res=0.0;
    for (int k=1; k<=n; k++) {
      double lambdak;
      res+=k*D[0][k]*D[n][k];
    }
    if (n==0) {
      res=1.0;
    }
    res*=exp(-lambda);


    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;
  }

  double fast_dcpoi(int n,double lambda,double theta, int alpha,double *p){
    
    //printf("fast_dcpoi call\n");
    //timer T;
    //T.start();

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation
    //printf("C=%e\n",C);

    double **D;
    double *CC;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(alpha+1));
      //printf("D=%p\n",D);
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // alpha+2 (and not alpha+1) for C
      D[0]=(double *)malloc(sizeof(double )*(alpha+1)*(n+1));
      //printf("D[0]=%p\n",D[0]);
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=alpha; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      CC=pos;
    } // end memory allocation

    { // initialization
      D[1][0]=exp(-lambda);
      double aux=p[0]*lambda;
      for (int j=1; j<=n; j++) {
	D[1][j]=D[1][j-1]*aux/(double)j;
      }
    } // end initialization
    
    { // precompute all necessary C[j]
      // allocate memory
      double *workspace=(double *)malloc(sizeof(double)*2*(n/(alpha+1)));
      if (workspace==NULL) {
	fprintf(stderr,"not enough memory in dcpoi. Aborting.\n");
	exit(EXIT_FAILURE);
      }
      double *X,*Y;
      X=workspace; // for (lambda*C)^m/m!*theta^(alpha+j-m*alpha)
      Y=workspace+(n/(alpha+1)); // for binomial(alpha+j-m*alpha-1,m-1)
      int mmax=1;
      X[0]=lambda*C;
      Y[0]=1.0;
      double aux=pow(theta,-(double)alpha);
      //printf("lambda*C=%e\ttheta=%e\talpha=%i\ttheta^{-alpha}=%e\n",X[0],theta,alpha,aux);
      // main loop
      CC[1]=lambda*C*theta;
      X[0]*=theta;
      for (int j=2; j<=n-alpha; j++) {
	//printf("j=%i\n",j);
	// initialization
	CC[j]=0.0;
	// update X and Y
	X[0]*=theta;
	CC[j]+=X[0]*Y[0];
	for (int m=1; m<mmax; m++) { // beware, replace m by m+1 in theoretical formulas
	  X[m]*=theta;
	  Y[m]*=( (double)(j-alpha*m-1)/(double)(j-alpha*m-m-1) );
	  CC[j]+=X[m]*Y[m];
	}
	// check mmax
	int newmmax=(alpha+j)/(alpha+1);
	if (newmmax>mmax) {
	  //printf("newmmax=%i\tmmax=%i\n",newmmax,mmax);
	  mmax=newmmax;
	  X[mmax-1]=X[mmax-2]*lambda*C/(double)(mmax)*aux;
	  Y[mmax-1]=1.0;
	  // update CC[j]
	  CC[j]+=X[mmax-1]*Y[mmax-1];
	}
	// print stuff
	//printf("j=%i 1<=m<=%i:\n",j,mmax);
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("X(%i)=%e\t",m+1,X[m]);
	//}
	//printf("\n");
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("Y(%i)=%e\t",m+1,Y[m]);
	//}
	//printf("\n");
	//printf("C(alpha,%i)=%e\n",j,CC[j]);
      } // end j loop;
      // free memory
      free(workspace);
    } // end precompute
    //printf("%.2fs to compute all C(a,j)\n",T.elapsed_time()); 

    { // main loop
      for (int i=2; i<=alpha; i++) {
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux=1.0;
	  for (int k=0; k<=kmax; k++) {
	    //D[i][j]+=pow(lambda*p[i-1],k)/factrl(k)*D[i-1][j-i*k];
	    D[i][j]+=aux*D[i-1][j-i*k];
	    aux*=lambda*p[i-1]/(k+1);
	  } // end k loop
	  //if (D[i][j]==0.0) {
	  //  warning=true;
	  //}
	} // end j loop
      } // end i loop
    } // end main loop
    //printf("%.2fs to compute all D(i,j)\n",T.elapsed_time()); 

    // get result
    double res=D[alpha][n];
//    printf("fast_dcpoi:\n");
//    for (int j=1; j<=n-alpha; j++) { 
//      printf("%e\t",D[alpha][n-alpha-j]);
//    }
//    printf("\n");
    for (int j=1; j<=n-alpha; j++) {
      res+=D[alpha][n-j-alpha]*CC[j];
      //printf("%e*%e\t",D[alpha][n-alpha-j],CC[j]);
    } // end j loop
    //printf("\n");
    //printf("%.2fs to compute the final result\n",T.elapsed_time()); 

    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    //printf("D[0]=%p\n",D[0]);
    //printf("D=%p\n",D);
    free(D[0]);
    free(D);


    //printf("%.2fs global computional time\n",T.elapsed_time()); 

    // return result
    return res;

  }

  double complete_log_dcpoi(int n,double lambda,double theta, int alpha,double *p){


    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation

    double **D;
    double *P;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(n+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // n+2 (and not n+1) for P
      D[0]=(double *)malloc(sizeof(double )*(n+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=n; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      //printf("test: %p = %p ?\n",pos,D[0]+(n+1)*(n+1));
      P=pos;
    } // end memory allocation

    // precompute all P[i]
    for (int i=1; i<=n; i++) {
      if (i<=alpha) {
	P[i]=p[i-1];
      } else {
	P[i]=C*pow(theta,i-alpha);
      }
    }

    { // initialization
      D[1][0]=-lambda;
      double aux=log(lambda*P[1]);
      for (int j=1; j<=n; j++) {
	D[1][j]=D[1][j-1]+aux-log((double)j);
	//printf("D[%i][%i]=%f\n",1,j,D[1][j]);
      }
    } // end initialization
    
    
    { // main loop
      for (int i=2; i<=n; i++) {
	double aux1=log(lambda*P[i]);
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux2=0.0;	 
	  for (int k=0; k<=kmax; k++) {
	    //printf("k=%i\texp=%f\n",k,aux2+D[i-1][j-i*k]-D[i][j-1]);
	    D[i][j]+=exp(aux2+D[i-1][j-i*k]-D[i][j-1]);
	    aux2+=aux1-log(k+1.0);
	  } // end k loop
	  D[i][j]=log(D[i][j])+D[i][j-1];
	  //printf("D[%i][%i]=%f\tD[%i][%i]=%f\n",i,j-1,D[i][j-1],i,j,D[i][j]);
	} // end j loop
      } // end i loop
    } // end main loop

    //printf("log D[alpha][n]=%f\n",D[alpha][n]);

    // get result
    double res=D[n][n];

    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;

  }

  double fast_log_dcpoi(int n,double lambda,double theta, int alpha,double *p){
    
    //printf("dcpoi call\n");
    //timer T;
    //T.start();

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation
    //printf("C=%e\n",C);

    double **D;
    double *CC;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(alpha+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // alpha+2 (and not alpha+1) for C
      D[0]=(double *)malloc(sizeof(double )*(alpha+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=alpha; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      CC=pos;
    } // end memory allocation

    { // initialization
      D[1][0]=-lambda;
      double aux=log(p[0]*lambda);
      for (int j=1; j<=n; j++) {
	D[0][j]=D[0][j-1]+aux-log((double)j);
      }
    } // end initialization
    
    { // precompute all necessary C[j]
      // allocate memory
      double *workspace=(double *)malloc(sizeof(double)*2*(n/(alpha+1)));
      if (workspace==NULL) {
	fprintf(stderr,"not enough memory in dcpoi. Aborting.\n");
	exit(EXIT_FAILURE);
      }
      double *X,*Y;
      X=workspace; // for log (lambda*C)^m/m!*theta^(alpha+j-m*alpha)
      Y=workspace+(n/(alpha+1)); // for log binomial(alpha+j-m*alpha-1,m-1)
      int mmax=1;
      X[0]=log(lambda*C);
      Y[0]=log(1.0);
      double aux=-(double)alpha*log(theta);
      //printf("lambda*C=%e\ttheta=%e\talpha=%i\ttheta^{-alpha}=%e\n",X[0],theta,alpha,aux);
      // main loop
      CC[1]=log(lambda*C*theta);
      X[0]+=log(theta);
      for (int j=2; j<=n-alpha; j++) {
	//printf("j=%i\n",j);
	CC[j]=0.0;
	// update X and Y
	X[0]+=log(theta);
	CC[j]=exp(X[0]+Y[0]);
	for (int m=1; m<mmax; m++) { // beware, replace m by m+1 in theoretical formulas
	  X[m]+=log(theta);
	  Y[m]+=log( (double)(j-alpha*m-1)/(double)(j-alpha*m-m-1) );
	  CC[j]+=exp(X[m]+Y[m]-CC[j-1]);
	}
	// check mmax
	int newmmax=(alpha+j)/(alpha+1);
	if (newmmax>mmax) {
	  //printf("newmmax=%i\tmmax=%i\n",newmmax,mmax);
	  mmax=newmmax;
	  X[mmax-1]=X[mmax-2]+log(lambda*C/(double)(mmax))+aux;
	  Y[mmax-1]=log(1.0);
	  // update CC[j]
	  CC[j]+=exp(X[mmax-1]+Y[mmax-1]-CC[j-1]);
	}
	CC[j]=log(CC[j])+CC[j-1];
	// print stuff
	//printf("j=%i 1<=m<=%i:\n",j,mmax);
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("X(%i)=%e\t",m+1,exp(X[m]));
	//}
	//printf("\n");
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("Y(%i)=%e\t",m+1,exp(Y[m]));
	//}
	//printf("\n");
	//printf("C(alpha,%i)=%e\n",j,exp(CC[j]));
	//printf("log C(alpha,%i)=%f\n",j,CC[j]);
      } // end j loop;
      // free memory
      free(workspace);
    } // end precompute
    //printf("%.2fs to compute all C(a,j)\n",T.elapsed_time()); 

    { // main loop
      for (int i=2; i<=alpha; i++) {
	double aux1=log(lambda*p[i-1]);
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux2=0.0;
	  for (int k=0; k<=kmax; k++) {
	    D[i][j]+=exp(aux2+D[i-1][j-i*k]-D[i][j-1]);
	    aux2+=aux1-log(k+1.0);
	  } // end k loop
	  D[i][j]=log(D[i][j])+D[i][j-1];
	} // end j loop
      } // end i loop
    } // end main loop
    //printf("%.2fs to compute all D(i,j)\n",T.elapsed_time()); 

    //printf("log D[alpha][n]=%f\n",D[alpha][n]);

    // get result    
    double res;
    {
      double max=D[alpha][n];
      double aux;
      for (int j=1; j<=n-alpha; j++) {
	aux=CC[j]+D[alpha][n-j-alpha];
	if (aux>max)
	  max=aux;
      } // end j loop
      //printf("max=%f\n",max);
      //printf("aux[0]=%f\n",D[alpha][n]-max);
      //for (int j=1; j<=n-alpha; j++) {
      //	aux=CC[j]+D[alpha][n-j-alpha];
      //	printf("aux[%i]=%f\n",j,aux-max);
      //}
      res=exp(D[alpha][n]-max);
      for (int j=1; j<=n-alpha; j++) {
	res+=exp(CC[j]+D[alpha][n-j-alpha]-max);
      } // end j loop
      res=log(res)+max;
    } // end get result

    //printf("\n");
    //printf("%.2fs to compute the final result\n",T.elapsed_time()); 

    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    free(D[0]);
    free(D);

    //printf("%.2fs global computional time\n",T.elapsed_time()); 

    // return result
    return res;

  }

  double barbour_pcpoi(int n,double lambda,double theta, int alpha,double *p) {

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation

    double **D;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(n+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      D[0]=(double *)malloc(sizeof(double )*(n+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=0; i<=n; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
    } // end memory allocation

    { // initialization
      for (int j=1; j<=n; j++) {
	D[j][j]=1.0/(double)j;
	// put lambdaj in D[0][j]
	D[0][j]=0.0;
	if (j<=alpha) {
	  D[0][j]=lambda*p[j-1];
	} else {
	  D[0][j]=lambda*C*pow(theta,j-alpha);
	}
      }
    } // end initialization
    
    { // main loop
      for (int i=1; i<=n; i++) {
	for (int j=1; j<i; j++) {
	  //printf("C_{%i,%i}:\n",i,j);
	  // compute c_{i,j}
	  D[i][j]=0.0;
	  for (int k=1; k<=i-j; k++) {
	    //printf("%i * lambda_%i / %i * C_{%i,%i}\n",k,k,i,i-k,j);
	    D[i][j]+=k*D[0][k]/(double)i*D[i-k][j];
	  } // end k loop
	} // end loop on j
      } // end loop on i
    } // end main loop

    // get result
    double res=0.0;
    if (n==0) {
      res=1.0;
    } else {
      double sum=1.0;
      for (int i=1; i<=n; i++) {
	res=0.0;
	for (int k=1; k<=i; k++) {
	  res+=k*D[0][k]*D[i][k];
	}
	sum+=res;
      }
      res=sum*exp(-lambda);
    }
    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;
  }

  double complete_pcpoi(int n,double lambda,double theta, int alpha,double *p){

    
    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation

    double **D;
    double *P;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(n+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // n+2 (and not n+1) for P
      D[0]=(double *)malloc(sizeof(double )*(n+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=n; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      P=pos;
    } // end memory allocation

    // precompute all P[i]
    for (int i=1; i<=n; i++) {
      if (i<=alpha) {
	P[i]=p[i-1];
      } else {
	P[i]=C*pow(theta,i-alpha);
      }
    }

    { // initialization
      double aux1=P[1]*lambda;
      double aux2=exp(-lambda);
      double aux3=1.0;
      D[1][0]=aux2;
      double m,e;
      for (int j=1; j<=n; j++) {
	aux3*=aux1/(double)j;
	D[1][j]=D[1][j-1]+aux2*aux3;
      }
    } // end initialization
    
    
    { // main loop
      //printf("lambda*P[1]=%e\n",lambda*P[1]);
      for (int i=2; i<=n; i++) {
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux=1.0;
	  for (int k=0; k<=kmax; k++) {
	    D[i][j]+=aux*D[i-1][j-i*k];
	    aux*=lambda*P[i]/(k+1.0);
	  } // end k loop
	} // end j loop
      } // end i loop
    } // end main loop

    // get result
    double res=D[n][n];

    //printf("complete_dcpoi:\n");
    //for (int j=1; j<=n-alpha; j++) { 
    //  printf("%e\t",D[alpha][n-alpha-j]);
    //}
    //printf("\n");
    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;

  }

  double fast_pcpoi(int n,double lambda,double theta, int alpha,double *p){
    
    //printf("fast_dcpoi call\n");
    //timer T;
    //T.start();

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation
    //printf("C=%e\n",C);

    double **D;
    double *CC;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(alpha+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // alpha+2 (and not alpha+1) for C
      D[0]=(double *)malloc(sizeof(double )*(alpha+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=alpha; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      CC=pos;
    } // end memory allocation

    { // initialization
      double aux1=p[0]*lambda;
      double aux2=exp(-lambda);
      double aux3=1.0;
      D[1][0]=aux2;
      double m,e;
      for (int j=1; j<=n; j++) {
	aux3*=aux1/(double)j;
	D[1][j]=D[1][j-1]+aux2*aux3;
      }
    } // end initialization
    
    { // precompute all necessary C[j]
      // allocate memory
      double *workspace=(double *)malloc(sizeof(double)*2*(n/(alpha+1)));
      if (workspace==NULL) {
	fprintf(stderr,"not enough memory in dcpoi. Aborting.\n");
	exit(EXIT_FAILURE);
      }
      double *X,*Y;
      X=workspace; // for (lambda*C)^m/m!*theta^(alpha+j-m*alpha)
      Y=workspace+(n/(alpha+1)); // for binomial(alpha+j-m*alpha-1,m-1)
      int mmax=1;
      X[0]=lambda*C;
      Y[0]=1.0;
      double aux=pow(theta,-(double)alpha);
      //printf("lambda*C=%e\ttheta=%e\talpha=%i\ttheta^{-alpha}=%e\n",X[0],theta,alpha,aux);
      // main loop
      CC[1]=lambda*C*theta;
      X[0]*=theta;
      for (int j=2; j<=n-alpha; j++) {
	//printf("j=%i\n",j);
	// initialization
	CC[j]=0.0;
	// update X and Y
	X[0]*=theta;
	CC[j]+=X[0]*Y[0];
	for (int m=1; m<mmax; m++) { // beware, replace m by m+1 in theoretical formulas
	  X[m]*=theta;
	  Y[m]*=( (double)(j-alpha*m-1)/(double)(j-alpha*m-m-1) );
	  CC[j]+=X[m]*Y[m];
	}
	// check mmax
	int newmmax=(alpha+j)/(alpha+1);
	if (newmmax>mmax) {
	  //printf("newmmax=%i\tmmax=%i\n",newmmax,mmax);
	  mmax=newmmax;
	  X[mmax-1]=X[mmax-2]*lambda*C/(double)(mmax)*aux;
	  Y[mmax-1]=1.0;
	  // update CC[j]
	  CC[j]+=X[mmax-1]*Y[mmax-1];
	}
	// print stuff
	//printf("j=%i 1<=m<=%i:\n",j,mmax);
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("X(%i)=%e\t",m+1,X[m]);
	//}
	//printf("\n");
	//printf("\t");
	//for (int m=0; m<mmax; m++) {
	//  printf("Y(%i)=%e\t",m+1,Y[m]);
	//}
	//printf("\n");
	//printf("C(alpha,%i)=%e\n",j,CC[j]);
      } // end j loop;
      // free memory
      free(workspace);
    } // end precompute
    //printf("%.2fs to compute all C(a,j)\n",T.elapsed_time()); 

    { // main loop
      for (int i=2; i<=alpha; i++) {
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux=1.0;
	  for (int k=0; k<=kmax; k++) {
	    //D[i][j]+=pow(lambda*p[i-1],k)/factrl(k)*D[i-1][j-i*k];
	    D[i][j]+=aux*D[i-1][j-i*k];
	    aux*=lambda*p[i-1]/(k+1);
	  } // end k loop
	  //if (D[i][j]==0.0) {
	  //  warning=true;
	  //}
	} // end j loop
      } // end i loop
    } // end main loop
    //printf("%.2fs to compute all D(i,j)\n",T.elapsed_time()); 

    // get result
    double res=D[alpha][n];
//    printf("fast_dcpoi:\n");
//    for (int j=1; j<=n-alpha; j++) { 
//      printf("%e\t",D[alpha][n-alpha-j]);
//    }
//    printf("\n");
    for (int j=1; j<=n-alpha; j++) {
      res+=D[alpha][n-j-alpha]*CC[j];
      //printf("%e*%e\t",D[alpha][n-alpha-j],CC[j]);
    } // end j loop
    //printf("\n");
    //printf("%.2fs to compute the final result\n",T.elapsed_time()); 

    //for (int i=0; i<=n; i++) {
    //  for (int j=0; j<=n; j++) {
    //	printf("D(%i,%i)=%f\t",i,j,D[i][j]);
    //  }
    //  printf("\n");
    //}

    // free memory
    free(D[0]);
    free(D);


    //printf("%.2fs global computional time\n",T.elapsed_time()); 

    // return result
    return res;

  }

  double complete_qcpoi(int n,double lambda,double theta, int alpha,double *p){
    if (n<1)
      return 1.0;
    else
      return 1.0-complete_pcpoi(n-1,lambda,theta,alpha,p);
  };
  
  double barbour_qcpoi(int n,double lambda,double theta, int alpha,double *p){
    if (n<1)
      return 1.0;
    else
      return 1.0-barbour_pcpoi(n-1,lambda,theta,alpha,p);
  };

  double fast_qcpoi(int n,double lambda,double theta, int alpha,double *p){

    if (n<0 || lambda<0.0 || theta>=1.0 || theta<0.0 || alpha<0) {
      fprintf(stderr,"slow_dcpoi: check parameters !\n");
      exit(EXIT_FAILURE);
    }

    double C=0.0;
    { // compute C so that sum lambda_k = 1.0
      // compute sum p_i
      double sum_p=0.0;
      for (int i=0; i<alpha; i++) {
	sum_p+=p[i];
      }
      if (sum_p>1 || sum_p<0) {
	fprintf(stderr,"slow_dcpoi: sum out of range, check parameters !\n");
	exit(EXIT_FAILURE);
      }
      if (theta>0) {
	C=(1.0-sum_p)*(1.0-theta)/theta;
      }
    } // end C computation
    //printf("C=%e\n",C);

    double **D;
    double *CC;
    { // memory allocation
      D=(double **)malloc(sizeof(double *)*(alpha+1));
      if (D==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // alpha+2 (and not alpha+1) for C
      D[0]=(double *)malloc(sizeof(double )*(alpha+1)*(n+1));
      if (D[0]==NULL) {
	fprintf(stderr,"slow_dcpoi: not enough memory !\n");
	exit(EXIT_FAILURE);
      }
      // connect
      double *pos=D[0];
      for (int i=1; i<=alpha; i++) {
	D[i]=pos;
	pos+=(n+1);
      }
      CC=pos;
    } // end memory allocation

    { // initialization
      double aux1=p[0]*lambda;
      double aux2=exp(aux1-lambda);
      D[1][0]=aux2;
      double m,e;
      for (int j=1; j<=n; j++) {
	D[1][j]=exp(aux1-lambda+qcdfpoi(j,aux1,&m,&e));
      }
    } // end initialization
    
    { // precompute all necessary C[j]
      // allocate memory
      double *workspace=(double *)malloc(sizeof(double)*2*(n/(alpha+1)));
      if (workspace==NULL) {
	fprintf(stderr,"not enough memory in dcpoi. Aborting.\n");
	exit(EXIT_FAILURE);
      }
      double *X,*Y;
      X=workspace; // for (lambda*C)^m/m!*theta^(alpha+j-m*alpha)
      Y=workspace+(n/(alpha+1)); // for binomial(alpha+j-m*alpha-1,m-1)
      int mmax=1;
      X[0]=lambda*C;
      Y[0]=1.0;
      double aux=pow(theta,-(double)alpha);
      // main loop
      CC[1]=lambda*C*theta;
      X[0]*=theta;
      for (int j=2; j<=n-alpha; j++) {
	//printf("j=%i\n",j);
	// initialization
	CC[j]=0.0;
	// update X and Y
	X[0]*=theta;
	CC[j]+=X[0]*Y[0];
	for (int m=1; m<mmax; m++) { // beware, replace m by m+1 in theoretical formulas
	  X[m]*=theta;
	  Y[m]*=( (double)(j-alpha*m-1)/(double)(j-alpha*m-m-1) );
	  CC[j]+=X[m]*Y[m];
	}
	// check mmax
	int newmmax=(alpha+j)/(alpha+1);
	if (newmmax>mmax) {
	  //printf("newmmax=%i\tmmax=%i\n",newmmax,mmax);
	  mmax=newmmax;
	  X[mmax-1]=X[mmax-2]*lambda*C/(double)(mmax)*aux;
	  Y[mmax-1]=1.0;
	  // update CC[j]
	  CC[j]+=X[mmax-1]*Y[mmax-1];
	}
      } // end j loop;
      // free memory
      free(workspace);
    } // end precompute

    { // main loop
      for (int i=2; i<=alpha; i++) {
	D[i][0]=D[1][0];
	for (int j=1; j<=n; j++) {
	  int kmax=j/i;
	  //printf("E(%i/%i)=%i\n",j,i,kmax);
	  D[i][j]=0.0;
	  double aux=1.0;
	  for (int k=0; k<=kmax; k++) {
	    D[i][j]+=aux*D[i-1][j-i*k];
	    aux*=lambda*p[i-1]/(k+1);
	  } // end k loop
	} // end j loop
      } // end i loop
    } // end main loop

    double Csum=0.0;
    { // compute sum_{j=n-alpha}^\infty C(alpha,j)
      int mmax=n/(alpha+1);
      { // the easy part
	double m,e;	
	Csum=exp(lambda*C*theta/(1.0-theta)+qcdfpoi(mmax+1,lambda*C*theta/(1.0-theta),&m,&e));
      } // end easy part
      // fixme: H should be computed on demand; same time complexity but O(1) in memory
      double **H;
      { // allocate memory for H
	H=(double **)malloc(sizeof(double *)*(mmax+1));
	if (H==NULL) {
	  fprintf(stderr,"Not enough memory in fast_qcpoi. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	H[0]=(double *)malloc(sizeof(double)*( mmax*(n-alpha) ));
	if (H[0]==NULL) {
	  fprintf(stderr,"Not enough memory in fast_qcpoi. Aborting.\n");
	  exit(EXIT_FAILURE);
	}
	{ // connecting
	  double *dpos=H[0];
	  for (int i=1; i<=mmax; i++) {
	    H[i]=dpos;
	    dpos+=(n-alpha);
	  } // end i loop
	} // end connect
      } // end memory allocation
      { // compute all H
	// simple terms first
	double aux1=1.0/(1.0-theta);
	double aux2=aux1*aux1;
	for (int q=0; q<=n-alpha-1; q++) {
	  H[1][q]=aux1;
	  H[2][q]=(q*(1.0-theta)+1.0)/(q+1.0)*aux2;
	} // end simple terms
	{ // recurrence part
	  for (int m=3; m<=mmax; m++) {
	    for (int q=0; q<=n-alpha-1; q++) {
	      H[m][q]=(q+2*m-3-(q+m-2)*theta)*H[m-1][q]+(2-m)*H[m-2][q];
	      H[m][q]/=(q+m-1)*(1.0-theta);
	    } // end q loop
	  } // end m loop
	} // end recurrence part
      } // end H computation
      // verif
      //printf("H(%i,%i)=%e\n",n-alpha-1,mmax,H[mmax][n-alpha-1]);
      //printf("Csum=%e\n",Csum);
      { // compute the sum
	double aux1=lambda*C*pow(theta,(double)(n-alpha));
	double aux2=log(lambda*C)-alpha*log(theta);
	Csum+=aux1*H[1][n-alpha-1];
	for (int m=2; m<=mmax; m++) {
	  int aux3=m*alpha;
	  aux1*=exp(aux2-log((m-1.0)*m)
		    +gammln((double)(n-aux3+alpha-m+2))
		    +gammln((double)(n-aux3))
		    -gammln((double)(n-aux3+alpha))
		    -gammln((double)(n-aux3-m+1)));
	  //printf("m=%i\taux1=%e\tH(%i,%i)=%e\n",m,aux1,n-aux3-m,m,H[m][n-aux3-m]);
	  Csum+=aux1*H[m][n-aux3-m];
	} // end m loop
      } // end sum
      // free H
      free(H[0]);
      free(H);
    } // end sum C(alpha,j) computation
    
    // get result
    //printf("\n***\nQ(alpha,n)=%e\tQ(alpha,0)=%e\tCsum=%e\n***\n",D[alpha][n],D[alpha][0],Csum);
    double res=D[alpha][n];
    for (int j=1; j<=n-alpha-1; j++) {
      res+=D[alpha][n-j-alpha]*CC[j];
    } // end j loop
    res+=D[alpha][0]*Csum;

    // free memory
    free(D[0]);
    free(D);

    // return result
    return res;
    
  };

};
