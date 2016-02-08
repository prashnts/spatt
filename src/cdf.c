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
/*  Provides support for Poisson and Gaussian cdf with special	     */
/*  care for very small p-value (log of p-value are given here)	     */
/*  								     */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "cdf.h"

double gammp(double a, double x){
  double gamser,gammcf,gln;

  /* test arguments */
  if (x<0.0 || a<=0.0) {
    fprintf(stderr,"Critical error in gammp ! argument out of range.\n");
    exit(EXIT_FAILURE);
  } else {    
    if (x < (a+1.0)) {
      /* use serie representation */
      gser(&gamser,a,x,&gln);
      return gamser;
    } else {
      /* use continued fraction */
      gcf(&gammcf,a,x,&gln);
      /* and returns complement */
      return 1.0-gammcf;
    }
  }
  
};

double gammq(double a, double x){
  double gamser,gammcf,gln;

  /* test arguments */
  if (x<0.0 || a<=0.0) {
    fprintf(stderr,"Critical error in gammq ! argument out of range.\n");
    exit(EXIT_FAILURE);
  } else {
    if (x < (a+1.0)) {
      /* use serie representation */
      gser(&gamser,a,x,&gln);
      /* and returns complement */
      return 1.0-gamser;
    } else {
      /* use continued fraction */
      gcf(&gammcf,a,x,&gln);
      return gammcf;
    }
  }
  
};

void gser(double *gamser, double a, double x, double *gln) {
  int n;
  double sum,del,ap;

  *gln=gammln(a);
  if (x<=0) {
    fprintf(stderr,"Critical error in gser ! argument out of range.\n");
    fprintf(stderr,"gser(%f,%x)\n",a,x);
    exit(EXIT_FAILURE);
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    fprintf(stderr,"Warning in gser ! argument a too large and/or ITMAX too small.\n");
    fprintf(stderr,"failed test: %e < %e ?\n",fabs(del),fabs(sum)*EPS);
    //*gamser=0.0;
    *gamser=sum*exp(-x+a*log(x)-(*gln));
    return;
  }

};

void gcf(double *gammcf, double a, double x, double *gln) {

  int i;
  double an,b,c,d,del,h;

  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1; i<=ITMAX; i++){
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN)
      d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN)
      c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i>ITMAX) {
    fprintf(stderr,"Warning in gcf ! argument a too large and/or ITMAX too small.\n");
    fprintf(stderr,"failed test: %e < %e ?\n",fabs(del-1.0),EPS);
    //*gammcf=0.0;
  }
  *gammcf=exp(-x+a*log(x)-(*gln))*h;

};

double gammln(double xx){
  double x,y,tmp,ser;
  static double cof[6]={
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  };
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp-= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0; j<6; j++)
    ser+= cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);

};

double pcdfpoi(double k,double x,double *m,double *e){
  double res;
  long le;
  int i;
  double add,num,del;

  //printf("pcdfpoi !\n");
  res=gammq(k+1.0,x);
  //printf("gammq(%e,%e)=%e\n",k+1,x,res);
  if (res==0.0) {
    /* log computations here */
    //fprintf(stderr,"Warning in pcdfpoi ! log computations not yet implemented.\n");
    if (k<=32)
      res=-x+k*log(x)-log(factrl((int)k));
    else
      res=-x+k*log(x)-gammln(k+1);
    add=1.0;
    del=1.0;
    i=0;
    while (i<k && del>EPS) {
      del*=(k-i)/x;
      //printf("i=%i: add=%e del=%e\n",i,add,del);
      add+=del;
      i++;
    }
    //printf("end loop: i=%i k=%f\n",i,k);
    res+=log(add);
  } else {
    res=log(res);
  }
  *e=res/log(10);
  le=(long)*e;
  if ( (le-*e) > 0) 
    le--;
  *e=le;
  *m=exp(res-(*e)*log(10));
  return res;
};

double qcdfpoi(double k,double x,double *m,double *e){
  double res;
  long le;
  double add,del,num,den;

  if (k==0.0) {
    res=0.0;
    *m=1.0;
    *e=1.0;
  } else {
    //printf("gammp(%e,%e)\n",k,x);
    res=gammp(k,x);
    if (res==0.0) {
      /* log computations here */
      //fprintf(stderr,"Warning in qcdfpoi ! log computations not yet implemented.\n");
      if (k<=32)
	res=-x+(k-1)*log(x)-log(factrl((int)(k-1)));
      else
	res=-x+(k-1)*log(x)-gammln(k);
      del=1.0;
      add=1;
      num=x;
      den=k;
      //printf("del=%f\n",res);
      while (del>EPS) {
	del=num/den;
	add+=del;
	num*=x;
	den*=den+1;	
      }
      //printf("add=%f\n",add);      
      res+=log(add);
    } else {
      res=log(res);
    }
  }
  *e=res/log(10);
  le=(long)*e;
  if ( (le-*e) > 0) 
    le--;
  *e=le;
  *m=exp(res-(*e)*log(10));
  return res;
};

double factrl(int n){
  static int ntop=4;
  static double a[33]={1.0,1.0,2.0,6.0,24.0};
  int j;

  if (n<0) {
    fprintf(stderr,"Critical error in factrl ! argument out of range.\n");
    exit(EXIT_FAILURE);
  } else {
    if (n>32)
      return exp(gammln(n+1.0));
    while (ntop<n) {
      j=ntop++;
      a[ntop]=a[j]*ntop;
    }
    return a[n];
  }
};

double betai(double a, double b, double x){
  double bt;

  if (x<0.0 || x>1.0) {
    fprintf(stderr,"Critical error in betai ! 3rd argument out of range.\n");
    fprintf(stderr,"betai(%f,%f,%f)\n",a,b,x);
    exit(EXIT_FAILURE);
  }

  if (x==0.0 || x==1.0)
    bt=0.0;
  else 
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x<(a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
};

double betacf(double a, double b, double x){

  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  
  //printf("betacf(%f,%f,%e)\n",a,b,x);
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  //printf("d=%e\n",d);
  if (fabs(d) < FPMIN)
    d=FPMIN;
  d=1.0/d;
  h=d;
  //printf("h=%f\n",h);
  for (m=1; m<MAXIT; m++){
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    //printf("aa=%e\n",aa);
    d=1.0+aa*d;
    if (fabs(d)<FPMIN)
      d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c)<FPMIN)
      c=FPMIN;
    d=1.0/d;
    h*=d*c;
    //printf("h=%f\n",h);
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    //printf("aa=%e\n",aa);
    d=1.0+aa*d;
    if (fabs(d)<FPMIN)
      d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c)<FPMIN)
      c=FPMIN;
    d=1.0/d;
    del=d*c;
    h*=del;
    //printf("h=%f\n",h);
    //printf("test(m=%i): %e < %e ?\n",m,fabs(del-1.0),EPS);
    if (fabs(del-1.0)<EPS) break;
  }
  if (m>MAXIT) {
    fprintf(stderr,"a or b too big, or MAXIT too small in betacf\n");
  }
  //printf("res=%f\n",h);
  return h;
};

double pcdfbin(double k,long n,double p,double *m,double *e){
  double res;
  long le;

  //printf("pcdfbin !\n");

  //printf("pcdfbin(%f,%li,%e)\n",k,n,p);
  //printf("betai=%f\n",betai(k+1,n-(k+1)+1,p));
  // res=1.0-betai(k+1,n-(k+1)+1,p);
  // thanks to betai(a,b,x)=1-betai(b,a,1.0-x)
  res=betai(n-(k+1)+1,k+1,1.0-p);
  //printf("1.0-betai(%f,%f,%e)=%e\n",k+1,n-(k+1)+1,p,res);
  //if (res<=0.0) {
  if (res==0.0) {
    /* use Poisson approximation */
    res=pcdfpoi(k,n*p,m,e);
  } else {
    res=log(res);
    *e=res/log(10);
    le=(long)*e;
    if ( (le-*e) > 0) 
      le--;
    *e=le;
    *m=exp(res-(*e)*log(10));    
  }
  return res;
};

double qcdfbin(double k,long n,double p,double *m,double *e){
  double res;
  long le;

  //printf("qcdfbin(%f,%li,%e)\n",k,n,p);
  res=betai(k,n-k+1,p);
  //printf("betai(%f,%f,%e)=%e\n",k,n-k+1,p,res);
  //if (res==0.0 || res>1.0) {
  if (res==0.0) {
    /* use Poisson approximation */
    res=qcdfpoi(k,n*p,m,e);
  } else {
    res=log(res);
    *e=res/log(10);
    le=(long)*e;
    if ( (le-*e) > 0) 
      le--;
    *e=le;
    *m=exp(res-(*e)*log(10));    
  }
  return res;
};

/* returns mantiss (in m) and exposant (in e) */
/* of the P(X<=x) with X Gaussian N(0,1) */
/* returns also log of that value */
double pcdfnor(double x,double *m,double *e){
  
  double p;
  double lp;

  if (x<-MAGNITUDE_LIMIT) {
    /* we use the tail approximation */
    lp=1/log(10)*(log(1/sqrt(2*PI)/(-x))-x*x/2);
    *e=floor(lp);
    *m=exp(log(1/sqrt(2*PI)/(-x))-x*x/2-*e*log(10));
    //lp=*e+log(*m)/log(10);
  } else {
    /* we use the GSL */
    p=gsl_cdf_gaussian_P(x,1.0);
    lp=log(p)/log(10);
    *e=floor(lp);
    *m=p/exp(*e * log(10) );
  }

  return lp;
};

/* returns mantiss (in m) and exposant (in e) */
/* of the P(X>=x) with X Gaussian N(0,1) */
/* returns also log of that value */
double qcdfnor(double x,double *m,double *e){

  double p;
  double lp;

  if (x>MAGNITUDE_LIMIT) {
    /* we use the tail approximation */
    lp=1/log(10)*(log(1/sqrt(2*PI)/x)-x*x/2);
    *e=floor(lp);
    *m=exp(log(1/sqrt(2*PI)/x)-x*x/2-*e*log(10));
    //lp=*e+log(*m)/log(10);
  } else {
    /* we use the GSL */
    p=gsl_cdf_gaussian_Q(x,1.0);
    lp=log(p)/log(10);
    *e=floor(lp);
    *m=p/exp(*e * log(10) );
  }

  return lp;
};
