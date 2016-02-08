#ifndef AMARKOV_H
#define AMARKOV_H

#include <cstdio>
#include <string>
#include <cmath>
#include <cstring>

#define ABUFFER_SIZE 100
#define ACOMMENT_CHAR '#'
#define DELIM " \t\n,:;"
#define MASS_KMAX 100
#define MASS_TOL 1e-100
#define WORK_MARGIN 20
#define EPSILON 1e-15

#include "fortran.h"
#include "cdf.h"
#include "cp.h"

class amarkov {

 public:
  // constructor from an automaton file
  amarkov(const char *sfile,const char *mfile,bool verbose=true,bool debug=true);

  // store in (previously allocated) res
  // x * ( P + t*Q )
  void leftprod(double *x,double t,double *res);

  // store in (previously allocated) res
  // x * ( P + t*Q )
  void rightprod(double *x,double t,double *res);

  // return stationary probability of pattern
  double nexp(long n);

  // return the binomial statistic
  double bstat(long n,long nobs, double nexp);

  // return the compound poisson statistic
  double cpstat(long n,long nobs, double nexp,bool verbose=true,bool debug=false);

  // return the exact statistic
  //  double xstat(long n,long nobs, double nexp,bool verbose=true,bool debug=false);

  // destructor
  ~amarkov();

 private:
  int _k; // alphabet size
  int _m; // Markov order
  int _nstart; // number of starting states
  int *_start; // index of starting states
  double *_dstart; // distribution of starting states
  double *_stationary; // stationary distribution
  double _secondmag; // second magnitude eigenvalue
  int _nstate; // number of states
  int _nfinal; // number of final states
  int *_final; // final states
  int *_h; // min word size for each final state
  int _min_h; // min of _h 
  int _nz; // number of non zero terms in transition matrix
  int _Pnz; // number of non zero terms in P
  int _Qnz; // number of non zero terms in Q
  long _data_size; // number of allocated octets
  int *_data; // allocated memory segment
  // main transition pointers
  int *_from;
  int *_to;
  double *_proba;
  // same pointer for P
  int *_Pfrom;
  int *_Pto;
  double *_Pproba;
  // same pointers for Q
  int *_Qfrom;
  int *_Qto;
  double *_Qproba;

  // normalize _dstart and the transition matrix
  // overwritting _stationary in the process 
  void normalize();

  // compute stationary distribution
  void compute_stationary(bool verbose=true);

  // compute _h
  void compute_h(bool verbose=true);

  // test if a vector is null
  bool is_null(double *x,int size);

}; 


#endif
