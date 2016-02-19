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
#include "xwaiting.h"

using namespace std;

namespace spatt {

xwaiting::xwaiting(dfa &D,pmc &Y,sequence &S,int rep,bool verbose,bool debug){

  if (verbose)
    printf(">>> call xstat::xstat\n");
  _rep=rep;
  if ((rep!=OVER)&&(rep!=UNDER)) {
    fprintf(stderr,"xwaiting::xwaiting : wrong value for rep. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _Dfa=&D;
  _Y=&Y;
  _Seq=&S;
  _m=_Y->_m;

  { // parse the sequences
    _Seq->reset();    
    int c;
    unsigned long token=_Dfa->reset();
    unsigned long Ym=_Y->recode(token);
    unsigned long i=0;
    _pos.push_back(_m);
    _final.push_back(false);
    while ((c=_Seq->next())!=SEQEMPTY) {
      //printf("c=%i\n",c);
      if (c>=0) { // valid char
	token=_Dfa->process(c);
	if (i==_m)
	  _obs.push_back(Ym);
	if (i++<_m) {
	  Ym=_Y->recode(token);
	}	
	if (_Dfa->is_final(token)) {
	  _pos.push_back(_Seq->nvalidchar());
	  _obs.push_back(_Y->recode(token));
	  _final.push_back(true);
	}
      } else { // sequence end
	i=0;
	token=_Dfa->reset();
	Ym=_Y->recode(token);
	_pos.push_back(_Seq->nvalidchar());
	_final.push_back(false);
      }
    }
    _obs.push_back(0);
  } // end parse the sequence

  if (verbose) {
    printf("observations(position;state;final)=");
    for (unsigned long i=0; i<_pos.size(); i++)
      printf("(%i;%i;%i)",_pos[i],_obs[i],(short)_final[i]);
    printf("\n");
  }
    

  precompute(debug);
  compute(debug);

};

void xwaiting::precompute(bool debug) {
  if (debug)
    printf(">>> call xwaiting::precompute()\n");

  // _u[0]=SQ
  _u.push_back(_Y->_SQ);
  bool all_positive=false;
  while (!all_positive) {
    // aux=P*last vector
    vector<double> aux(_Y->_nstates);
    _Y->zero(aux);
    _Y->add_Mx(_Y->_Pt,_u.back(),aux);
    { // check if aux >> 0
      //printf("aux=[ ");
      //for (unsigned long i=0; i<_Y->_nstates; i++)
      //	printf("%e ",aux[i]);
      //printf("]\n");      
      unsigned long i=0;
      while (i<_Y->_nstates) {
	if (!(aux[i]>0.0))
	  break;
	i++;
      }
      if (i==_Y->_nstates)
	all_positive=true;      
    } // end check
    // add aux to _u
    _u.push_back(aux);
  }
  //printf("ok\n");
  bool converged=false;
  while (!converged) {
    // aux=P*last vector
    vector<double> aux(_Y->_nstates);
    _Y->zero(aux);
    _Y->add_Mx(_Y->_Pt,_u.back(),aux);
    { // check if aux = lambda * last vector
      //printf("aux/last=[ ");
      //for (unsigned long i=0; i<_Y->_nstates; i++)
      //	printf("%e/%e ",aux[i],_u.back()[i]);
      //printf("]\n");      
      _lambda=aux[0]/_u.back()[0];
      unsigned long i=1;      
      while (i<_Y->_nstates) {
	if (fabs((_lambda-aux[i]/_u.back()[i])/_lambda)>TOL)
	  break;
	i++;
      }
      if (i==_Y->_nstates)
	converged=true;
    } // end check
    // add aux to _u
    _u.push_back(aux);
  } // end while (!converged)

  _t0=_u.size()-1;
  
  if (debug)
    printf("_t0=%i\t_lambda=%e\n",_t0,_lambda);

};

void xwaiting::compute(bool debug) {
  if (debug)
    printf(">>> call xwaiting::compute()\n");

  _pvalue.push_back(-1.0);
  for (unsigned long i=1; i<_obs.size(); i++) {
    if (_final[i]) {
      if (_rep==UNDER)
	_pvalue.push_back(tail_tau(_obs[i-1],_pos[i]-_pos[i-1]-1));
      else 
	_pvalue.push_back(cumsum_tau(_obs[i-1],_pos[i]-_pos[i-1]));
    } else {
      _pvalue.push_back(tail_tau(_obs[i-1],_pos[i]-_pos[i-1]));
    }
  }
  if (debug) {
    printf("pvalue=[ ");
    for (vector<double>::iterator it=_pvalue.begin(); it!=_pvalue.end(); it++)
      printf("%e ",*it);
    printf("]\n");
  }
  
};

double xwaiting::cumsum_tau(unsigned long from,unsigned long t){
  
  double res=0.0;
  if (t<_t0) {
    for (unsigned long i=0; i<t; i++)
      res+=_u[i][from];
  } else {
    for (unsigned long i=0; i<_t0; i++)
      res+=_u[i][from];
    res+=(1.0-pow(_lambda,(double)(t-_t0)))*_u[_t0][from]/(1.0-_lambda);
  }
  return res;
};


double xwaiting::tail_tau(unsigned long from,unsigned long t) {
  if (t>=_t0) {
    // add geometric part
    return pow(_lambda,(double)(t-_t0))*_u[_t0][from]/(1.0-_lambda);
  } else {
    double res=0.0;
    // add non geometric part
    for (unsigned long i=t; i<_t0; i++)
      res+=_u[i][from];
    // add geometric part
    res+=_u[_t0][from]/(1.0-_lambda);
    return res;
  }
};

void xwaiting::print() {
  for (unsigned long i=1; i<_obs.size(); i++)
    printf("%i\t%e\t%i\n",_pos[i],_pvalue[i],(short)_final[i]);
};

};
