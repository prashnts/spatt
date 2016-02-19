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
#include "pmc.h"

using namespace std;

namespace spatt {

pmc::pmc(dfa &A,vector<double> &markov_param,bool verbose,bool debug) {

  if (verbose)
    printf(">>> call pmc::pmc(dfa&,vector<double>)\n");


  // get statistics 
  _alphabet_size=A._alphabet_size;
  _m=A._m;
  _nstates=0;

  { // check markov_param size
    unsigned long n=1;
    for (unsigned short i=0; i<=_m; i++)
      n*=_alphabet_size;
    if (markov_param.size()!=n) {
      fprintf(stderr,"pmc::pmc : size of markov_param does not match. Aborting._n");
      exit(EXIT_FAILURE);
    }
  } // end check markov_param

  // set _recode table 
  if (_m==0) {
    // get all states
    for (long p=0; p<A._nstates; p++) {
      _recode.push_back(p+1);
      _nstates++;
    }
  } else {
    // get rid of states with no past
    long number=1;
    for (long p=0; p<A._nstates; p++) {
      string s=A.inv_delta(_m,p);
      if (debug) {
	printf("inv_delta(%i,%i)=",_m,p);
	A.print_string(s);
	printf("\n");
      }
      if (!s.empty()) {
	_recode.push_back(number++);
	_nstates++;
      } else {
	_recode.push_back(0);
	if (debug)
	  printf("state %i with no %i-past has been removed\n",p,_m);
      }
    }    
  }

  { // fill _Pt and _Qt
    unsigned long dest;
    unsigned long c;
    string sc;
    term t;
    for (unsigned long p=0; p<A._nstates; p++)
      if (_recode[p]!=0) {
	if (_m>0) {
	  sc=(A._inv_delta[p].begin())->first;
	  c=A.code(sc);
	} else {
	  c=0;
	}
	for (unsigned short a=0; a<_alphabet_size; a++) {
	  dest=A._delta[a][p];
	  t.from=_recode[p]-1;
	  t.to=_recode[dest]-1;
	  //t.proba=c*_alphabet_size+a;
	  t.index=c*_alphabet_size+a;
	  t.proba=markov_param[c*_alphabet_size+a];
	  if (A._is_final[dest]) 
	    _Qt.push_back(t);
	  else
	    _Pt.push_back(t);	    
	}
      }	
  } // end fill

  { // compute _SQ
    for (unsigned long i=0; i<_nstates; i++)
      _SQ.push_back(0.0);
    for (vector<term>::iterator it=_Qt.begin(); it!=_Qt.end(); it++)
      _SQ[it->from]+=it->proba;
  } // end compute _SQ


  { // fill _start
    _start=A.starts(debug);
    for (vector<unsigned long>::iterator it=_start.begin(); it!=_start.end(); it++)
      *it=_recode[*it]-1;
    if (_start.size()==_nstates)
      fprintf(stderr,"pmc::pmc(dfa&) : Warning, degenerated case detected ! Try to lower the Markov order or to consider a longer pattern.\n");
    // fill _final
    _final=A._final;
    for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++)
      *it=_recode[*it]-1;    
  } // end fill

  if (verbose)
    print();
  
};

pmc::~pmc() {
  
};

void pmc::print() {

  printf("m=%i\tnstates=%i\talphabet_size=%i\n",_m,_nstates,_alphabet_size);
  printf("starts = [ ");
  if (_start.size()<=20)
    for (vector<unsigned long>::iterator it=_start.begin(); it!=_start.end(); it++)
      printf("%i ",*it);
  else
    printf("more than 20 states");
  printf("]\n");
  printf("P (nz=%i) = ",_Pt.size());
  if (_Pt.size()<=20)
    for (vector<term>::iterator it=_Pt.begin(); it!=_Pt.end(); it++)
      printf("(%i,%i;%.2f) ",it->from,it->to,it->proba);
  else
    printf("more than 20 transitions");
  printf("\n");
  printf("Q (nz=%i) = ",_Qt.size());
  if (_Qt.size()<=20)
    for (vector<term>::iterator it=_Qt.begin(); it!=_Qt.end(); it++)
      printf("(%i,%i;%.2f) ",it->from,it->to,it->proba);
   else
     printf("more than 20 transitions");
  printf("\n");  
  printf("SQ = [ ");
  if (_SQ.size()<=20)
    for (vector<double>::iterator it=_SQ.begin(); it!=_SQ.end(); it++)
      printf("%.2f ",*it);
  else
    printf("more than 20 terms");
  printf("]\n");

};

void pmc::sci_export(const char *file) {
  
  FILE *out=NULL;
  out=fopen(file,"w");
  if (out==NULL) {
    fprintf(stderr,"pmc::sci_export(): Cannot write on \"%s\". Aborting.\n",file);
    exit(EXIT_FAILURE);
  }
  fprintf(out,"m=%i;\nL=%i;\nk=%i;\n",_m,_nstates,_alphabet_size);
  fprintf(out,"starts = [ ");
  for (vector<unsigned long>::iterator it=_start.begin(); it!=_start.end(); it++)
      fprintf(out,"%i ",*it+1);
  fprintf(out,"];\n");
  fprintf(out,"tmp = [\n");
  for (vector<term>::iterator it=_Pt.begin(); it!=_Pt.end(); it++)
    fprintf(out,"%i,%i,%e;\n",it->from+1,it->to+1,it->proba);
  fprintf(out,"];\n");
  fprintf(out,"P=sparse(tmp(:,1:2),tmp(:,3),[%i %i]);\n",_nstates,_nstates);
  fprintf(out,"tmp = [\n");
  for (vector<term>::iterator it=_Qt.begin(); it!=_Qt.end(); it++)
    fprintf(out,"%i,%i,%e;\n",it->from+1,it->to+1,it->proba);
  fprintf(out,"];\n");  
  fprintf(out,"Q=sparse(tmp(:,1:2),tmp(:,3),[%i %i]);\n",_nstates,_nstates);
  fprintf(out,"index=1:L; final=index(full(sum(Q,1))>0);\n");
  fclose(out);
};

  void pmc::indexed_export(const char *file) {
  
    FILE *out=NULL;
    out=fopen(file,"w");
    if (out==NULL) {
      fprintf(stderr,"pmc::indexed_export(): Cannot write on \"%s\". Aborting.\n",file);
      exit(EXIT_FAILURE);
    }
    fprintf(out,"# Markov indexed file generated by SPatt\n");
    fprintf(out,"# alphabet_size model_order\n");
    fprintf(out,"%i\t%i\n",_alphabet_size,_m);
    fprintf(out,"# starting state(s)\n");
    for (vector<unsigned long>::iterator it=_start.begin(); it!=_start.end(); it++)
      fprintf(out,"%i\n",*it);
    fprintf(out,"# number of final states\n");
    fprintf(out,"%i\n",_final.size());
    fprintf(out,"# final state(s)\n");
    for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++)
      fprintf(out,"%i\n",*it);
    fprintf(out,"# nstates nz Pnz Qnz\n");
    fprintf(out,"%i\t%i\t%i\t%i\n",_nstates,_Pt.size()+_Qt.size(),_Pt.size(),_Qt.size());
    fprintf(out,"# P\n");
    for (vector<term>::iterator it=_Pt.begin(); it!=_Pt.end(); it++)
      fprintf(out,"%i\t%i\t%i\n",it->from,it->to,it->index);
    fprintf(out,"# Q\n");
    for (vector<term>::iterator it=_Qt.begin(); it!=_Qt.end(); it++)
      fprintf(out,"%i\t%i\t%i\n",it->from,it->to,it->index);
  };
  
};
