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
#include "stat.h"

using namespace std;

namespace spatt {

stat::stat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool verbose){

  if (verbose)
    printf(">>> call stat::stat\n");
  _rep=rep;
  if ((rep!=OVER)&&(rep!=UNDER)) {
    fprintf(stderr,"stat::stat : wrong value for rep. Aborting.\n");
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
    seq s;
    _nobs=0;
    _npresence=0;
    bool seen=false;
    while ((c=_Seq->next())!=SEQEMPTY) {
      //printf("c=%i\n",c);
      if (c>=0) { // valid char
	token=_Dfa->process(c);
	if (i++<_m) {
	  Ym=_Y->recode(token);
	}
	if (_Dfa->is_final(token)) {
	  _nobs++;
	  seen=true;
	}
      } else { // sequence end
	if (seen) {
	  _npresence++;
	  seen=false;
	}
	s.start=Ym;
	s.length=_Seq->lastnvalidchar();
	_seq.push_back(s);
	i=0;
	token=_Dfa->reset();
	Ym=_Y->recode(token);
      }
    }
  } // end parse the sequence

  if (nobs>=0) {
    _nobs=nobs;
    _npresence=nobs;
  }

  if (verbose) {
    printf("nobs=%i\n",_nobs);
    printf("npresence=%i\n",_npresence);
    printf("m=%i\n",_m);
    printf("seq(start;length)=");
    for (vector<seq>::iterator it=_seq.begin(); it!=_seq.end(); it++)
      printf("(%i;%i)",it->start,it->length);
    printf("\n");
  }

  if (verbose)
    fprintf(stderr,"end of stat::stat()\n");

};

stat::stat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool verbose){

  if (verbose)
    printf(">>> call stat::stat\n");
  _rep=rep;
  if ((rep!=OVER)&&(rep!=UNDER)) {
    fprintf(stderr,"stat::stat : wrong value for rep. Aborting.\n");
    exit(EXIT_FAILURE);
  }
  _Dfa=&D;
  _Y=&Y;
  _m=_Y->_m;

    {
        seq s;
        s.start=1;
        s.length=sequence_length;
        _seq.push_back(s);
    }

  if (nobs>=0) {
    _nobs=nobs;
    _npresence=nobs;
  }

  if (verbose) {
    printf("nobs=%i\n",_nobs);
    printf("npresence=%i\n",_npresence);
    printf("m=%i\n",_m);
    printf("seq(start;length)=");
    for (vector<seq>::iterator it=_seq.begin(); it!=_seq.end(); it++)
      printf("(%i;%i)",it->start,it->length);
    printf("\n");
  }

  if (verbose)
    fprintf(stderr,"end of stat::stat()\n");

};

void stat::compute(bool interrupt,bool verbose) {
  fprintf(stderr,"stat::compute() : nothing to do\n");
};

void stat::compute(vector<double> &mu0,bool interrupt,bool verbose) {
  fprintf(stderr,"stat::compute() : nothing to do\n");
};

void stat::print(const char *label) {
  fprintf(stderr,"stat::print() : nothing to do\n");
};

};
