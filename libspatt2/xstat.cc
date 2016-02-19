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
#include "xstat.h"

using namespace std;

namespace spatt {

  xstat::xstat(dfa &D,pmc &Y,sequence &S,int rep,long nobs,bool presence,bool verbose)
    : stat(D,Y,S,rep,nobs,verbose) {

    _presence=presence;

//    if (_presence) {
//      fprintf(stderr,"presence option not implemented in xstat is ignored\n");
//      _presence=false;
//    }

    if (_presence) {
      _c=1;
    } else {
      if (_rep==OVER)
	_c=_nobs;
      else
	_c=_nobs+1;
    }

    // set vector _x
    for (unsigned long j=0; j<_c; j++)
      _x.push_back(vector<double>(_Y->_nstates));
    _x.push_back(vector<double>(1));

  };

  xstat::xstat(dfa &D,pmc &Y,unsigned long sequence_length,int rep,long nobs,bool presence,bool verbose)
    : stat(D,Y,sequence_length,rep,nobs,verbose) {

    _presence=presence;

//    if (_presence) {
//      fprintf(stderr,"presence option not implemented in xstat is ignored\n");
//      _presence=false;
//    }

    if (_presence) {
      _c=1;
    } else {
      if (_rep==OVER)
    _c=_nobs;
      else
    _c=_nobs+1;
    }

    // set vector _x
    for (unsigned long j=0; j<_c; j++)
      _x.push_back(vector<double>(_Y->_nstates));
    _x.push_back(vector<double>(1));

  };

  void xstat::compute(bool interrupt,int offset,bool verbose) {
    vector<double> mu0;
    compute(mu0,interrupt,offset,verbose);
  };

  void xstat::compute(vector<double> &mu0,bool interrupt,int offset,bool verbose) {

    if (verbose)
      printf(">>> call xstat::compute()\n");

    if (_presence) { // presence case

      // gather information
      map<unsigned long,list<unsigned long> > starts;
      for (vector<seq>::iterator it=_seq.begin(); it!=_seq.end(); it++) {
	starts[it->start].push_back(it->length);
      }

      // initialization
      for (unsigned long j=0; j<=_c; j++)
	for (vector<double>::iterator it=_x[j].begin(); it!=_x[j].end(); it++)
	  *it=0.0;
      _x[0][0]=1.0;

      if (!interrupt) { // with ignore-interrupt
	fprintf(stderr,"option ignore-interrupt not implemented with presence in xstat. Aborting\n");
	exit(EXIT_FAILURE);
      } else {  // with interruptions
	unsigned long presence_size=_seq.size()+1;
	double *ppresence= new double[presence_size];
	double *ppresence_bis= new double[presence_size];
	double *ppresence_aux= NULL;
	double p,q;
	for (unsigned long i=0; i<presence_size; i++) {
	  ppresence[i]=0.0;
	  ppresence_bis[i]=0.0;
	}

	vector<double> aux(_Y->_nstates);
	bool first=true;
	unsigned long nseq=0;

	for (map<unsigned long,list<unsigned long> >::iterator it=starts.begin(); it!=starts.end(); it++) {
	  // sort the list (increasing order)
	  it->second.sort();
	  if (verbose) {
	    printf("start = %i, ell = [ ",it->first);
	    for (list<unsigned long>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++)
	      printf("%i ",*it2);
	    printf("]\n");
	  }
	  // initialization
	  for (unsigned long j=0; j<=_c; j++)
	    for (vector<double>::iterator it2=_x[j].begin(); it2!=_x[j].end(); it2++)
	      *it2=0.0;
	  _x[0][0]=0.0;
	  if (mu0.empty()) {
	    _x[0][it->first]=1.0;
	  } else {
	    unsigned long imax=_Y->_start.size();
	    for (unsigned long i=0; i<imax; i++)
	      _x[0][_Y->_start[i]]=mu0[i];
	  }
	  // initialize position in the sequence to _m
	  unsigned long pos=_m;
	  for (list<unsigned long>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++) {
	    unsigned long new_pos=*it2;
	    nseq++;
	    // loop on position in the sequence
	    for (unsigned long i=pos; i<new_pos; i++) {
	      // _x[_c] = _x[_c] + _x[_c-1]*S_Qt
	      _x[_c][0]+=_Y->scalprod(_x[_c-1],_Y->_SQ);
	      // for j=(_c-1)...1 : _x[j]=_x[j]*P+_x[j-1]*Q
	      for (unsigned long j=(_c-1); j>=1; j--) {
		_Y->zero(aux);
		_Y->add_xM(_x[j],_Y->_Pt,aux);
		_Y->add_xM(_x[j-1],_Y->_Qt,aux);
		_x[j]=aux;
	      }
	      // _x[0]=_x[0]*P
	      _Y->zero(aux);
	      _Y->add_xM(_x[0],_Y->_Pt,aux);
	      _x[0]=aux;
	    } // end loop on positions in the sequence

	    // get p=P(N>=1) and q=1-p=P(N=0)
	    p=_x[1][0];
	    q=0.0;
	    for (vector<double>::iterator it3=_x[0].begin(); it3!=_x[0].end(); it3++)
	      q+=*it3;

	    // update pos
	    pos=new_pos;
	    // update ppresence
	    if (first) {
	      ppresence[0]=q;
	      ppresence[1]=p;
	      first=false;
	    } else {
	      // update position 0
	      ppresence_bis[0]=ppresence[0]*q;
	      // update positions 1 to nseq
	      for (unsigned long i=1; i<=nseq; i++)
		ppresence_bis[i]=ppresence[i-1]*p+ppresence[i]*q;
	      // switch ppresence_bis and ppresence
	      ppresence_aux=ppresence;
	      ppresence=ppresence_bis;
	      ppresence_bis=ppresence_aux;
	    }

	  } // end for it2

	} // it map loop

	// get the pvalue
	_pvalue=0.0;
	if (_rep==OVER) {
	  for (unsigned long i=_npresence; i<=presence_size; i++)
	    _pvalue+=ppresence[i];
	} else {
	  for (unsigned long i=0; i<=_npresence; i++)
	    _pvalue+=ppresence[i];
	}
	// free allocated memory
	delete[] ppresence;
	delete[] ppresence_bis;
      }

    } else { // not presence case

      if ((_nobs==0)&&(_rep==OVER)) {
	_x[0][0]=1.0;
	_pvalue=1.0;
      } else {
	// initialization
	for (unsigned long j=0; j<=_c; j++)
	  for (vector<double>::iterator it=_x[j].begin(); it!=_x[j].end(); it++)
	    *it=0.0;
	_x[0][0]=1.0;

	// loop on sequences
	vector<double> aux(_Y->_nstates);
	bool first=true;
	for (vector<seq>::iterator it=_seq.begin(); it!=_seq.end(); it++) {
	  if (!interrupt) {
	    //fprintf(stderr,"option --ignore-interrupt not (yet) implemented in xstat\n");
	    // find the first sequence which length is longuer than _m

	    if (first) {
	      if (it->length>=_m) {
		first=false;
		_x[0][0]=0.0;
		if (mu0.empty()) {
		  _x[0][it->start]=1.0;
		} else {
		  unsigned long imax=_Y->_start.size();
		  for (unsigned long i=0; i<imax; i++)
		    _x[0][_Y->_start[i]]=mu0[i];
		}
		// loop on position in the sequence
		for (unsigned long i=_m; i<it->length; i++) {
		  // _x[_c] = _x[_c] + _x[_c-1]*S_Qt
		  _x[_c][0]+=_Y->scalprod(_x[_c-1],_Y->_SQ);
		  // for j=(_c-1)...1 : _x[j]=_x[j]*P+_x[j-1]*Q
		  for (unsigned long j=(_c-1); j>=1; j--) {
		    _Y->zero(aux);
		    _Y->add_xM(_x[j],_Y->_Pt,aux);
		    _Y->add_xM(_x[j-1],_Y->_Qt,aux);
		    _x[j]=aux;
		  }
		  // _x[0]=_x[0]*P
		  _Y->zero(aux);
		  _Y->add_xM(_x[0],_Y->_Pt,aux);
		  _x[0]=aux;
		} // end loop on positions in the sequence
	      }
	    } else {
	      // loop on position in the sequence
	      for (unsigned long i=0; i<it->length-offset; i++) {
		// _x[_c] = _x[_c] + _x[_c-1]*S_Qt
		_x[_c][0]+=_Y->scalprod(_x[_c-1],_Y->_SQ);
		// for j=(_c-1)...1 : _x[j]=_x[j]*P+_x[j-1]*Q
		for (unsigned long j=(_c-1); j>=1; j--) {
		  _Y->zero(aux);
		  _Y->add_xM(_x[j],_Y->_Pt,aux);
		  _Y->add_xM(_x[j-1],_Y->_Qt,aux);
		  _x[j]=aux;
		}
		// _x[0]=_x[0]*P
		_Y->zero(aux);
		_Y->add_xM(_x[0],_Y->_Pt,aux);
		_x[0]=aux;
	      } // end loop on positions in the sequence
	    }
	    // interrupt part now
	  } else {
	    // treat only sequences longuer than _m
	    if (it->length>=_m) {
	      // sum and reset
	      for (unsigned long j=0; j<_c; j++) {
		double sum=0.0;
		for (vector<double>::iterator it2=_x[j].begin(); it2!=_x[j].end(); it2++)
		  sum+=*it2;
		for (vector<double>::iterator it2=_x[j].begin(); it2!=_x[j].end(); it2++)
		  *it2=0.0;
		if (mu0.empty()) {
		  _x[j][it->start]=sum;
		} else {
		  unsigned long imax=_Y->_start.size();
		  for (unsigned long i=0; i<imax; i++)
		    _x[j][_Y->_start[i]]=mu0[i]*sum;
		}
	      } // end sum and reset
	      // loop on positions in the sequence
	      for (unsigned long i=_m; i<it->length; i++) {
		// _x[_c] = _x[_c] + _x[_c-1]*S_Qt
		_x[_c][0]+=_Y->scalprod(_x[_c-1],_Y->_SQ);
		// for j=(_c-1)...1 : _x[j]=_x[j]*P+_x[j-1]*Q
		for (unsigned long j=(_c-1); j>=1; j--) {
		  _Y->zero(aux);
		  _Y->add_xM(_x[j],_Y->_Pt,aux);
		  _Y->add_xM(_x[j-1],_Y->_Qt,aux);
		  _x[j]=aux;
		}
		// _x[0]=_x[0]*P
		_Y->zero(aux);
		_Y->add_xM(_x[0],_Y->_Pt,aux);
		_x[0]=aux;
	      } // end loop on positions in the sequence
	    } // end if length >= _m
	  } // end if (!interrupt)
	} // end loop on sequences

	// compute the pvalue
	_pvalue=0.0;
	if (_rep==OVER)
	  _pvalue=_x[_c][0];
	else
	  for (unsigned long j=0; j<_c; j++)
	    for (vector<double>::iterator it2=_x[j].begin(); it2!=_x[j].end(); it2++)
	      _pvalue+=*it2;
      }
    }
  };

  void xstat::print_distribution(std::string &format_float) {

    if (!_presence) {
      printf("distribution:\n");
      for (unsigned long j=0; j<_c; j++) {
	double sum=0.0;
	for (vector<double>::iterator it2=_x[j].begin(); it2!=_x[j].end(); it2++)
	  sum+=*it2;
        printf("P(N=%i)=",j);
        printf(format_float.c_str(),sum);
        printf("\n");
      }
      printf("P(N>=%i)=",_c);
      printf(format_float.c_str(),_x[_c][0]);
      printf("\n");
    }

  };

  void xstat::print(const char *label,std::string &format_float) {
    if (!_presence)
      printf("pattern=%s\tNobs=%i\t",label,_nobs);
    else
      printf("pattern=%s\tNpresence=%i\t",label,_npresence);
    if (_rep==OVER)
      printf("P(N>=Nobs)=");
    else
      printf("P(N<=Nobs)=");
    printf(format_float.c_str(),_pvalue);
    printf("\n");
  };

};
