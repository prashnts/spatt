/* $Id: count.cc 547 2005-11-22 10:30:51Z gnuel $ */
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

#include "count.h"
#include <iostream>

using namespace std;

namespace spatt {

/* count from a countfile */
count::count(const spattparameters &params,
	     const alphabet &alpha): 
  _params(&params), _alpha(&alpha),_seq(NULL),_n(0)
{

  
  _length=_params->get_length();
  if (_length>0) {
    sprintf(_filename,_params->get_count_filename());

    /* opens stream */
    _stream=fopen(_filename,"r");
    if (_stream==NULL) {
      cerr << "Cannot open count file \"" << _filename << "\"" << endl;
      exit(EXIT_FAILURE);
    }
    
    this->init_occ();
    
    /* reading occ from stream */
    cerr << "count constructor from a count file not yet implemented !" << endl;
    exit(EXIT_FAILURE);
  } else {
    _occ=NULL;
  }
  
};

/* count from a sequence */
  count::count(const spattparameters &params,
	       const alphabet &alpha,
	       sequence &seq) : 
    _params(&params),_alpha(&alpha),_seq(&seq),_stream(NULL),_n(0),
    _occ(NULL), _dim(NULL)
{

  _length=_params->get_length();

  if (_length>0) {
    this->init_occ();
    
    //printf("ok\n");
    
    /* reading from sequence */
    {
      int ntoken=0;
      long token[MAX_LENGTH];
      for (int i=0; i<MAX_LENGTH; i++)
	token[i]=0;
      
      /* restarting sequence (in any case) */
      _seq->restart();
      
      /* main loop */
      int c;
      c=_seq->next();
      while ( c!= EOF_CODE) {
	switch(c){
	case(IGNORE_CODE):
	  /* nothing to do */
	  break;
	case(COMMENT_CODE):
	  /* nothing to do */
	  break;
	case(INVALID_CODE):
	  /* reset all tokens */
	  ntoken=0;
	  for (int i=0; i<MAX_LENGTH; i++)
	    token[i]=0;
	  break;
	case(SEQUENCE_CODE):
	  /* reset all tokens */
	  ntoken=0;
	  for (int i=0; i<MAX_LENGTH; i++)
	    token[i]=0;
	  break;
	default:
	  /* c is a valid char code */
	  _n++;
	  if (ntoken<_length)
	    ntoken++;
	  for (int i=ntoken-1; i>0; i--) {
	    token[i]=_alpha->get_size()*token[i-1]+c;
	    _occ[i][token[i]]++;
	  }
	  token[0]=c;
	  _occ[0][c]++;
	  break;
	}
	c=_seq->next();
      }
      
    }
  } else {
    _occ=NULL;
  }


  if (_params->get_debug_level()>=1)
    cout << "sequence length = " << _n << endl;

};


/* returns count for code word of length h */
long count::get(int h,long code) const {
  if (_occ!=NULL) {
    //printf("get begin\n");
    if (h>=1 && h<=_length) {
      //printf("_dim[%i]=%li\n",h-1,_dim[h-1]);
      if (code<_dim[h-1])
	return _occ[h-1][code];
      else {
	cerr << "Call of count::get(" << h << "," << code <<"). Second parameter out of range [0;" << _dim[h-1]-1 << "]!" << endl;
	exit(EXIT_FAILURE);      
      }
    } else {
	cerr << "Call of count::get(" << h << "," << code <<"). First parameter out of range [1;" << _length << "]!" << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    return -1;
  }
};

/* alloc and intitalize _occ */
void count::init_occ(){

  _memory=0;

  _dim=(long *)malloc(sizeof(long)*_length);

  if (!_dim) {
    cerr << "Not enough memory for alloctation in count" << endl;
    exit(EXIT_FAILURE);
  }
  /* memory allocation */
  {
    long current_dim=1;
    _occ=(long **)malloc(sizeof(long *)*_length);
    if (!_occ) {
      cerr << "Not enough memory for alloctation in count" << endl;
      cerr << "length specified (l=" << _length 
	   << ") is too long, try a shorter one (option -l)" << endl;
      exit(EXIT_FAILURE);
    }
    for (int i=0; i<_length; i++) {
      current_dim*=_alpha->get_size();
      _dim[i]=current_dim;
      //printf("dim[%i]=%li\n",i,_dim[i]);
      _occ[i]=(long *)malloc(sizeof(long)*current_dim);
      if (!_occ[i]) {
	cerr << "Not enough memory for alloctation in count" << endl;
	cerr << "length specified (l=" << _length 
	     << ") is too long, try a shorter one (option -l)" << endl;
	exit(EXIT_FAILURE);
      }
      _memory+=sizeof(long)*current_dim;
      for (int j=0; j<current_dim; j++)
        _occ[i][j]=0;
    }
  }

  if (_params->get_debug_level()>=1)
    cout <<"count memory usage (bytes) = " << _memory << endl;
};

/* close stream */
count::~count(){
  if (_stream)
    fclose(_stream);
  for (int i=0; i<_length; i++) {
    free(_occ[i]);
  }
  free(_occ);

  free(_dim);
};


};
