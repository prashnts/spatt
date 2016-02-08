/* $Id: word.cc 504 2005-10-19 11:53:32Z mhoebeke $ */
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
/*  Class word:                                                      */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "word.h"

using namespace std;

namespace spatt {

/* constructor from a label */
  word::word(const spattparameters &params,
	     const alphabet &alpha,
	     sequence &seq,
	     const count &occ,
	     char *label): 
    _params(&params), _alpha(&alpha), _seq(&seq), _occ(&occ),
    _dfa(NULL), _states(NULL), _token(NULL)
{

  sprintf(_label,label);
  //if (_alpha->_is_case_sensitive==0)
  //  for (int i=0; i<strlen(_label); i++)
  //    _label[i]=(char)tolower(_label[i]);
  for (unsigned int i=0; i<strlen(_label); i++)
    _label[i]=_alpha->code2char(_alpha->char2code(label[i]));
  _count=-1;
  _code=-1;
  _size=strlen(_label);
  if (_size<1) {
    fprintf(stderr,"Critical error in word::word(%p,%p,%p,\"%s\") ! Label is too short\n",&params,&alpha,&occ,label);
    exit(EXIT_FAILURE); 
  }

  //printf("word:word(%s)\n",label);
  /* if possible, update _count with _occ */
  if (_size<=_params->get_length()) {
    /* code word */
    {
      _code=0;
      int power=1;
      int c;
      for (int i=_size-1; i>=0; i--) {
	c=_alpha->char2code(_label[i]);
	if (c<0 || c>=_alpha->get_size()) {
	  fprintf(stderr,"Critical error in word:word(%p,%p,%p,\"%s\") ! \"%s\" uses invalid char\n",&params,&alpha,&occ,label,_label);
	  exit(EXIT_FAILURE);
	}
	_code+=power*c;
	power*=_alpha->get_size();
      }
    }
    /* get count */
    _count=_occ->get(_size,_code);
    if (_params->get_debug_level()>=3) 
      printf("word=\"%s\" size=%i code=%li count=%li\n",_label,_size,_code,_count);
  }
};

/* build dfa if: */
/* 1) word not already counted */
/* 2) a seq is available */
/* if any of these condition is not fullfilled */
/* function will do nothing */
void word::build_dfa(){
  //printf("Begin build_dfa\n");
  /* check conditions */
  //printf("_count=%li\n",_count);
  //printf("_seq=%p\n",_seq);
  if (_seq!=NULL && _count<0) {
    /* creating an empty dfa */
    //printf("real start\n");
    /* allocating memory */
    /* dfa alloc */
    {
      struct state **current;
      // line alloc
      _states=(struct state ***)malloc(sizeof(struct state **)*(_size+1));
      if (_states==NULL){
	fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	exit(EXIT_FAILURE);
      }
      // table alloc
      _states[0]=(struct state **)malloc(sizeof(struct state *)*(_size+1)*_alpha->get_size());
      if (_states[0]==NULL){
	fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	exit(EXIT_FAILURE);
      }
      // connexion
      current=_states[0];
      for (int i=0; i<(_size+1); i++) {
	_states[i]=current;
	current+=_alpha->get_size();
      }
    } // end states alloc
    /* dfa alloc */
    _dfa=(struct state *)malloc(sizeof(struct state)*(_size+1));
    if (_dfa==NULL){
      fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
      exit(EXIT_FAILURE);
    }
    /* end memory allocation */

    /* initialization */
    {
      _dfa[0].num=0;
      _dfa[0].tag=INITIAL_TAG;
      _dfa[0].link=_states[0];
      for (int i=0; i<_size; i++){
	_dfa[i].num=i;
	_dfa[i].tag=REGULAR_TAG;
	_dfa[i].link=_states[i];
      }
      _dfa[_size].num=_size;
      _dfa[_size].tag=COUNT_TAG;
      _dfa[_size].link=_states[_size];
      int imax=(_size+1)*_alpha->get_size();
      for (int i=0; i<imax; i++)
	_states[0][i]=&_dfa[0];
    }
    /* end initialization */
      
    //printf("end allocation and initialization\n");
    
    /* building dfa */
    {
      /* special treatment for initial state */
      int code=_alpha->char2code(_label[0]);
      _dfa[0].link[code]=&_dfa[1];
      /* main loop */
      int pos;
      char c;
      for (int i=1; i<=_size; i++) {
	/* get init pos for longest suffix */
	c=_label[i];
	//printf("working with \'%c\'\n",c);
	_label[i]='\0';
	//printf("get_state(\"%s\")\n",&_label[1]);
	pos=get_state(&_label[1]);	
	//printf("return=%i\n",pos);
	_label[i]=c;
	/* loop on code */
	for (int j=0; j<_alpha->get_size(); j++) {
	  if (_alpha->code2char(j)==c && i!=_size) {
	    _dfa[i].link[j]=&_dfa[i+1];
	    //_states[i][j]=&_dfa[i+1];
	  } else {
	    _dfa[i].link[j]=_dfa[pos].link[j];
	    //_states[i][j]=_dfa[pos].link[j];
	  }
	}
      }
    }

    if (_params->get_debug_level()>=2)
      display_dfa();

  } // end _seq!=NULL and _count<0
};

/* initialize count */
/* no effect without dfa */
void word::dfa_initialize_count(){
  _count=0;
  if (_dfa==NULL) {
    fprintf(stderr,"Critical error in word:dfa_initialize_count ! no dfa found\n");
    exit(EXIT_FAILURE);
  }
  _token=&_dfa[0];
};

/* update count processing dfa */
/* and return current count */
/* no effect with no dfa */
long word::dfa_update_count(int code){
  //printf("update(%i)->",code);
  if (_dfa==NULL) {
    fprintf(stderr,"Critical error in word:dfa_update_count ! no dfa found\n");
    exit(EXIT_FAILURE);
  }
  switch(code){
  case(IGNORE_CODE):
    break;
  case(COMMENT_CODE):
    break;
  case(INVALID_CODE):
    _token=&_dfa[0];
    break;
  case(SEQUENCE_CODE):
    _token=&_dfa[0];
    break;
  default:
    _token=_token->link[code];
    //printf("state(%i)\n",token[1]);
    if (_token->tag==COUNT_TAG) {
      _count++;
    }
  }
  return _count;
};


/* displays automaton */
void word::display_dfa(){
  printf("dfa(%s):\n",_label);
  for (int i=0; i<=_size; i++) {
    printf("  state(%i):",_dfa[i].num);
    switch(_dfa[i].tag) {
    case(INITIAL_TAG):
      printf("initial\n    |");
      break;
    case(COUNT_TAG):
      printf("count\n    |");
      break;
    default:
      printf("\n    |");
      break;
    }
    for (int j=0; j<_alpha->get_size(); j++) {
      if (_dfa[i].link[j])
	printf("\'%c\'->%i|",_alpha->code2char(j),_dfa[i].link[j]->num);
      else
	printf("\'%c\'->na|",_alpha->code2char(j));
    }
    printf("\n");
  }
}

/* return state after processing seq through dfa */
int word::get_state(char *seq){
  int i=0;
  char c;
  state *token=&_dfa[0];
  while ((c=seq[i])!='\0') {
    token=token->link[_alpha->char2code(c)];
    i++;
  }
  return token->num;
}

/* empty constructor */
  word::word() : _size(-1), _alpha(NULL), _seq(NULL), _occ(NULL), 
		 _count(-1), _dfa(NULL), _token(NULL),
		 _code(-1), _states(NULL)
{
  sprintf(_label,"");
};

/* copy constructor */
word::word(const word &source){
  //fprintf(stderr,"Copy constructor call for object word\n");
  //exit(EXIT_FAILURE);
  /* basic copy */
  {
    _params=source._params;
    _alpha=source._alpha;
    _seq=source._seq;
    _occ=source._occ;
    sprintf(_label,source._label);
    _size=source._size;
    _code=source._code;
    _count=source._count;
  }
  /* dfa and states copy */
  {
    if (source._dfa!=NULL && source._states!=NULL) {
      /* alloc */
      {
	struct state **current;
	// line alloc
	_states=(struct state ***)malloc(sizeof(struct state **)*(_size+1));
	if (_states==NULL){
	  fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	  exit(EXIT_FAILURE);
	}
	// table alloc
	_states[0]=(struct state **)malloc(sizeof(struct state *)*(_size+1)*_alpha->get_size());
	if (_states[0]==NULL){
	  fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	  exit(EXIT_FAILURE);
	}
	// connexion
	current=_states[0];
	for (int i=0; i<(_size+1); i++) {
	  _states[i]=current;
	  current+=_alpha->get_size();
	}
      } // end states alloc
      /* dfa alloc */
      _dfa=(struct state *)malloc(sizeof(struct state)*(_size+1));
      if (_dfa==NULL){
	fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	exit(EXIT_FAILURE);
      }      
      /* end memory allocation */

      /* copy */
      {
	for (int i=0; i<=_size; i++){
	  _dfa[i].num=i;
	  _dfa[i].tag=source._dfa[i].tag;
	  _dfa[i].link=_states[i];
	}
	int imax=(_size+1)*_alpha->get_size();
	for (int i=0; i<imax; i++)
	  _states[0][i]=&_dfa[source._states[0][i]->num];
      }
      /* copy */
    } else {
      _dfa=NULL;
      _states=NULL;
    }
  }
};

/* affectation operator */
word & word::operator=(const word &source){
  //fprintf(stderr,"affectation operator call for object word\n");
  //printf("source=\"%s\" _dfa=%p _states=%p\n",source._label,source._dfa,source._states);
  //exit(EXIT_FAILURE);    
  /* basic affectation */
	if (this != &source) {
  {
    _params=source._params;
    _alpha=source._alpha;
    _seq=source._seq;
    _occ=source._occ;
    sprintf(_label,source._label);
    _size=source._size;
    _code=source._code;
    _count=source._count;
  }
  /* dfa and states copy */
  {
    if (source._dfa!=NULL && source._states!=NULL) {
      
      /* free existing _dfa and _states */
      if (_dfa!=NULL)
	free(_dfa);
      if (_states!=NULL) {
	if (_states[0]!=NULL){
	  if (_states[0]!=NULL)
	    free(_states[0]);
	  free(_states);
	}
      }
      /* then alloc and copy source */
      /* alloc */
      {
	struct state **current;
	// line alloc
	_states=(struct state ***)malloc(sizeof(struct state **)*(_size+1));
	if (_states==NULL){
	  fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	  exit(EXIT_FAILURE);
	}
	// table alloc
	_states[0]=(struct state **)malloc(sizeof(struct state *)*(_size+1)*_alpha->get_size());
	if (_states[0]==NULL){
	  fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	  exit(EXIT_FAILURE);
	}
	// connexion
	current=_states[0];
	for (int i=0; i<(_size+1); i++) {
	  _states[i]=current;
	  current+=_alpha->get_size();
	}
      } // end states alloc
      /* dfa alloc */
      _dfa=(struct state *)malloc(sizeof(struct state)*(_size+1));
      if (_dfa==NULL){
	fprintf(stderr,"Critical error in word:build_dfa ! Not enough memory for allocation.\n");
	exit(EXIT_FAILURE);
      }      
      /* end memory allocation */

      /* copy */
      {
	for (int i=0; i<=_size; i++){
	  _dfa[i].num=i;
	  _dfa[i].tag=source._dfa[i].tag;
	  _dfa[i].link=_states[i];
	}
	int imax=(_size+1)*_alpha->get_size();
	for (int i=0; i<imax; i++)
	  _states[0][i]=&_dfa[source._states[0][i]->num];
      }
      /* copy */
    } else {
      _dfa=NULL;
      _states=NULL;
    }
  }
	}
	return *this;
};

/* destructor */
word::~word(){
  if (_states!=NULL) {
    if (_states[0]!=NULL)
      free(_states[0]);
    free(_states);
  }
  if (_dfa!=NULL)
    free(_dfa);
};

/* return word code */
/* and compute it if necessary */
long word::code(){
  if (_code==-1) {
    _code=0;
    int power=1;
    int c;
    for (int i=_size-1; i>=0; i--) {
      c=_alpha->char2code(_label[i]);
      if (c<0 || c>=_alpha->get_size()) {
	fprintf(stderr,"Critical error in word:word \"%s\" ! word use invalid char\n",_label);
	exit(EXIT_FAILURE);
      }
      _code+=power*c;
      power*=_alpha->get_size();
    }
  }
  return _code;
};
};
