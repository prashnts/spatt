/* $Id: alphabet.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "alphabet.h"

namespace spatt {
/* constructor using a label */ 
  alphabet::alphabet(const spattparameters &params) : 
    _is_case_sensitive(false),
    _has_complementary(false),
    _params(&params)
{

  _label=_params->get_alphabet_label();
  if (_params->get_debug_level()>=3)
    std::cout << "Begin: " << this << "=alphabet(" << _label << ")" << std::endl;

  /* code initialization */
  for (int i=0; i<CHAR_RANGE; i++)
    _code[i]=IGNORE_CODE;
  /* set alpha char to invalid */
  {
    char alpha[MAX_STRING_LENGTH];
    sprintf(alpha,ALPHA_CHAR);
    int imax=strlen(alpha);
    for (int i=0; i<imax; i++)
      _code[(unsigned char)alpha[i]]=INVALID_CODE;
  }
  /* set comment and sequence codes */
  _code[(unsigned char)COMMENT_CHAR]=COMMENT_CODE;
  _code[(unsigned char)SEQUENCE_CHAR]=SEQUENCE_CODE;
  /* ignore whitespaces */
  /* no more necessary */
  //_code[(unsigned char)' ']=IGNORE_CODE;
  //_code[(unsigned char)'\t']=IGNORE_CODE;
  //_code[(unsigned char)'\n']=IGNORE_CODE;
  //_code[(unsigned char)'\r']=IGNORE_CODE;
  /* check if alphabet is case sensitive */
  
  int lower=0,upper=0;
  int i=0;
  char current=_label[0];
  while (current!='\0') {
    if (current!='[' && current!=']' && current!=':') {
      if (current==tolower(current))
	lower=1;
      if (current==toupper(current))
	upper=1;
    }
    current=_label[++i];
  }
  if (lower*upper==1) 
      _is_case_sensitive=true;

  /* parsing label */
  {
    int i=0,c=0;
    int is_bracket=0;
    char current=_label[0];
    while (current!='\0' && current!=':') {
      switch(current) {
      case('['):
	/* begin a char definition */
	is_bracket=1;
	break;
      case(']'):
	/* end char definition */
	is_bracket=0;
	c++;
	break;
      default:
	if (_is_case_sensitive) {
	  if (_params->get_debug_level()>=3)	
	    cout << current << "->" << (unsigned char)current << "->" << c << endl;
	  _code[(unsigned char)current]=c;
	  _decode[c]=current;
	} else {
	  if (_params->get_debug_level()>=3)
	    cout << "[" << tolower(current) << toupper(current) << "]->" 
		 << (unsigned char)tolower(current) << " " 
		 << (unsigned char)toupper(current) << "->" << c << endl;
	  _code[(unsigned char)tolower(current)]=c;
	  _code[(unsigned char)toupper(current)]=c;
	  _decode[c]=tolower(current);
	}
	if (!is_bracket)
	  c++;
      }
      current=_label[++i];
    }
    /* ending decode table */
    _decode[c]='\0';
    if (_params->get_debug_level()>=3)
      cout << "decode table \"" << _decode << "\"" << endl;
    /* alphabet size */
    _size=c;
    if (_params->get_debug_level()>=3)
      cout << "alphabet size " << _size << endl;

    /* complementary part */
    if (current!='\0') {
      i++;
      if (strlen(&_label[i])!=_size) {
	if (_params->get_debug_level()>=0)
	  cerr << "Warning: complementary string does not match alphabet size and will be ignored!" << endl;
      } else {
	/* testing validity of complementary part */
	int invalid=0;
	for (int j=0;j<_size; j++)
	  if (_code[(unsigned char)_label[i+j]]<0)
	    invalid=1;
	if (invalid) {
	  if (_params->get_debug_level()>=0)
	  cerr << "Warning: complementary string does not match alphabet definition and will be ignored!" << endl;
	} else {
	  sprintf(_idecode,&_label[i]);
	  _has_complementary=1;
	  if (_params->get_debug_level()>=3)
	    cout << "idecode table \"" << _idecode << "\"" << endl;
	}
      }
    }

  }
  /* end parsing */
  
  /* Ignore numbers when no one belongs to alphabet */
  {
    int is_number=0;
    /* constructing a string for numbers */
    char number[2];
    /* testing */
    for (int k=0; k<10; k++) {
      sprintf(number,"%i",k);
      if (_code[(unsigned char)number[0]]>=0)
	is_number=1;
    }
    /* ignoring */
    if (!is_number) {
      for (int k=0; k<10; k++) {
	sprintf(number,"%i",k);
	_code[(unsigned char)number[0]]=IGNORE_CODE;
      }
    }      
  }
  
  /* displays code */
  if (_params->get_debug_level()>=4) {
    printf("code:\n");
    for (int i=0; i<CHAR_RANGE; i++) {
      switch(_code[(unsigned char)i]) {
      case(IGNORE_CODE):
	/* no output */
	break;
      case(INVALID_CODE):
	cout << "'" << (unsigned char)i <<"'-> ignored" << endl;
	break;
      case(COMMENT_CODE):
	cout << "'" << (unsigned char)i <<"'-> comment" << endl;
	break;
      case(SEQUENCE_CODE):
	cout << "'" << (unsigned char)i <<"'-> sequence" << endl;
	break;
      default:
	cout << "'" << (unsigned char)i << "'->" << _code[(unsigned char)i] << endl;
      }
    }
  } 

  //if (_is_case_sensitive==0)
  //  for (int i=0; i<strlen(_label); i++)
  //    _label[i]=(char)tolower(_label[i]);

  if (_params->get_debug_level()>=3)
    cout << "End: " << this << "=alphabet(\"" << _label << "\")" << endl;
};

alphabet::~alphabet() {
};

string* alphabet::code2string(long code,int size) const {
  string *res=new string(size,'*');
  // code to add
  double tempcode=code;
  long div=(long)pow((double)_size,(double)(size-1));
  int charcode;
  for (int i=0; i<size; i++) {
    charcode=(int)(tempcode/div);
    (*res)[i]=code2char(charcode);
    tempcode-=charcode*div;
    div/=_size;
  }
  return res;
};

};
