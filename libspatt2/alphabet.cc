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
#include "alphabet.h"

using namespace std;

namespace spatt {

alphabet::alphabet (string &alphabet_label) {

  _alphabet_label=alphabet_label;
  _alphabet_size=_alphabet_label.size();
  
  bool case_sensitive=false;

  { // check if the alphabet is case sensitive
    int lower=0,upper=0;
    int i=0;
    unsigned char current=alphabet_label[0];
    while (current!='\0') {
      if (current==tolower(current))
	lower=1;
      if (current==toupper(current))
	upper=1;
      current=alphabet_label[++i];
    }
    if (lower*upper==1)
      case_sensitive=true;
  }


  { // filling code table
    for (unsigned short i=0; i<256; i++)
      _alphabet_code[i]=UNDEF;
    if (case_sensitive) {
    for (unsigned short i=0; i<_alphabet_size; i++) 
      _alphabet_code[(unsigned short)_alphabet_label[i]]=i;
    } else {
      for (unsigned short i=0; i<_alphabet_size; i++) {
	_alphabet_code[(unsigned short)tolower(_alphabet_label[i])]=i;
	_alphabet_code[(unsigned short)toupper(_alphabet_label[i])]=i;
      }	  
    }
    /* separators */
    _alphabet_code[(unsigned short)'-']=SEPARATOR;
    _alphabet_code[(unsigned short)',']=SEPARATOR;
    _alphabet_code[(unsigned short)':']=SEPARATOR;
    _alphabet_code[(unsigned short)';']=SEPARATOR;
    _alphabet_code[(unsigned short)'|']=SEPARATOR;
    /* end block */
    _alphabet_code[(unsigned short)')']=END_BLOCK;
    _alphabet_code[(unsigned short)'}']=END_BLOCK;
    _alphabet_code[(unsigned short)']']=END_BLOCK;
    /* repeat */
    _alphabet_code[(unsigned short)'(']=REPEAT;
    /* choice */
    _alphabet_code[(unsigned short)'[']=CHOICE;
    /* nchoice */
    _alphabet_code[(unsigned short)'{']=NCHOICE;
    /* any char */
    _alphabet_code[(unsigned short)'_']=ANY_CHAR;
    _alphabet_code[(unsigned short)'.']=ANY_CHAR;
    /* sequence */
    _alphabet_code[(unsigned short)'>']=SEQUENCE;
    /* comment */
    _alphabet_code[(unsigned short)'#']=COMMENT;
    _alphabet_code[(unsigned short)'%']=COMMENT;
    /* ignore */
    _alphabet_code[(unsigned short)' ']=IGNORE;
    _alphabet_code[(unsigned short)'/']=IGNORE;
    _alphabet_code[(unsigned short)'\t']=IGNORE;
    _alphabet_code[(unsigned short)'\r']=IGNORE;
    _alphabet_code[(unsigned short)'\n']=IGNORE;
    _alphabet_code[(unsigned short)'0']=IGNORE;
    _alphabet_code[(unsigned short)'1']=IGNORE;
    _alphabet_code[(unsigned short)'2']=IGNORE;
    _alphabet_code[(unsigned short)'3']=IGNORE;
    _alphabet_code[(unsigned short)'4']=IGNORE;
    _alphabet_code[(unsigned short)'5']=IGNORE;
    _alphabet_code[(unsigned short)'6']=IGNORE;
    _alphabet_code[(unsigned short)'7']=IGNORE;
    _alphabet_code[(unsigned short)'8']=IGNORE;
    _alphabet_code[(unsigned short)'9']=IGNORE;
  } // end filling code table

};

};
