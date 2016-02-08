/* $Id: alphabet.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class alphabet                                                   */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef ALPHABET_H
#define ALPHABET_H

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <string>
#include <cmath>

#include "spattglobals.h"
#include "spattparameters.h"

using namespace std;

namespace spatt {

  const int CHAR_RANGE     = 256;
  const int INVALID_CODE   =  -1;
  const int COMMENT_CODE   =  -2;
  const int SEQUENCE_CODE  =  -3;
  const int IGNORE_CODE    =  -4;
  const int EOF_CODE       =  -5;
  const char SEQUENCE_CHAR = '>';
  const char ALPHA_CHAR[]  = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  
  class alphabet {
  

  public:
    
    /* constructor using main arguments */
    alphabet(const spattparameters &);
    
    /* destructor */
    ~alphabet();


    /* return the string corresponding to a word */
    string* code2string(long code,int size) const ;

    /* returns size */
    inline int get_size() const {
      return _size;
    };
    
    /* gives code for a char */
    inline int char2code(char c) const {
      return _code[(unsigned char)c];
    };
    
    /* gives char for a code */
    /* '\0' if not available */
    inline char code2char(int c) const {
      if (c>=0 && c<_size)
	return _decode[c];
      else
	return '\0';
    };
    
    /* gives code of complementary */
    /* -1 if not available */
    inline int complement_code(int c) const {
      if (c>=0 && c<_size)
	return char2code(_idecode[c]);
      else
	return -1;
    };
    
    
    /* gives char for complementary */
    /* '\0' if not available */
    inline char complement_char(char c) const {
      //printf("compl: idecode=\"%s\"\n",_idecode);
      if (char2code(c)>=0)
	return _idecode[char2code(c)];
      else
	return '\0';
    };
    
    inline string get_label() const {
      return _label;
    };

    inline bool is_case_sensitive() const {
      return _is_case_sensitive;
    }

    inline bool has_complementary() const {
      return _has_complementary;
    }
					   
    
  private:
    /* empty constructor */
    alphabet();
    
    /* copy constructor */
    alphabet(const alphabet &source);
    
    /* affectation operator */
    alphabet & operator=(const alphabet &);
    
    int _code[CHAR_RANGE];
    int _size;
    char _decode[MAX_STRING_LENGTH];
    char _idecode[MAX_STRING_LENGTH];
    string _label;
    bool _is_case_sensitive;
    bool _has_complementary;
    const spattparameters *_params;
  };
  
};

#endif
