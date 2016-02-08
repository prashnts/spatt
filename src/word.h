/* $Id: word.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class word                                                       */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef WORD_H
#define WORD_H

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "count.h"

using namespace std;


namespace spatt {

  const int INITIAL_TAG = -1;
  const int REGULAR_TAG =  0;
  const int COUNT_TAG   =  1;

  struct state {
    int num;
    int tag;
    struct state **link;
  };

  class word {

  public:
    
    /* constructor from a label */
    /* seq could be NULL but then no dfa stuff */
    /* will be done */
    word(const spattparameters &, const alphabet &, sequence &, const count &, 
	 char *);
    
    /* empty constructor */
    word();
    
    /* copy constructor */
    word(const word &source);
    
    /* affectation operator */
    word & operator=(const word &source);
    
    /* destructor */
    ~word();
    
    /* build dfa if: */
    /* 1) word not already counted */
    /* 2) a seq is available */
    /* if any of these condition is not fullfilled */
    /* function will do nothing */
    void build_dfa();
    
    /* initialize count */
    /* no effect without dfa */
    void dfa_initialize_count();
    
    /* update count processing dfa */
    /* and return current count */
    /* no effect with no dfa */
    long dfa_update_count(int code);
    
    /* return word code */
    /* and compute it if necessary */
    long code();

    inline double get_expected() const { return _expected; }
    inline void set_expected(const double expected) { _expected=expected; }

    inline char * get_label() { return _label; }
    
    inline int get_size() const { return _size;  }

    inline long get_code() const { return _code; }
    inline void set_code(long code) { _code=code; }

    inline const long get_count() const { return _count; }
    inline void set_count(int count) { _count=count; }


    inline const alphabet *get_alpha() const { return _alpha; }

    inline const struct state *get_dfa() const { return _dfa; }

  private:
    
    const spattparameters *_params;
    const alphabet *_alpha;
    sequence *_seq;
    const count *_occ;
    
    char _label[MAX_STRING_LENGTH];
    int _size;
    long _code;
    long _count;
    double _expected;
    
    struct state *_dfa;
    struct state ***_states;
    struct state *_token;
    
    /* displays automaton */
    void display_dfa();
    
    /* return state after processing seq through dfa */
    int get_state(char *seq);
    
  };

};

#endif
