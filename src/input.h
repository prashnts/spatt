/* $Id: input.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class input: for core program inputs                             */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef INPUT_H
#define INPUT_H

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>

#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "word.h"
#include "count.h"
#include "pattern.h"

using namespace std;

namespace spatt {
  const int INVALID_WORD = -1;
  const int COUNT_WORD   =  1;
  const int DFA_WORD     =  2;

class input {

 public:

  /* main constructor */
  input(spattparameters &, const alphabet &, sequence &, spatt::count &);

  /* destructor */
  ~input();

  /* test a word descriptor */
  /* returns: INVALID_WORD if invalid */
  /* COUNT_WORD if size included in count */
  /* DFA_WORD else */
  int test_word(string word);


  inline const int get_optimized_dfa_word_size() const
    { return _optimized_dfa_word_size; }

  inline word ** get_optimized_dfa_word() const 
    { return _optimized_dfa_word; }

  inline const int get_optimized_all_word_size() const 
    { return _optimized_all_word_size; }

  inline word ** get_optimized_all_word() const
    { return _optimized_all_word; }


  inline vector<pattern> & get_pattern_list() 
    { return _pattern_list; }

 private:
  input();
  input(const input &);
  input & operator=(const input &);
  
  spattparameters *_params;
  const alphabet *_alpha;
  sequence *_seq;
  spatt::count *_count;


  FILE *_stream;
  
  vector<pattern> _pattern_list;
  map<string,word> _all_word_table;
  vector<string> _all_word_list;
  vector<string> _count_word_list;
  vector<string> _dfa_word_list;

  int _optimized_all_word_size;
  word ** _optimized_all_word;
  int _optimized_count_word_size;
  word ** _optimized_count_word;
  int _optimized_dfa_word_size;
  word ** _optimized_dfa_word;

  //int _optimized_pattern_size;
  //pattern ** _optimized_pattern;
  
};

};

#endif
