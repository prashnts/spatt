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
#ifndef NFA_H
#define NFA_H 1

#include <iterator>
#include <map>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "alphabet.h"

//#define UNDEF 0
//#define SEPARATOR -1 
//#define REPEAT -2
//#define CHOICE -3
//#define NCHOICE -4
//#define END_BLOCK -5
//#define ANY_CHAR -6

#define MAX_DOT 50

namespace spatt {

typedef struct {
  unsigned int from;
  unsigned int to;
} trans;

class nfa {
 public:

  /* alphabet part */
  alphabet *_Alpha;
  unsigned short _alphabet_size;

//  std::string  _alphabet_label;
//  int _alphabet_code[256]; /* coding table */

  /* pattern part */
  std::string  _pattern_label;

  /* statistics */
  unsigned int _nstates;
  unsigned int *_ntrans; /* table of dim _alphabet_size+1; _ntrans[0]
			    give total number of epsilon-trans,
			    _ntrans[1] 1-trans, etc ... */
  unsigned int _total_ntrans; /* sum of _ntrans */

  /* data */
  unsigned int _start;
  unsigned int _final;
  unsigned int _data_size;
  trans *_data; /* table of dim _total_ntrans with _trans[0]
		   (_ntrans[0] terms) then _trans[1] (_ntrans[1]
		   terms), etc ... */
  trans **_trans; /* table of dim _alphabet_size+1 pointing to
		     portions of _data */ 

  /* functions */
  
  /* constructor from a pattern label */
  nfa(alphabet &A,
      //nfa(std::string &alphabet_label,
      std::string &pattern_label,
      bool verbose=true,
      bool debug=false);

  /* destructor */
  ~nfa();

  /* from and to are two array of size _nstates, 0 <= a <=
     _alphabet_size, computes in to the image of from by a-transition,
     returns true if to is modified from init (all false for
     alpha-transitions and a copy of from for epsilon-transition);
     content of from unchanged */
  bool delta(bool *from,
	     unsigned int a,
	     bool *to); 

  /* compute the epsilon_closure of from in to (array previously
     allocated); returns true if computations done. from is changed in
     the process */
  bool epsilon_closure(bool *from,
		       bool *to);

  /* export the nfa to a file in dot format */
  void dot(std::string &file);

  /* add to the nfa all patterns with at most nerr */
  /* not yet implemented */
  // void hamming(unsigned short nerr);

};
};
#endif /* nfa.h */
