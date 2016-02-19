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
#ifndef DFA_H
#define DFA_H 1

#include <map>
#include <vector>
#include <iterator>
#include "nfa.h"

#define MAX_DISPLAY 10

namespace spatt {

class dfa {

 public:
  
  /* alphabet part */
  unsigned short _alphabet_size;
  alphabet *_Alpha;
  //std::string _alphabet_label;
  unsigned long _token;

  unsigned long _nstates;
  unsigned long _start;
  std::vector<unsigned long> _final;
  std::vector<unsigned long> _state;
  std::vector<bool> _is_final;
  
  //std::vector< std::map<unsigned long,unsigned long> > _delta; 
  std::vector< std::vector<unsigned long> > _delta; 

  std::vector< std::vector<unsigned long> > _partition;

  /* ambiguity stuff */
  unsigned short _m; // unambiguity level (default: 0)
  bool _renewal; // indicate if the dfa count renewal ocurrence or not (default: false)
  std::vector< std::map<std::string,bool> > _inv_delta; // _inv_delta[q] = a in A^m, exists p, delta(p,a)=q
  std::vector< std::map<unsigned long,bool> > _inv_Delta; // _inv_Delta[q] = p in Q, exists a in A, delta(p,a)=q
 
  /* functions */

  /* build the dfa from a nfa through subset construction */
  dfa(nfa &N,bool verbose=true,bool debug=false);

  ~dfa();

  /* providing x a bool[size], he function returns the vector of
     positions i where x[i] is true */
  std::vector<unsigned int> bool2vector(bool *x,unsigned int size);

  /* printf a vector of unsigned int */
  void print_vector(std::vector<unsigned int> V);
  void print_lvector(std::vector<unsigned long> V);
  
  /* print the dfa */
  void print();

  /* compute the minimum partition */
  void  minimize(bool verbose=true,bool debug=false);

  /* rebuild the dfa according to the minimum partition */
  void rebuild(bool verbose=true);

  /* returns the intersection of two (ordered) vectors */
  std::vector<unsigned long>
    intersection(std::vector<unsigned long> &v1,
		 std::vector<unsigned long> &v2);

  /* returns the difference of two (ordered) vectors (v2 included in
     v1) */
  std::vector<unsigned long>
    difference(std::vector<unsigned long> &v1,
		 std::vector<unsigned long> &v2);

  /* returns the union of v1 and v2 (both ordered vectors) */
  std::vector<unsigned long>
    reunion(std::vector<unsigned long> &v1,
	       std::vector<unsigned long> &v2);

  /* export the dfa to a file in dot format */
  void dot(std::string file="dfa.dot");

  /* locate and count occurrences in a sequence. Returns a vector with
     ending positions of the occurrences */
  std::vector<unsigned long> locate_occ(std::string &file,bool verbose=true);

  // initialize non ambiguity vectors _inv_delta and _inv_Delta
  void initialize_ambiguity_vectors(bool verbose=true,bool debug=false);

  // print non ambiguity vectors _inv_delta and _inv_Delta
  void print_ambiguity_vectors();

  // remove local ambiguity at order m (do nothing if _m != m-1)
  void remove_local_ambiguity(unsigned short m,unsigned long state,bool verbose=false,bool debug=false);
  
  // return the size m string corresponding to delta^(-m)(state)
  // abort if m>_m
  std::string inv_delta(unsigned short m,unsigned long state);

  // print a string of unsigned short using the alphabet 
  void print_string(std::string w);
  
  // remove order m>=1 ambiguity
  void remove_ambiguity(unsigned short m,bool verbose=true,bool debug=false);

  // switch automaton to renewal counting
  // not reversible
  void renewal(bool verbose);

  // return a vector with delta(_start,A^_m)
  std::vector<unsigned long> starts(bool verbose);
  void recursive_starts(std::string &path, std::vector<unsigned long> &res, unsigned long pos, unsigned short m, unsigned short i);

  // return the code corresponding to a string of short
  unsigned long code(std::string &word);

  inline unsigned long reset() {
    _token=_start;
    return _token;
  };

  inline unsigned long process(unsigned short a){
    _token=_delta[a][_token];
    return _token;
  };

  inline bool is_final(unsigned long p) {
    return _is_final[p];
  };

};

};

#endif 
