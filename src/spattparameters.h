/* $Id: spattparameters.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class parameter allows reads the command line and returns the    */
/*  the parameters.                                                  */
/*  								     */
/*  This class uses argtable2                 			     */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef SPATTPARAMETERS_H
#define SPATTPARAMETERS_H

#include <string>
#include <vector>

#include "spattglobals.h"

using namespace std;

namespace spatt {

class spattparameters  {

public:

  virtual ~spattparameters() {};

  inline bool get_all_words() const {
    return _all_words;
  }

  inline const char *get_alphabet_label() const {
    return _alphabet_label.c_str();
  }

  inline bool get_both_strands() const {
    return _both_strands;
  }
  inline void set_both_strands(bool bs) {
    _both_strands=bs;
  }
  

  inline const char * get_count_filename() const {
    return _count_filename.c_str();
  }

  inline int get_debug_level() const {
    return _debug_level;
  }

  inline int get_length() const {
    return _length;
  }

  inline const char * get_markov_filename() const {
    return _markov_filename.c_str();
  }

  inline int get_markov_order() const {
    return _markov_order;
  }

  inline const char * get_model_filename() const {
    return _model_filename.c_str();
  }

  inline bool get_normalize() const {
    return _normalize;
  }
			
  inline const char * get_pattern_filename() const {
    return _pattern_filename.c_str();
  }

  inline const char * get_sequence_filename() const {
    return _sequence_filename.c_str();
  }

  inline const char * get_stationary_filename() const {
    return _stationary_filename.c_str();
  }


  inline bool ignore_case() const {
    return _ignore_case;
  }

  inline bool use_count_file() const {
    return _use_count_file;
  }

  inline bool use_markov_file() const {
    return _use_markov_file;
  }

  inline bool use_model_file() const {
    return _use_model_file;
  }

  inline bool use_pattern_file() const {
    return _use_pattern_file;
  }

  inline bool use_sequence_file() const {
    return _use_sequence_file;
  }

  inline bool use_stationary_file() const {
    return _use_stationary_file;
  }


  inline vector<string> & get_pattern_label_list() {
    return _pattern_label_list;
  }
    
  inline double get_max_pvalue() const {
    return _max_pvalue;
  }
  
  inline long get_nobs() const {
    return _nobs;
  }


protected:
  spattparameters();

  bool _all_words;
  string _alphabet_label;
  bool _both_strands;
  string _count_filename;
  int _debug_level;
  bool _ignore_case;
  int _length;
  string _markov_filename;
  int _markov_order;
  string _model_filename;
  bool _normalize;
  string _pattern_filename;
  vector<string> _pattern_label_list;
  string _sequence_filename;
  string _stationary_filename;
  bool _use_count_file;
  bool _use_markov_file;
  bool _use_model_file;
  bool _use_pattern_file;
  bool _use_sequence_file;
  bool _use_stationary_file;
  long _nobs;
  double _max_pvalue;

  friend class spattparameterparser;

 private:
  spattparameters(const spattparameters &) {};
  spattparameters &operator=(const spattparameters &) { return *this; };
};

}
#endif
