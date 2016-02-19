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
#ifndef SEQUENCE_H
#define SEQUENCE_H 1

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <iterator>
#include "alphabet.h"

#define BUFFER_SIZE 300
#define SEQEND -1
#define SEQEMPTY -2

namespace spatt {

class sequence {

 private:
  
  alphabet *_Alpha;
  std::vector<std::string> _sequence_file;
  std::vector<std::string>::iterator _current;
  FILE *_in;
  char _buffer[BUFFER_SIZE];
  int _buffer_size;
  int _pos;
  unsigned long _nvalidchar;
  unsigned long _lastnvalidchar;
  unsigned long _nvalidseq;
  bool _seqend;
  bool _seqinterrupt;
    
  // get to the next line
  // returns false if there is no more line
  bool next_line();

  // get to the next valid file
  // returns false if there is no more files
  bool next_file();

 public:

  sequence(alphabet &A,std::vector<std::string> &sequence_file);
  
  sequence(alphabet &A,std::string &sequence_file);

  void reset();
  
  // returns the next caracter code
  // SEQEND if a sequence end
  // SEQEMPTY their is no more char
  int next();
  
  inline unsigned long nvalidchar() {
    return _nvalidchar;
  };

  inline unsigned long lastnvalidchar() {
    return _lastnvalidchar;
  };

  inline unsigned long nvalidseq() {
    return _nvalidseq;
  };

  ~sequence();

};

};
#endif
