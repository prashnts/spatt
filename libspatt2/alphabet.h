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
#ifndef ALPHABET_H
#define ALPHABET_H 1

#include <cstdio>
#include <string>

#define UNDEF -1
#define SEPARATOR -2 
#define REPEAT -3
#define CHOICE -4
#define NCHOICE -5
#define END_BLOCK -6
#define ANY_CHAR -7
#define SEQUENCE -8
#define COMMENT -9
#define IGNORE -10

namespace spatt {

class alphabet {

 private:
  
  unsigned short _alphabet_size;
  std::string _alphabet_label;
  int _alphabet_code[256];

 public:

  alphabet(std::string &alphabet_label);

  inline int code(char c) {
    return  _alphabet_code[(unsigned short)c];
  };

  inline char decode(unsigned short i) {
    return _alphabet_label[i];
  }

  inline unsigned short size() {
    return _alphabet_size;
  };

  inline std::string & label() {
    return _alphabet_label;
  };

};

};
#endif
