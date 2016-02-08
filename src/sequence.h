/* $Id: sequence.h 676 2006-02-22 14:25:10Z gnuel $ */
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
/*  Class sequence                                                   */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <cstdio>
#include <cstdlib>
#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"

using namespace std;

namespace spatt {

  const int BUFFER_DONE=-1;

class sequence {

 private:
  char _filename[MAX_STRING_LENGTH];
  FILE *_stream;
  long _length;
  char _buffer[BUFFER_SIZE];
  int _position;

  const alphabet *_alpha;
  const spattparameters *_params;

  sequence();
  sequence(const sequence &);
  sequence &operator=(const sequence &);

 public:
  long _nseq;
  long _valid_char;

  long inline get_length() const {
    return _length;
  }

  /* opens the stream and computes length */
  sequence(const spattparameters &,const alphabet &);

  /* returns next code */
  int next();

  /* return filename */
  inline char * get_filename() {
    return _filename;
  }
  
  /* return valid char */
  inline long get_valid_char() const {
    return _valid_char;
  }
  
  /* restart sequence */
  void restart();

  /* close stream */
  ~sequence();

};
};

#endif
