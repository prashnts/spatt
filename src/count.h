/* $Id: count.h 676 2006-02-22 14:25:10Z gnuel $ */
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
/*  Class count: counting occurrences of words                       */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef COUNT_H
#define COUNT_H

#include <cstdio>
#include <cstdlib>
#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "sequence.h"

using namespace std;

namespace spatt {

class count {

 public:

  /* count from a countfile */
  count(const spattparameters &, const alphabet &);

  /* count from a sequence */
  count(const spattparameters &, const alphabet &, sequence &);

  /* returns count for code word of length h */
  long get(int h,long code) const;

  /* close stream */
  ~count();

  inline const long get_n() const { return _n; };

  inline const long get_length() const { return _length; };

  inline  sequence * get_seq() const {return _seq; };

 private:
  const spattparameters *_params;
  const alphabet *_alpha;
  sequence *_seq;
  FILE *_stream;
  char _filename[MAX_STRING_LENGTH];
  long **_occ;
  long *_dim;
  long _memory;
  long _length;
  long _n;


  /* alloc and intitalize _occ */
  void init_occ();

  /* empty constructor */
  count();

  /* copy constructor */
  count(count &source);

  /* affectation operator */
  count & operator=(count &source);

};
};
#endif
