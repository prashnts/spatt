/* $Id: markov.h 705 2006-03-09 07:49:19Z gnuel $ */
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

/*********************************************************************/   /*  								     */
/*  Class Markov: for Markov chains                                  */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

/* todo: finish the zgeev_par class and move all *_par classes in another source files */

#ifndef MARKOV_H
#define MARKOV_H

#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "sequence.h"
#include "count.h"
#include "word.h"
#include "fortran.h"

using namespace std;

namespace spatt {

class markov {

 public:

  /* markov from a markovfile */
  markov(const spattparameters &param, const alphabet &);

  /* markov from a count */
  markov(const spattparameters &param, const alphabet &alpha, const count &);

  /* empty constructor */
  markov();

  /* affectation operator */
  markov & operator=(const markov &source);

  /* close stream */
  ~markov();

  /* compute expected count for a word */
  double expect(word &w) const;

  /* compute expected count for a label */
  double expect(char *label) const;

  /* compute P(fromto|from) */
  double trans(const char *from,const char *to) const;

  /* compute augmented model */
  void compute_augmented(int order);

  inline const alphabet *get_alpha() const { return _alpha; }

  inline const long get_dim() const { return _dim; }

  inline const long get_DIM() const { return _DIM; }

  inline double **get_augmented_model() const 
    { return _augmented_model; }

  inline int get_augmented_order() const 
    { return _augmented_order; }

  inline double *get_stationary() const { return _stationary; }

  inline const int get_order() const { return _order; }

  inline const long get_stationary_size() const {return _n; }

  inline void set_occ(count &occ) { _occ=&occ; }

  inline void set_seq(sequence &seq) { _seq=&seq; }

  inline const count * get_occ() const { return _occ; }

  inline double **get_model() const { return _model; }

  inline const spattparameters *get_params() const { return _params;}

  inline const long get_n() const { return _n; }

  inline const double get_secondmag() const { return _secondmag; }

  inline const double get_trans(long i,int j) const {
    if (_order>=1)
      return _model[i][j];
    else
      return _stationary[j];
  }

  /* compute probability of one occurrence for a label */
  double mu(const char *label) const;
  double mu(string *s) const;
  double mu(long code,int size) const;

 private:
  const spattparameters *_params;
  const alphabet *_alpha;

  sequence *_seq;
  const count *_occ;

  FILE *_stream;
  double **_model;
  double **_augmented_model;
  int _order;
  int _augmented_order;
  double *_stationary;
  long _n;
  long _dim;
  int _k;
  char _filename[MAX_STRING_LENGTH];
  double _secondmag;

  /* augmented model dimension */
  long _DIM;

  /* copy constructor */
  markov(const markov &) {};

  /* compute stationary distribution through Arnoldi */
  void compute_stationary();

  /* dump _model in a file */
  void dump_model(const char* filename);

  /* dump stationary in a file */
  void dump_stationary(const char* filename);

};

};
#endif
