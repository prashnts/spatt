/* $Id: stat.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class stat: class for statistics computations                    */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef STAT_H
#define STAT_H

#include <cstdio>
#include <cstring>
#include <cmath>
#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "sequence.h"
#include "count.h"
#include "word.h"
#include "markov.h"
#include "pattern.h"
#include "input.h"

namespace spatt {

  const int MAX_LABEL_SIZE = 500;

  class stat {
    
  public:
    stat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    virtual ~stat();
    /* compute the stat (here nobs) for simple word, word and ic or pattern */
    virtual void compute(word *w);
    virtual void compute(word *w1,word *w2);
    virtual void compute(pattern *patt);
    /* if called with no previous call to compute, displays the output format */
    virtual void print(FILE *fout,double threshold);
    virtual void print_format(FILE *fout);
    virtual void print_regular(FILE *fout);
    virtual void print_normalized(FILE *fout);
    virtual void post_comp();
    
  protected:
    const spattparameters *_params;
    const alphabet *_alpha;
    const sequence *_seq;
    const count *_occ;
    const markov *_M;
    const input *_In;
    long _nobs;
    double _natt;
    double _stat;
    double _pvalue;
    char _label[MAX_LABEL_SIZE];
    long _nmax;
    
  };
};
#endif
