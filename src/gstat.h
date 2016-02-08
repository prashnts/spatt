/* $Id: gstat.h 455 2005-09-01 14:27:20Z gnuel $ */
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

#ifndef GSTAT_H
#define GSTAT_H

#include "gspattparameters.h"
#include "stat.h"
#include "cdf.h"

namespace spatt {

  const double GSTAT_TARGET_RELATIVE_PRECISION = 1e-5;
  const double GSTAT_MACHINE_RELATIVE_PRECISION = 1e-16;
  const int SECURITY_MARGIN = 0;

  class gstat: public stat {
    
  public:
    gstat(const gspattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    virtual ~gstat();
    /* compute the stat (here nobs) for simple word, word and ic or pattern */
    virtual void compute(word *w);
    virtual void compute(word *w1,word *w2);
    virtual void compute(pattern *patt);
    ///* if called with no previous call to compute, displays the output format */
    virtual void print_format(FILE *fout);
    virtual void print_regular(FILE *fout);
    virtual void print_normalized(FILE *fout);

  private:
    double _sd; //standard deviation
    //long _n;
    //double _Vno,_Vo; // non overlapping and overlapping variance
    void comp();
    void precomp(word *w);

    double _overlap; // O(v,w)
    double _main; // M(v,w)
    int _g,_h; // length of v and w
    int _ell; // length of the sequence
    double comp_overlap(word *w);
    double comp_overlap(word *v,word *w);
    double overlap(word *v,word *w,int d);
    double fit(word *v,word *w,int i);
    double comp_main(word *w);
    double comp_main(word *v,word *w);
    long _rank;
    double *_power;
    long _size_power;
    long _n_power;
    long _one_size;
    long _order;
    double _zscore;

  protected:
    /* added here new members */
    
  };
};
#endif
