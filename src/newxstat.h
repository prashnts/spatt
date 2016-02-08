/* $Id: xstat.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#ifndef NEWXPSTAT_H
#define NEWXPSTAT_H

#include "sstat.h"
#include "transition.h"

#define LIMIT_LOG_COMP 290.0

namespace spatt {

  class newxstat : public sstat {
    
  public:
    newxstat(spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    virtual ~newxstat();
    /* compute the stat (here nobs) for simple word, word and ic or pattern */
    virtual void compute(word *w);
    virtual void compute(word *w1,word *w2);
    virtual void compute(pattern *patt);
    /* if called with no previous call to compute, displays the output format */
    virtual void print_format(FILE *fout);
    virtual void print_regular(FILE *fout);
    virtual void print_normalized(FILE *fout);
      
  private:
    transition *_T;
    long _order; // order of the transition matrix
    long _nblock; // number of block in the fmci matrix
    long _dim; // _dim = _order * _nblock
    long _n; // sequence length
    double *_workspace; // workspace of size 2*_noblock*_order+1*_order
    double *_current_u,*_last_u; // vectors of size _nblock*_order (in _workspace)
    double *_sum; // aux vector of size _order (in _workspace)
    unsigned long _workspace_size; // size of workspace in sizeof(double)
    /* main routines */
    void comp(); // computations
    void log_comp(); // log computations
  };
};
#endif
