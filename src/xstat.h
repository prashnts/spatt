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

#ifndef XPSTAT_H
#define XPSTAT_H

#include "x_WordFam.h"
#include "x_PSucceed.h"
#include "x_PAppearFast.h"

#include "cpstat.h"


//#include <string>
//using namespace std;


namespace spatt {

  const int XSTAT_SECURITY_MARGIN = 0;

  const double XSTAT_TARGET_RELATIVE_PRECISION = 1e-5;
  const double  XSTAT_MACHINE_RELATIVE_PRECISION = 1e-16;

  class xstat : public cpstat {
    
  public:
    xstat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I);
    virtual ~xstat();
    /* compute the stat (here nobs) for simple word, word and ic or pattern */
    virtual void compute(word *w);
    virtual void compute(word *w1,word *w2);
    virtual void compute(pattern *patt);
    ///* if called with no previous call to compute, displays the output format */
    virtual void print_format(FILE *fout);
    virtual void print_regular(FILE *fout);
    virtual void print_normalized(FILE *fout);
      
    inline long int get_rank() const { return _rank; }
    inline double *get_power() const { return _power; }

  private:
    double *_data; // workspace
    long _size_data; // workspace size 
    long _ell; // length of the sequence
    long _ell0; // limit for tail sums
    double *_ac; // auto-correlation (connected to _data)
    double *_current; // connected to _data
    double *_old; // connected to _data
    int _k; //size of the pattern
    bool _tail_sum; // true if tail sum is required
    void init(); // check tail sum or no, check memory requirements and connect
    void comp(); // compute _ac, process recurrence and get the result
    void comp_ell();
    void x_comp();

    long _order;
    long _rank;
    double *_power;
    long _size_power;
    long _n_power;
    long _one_size;

    void comp_rank();
    void update_power();
    
    WordFam *_W;
    PAppearFast *_P_W;

  protected:
    /* added here new members */    
  };
};
#endif
