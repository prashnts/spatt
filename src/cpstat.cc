/* $Id: cpstat.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "cpstat.h"

namespace spatt {

  cpstat::cpstat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
		const markov &M,const input &I) : sstat(params,alpha,seq,occ,M,I)//,_toto(42) 
  {
    /* add here specific stuff */
  };
  
  void cpstat::compute(word *w){
    this->sstat::compute(w);
    _a=cpA(_label,_M);
    if (_a>0) {
      /* real cp case */
      _stat=newcpstat(_natt,_a,_stat,_nobs);
    }
    post_comp();
  };
  
  void cpstat::compute(word *w1,word *w2){
    this->sstat::compute(w1,w2);
    fprintf(stderr,"Warning ! Compound Poisson statistic for degenerate pattern not yet implemented using binomial statistic instead\n");
  };
  
  void cpstat::compute(pattern *patt){
    this->sstat::compute(patt);
    if (patt->get_optimized_word_list_size()>1 || _params->get_both_strands()) {
      /* case of a degenerate pattern */
      fprintf(stderr,"Warning ! Compound Poisson statistic for degenerate pattern not yet implemented using binomial statistic instead\n");
    } else {
      /* case of a simple word */
      _a=cpA(_label,_M);
      if (_a>0) {
	_stat=newcpstat(_natt,_a,_stat,_nobs);
      }
    }
    post_comp();
  };  

  cpstat::~cpstat() {
  };

  void cpstat::print_format(FILE *fout) {
    fprintf(fout,"pattern\tnobs\tnatt\ta\tstat\n");
  };

  void cpstat::print_regular(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%.2e\t%f\n",_label,_nobs,_natt,_a,_stat);
  };

  void cpstat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%.2e\t%e\n",_label,_nobs,_natt,_a,_stat);
  }

};
