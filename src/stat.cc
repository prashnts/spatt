/* $Id: stat.cc 455 2005-09-01 14:27:20Z gnuel $ */
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

#include "stat.h"

namespace spatt {

  stat::stat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
	     const markov &M,const input &I) : _params(&params), _alpha(&alpha), _seq(&seq), 
					       _occ(&occ), _M(&M), _In(&I), _nobs(-1), _natt(-1.0),
					       _pvalue(-1.0), _stat(0.0)
  {
    //printf("_nobs=%i\n",_nobs);
  };

  stat::~stat() {
  };

  void stat::compute(word *w){
    if (_params->get_nobs()>-1) {
      _nobs=_params->get_nobs();
    } else {
      // does not work if w from a dfa pattern
      _nobs=_occ->get(w->get_size(),w->get_code());
    }
    strcpy(_label,w->get_label());
  };
  
  void stat::compute(word *w1,word *w2){
    if (_params->get_nobs()>-1) {
      _nobs=_params->get_nobs();
    } else {
      _nobs=_occ->get(w1->get_size(),w1->get_code());
      _nobs+=_occ->get(w2->get_size(),w2->get_code());
    }
    strcpy(_label,w1->get_label());
    strcat(_label,"|");
    strcat(_label,w2->get_label());
  };

  void stat::compute(pattern *patt){
    
    if (_params->get_nobs()>-1) {
      _nobs=_params->get_nobs();
    } else {
      word *w=NULL;
      if (patt->get_count()<0) {
	patt->set_count(0);
	int h=0;
	for (int j=0; j<patt->get_optimized_word_list_size(); j++) {
	  w=patt->get_optimized_word_list()[j];
	  if (h<w->get_size())
	    h=w->get_size();
	  patt->set_count(patt->get_count()+w->get_count());
	}
      }
      _nobs=patt->get_count();
    }
    strcpy(_label,patt->get_label());
    
  };

  void stat::print(FILE *fout,double threshold){
    if (_nobs<0) {
      print_format(fout);
    } else {
      if (_pvalue<=threshold) {
	if (_params->get_normalize())
	  print_normalized(fout);
	else
	  print_regular(fout);
      }
    }
  };
  
  void stat::post_comp() {

    _pvalue=exp(-fabs(_stat)*log(10.0));
    /* check normalize */
    if (_params->get_normalize()) {
      _stat/=(double)_occ->get_n();
    }

  };

  void stat::print_format(FILE *fout) {
    fprintf(fout,"pattern\tnobs\n");
  };

  void stat::print_regular(FILE *fout) {
    fprintf(fout,"%s\t%li\n",_label,_nobs);
  };

  void stat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\n",_label,_nobs);
  }

};
