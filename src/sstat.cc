/* $Id: sstat.cc 705 2006-03-09 07:49:19Z gnuel $ */
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

#include "sstat.h"

namespace spatt {

  sstat::sstat(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,
	       const markov &M,const input &I) : stat(params,alpha,seq,occ,M,I) //,_toto(42) 
  {
    /* add here specific stuff */
    //printf("sstat object created: _seq=%p\n",_seq);
    //printf("sstat object created: _seq->_nseq=%i\t_seq->_valid_char=%i\n",_seq->_nseq,_seq->_valid_char);
  };

  void sstat::compute(word *w){
    this->stat::compute(w);
    _natt=_M->expect(*w);
    _h=w->get_size();
    comp();
  };
  
  void sstat::compute(word *w1,word *w2){
    this->stat::compute(w1,w2);
    _natt=_M->expect(*w1);
    _natt+=_M->expect(*w2);
    _h=w1->get_size();
    {
      int h=w1->get_size();
      if (h>_h)
	_h=h;
    }
    comp();
  };
  
  void sstat::compute(pattern *patt){
    word *w=NULL;
    _h=0;
    if (patt->get_count()<0) {
      patt->set_count(0);
      patt->set_expected(0.0);
      for (int j=0; j<patt->get_optimized_word_list_size(); j++) {
	w=patt->get_optimized_word_list()[j];
	if (w->get_size()>_h)
	  _h=w->get_size();
	patt->set_count(patt->get_count()+w->get_count());
	patt->set_expected(patt->get_expected()+w->get_expected());
      }
    }
    if (_params->get_nobs()>-1) {
      _nobs=_params->get_nobs();
    } else {
      _nobs=patt->get_count();
    }
    _natt=patt->get_expected();
    strcpy(_label,patt->get_label());
    comp();
  };

  void sstat::comp(){
    //printf("_h=%i\n",_h);
    //printf("nexp=%f\n",_natt);
    double m,e,r,s;
    //long n=_occ->get_n()-_h+1;
    //printf("nseq=%i\n",_seq->_nseq);
    //printf("sstat::comp(): _seq=%p\n",_seq);
    //printf("sstat::comp(): _seq->_nseq=%i\t_seq->_valid_char=%i\n",_seq->_nseq,_seq->_valid_char);
    long n=_seq->_valid_char-_seq->_nseq*(_h-1);
    //printf("ok\n");
    if (_natt!=0.0) {
      if (_nobs>_natt) {
	/* pattern is over-represented */
	//printf("nobs=%li\tnatt=%f\tn=%li\n",_nobs,_natt,n);
	r=qcdfbin(_nobs,
		  n,
		  _natt/n,&m,&e);
	s=-e-log(m)/log(10.0);
      } else {
	/* pattern is under-represented */
	//printf("nobs=%li\tnatt=%f\tn=%li\n",_nobs,_natt,n);
	r=pcdfbin(_nobs,n,
		  _natt/n,&m,&e);
	s=e+log(m)/log(10.0);
      }
      //_pvalue=exp(-fabs(s)*log(10));
      _stat=s;
    } else {
      fprintf(stderr,"pattern %s expected 0 times skipped\n",_label);
      _pvalue=0.0;
      _stat=0.0;
    }
    ///* check normalize */
    //if (_params->get_normalize()) {
    //  _stat/=(double)n;
    //}
    post_comp();
  };

  sstat::~sstat() {
  };
  
  void sstat::print_format(FILE *fout) {
    fprintf(fout,"pattern\tnobs\tnatt\tstat\n");
  };

  void sstat::print_regular(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%f\n",_label,_nobs,_natt,_stat);
  };

  void sstat::print_normalized(FILE *fout) {
    fprintf(fout,"%s\t%li\t%.2f\t%e\n",_label,_nobs,_natt,_stat);
  }

};
