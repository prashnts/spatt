/* $Id: process.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "process.h"

namespace spatt {

  process::process(const spattparameters &params,const alphabet &alpha,sequence &seq,const count &occ,const markov &M,input &I,stat &S) : _params(&params), _alpha(&alpha), _seq(&seq), _occ(&occ), _M(&M), _In(&I), _Sn(&S){

  };

  process::~process() {
  };

  void process::count_all_dfa(){
    
    /* getting dfa_word_list size */
    int imax=_In->get_optimized_dfa_word_size();
    // read sequence only if necessary
    if (imax>0) {
      /* restarting sequence */
      _seq->restart();
      //std::string w;
      word *current;
      /* dfa initialization */
      for (int i=0; i<imax; i++) {
	//w=_In->_dfa_word_list[i];
	//current=&_In->_all_word_table[w];
	current=_In->get_optimized_dfa_word()[i];
	//printf("current(%p)=\"%s\"\n",current,current->_label);
	current->dfa_initialize_count();
      }
      /* processing sequence */
      int c = _seq->next();
      while (c!=EOF_CODE) {
	/* updating all dfa */
	for (int i=0; i<imax; i++) {
	  //w=_In->_dfa_word_list[i];
	  //current=&_In->_all_word_table[w];
	  current=_In->get_optimized_dfa_word()[i];
	  current->dfa_update_count(c);
	}      
	c = _seq->next();
      }    
    } // end if (imax>0)
    
    /* if a model is given */
    /* computing expectation for all dfa words */
    if ( (_M->get_model()!=NULL || _M->get_order()<=0) &&  _M->get_stationary()!=NULL) {
      word *current;
      int imax=_In->get_optimized_all_word_size();
      for (int i=0; i<imax; i++) {
	current=_In->get_optimized_all_word()[i];
	_M->expect(*current);
	//printf("%s\t%f\n",current->_label,current->_expected);
      }
    }
    
  };

  void process::select_case(){
    //printf("process: debug_level=%i\n",_params->get_debug_level());
    if (_params->get_debug_level()>0) {
      //printf("process: format\n");
      _Sn->print(stdout,_params->get_max_pvalue());
    }
    /* general case */
    if (_params->get_all_words() == 0) {
      all_patterns_stat();
    }    
    /* --all-words case */
    if (_params->get_all_words() == 1) {
      all_words_stat();
    }
  }
  

  void process::all_words_stat(){
    if (_params->get_both_strands()==0) {
      /* easy case, without both-strands option */
      char init_label[MAX_STRING_LENGTH];
      char *current_label;
      char c;
      double m,e,r,s;
      int h=_params->get_length();
      /* initialization */
      c=_alpha->code2char(0);
      for (int i=0; i<h; i++) {
	init_label[i]=c;
      }
      init_label[h]='\0';
      word current_word(*_params,*_alpha,*_seq,*_occ,init_label);
      current_label=current_word.get_label();
      /* loop on all words */
      int k=_alpha->get_size();
      long imax=(long)pow((float)k,h);
      int cc;
      for (int i=0; i<imax; i++) {
	/* update word */
	current_word.set_code(i);
	//current_word.set_count(_occ->get(h,i));
	/* compute stat */
	_Sn->compute(&current_word);
	/* output result */
	_Sn->print(stdout,_params->get_max_pvalue());
	//_Sn->print(stdout,1.0);
	/* output result */
	//printf("%s\t%li\n",current_word.get_label(),current_word.get_count());
	/* next word */
	int j=1;
	while (j<=h) {
	  //printf("j=%i\n",j);
	  cc=_alpha->char2code(current_label[h-j]);
	  cc++;
	  if (cc==k) {
	    cc=0;
	    current_label[h-j]=_alpha->code2char(cc);
	    j=j+1;
	  } else {
	    current_label[h-j]=_alpha->code2char(cc);
	    j=h+1;
	  }
	} // end while
      } // end for

    } else {
      /* tough case, with both-strands option */
      char init_label[MAX_STRING_LENGTH];
      char init_ic_label[MAX_STRING_LENGTH];
      char *current_label;
      char *current_ic_label;
      char c;
      double m,e,r,s;
      long total_count;
      int h=_params->get_length();
      /* initialization */
      c=_alpha->code2char(0);
      for (int i=0; i<h; i++) {
	init_label[i]=c;
	init_ic_label[h-1-i]=_alpha->complement_char(c);
      }
      init_label[h]='\0';
      init_ic_label[h]='\0';
      word current_word(*_params,*_alpha,*_seq,*_occ,init_label);
      word current_ic_word(*_params,*_alpha,*_seq,*_occ,init_ic_label);
      current_label=current_word.get_label();
      current_ic_label=current_ic_word.get_label();
      /* loop on all words */
      int k=_alpha->get_size();
      long imax=(long)pow((float)k,h);
      int cc;
      for (int i=0; i<imax; i++) {
	//printf("code=%li ic_code=%li\n",current_word._code,current_ic_word._code);
	if (current_word.get_code()<=current_ic_word.get_code()) {
	  /* update word */
	  current_word.set_code(i);
	  //current_word.set_count(_occ->get(h,i));
	  /* update ic_word */
	  //current_ic_word.set_count(_occ->get(h,current_ic_word.get_code()));
	  /* compute stat */
	  _Sn->compute(&current_word,&current_ic_word);
	  // _Sn->compute(h,i,current_ic_word.get_code(),current_word.get_label(),current_ic_word.get_label());
	  /* output result */
	  _Sn->print(stdout,_params->get_max_pvalue());
	  //_Sn->print(stdout,1.0);
	  /* update total */
	  //total_count=current_word.get_count()+current_ic_word.get_count();
	  /* output result */
	  //printf("%s|%s\t%li\n",current_word.get_label(),
	  // current_ic_word.get_label(),total_count);      	
	} // end _code <= ic_code
	/* next word */
	int j=1;
	int diff;
	long power=(long)pow((float)k,h);
	while (j<=h) {
	  //printf("j=%i\n",j);
	  power/=k;
	  cc=_alpha->char2code(current_label[h-j]);
	  diff=_alpha->complement_code(cc);
	  cc++;
	  if (cc==k) {
	    cc=0;
	    current_label[h-j]=_alpha->code2char(cc);
	    current_ic_label[j-1]=_alpha->complement_char(current_label[h-j]);
	    diff-=_alpha->complement_code(cc);
	    current_ic_word.set_code(current_ic_word.get_code()-(diff*power));
	    j=j+1;
	  } else {
	    current_label[h-j]=_alpha->code2char(cc);
	    current_ic_label[j-1]=_alpha->complement_char(current_label[h-j]);
	    diff-=_alpha->complement_code(cc);
	    current_ic_word.set_code(current_ic_word.get_code()-(diff*power));
	    j=h+1;
	  }
	} // end while
      } // end for    
    }
  };
  

  void process::all_patterns_stat(){
    pattern *current;
    for (unsigned int i=0; i<_In->get_pattern_list().size(); i++) {
      current=&_In->get_pattern_list()[i];
      /* compute stat */
      _Sn->compute(current);
      /* output result */
      _Sn->print(stdout,_params->get_max_pvalue());
      //_Sn->print(stdout,1.0);      
    };
  }


};
