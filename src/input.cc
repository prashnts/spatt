/* $Id: input.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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
/*  Class input: for core program inputs                             */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "input.h"

namespace spatt {

/* main constructor */
input::input(spattparameters &params,
	     const alphabet &alpha,
	     sequence &seq,
	     spatt::count &count) :
  _params(&params), _alpha(&alpha), _seq(&seq), _count(&count)
{
  
  /* determines mode */
  if (_params->use_pattern_file()==1) {
    /* open the file and read patterns */
    //printf("todo: read pattern file \"%s\"\n",_params->_pattern_filename);
    /* open the stream */
    _stream=fopen(_params->get_pattern_filename(),"r");
    if (_stream==NULL) {
      fprintf(stderr,"Cannot open sequence file \"%s\"\n",_params->get_pattern_filename());
      exit(EXIT_FAILURE);
    }
    /* read patterns and add them to _pattern_label_list */
    char *test=NULL;
    char buffer[BUFFER_SIZE];
    test=fgets(buffer,BUFFER_SIZE,_stream);
    std::string pattern;
    while (test!=NULL) {
      /* do nothing if comment */
      if (buffer[0]!=COMMENT_CHAR) {
	pattern=std::string(buffer);
	/* remove last char */
	pattern=pattern.substr(0,pattern.length()-1);
	/* add pattern to list */
	_params->get_pattern_label_list().push_back(pattern);
      }
      /* gets new line */
      test=fgets(buffer,BUFFER_SIZE,_stream);
    }
  }

  /* creating _pattern_list */  
  if (_params->get_pattern_label_list().size()>0) {
    for (int i=0; i<_params->get_pattern_label_list().size()>0; i++) {
      _pattern_list.push_back(pattern(*_params,*_alpha,_params->get_pattern_label_list()[i].c_str()));
    }
  }

  /* verification */
//  if (_pattern_list.size()>0) {
//    for (int i=0; i<_pattern_list.size()>0; i++) {
//      printf("pattern(%i)=\"%s\"\n",i,_pattern_list[i]._label);
//      printf("pattern(%i)._is_counted=%i\n",i,_pattern_list[i]._is_counted);
//      _pattern_list[i].display();
//    }
//  }
  
  //printf("creating words\n");

  if (_params->get_both_strands()==1 && _alpha->has_complementary()==1) {
    /* reading all patterns */
    pattern *current=NULL;
    for (int i=0; i<_pattern_list.size(); i++) {
      current=&_pattern_list[i];
      //printf("pattern \"%s\"\n",current->_label);
      int jmax=current->get_label_list().size();
      std::string w;
      char w1[MAX_STRING_LENGTH];
      char w2[MAX_STRING_LENGTH];
      for (int j=0; j<jmax; j++) {
	w=current->get_label_list()[j];
	//printf("  w=\"%s\"\n",w.c_str());
	/* rewrite w with non ambigous char */
	int code;
	for (int k=0; k<w.size(); k++) {
	  code=_alpha->char2code(w[k]);
	  if (code<0) {
	    w1[k]='\0';
	    k=w.size();
	  } else {
	    w1[k]=_alpha->code2char(code);
	  }
	}
	w1[w.size()]='\0';
	//printf("  w1=\"%s\"\n",w1);
	/* create w2 inverse complementary of w1 */
	{ 
	  int kmax=strlen(w1);
	  //printf("kmax=%i\n",kmax);
	  for (int k=0; k<kmax;k++) {
	    //printf("w1[kmax-k-1]=\'%c\'\n",w1[kmax-k-1]);
	    w2[k]=_alpha->complement_char(w1[kmax-k-1]);
	    //printf("w2[k]=\'%c\'\n",w2[k]);
	  }
	  w2[kmax]='\0';
	}
	//printf("  w2=\"%s\"\n",w2);
	/* add word if necessary */
	if (strcmp(w1,w2)!=0) {
	  //printf("added\n");
	  current->get_label_list().push_back(w2);
	}
      } // end word loop
    } // end pattern loop
  } // end both-strands

  /* reading _pattern_list creating words */
  if (_pattern_list.size()>0) {
    pattern *current=NULL;
    int code;
    for (int i=0; i<_pattern_list.size(); i++) {
      current=&_pattern_list[i];
      //printf("working with pattern(%i)=\"%s\"\n",i,current->_label);
      /* reading _label_list */
      std::string label;
      word *current_word;
      word *test;
      for (int j=0; j<current->get_label_list().size(); j++) {
	label=current->get_label_list()[j];
	//printf("  working with label(%i)=\"%s\"\n",j,label.c_str());
	code=test_word(label);
	//printf("    code=%i\n",code);
	switch(code) {
	case(INVALID_WORD):
	  if (_params->get_debug_level()>0)
	    fprintf(stderr,"Warning ! word \"%s\" is invalid and will not be treated\n",label.c_str());
	  break;
	case(COUNT_WORD):
	  if (_params->get_debug_level()>1)
	    printf("\"%s\" is a valid count pattern\n",label.c_str());
	  current_word = new word(*_params,*_alpha,*_seq,*_count,(char *)label.c_str());
	  /* two cases */
	  test=&_all_word_table[current_word->get_label()];
	  if (test->get_alpha()==NULL) {
	    /* it is a new word */
	    _all_word_table[current_word->get_label()]=*current_word;
	    current->get_word_list().push_back(current_word->get_label());
	    _all_word_list.push_back(current_word->get_label());
	    _count_word_list.push_back(current_word->get_label());
	  } else {
	    /* it is an old word */
	    if (_params->get_debug_level()>1)
	      printf("its an old word, deletion\n");
	    delete current_word;
	    current->get_word_list().push_back(test->get_label());
	    //_all_word_list.push_back(test);
	    //_count_word_list.push_back(test);
	  }
	  current_word=NULL;
	  break;
	case(DFA_WORD):
	  if (_params->get_debug_level()>1)
	    printf("\"%s\" is a valid dfa pattern\n",label.c_str());
	  current_word = new word(*_params,*_alpha,*_seq,*_count,(char *)label.c_str());
	  //printf("ok\n");
	  /* two cases */
	  test=&_all_word_table[current_word->get_label()];
	  //test2=_all_word_table[current_word->_label];
	  //printf("test(%p)=\"%s\"\n",test,test->_label);
	  if (test->get_alpha()==NULL) {
	    /* it is a new word */
	    current_word->build_dfa();	    
	    _all_word_table[current_word->get_label()]=*current_word;
	    current->get_word_list().push_back(current_word->get_label());
	    _all_word_list.push_back(current_word->get_label());
	    _dfa_word_list.push_back(current_word->get_label());
	  } else {
	    /* it is an old word */
	    delete current_word;
	    current->get_word_list().push_back(test->get_label());
	    //_all_word_list.push_back(test);
	    //_dfa_word_list.push_back(test);
	  }
	  current_word=NULL;
	  break;
	}	
      }
    }
  }

  /* verification at the end */
  //printf("Verification at end of input: _all_word_list\n");
  //{
  //  int imax=_all_word_list.size();
  //  std::string current;
  //  for (int i=0; i<imax; i++) {
  //    current=_all_word_list[i];
  //    printf("current(%p)=\"%s\"\n",&_all_word_table[current],_all_word_table[current]._label);
  //  }
  //}
  //{
  //  printf("_pattern_list\n");
  //  int imax=_pattern_list.size();
  //  pattern *current;
  //  for (int i=0; i<imax; i++) {
  //    current=&_pattern_list[i];
  //    printf("pattern \"%s\"\n",current->_label);
  //    for (int j=0; j<current->_word_list.size(); j++)
  //      printf("  word(%p)=\"%s\"\n",&_all_word_table[current->_word_list[j]],_all_word_table[current->_word_list[j]]._label);
  //  }
  //}

  /* creating _optimized_*_word pointing on map */
  _optimized_all_word_size=_all_word_list.size();
  _optimized_all_word=NULL;
  if (_optimized_all_word_size>0) {
    _optimized_all_word=(word **)malloc(sizeof(word *)*_optimized_all_word_size);
    if (_optimized_all_word==NULL) {
      fprintf(stderr,"Critical error in input::input ! Not enough memory\n");
      exit(EXIT_FAILURE);
    }
  }
  _optimized_count_word_size=_count_word_list.size();
  _optimized_count_word=NULL;
  if (_optimized_count_word_size>0) {
    _optimized_count_word=(word **)malloc(sizeof(word *)*_optimized_count_word_size);
    if (_optimized_count_word==NULL) {
      fprintf(stderr,"Critical error in input::input ! Not enough memory\n");
      exit(EXIT_FAILURE);
    }
  }
  _optimized_dfa_word_size=_dfa_word_list.size();
  _optimized_dfa_word=NULL;
  if (_optimized_dfa_word_size>0) {
    _optimized_dfa_word=(word **)malloc(sizeof(word *)*_optimized_dfa_word_size);
    if (_optimized_dfa_word==NULL) {
      fprintf(stderr,"Critical error in input::input ! Not enough memory\n");
      exit(EXIT_FAILURE);
    }
  }
  {
    int count_pos=0;
    int dfa_pos=0;
    for (int i=0; i<_optimized_all_word_size; i++) {
      _optimized_all_word[i]=&_all_word_table[_all_word_list[i]];
      if (_optimized_all_word[i]->get_dfa()==NULL) {
	_optimized_count_word[count_pos]=_optimized_all_word[i];
	count_pos++;
      } else {
	_optimized_dfa_word[dfa_pos]=_optimized_all_word[i];
	dfa_pos++;
      }
    }
  }

  /* same thing inside each pattern */
  {
    pattern *current=NULL;
    int imax=_pattern_list.size();
    /* loop on patterns */
    for (int i=0; i<imax; i++) {
      current=&_pattern_list[i];
      /* alloc */
      current->set_optimized_word_list_size(current->get_word_list().size());
      current->set_optimized_word_list((word **)malloc(sizeof(word *)*current->get_optimized_word_list_size()));
      if (current->get_optimized_word_list()==NULL) {
	fprintf(stderr,"Critical error in input::input ! Not enough memory\n");
	exit(EXIT_FAILURE);
      }
      for (int j=0; j<current->get_optimized_word_list_size(); j++) {
	current->get_optimized_word_list()[j]=&_all_word_table[current->get_word_list()[j]];
      }
    }
  }
};

/* test a pattern descriptor */
/* returns: INVALID_WORD if invalid */
/* COUNT_WORD if size included in count */
/* DFA_WORD else */
int input::test_word(std::string word){
  
  /* testing word contain only valid char */
  {
    int test=1;
    for (int i=0; i<word.length(); i++) {
      if (_alpha->char2code(word[i])<0)
	test=0;
    }
    if (test==0)
      return INVALID_WORD;
  }
  /* sorting with word length */
  if (word.length()<=_count->get_length()) {
    /* it is a count word */
    return COUNT_WORD;
  } else {
    /* it is a dfa word */
    return DFA_WORD;
  }
}


/* destructor */
input::~input(){
//  for (int i=0; i<_count_word_list.size(); i++) {
//    //printf("count_list(%i)=%p\n",i,_count_word_list[i]);
//    delete _count_word_list[i];
//  }
//  for (int i=0; i<_dfa_word_list.size(); i++) {
//    //printf("dfa_list(%i)=%p\n",i,_dfa_word_list[i]);
//    delete _dfa_word_list[i];
//  }
  if (_optimized_all_word!=NULL)
    free(_optimized_all_word);
  if (_optimized_count_word!=NULL)
    free(_optimized_count_word);
  if (_optimized_dfa_word!=NULL)
    free(_optimized_dfa_word);
};

};
