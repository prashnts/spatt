/* $Id: pattern.cc 504 2005-10-19 11:53:32Z mhoebeke $ */
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
/*  Class pattern                                                    */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "pattern.h"

using namespace std;

namespace spatt {
  
  node::node(char *word,int size,int deep) : _size(size), _deep(deep) {
    sprintf(_word,word);
    _next=(struct node **)malloc(sizeof(struct node *)*_size);
    if (_next==NULL) {
      fprintf(stderr,"Critical error in node::node ! Not enough memory\n");
      exit(EXIT_FAILURE);
    }
    for (int i=0; i<_size; i++)
      _next[i]=NULL;
  };
  
  node::~node(){
    if (_next!=NULL)
      free(_next);
  };
  
/* constructor from a label */
pattern::pattern(const spattparameters &params, const alphabet &alpha,
		 const char *label) :
  _params(&params), _alpha(&alpha), _count(-1), 
  _expected(-1.0), _optimized_word_list(NULL),
	_optimized_word_list_size(0)
{
  //sprintf(_label,label);
  strcpy(_label,label);

  /* splitting */
  char *split;
  char label_copy[MAX_STRING_LENGTH];
  //sprintf(label_copy,_label);
  strcpy(label_copy,_label);
  split=strtok(label_copy,"|");

  /* checking for --both-strands option */
  if (params.get_both_strands()){
    strcat(_label,"|{rc}");
  }

  /* parsing first part */
  _root=new struct node("",_alpha->get_size(),0);
  //printf("parse(%s)\n",split);
  parse(_root,split);
  read_tree(_root);
  _root=NULL;

  /* parsing remaining */
  while (split=strtok(NULL,"|")){
    _root=new struct node("",_alpha->get_size(),0);
    //printf("parse(%s)\n",split);
    parse(_root,split);
    read_tree(_root);
    _root=NULL;
  }

  /* display for verif */
  if (_params->get_debug_level()>2) {
    printf("pattern \"%s\": _count=%li\n",_label,_count);
    for (unsigned int i=0; i<_label_list.size(); i++)
      printf("  %s\n",_label_list[i].c_str());
  }
};

/* display _label_list */
void pattern::display(){
   printf("pattern \"%s\":\n",_label);
   for (unsigned int i=0; i<_label_list.size(); i++)
      printf("  %s\n",_label_list[i].c_str());
};

/* recursive parser */
void pattern::parse(struct node* root,char *label){
  //printf("parse(%p,%s)\n",root,label);
  /* if there is still something to parse */
  if (label[0]!='\0') {
    if (label[0]=='.') {
      for (int j=0; j<_alpha->get_size(); j++)
	update_parse(_alpha->code2char(j),root,&label[1]);
    } else if (label[0]=='[') {
      /* looking for next ']' */
      int pos=0;
      while (label[pos]!='\0' && label[pos]!=']') {
	//printf("label[pos]=\'%c\'\n",label[pos]);
	pos++;
      }
      //printf("\']\' pos=%i\n",pos);
      if (label[pos]=='\0') {
	fprintf(stderr,"Warning in pattern:parse ! no ']' matching '['\n");
      } else {
	for (int j=1; j<pos; j++) {
	  //printf("update_parse(\'%c\')\n",label[j]);
	  update_parse(label[j],root,&label[pos+1]);	  
	}
      }
    } else {
      update_parse(label[0],root,&label[1]);
    }
  }
  /* display for verif */
  //printf("pattern \"%s\":\n",_label);
  //for (int i=0; i<_nword; i++)
  //  printf("%s\n",_label_list[i].c_str());
};

/* updating parsing with c and continue with label */
void pattern::update_parse(char c,struct node *root,char *label){
  //printf("update_parse(%c,%p,%s)\n",c,root,label);
  struct node *next;
  char new_word[MAX_STRING_LENGTH];
  sprintf(new_word,root->_word);
  new_word[root->_deep]=c;
  new_word[root->_deep+1]='\0';
  next=new node(new_word,_alpha->get_size(),root->_deep+1);
  int code=_alpha->char2code(c);
  if (code<0) {
    fprintf(stderr,"Warning in pattern:update_parse ! \'%c\' is an invalid char in \"%s\"\n",c,_label);
    delete next;
  } else {
    root->_next[code]=next;
    parse(next,label);
  }
};

  /* recursive display function */
void pattern::read_tree(struct node *current){
  //printf("node(%p):%s\n",current,current->_word);
  int test=0;
  for (int i=0;i<current->_size;i++) {
    if(current->_next[i]!=NULL) {
      test=1;
      read_tree(current->_next[i]);      
    }
  }
  if (test==0) {
    //printf("%s\n",current->_word);
    string tmp(current->_word);
    _label_list.push_back(tmp);
  }
  //printf("delete(%p)\n",current);
  delete current;
};

/* copy constructor */
pattern::pattern(const pattern &source){
  //fprintf(stderr,"Copy constructor call for object pattern\n");
  //exit(EXIT_FAILURE);
  _params=source._params;
  _alpha=source._alpha;

  sprintf(_label,source._label);
  _count=source._count;
  _expected=source._expected;
  _label_list=source._label_list;
  _optimized_word_list_size=source._optimized_word_list_size;
  if (_optimized_word_list_size) {
  	_optimized_word_list=(word **)malloc(sizeof(word *)*_optimized_word_list_size);
	memcpy(_optimized_word_list,source._optimized_word_list,sizeof(word *)*_optimized_word_list_size);
  } else
	  _optimized_word_list=NULL;
  _root=NULL;
};


pattern &
pattern::operator=(const pattern &p) 
{
  if (this != &p) {
    _params=p._params;
    _alpha=p._alpha;
    
    sprintf(_label,p._label);
    _count=p._count;
    _expected=p._expected;
    if (_optimized_word_list_size)
	  free(_optimized_word_list);
    _optimized_word_list_size=p._optimized_word_list_size;
    if (_optimized_word_list_size) {
  	_optimized_word_list=(word **)malloc(sizeof(word *)*_optimized_word_list_size);
	memcpy(_optimized_word_list,p._optimized_word_list,sizeof(word *)*_optimized_word_list_size);
    }
    _label_list=p._label_list;
    _root=NULL;
  }

  return *this;
}

/* destructor */
pattern::~pattern(){
  if (_optimized_word_list_size)
    free(_optimized_word_list);
};

};
