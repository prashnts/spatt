/* $Id: pattern.h 504 2005-10-19 11:53:32Z mhoebeke $ */
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
/*  Class pattern parsing pattern descriptors                        */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef PATTERN_H
#define PATTERN_H

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>

#include "spattglobals.h"
#include "spattparameters.h"
#include "alphabet.h"
#include "count.h"
#include "word.h"

using namespace std;

namespace spatt {

  struct node {
    char _word[MAX_STRING_LENGTH];
    int _size;
    int _deep;
    struct node **_next;
  
    node(char *word,int size,int deep);
    ~node();
  };

  class pattern {
    
  public:
    /* constructor from a label */
    pattern(const spattparameters &, const alphabet &,const char *);
    
    /* copy constructor */
    pattern(const pattern &source);
    
    pattern & operator=(const pattern &);

    /* destructor */
    ~pattern();
    
    /* display _label_list */
    void display();
    
    inline const long get_count() const { return _count; };
    inline void set_count(long count) { _count=count; };


    inline const char* get_label() const { return _label; };

    inline const double get_expected() const { return _expected; };
    inline void set_expected(double expected) { _expected=expected; };

    inline int get_optimized_word_list_size() const 
      { return _optimized_word_list_size; }
    inline void set_optimized_word_list_size(const int size) 
      { _optimized_word_list_size = size; }


    inline word** get_optimized_word_list() const 
      { return _optimized_word_list; }
    inline void set_optimized_word_list(word **list) 
      { if (_optimized_word_list)
	free(_optimized_word_list);
      _optimized_word_list=list; }

    inline const alphabet* get_alpha() const { return _alpha; }

    inline vector<string> & get_label_list() { return _label_list; }

    inline vector<string> & get_word_list() { return _word_list; }

  private:
    const spattparameters *_params;
    const alphabet *_alpha;

    char _label[MAX_STRING_LENGTH];
    long _count;
    double _expected;
    
    vector<string> _label_list;
    vector<string> _word_list; // clean label
    struct node *_root;
    
    int _optimized_word_list_size;
    word **_optimized_word_list;

    /* empty constructor */
    pattern() {};
    
    /* recursive parser */
    void parse(struct node* root,char *label);
    
    /* updating parsing with c and continue with label */
    void update_parse(char c,struct node *root,char *label);
    
    /* recursive tree reading with deletion function */
    void read_tree(struct node *current);
    

};

};

#endif
