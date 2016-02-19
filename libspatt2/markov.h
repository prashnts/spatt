/**
    Copyright 2006 Mark Hoebeke, Vincent Miele & Gregory Nuel.

    This file is part of SPatt 2

    SPatt 2 is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SPatt 2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SPatt 2; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**/
#ifndef MARKOV_H
#define MARKOV_H 1

#define MAX_ITER 1000
#define BUFFER_SIZE 300
#define COMMENT_CHAR '#'
#define SEPARATORS "\n\t ,;"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>

namespace spatt {

class markov {

 public:

  /* statistics */
  unsigned short _alphabet_size;
  unsigned short _m; // order >=0

  /* boolean */
  bool _stationary; /* true if Markov chains are stationary */

  /* data */
  std::vector<double> _param;
  std::vector<double> _mu0; /* stationary distribution */

  markov(unsigned short alphabet_size,unsigned short m,std::string markov_file,bool stationary=false,bool verbose=false);
  markov(unsigned short alphabet_size,unsigned short m,double* transitions,bool stationary,bool verbose);
  markov(unsigned short alphabet_size,unsigned short m,bool stationary=false,bool verbose=false);

  // case uniform m=0
  markov(unsigned short alphabet_size,bool stationary=false,bool verbose=false);

  ~markov();

  void print();

  void normalize();

  void compute_mu0(bool verbose=false);

};

};
#endif
