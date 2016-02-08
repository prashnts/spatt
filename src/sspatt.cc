/* $Id: sspatt.cc 864 2006-09-13 13:21:10Z gnuel $ */
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
/*  sspatt main program                                              */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "sspattparameters.h"
#include "sspattparameterparser.h"

#include "alphabet.h"
#include "sequence.h"
#include "count.h"
#include "word.h"
#include "input.h"
#include "pattern.h"
#include "markov.h"
#include "process.h"
#include "stat.h"
#include "sstat.h"

using namespace std;
using namespace spatt;

int main (int argc, char **argv) {
  
  sspattparameters params;
  sspattparameterparser parser;
  parser.parse(argc,argv,&params);
  alphabet alpha(params);
  sequence seq(params,alpha);
  spatt::count occ(params,alpha,seq);
  input I(params,alpha,seq,occ);

  markov *M=NULL;
  if (params.use_markov_file()==true) {
    M = new markov(params,alpha);
  } else {
    M = new markov(params,alpha,occ);
  }
  // as _occ is needed even in the markov_file case
  M->set_occ(occ);
  M->set_seq(seq);

  /* check if alphabet compatible with both strand option */
  if (alpha.has_complementary()==0 && params.get_both_strands()==true) {
    fprintf(stderr,"Warning ! --both-strands option specified but not complementary alphabet supplied. Option will be ignored.\n");
    params.set_both_strands(false);
  }

  stat *S=NULL;
  if ((M->get_model()!=NULL || M->get_order()<=0) &&  M->get_stationary()!=NULL) {
    S= new sstat(params,alpha,seq,occ,*M,I);
  } else {
    S = new stat(params,alpha,seq,occ,*M,I);
  }
  process P(params,alpha,seq,occ,*M,I,*S);

  P.count_all_dfa();
  P.select_case();

  delete M;
  delete S;

}
