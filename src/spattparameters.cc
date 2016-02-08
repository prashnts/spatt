/* $Id: spattparameters.cc 708 2006-03-09 08:14:04Z gnuel $ */
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
/*  Class parameter allows reads the command line and returns the    */
/*  the parameters.                                                  */
/*  								     */
/*  This class uses argtable2                 			     */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "spattparameters.h"
#include "spattglobals.h"

#include <cstring>
#include <iostream>

using namespace std;

namespace spatt {

  spattparameters::spattparameters() :
  _debug_level(0),
  _ignore_case(true),
  _alphabet_label("acg[ut]:tgca"),
  _sequence_filename(""),
  _count_filename(""),
  _stationary_filename(""),
  _markov_filename(""),
  _model_filename(""),
  _length(DEFAULT_COUNT_LENGTH),
  _use_sequence_file(false),
  _use_count_file(false),
  _use_markov_file(false),
  _use_pattern_file(false),
  _use_model_file(false),
  _use_stationary_file(false),
  _both_strands(false),
  _markov_order(-2),
  _all_words(false),
  _nobs(-1),
  _max_pvalue(1.0)
{
}

};
