/* $Id: spattparameterparser.h 547 2005-11-22 10:30:51Z gnuel $ */
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

#ifndef SPATTPARAMETERPARSER_H
#define SPATTPARAMETERPARSER_H

#include <string>
#include <vector>

/**
 * Ugly fix (part I) to avoid macro redefinition warning messages.
 **/ 
#undef PACKAGE
#undef VERSION
#include "../tclap-1.0.5/config/config.h"
/**
 * Ugly fix (part II) to avoid macro redefinition warning messages.
 **/ 
#undef PACKAGE
#undef VERSION
#include <tclap/CmdLine.h>


#include "spattglobals.h"
#include "spattparameters.h"

using namespace std;
using namespace TCLAP;

namespace spatt {

class spattparameterparser  {

public:

  virtual ~spattparameterparser();


  virtual void parse(int, char **, spattparameters *);


protected:
  spattparameterparser(const string &, const string &);

  CmdLine *_cmdline;
  SwitchArg *_quietArg;
  SwitchArg *_verboseArg;
  SwitchArg *_debugArg;
  ValueArg<int> *_debuglevelArg;
  ValueArg<string> *_alphabetArg;
  UnlabeledValueArg<string> *_sequenceArg;
  ValueArg<string> *_modelArg;
  ValueArg<string> *_stationaryArg;
  ValueArg<string> *_markovArg;
  ValueArg<string> *_countArg;
  ValueArg<string> *_patternFileArg;
  ValueArg<int> *_lengthArg;
  MultiArg<string> *_patternListArg;
  SwitchArg *_strandsArg;
  ValueArg<int> *_orderArg;
  SwitchArg *_allwordsArg;
  SwitchArg *_normalizeArg;
  ValueArg<double> *_maxpvalueArg;
  ValueArg<int> *_nobsArg;
};

}
#endif
