/* $Id: gspattparameterparser.cc 455 2005-09-01 14:27:20Z gnuel $ */
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
#include "gspattparameterparser.h"

namespace spatt {

  gspattparameterparser::gspattparameterparser() : 
    spattparameterparser("gspatt",VERSION)

  {
    _output_zscoreArg=new SwitchArg("","zscore","output zscore instead of classical logscale statistics.",false);
    _cmdline->add(_output_zscoreArg);
  }

  void 
  gspattparameterparser::parse(int argc, char **argv,spattparameters *params)
 {

    gspattparameters *sparams=dynamic_cast<gspattparameters *>(params);      

    this->spattparameterparser::parse(argc,argv,sparams);

    if (_output_zscoreArg->getValue() == true)
      sparams->_output_zscore=true;

  }

};


