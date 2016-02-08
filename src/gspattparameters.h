/* $Id: gspattparameters.h 455 2005-09-01 14:27:20Z gnuel $ */
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
#ifndef _GSPATTPARAMETERS_H
#define _GSPATTPARAMETERS_H

#include "spattparameters.h"

namespace spatt {
  class gspattparameters : public spattparameters {
  public:
    gspattparameters();

    inline const bool get_output_zscore() const { return _output_zscore; }
    
    // fixme: _output_zscore really need to be public ?
    bool _output_zscore;

  };

};

#endif
