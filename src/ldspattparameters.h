/* $Id: ldspattparameters.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
#ifndef _LDSPATTPARAMETERS_H
#define _LDSPATTPARAMETERS_H

#include "spattparameters.h"

namespace spatt {
  class ldspattparameters : public spattparameters {

  public:
    ldspattparameters();

  inline const bool get_precise() const { return _precise; }

  // fixme: _precise really need to be public ?
  bool _precise;
  
  };

};

#endif
