/* $Id: sspattparameters.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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
#ifndef _SSPATTPARAMETERS_H
#define _SSPATTPARAMETERS_H

#include "spattparameters.h"

namespace spatt {
  class sspattparameters : public spattparameters {
  public:
    sspattparameters();
    virtual ~sspattparameters() {};

    inline double get_min_pvalue() const {
      return _min_pvalue;
    }
	
    inline int get_nobs() const {
      return _nobs;
    }


  protected:
    int _nobs;
    double _min_pvalue;

    friend class sspattparameterparser;

  private:
    sspattparameters(const sspattparameters &) {};
    sspattparameters &operator=(const sspattparameters &) { return *this;};
  };


};

#endif
