/* $Id: spattglobals.h 675 2006-02-22 13:33:34Z gnuel $ */
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
/*  header defining global values                                    */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#ifndef SPATTGLOBALS_H
#define SPATTGLOBALS_H

#include <config.h>

namespace spatt {
  const int BUFFER_SIZE          = 200;
  const int DEFAULT_COUNT_LENGTH =   8;
  const int MAX_LENGTH           =  30;
  const int MAX_STRING_LENGTH    = 2000;
  const int COMMENT_CHAR         = '#';

};

#endif
