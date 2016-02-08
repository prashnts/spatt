/* $Id: sequence.cc 676 2006-02-22 14:25:10Z gnuel $ */
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
/*  Class sequence                                                   */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "sequence.h"

using namespace std;

namespace spatt {
/* opens the stream and computes length */
  sequence::sequence(const spattparameters &params,
		     const alphabet &alpha) :
  _params(&params), _alpha(&alpha)
{

  /* gets filename */
  sprintf(_filename,_params->get_sequence_filename()); 
  /* opens stream */
  _stream=fopen(_filename,"r");
  if (_stream==NULL) {
    fprintf(stderr,"Cannot open sequence file \"%s\"\n",_filename);
    exit(EXIT_FAILURE);
  }
  /* computes length */
  fseek(_stream,0,SEEK_END);
  _length=ftell(_stream);
  rewind(_stream);
  if (_params->get_debug_level()>0)
    printf("length(\"%s\")=%li\n",_filename,_length);

  /* initialize position to BUFFER_DONE */
  _position=BUFFER_DONE;

  /* initialize valid_char */
  _valid_char=0;

  /* initialize _nseq */
  _nseq=0;
};

/* returns next code */
int sequence::next(){

  int c;
  char *test=NULL;

  /* gets new buffer if necessary */
  if (_position==BUFFER_DONE) {
    //printf("fgets call\n");
    test=fgets(_buffer,BUFFER_SIZE,_stream);
    //printf("test=%p\n",test);
    if (test && _params->get_debug_level()>=4)
      printf("%s",_buffer);
    _position=0;
  } else {
    test=_buffer;
    _position++;
  }
  /* tests buffer still validity */
  if (test==NULL) {
    //printf("valid_char=%i\n",_valid_char);
    return EOF_CODE;
  }

  c=_alpha->char2code(_buffer[_position]);
  //if (c==INVALID_CODE) 
  //  printf("\'%c\'=code(%i) is invalid\n",_buffer[_position],(unsigned char)_buffer[_position]);
  //printf("buffer[%i]=\'%c\'(%i)\n",_position,_buffer[_position],c);

  /* end buffer if comment or sequence */
  if (c==COMMENT_CODE || c==SEQUENCE_CODE)
    _position=BUFFER_DONE;

  /* update sequence count */
  if (c==SEQUENCE_CODE) 
    _nseq++;

  /* test end buffer */
  if (_buffer[_position]=='\0')
    _position=BUFFER_DONE;

  /* if ignored, get next one */
  if (c==IGNORE_CODE)
    return next();

  /* count valid char */
  if (c>=0)
    _valid_char++;
  
  /* return code */
  return c;
};

/* restart sequence */
void sequence::restart(){
  rewind(_stream);
  _position=BUFFER_DONE;
  _valid_char=0;
};

/* close stream */
sequence::~sequence(){
  fclose(_stream);
};

};

