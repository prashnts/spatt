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
#include "sequence.h"

using namespace std;

namespace spatt {

  sequence::~sequence() {
    if (_in)
      fclose(_in);
  };

sequence::sequence(alphabet &A,vector<string> &sequence_file) {

  _Alpha=&A;
  _sequence_file=sequence_file;
  reset();

};

sequence::sequence(alphabet &A,string &sequence_file) {

  _Alpha=&A;
  _sequence_file.push_back(sequence_file);
  reset();
  _in=NULL;
};

void sequence::reset() {

  _current=_sequence_file.begin();
  // fixme: there is still a missing free related to _in
  //if (_in)
  //  fclose(_in);
  _in=NULL;
  _nvalidchar=0;
  _lastnvalidchar=0;
  _nvalidseq=0;
  _seqend=false;
  _seqinterrupt=false;
  //printf("ok before next_line\n");
  next_line();
  //printf("ok after next_line\n");
  if (_in==NULL) {
    fprintf(stderr,"sequence::reset() : warning, no valid sequence file. Aborting.\n");
    exit(EXIT_FAILURE);
  }
};

bool sequence::next_file() {
  if (_in)
    fclose(_in);
  _in=NULL;
  while (_current!=_sequence_file.end()) {
    // fixme: next_file appear to open twice a single file (uncomment the line below to see the problem)
    //printf("next_file:: opening \"%s\"\n",_current->c_str()); 
    _in=fopen(_current->c_str(),"r");
    if (_in==NULL) {
      fprintf(stderr,"sequence::next_file() : warning, cannot read file \"%s\" (skipping it)\n", _current->c_str());
      _current++;
      return next_file();
    } else {
      _current++;
      break;
    }
  }
  if (!_in)
    return false;
  else
    return true;
};      

bool sequence::next_line() {
  //printf("next_line:: _in=%p\n",_in);
  if (!_in) {
    if (next_file())
      return next_line();
    else 
      return false;
  } else {
    if (!fgets(_buffer,BUFFER_SIZE,_in)) {
      if (next_file())
	return next_line();
      else
	return false;
    } else {
      _buffer_size=strlen(_buffer);
      _pos=0;
      return true;
    }
  }
};

int sequence::next() {
  //printf("(next:_pos=%i,_buffer_size=%i,_seqend=%i)",_pos,_buffer_size,(int)_seqend);
  if (_seqend) {
    _nvalidseq++;
    _lastnvalidchar=0;
    _seqend=false;
    if (_seqinterrupt) {
      _seqinterrupt=false;
      return next();
    } else {
      if (next_line())
	return next();
      else
	return SEQEMPTY;
    }
  } else {
    int code;
    if (_pos<_buffer_size) {
      // simple case
      code=_Alpha->code(_buffer[_pos++]);
      if (code==COMMENT) {
	if (next_line())
	  return next();
	else {
	  _seqend=true;
	  return SEQEND;
	}
      } else if (code==SEQUENCE) {
	if (_nvalidchar>0) {
	  _seqend=true;
	  return SEQEND;
	} else {
	  if (next_line())
	    return next();
	  else
	    return SEQEMPTY;
	}
      } else {
	if (code>UNDEF) {
	  _nvalidchar++;
	  _lastnvalidchar++;
	  return code;
	} else if (code!=IGNORE) {
	  fprintf(stderr,"sequence::next(): warning, unexpected sequence interruption with character '%c'\n",_buffer[_pos-1]);
	  _seqend=true;
	  _seqinterrupt=true;
	  return SEQEND;
	} else {
	  return next();
	}
      } 
    } else {
      if (next_line())
	return next();
      else {
	_seqend=true;
	return SEQEND;
      } // end if (next_line())
    } // end if (_pos<_buffer_size)
  } // end if (_seqend)
};
};
