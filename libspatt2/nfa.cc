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
#include "nfa.h"

using namespace std;

namespace spatt {

//nfa::nfa(string &alphabet_label,
nfa::nfa(alphabet &A,
	 string &pattern_label,
	 bool verbose,
	 bool debug) {
  
  if (verbose) {
    fprintf(stdout,">>> call nfa::nfa(alphabet_label=\"%s\",pattern_label=\"%s\")\n",A.label().c_str(),pattern_label.c_str()); 
  }

  _pattern_label=pattern_label;
  //_alphabet_label=alphabet_label;
  _Alpha=&A;
  _data=NULL;

  _alphabet_size=_Alpha->size();
  
//  { /* alphabet part */
//    _alphabet_size=alphabet_label.size();
//    /* filling code table */
//    for (unsigned short i=0; i<256; i++)
//      _alphabet_code[i]=UNDEF;
//    for (unsigned short i=0; i<_alphabet_size; i++)
//      _alphabet_code[(unsigned short)_alphabet_label[i]]=i+1;
//    /* separators */
//    _alphabet_code[(unsigned short)'-']=SEPARATOR;
//    _alphabet_code[(unsigned short)',']=SEPARATOR;
//    _alphabet_code[(unsigned short)':']=SEPARATOR;
//    _alphabet_code[(unsigned short)';']=SEPARATOR;
//    _alphabet_code[(unsigned short)'|']=SEPARATOR;
//    /* end block */
//    _alphabet_code[(unsigned short)')']=END_BLOCK;
//    _alphabet_code[(unsigned short)'}']=END_BLOCK;
//    _alphabet_code[(unsigned short)']']=END_BLOCK;
//    /* repeat */
//    _alphabet_code[(unsigned short)'(']=REPEAT;
//    /* choice */
//    _alphabet_code[(unsigned short)'[']=CHOICE;
//    /* nchoice */
//    _alphabet_code[(unsigned short)'{']=NCHOICE;
//    /* any char */
//    _alphabet_code[(unsigned short)'_']=ANY_CHAR;
//    _alphabet_code[(unsigned short)'.']=ANY_CHAR;
//  } /* end of alphabet part */

  /* allocating _ntrans */
  _ntrans= new unsigned int [_alphabet_size+1];
  
  { /* first parsing: filling statistics */
    /* we start with A* */
    _ntrans[0]=0;
    for (unsigned short i=1; i<=_alphabet_size; i++)
      _ntrans[i]=1;
    /* parsing */
    const char *buffer=_pattern_label.c_str();
    unsigned int buffer_size=_pattern_label.size();
    if (debug)
      printf("first parsing of pattern:\n");
    bool valid=false;
    bool repeat=false;
    unsigned int *which = new unsigned int[_alphabet_size];
    int repeat_m=-1;
    int repeat_n=-1;
    for (unsigned int pos=0; pos<buffer_size; pos++) {
      char c=buffer[pos];
      //int code=_alphabet_code[(unsigned short)c];
      int code=_Alpha->code(c);
      //if (debug)
      //	printf("'%c'>%i\n",buffer[pos],code);
      switch (code) {
      case ANY_CHAR:
	if (debug)
	  printf("any char\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=1;
	valid=true;
	break;
      case CHOICE:
	if (debug)
	  printf("choice\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=0;
	// read until end block
	{
	  unsigned int endpos=pos;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    //if (local_code<=0) {
	    if (local_code<=UNDEF) {	    
	      pos=endpos;
	      break;
	    } else {
	      //which[local_code-1]=1;
	      which[local_code]=1;
	    }
	  } 
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"syntax error: invalid choice block\nAborting !\n");
	    exit(EXIT_FAILURE);
	  }
	}
	valid=true;
	break;
      case NCHOICE:
	if (debug)
	  printf("nchoice\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=1;
	// read until end block
	{
	  unsigned int endpos=pos;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    //if (local_code<=0) {
	    if (local_code<=UNDEF) {
	      pos=endpos;
	      break;
	    } else {
	      //which[local_code-1]=0;
	      which[local_code]=0;
	    }
	  } 
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"syntax error: invalid nchoice block\nAborting !\n");
	    exit(EXIT_FAILURE);
	  }
	}
	valid=true;
	break;
      case REPEAT:
	if (!valid) {
	  fprintf(stderr,"invalid syntax: unary repeat operator '()' without left-hand term\nAborting !\n");
	  exit(EXIT_FAILURE);
	}
	if (debug)
	  printf("repeat\n");
	{
	  unsigned int endpos=pos;
	  unsigned int seppos=0;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    if (local_code==SEPARATOR) {
	      seppos=endpos;
	    }
	    if (local_code==END_BLOCK) {
	      break;
	    }
	  } // end endpos loop
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"invalid syntax: repeat block without end\nAborting !\n");
	     exit(EXIT_FAILURE);
	  }
	  if (seppos==0) {
	    repeat_m=atoi(&buffer[pos+1]);
	    repeat_n=repeat_m;
	    if (debug) {
	      printf("m=%i\n",repeat_m);
	    }
	    if (repeat_m<0) {
	      fprintf(stderr,"invalid repeat argument\nAborting !\n");
	      exit(EXIT_FAILURE);
	    }
	  } else {
	    repeat_m=atoi(&buffer[pos+1]);
	    repeat_n=atoi(&buffer[seppos+1]);
	    if (debug) {
	      printf("m=%i n=%i\n",repeat_m,repeat_n);
	    }
	    if ((repeat_m<0)||(repeat_n<repeat_m)||(repeat_n==0)) {
	      fprintf(stderr,"invalid repeat argument\nAborting !\n");
	      exit(EXIT_FAILURE);
	    }
	  }
	  pos=endpos;
	  repeat=true;
	}
	valid=false;
	break;
      case SEPARATOR:
	if (debug)
	  printf("pattern separator\n");
	valid=false;
	break;
      case UNDEF:
	if (debug)
	  printf("undef char ignored\n");
	 fprintf(stderr,"invalid syntax: '%c' is not allowed\nAborting !\n",c);
	 exit(EXIT_FAILURE);
	break;
      default:
	if (code>UNDEF) {
	  if (debug)
	    printf("simple letter\n");
	  for (unsigned short i=0; i<_alphabet_size; i++)
	    which[i]=0;
	  //which[code-1]=1;
	  which[code]=1;
	  valid=true;
	  break;
	}
      }
      if (valid) {
	for (unsigned short i=0; i<_alphabet_size; i++)
	  _ntrans[i+1]+=which[i];
	if (debug) {
	  for (unsigned short i=0; i<_alphabet_size; i++)
	    //printf("%c=%i ",_alphabet_label[i],which[i]);
	    printf("%c=%i ",_Alpha->decode(i),which[i]);
	  printf("\n");
	}
      } else {
	if (repeat) {
	  // its a repeat
	  for (unsigned short i=0; i<_alphabet_size; i++)
	    _ntrans[i+1]+=(repeat_n-1)*which[i];
	  _ntrans[0]+=repeat_n-repeat_m;
	  repeat=false;
	} else {
	  _ntrans[0]++;
	  // new pattern
	}
      }
    } // end loop on pos
    delete[] which;
    if (debug) {
      printf("statistics:\n");
      printf("\t%i epsilon-trans\n",_ntrans[0]);
      for (unsigned short i=0; i<_alphabet_size; i++)
	//printf("\t%i %c-trans\n",_ntrans[i+1],_alphabet_label[i]);
	printf("\t%i %c-trans\n",_ntrans[i+1],_Alpha->decode(i));
    }
  } /* end first parsing */

  { /* allocating memory */
    _total_ntrans=0;
    for (unsigned short i=0; i<=_alphabet_size; i++)
      _total_ntrans+=_ntrans[i];
    _data=new trans [_total_ntrans];
    _trans=new trans* [_alphabet_size+1];
    trans *pos=_data;
    for (unsigned short i=0; i<=_alphabet_size; i++) {
      _trans[i]=pos;
      pos+=_ntrans[i];
    }
  } /* end allocating memory */

  { /* second parsing: filling data */
    /* we start with A* */
    _start=0;
    _final=0;
    _nstates=1;
    unsigned int last=0;
    unsigned int *current_trans = new unsigned int [_alphabet_size+1];
    current_trans[0]=0;
    for (unsigned short i=1; i<=_alphabet_size; i++) {
      _trans[i][0].from=0;
      _trans[i][0].to=0;      
      current_trans[i]=1;
    }
    /* parsing */
    const char *buffer=_pattern_label.c_str();
    unsigned int buffer_size=_pattern_label.size();
    if (debug)
      printf("second parsing of pattern:\n");
    bool valid=false;
    bool repeat=false;
    bool first_pattern=true;
    unsigned int *which = new unsigned int[_alphabet_size];
    int repeat_m=-1;
    int repeat_n=-1;
    for (unsigned int pos=0; pos<buffer_size; pos++) {
      char c=buffer[pos];
      //int code=_alphabet_code[(unsigned short)c];
      int code=_Alpha->code(c);
      //if (debug)
      //	printf("'%c'>%i\n",buffer[pos],code);
      switch (code) {
      case ANY_CHAR:
	if (debug)
	  printf("any char\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=1;
	valid=true;
	break;
      case CHOICE:
	if (debug)
	  printf("choice\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=0;
	// read until end block
	{
	  unsigned int endpos=pos;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    //if (local_code<=0) {
	    if (local_code<=UNDEF) {
	      pos=endpos;
	      break;
	    } else {
	      //which[local_code-1]=1;
	      which[local_code]=1;
	    }
	  } 
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"syntax error: invalid choice block\nAborting !\n");
	    exit(EXIT_FAILURE);
	  }
	}
	valid=true;
	break;
      case NCHOICE:
	if (debug)
	  printf("nchoice\n");
	for (unsigned short i=0; i<_alphabet_size; i++)
	  which[i]=1;
	// read until end block
	{
	  unsigned int endpos=pos;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    //if (local_code<=0) {
	    if (local_code<=UNDEF) {
	      pos=endpos;
	      break;
	    } else {
	      //which[local_code-1]=0;
	      which[local_code]=0;
	    }
	  } 
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"syntax error: invalid nchoice block\nAborting !\n");
	    exit(EXIT_FAILURE);
	  }
	}
	valid=true;
	break;
      case REPEAT:
	if (!valid) {
	  fprintf(stderr,"invalid syntax: unary repeat operator '()' without left-hand term\nAborting !\n");
	  exit(EXIT_FAILURE);
	}
	if (debug)
	  printf("repeat\n");
	{
	  unsigned int endpos=pos;
	  unsigned int seppos=0;
	  int local_code;
	  for (endpos=pos+1; endpos<buffer_size; endpos++) {
	    //local_code=_alphabet_code[(unsigned short)buffer[endpos]];
	    local_code=_Alpha->code(buffer[endpos]);
	    if (local_code==SEPARATOR) {
	      seppos=endpos;
	    }
	    if (local_code==END_BLOCK) {
	      break;
	    }
	  } // end endpos loop
	  if (local_code!=END_BLOCK) {
	    fprintf(stderr,"invalid syntax: repeat block without end\nAborting !\n");
	     exit(EXIT_FAILURE);
	  }
	  if (seppos==0) {
	    repeat_m=atoi(&buffer[pos+1]);
	    repeat_n=repeat_m;
	    if (debug) {
	      printf("m=%i\n",repeat_m);
	    }
	    if (repeat_m<0) {
	      fprintf(stderr,"invalid repeat argument\nAborting !\n");
	      exit(EXIT_FAILURE);
	    }
	  } else {
	    repeat_m=atoi(&buffer[pos+1]);
	    repeat_n=atoi(&buffer[seppos+1]);
	    if (debug) {
	      printf("m=%i n=%i\n",repeat_m,repeat_n);
	    }
	    if ((repeat_m<0)||(repeat_n<repeat_m)||(repeat_n==0)) {
	      fprintf(stderr,"invalid repeat argument\nAborting !\n");
	      exit(EXIT_FAILURE);
	    }
	  }
	  pos=endpos;
	  repeat=true;
	}
	valid=false;
	break;
      case SEPARATOR:
	if (debug)
	  printf("pattern separator\n");
	valid=false;
	if (!first_pattern) {
	  _trans[0][current_trans[0]].from=last;
	  _trans[0][current_trans[0]].to=_final;
	  current_trans[0]++;
	}
	last=_start;
	first_pattern=false;
	break;
      case UNDEF:
	if (debug)
	  printf("undef char ignored\n");
	 fprintf(stderr,"invalid syntax: '%c' is not allowed\nAborting !\n",c);
	 exit(EXIT_FAILURE);
	break;
      default:
	if (code>UNDEF) {
	  if (debug)
	    printf("simple letter\n");
	  for (unsigned short i=0; i<_alphabet_size; i++)
	    which[i]=0;
	  //which[code-1]=1;
	  which[code]=1;
	  valid=true;
	  break;
	}
      }
      if (valid) {
	for (unsigned short i=0; i<_alphabet_size; i++) {
	  if (which[i]==1) {
	    _trans[i+1][current_trans[i+1]].from=last;
	    _trans[i+1][current_trans[i+1]].to=_nstates;
	    current_trans[i+1]++;
	  }
	}
	last=_nstates;
	if (first_pattern)
	  _final=last;
	_nstates++;
      } else {
	if (repeat) {
	  for (unsigned int state=_nstates; state<_nstates+repeat_n-1; state++) {
	    for (unsigned short i=0; i<_alphabet_size; i++) {
	      if (which[i]==1) {
		//printf("*** letter %i, trans %i, from=%i, to=%i\n",i+1,current_trans[i+1],last,state);
		_trans[i+1][current_trans[i+1]].from=last;
		_trans[i+1][current_trans[i+1]].to=state;
		current_trans[i+1]++;
	      }
	    }
	    last=state;
	  } // end state loop
	  // adding the epsilon-trans
	  for (unsigned int state=_nstates+repeat_m-2; state<_nstates+repeat_n-2; state++) {
	    //	    printf(">>>>>>>state=%i last=%i current_trans[0]=%i\n",state,last,current_trans[0]);
	    _trans[0][current_trans[0]].from=state;
	    _trans[0][current_trans[0]].to=last;
	     current_trans[0]++;
	  }
	  _nstates=_nstates+repeat_n-1;
	  //printf(">>>>>>_nstates=%i\n",_nstates);
	  repeat=false;
	} else {
	  // nothing
	}
      }
    } // end loop on pos
    if (!first_pattern) {
      _trans[0][current_trans[0]].from=last;
      _trans[0][current_trans[0]].to=_final;
      current_trans[0]++;
    }
    delete[] which;
    delete[] current_trans;
    if (verbose) {
      printf("transitions (%i states, start=%i, final=%i):\n",_nstates,_start,_final);
      printf("\t%i epsilon-trans: ",_ntrans[0]);
      for (unsigned int j=0; j<_ntrans[0]; j++)
	printf("(%i,%i) ",_trans[0][j].from,_trans[0][j].to);
      printf("\n");
      for (unsigned short i=0; i<_alphabet_size; i++) {
	//printf("\t%i %c-trans: ",_ntrans[i+1],_alphabet_label[i]);
	printf("\t%i %c-trans: ",_ntrans[i+1],_Alpha->decode(i));
	for (unsigned int j=0; j<_ntrans[i+1]; j++)
	  printf("(%i,%i) ",_trans[i+1][j].from,_trans[i+1][j].to);
	printf("\n");
      }
    }
  } /* end second parsing */

};

nfa::~nfa() {
  delete[] _ntrans;
  delete[] _trans;
  delete[] _data; 
}

bool nfa::delta(bool *from,
		unsigned int a,
		bool *to) {

  bool changed=false;

  if ((a<0)||(a>_alphabet_size)) {
    fprintf(stderr,"a out of range in nfa::delta\nAborting !");
    exit(EXIT_FAILURE);
  }
  
  if (a==0) {
    for (unsigned int i=0; i<_nstates; i++)
      to[i]=from[i];
  } else {
    for (unsigned int i=0; i<_nstates; i++)
      to[i]=false;
  }

  for (unsigned int i=0; i<_ntrans[a]; i++) {
    //printf("trans[%i][%i].from=%i,trans[%i][%i].to=%i\n",a,i,_trans[a][i].from,a,i,_trans[a][i].to);
    if (from[_trans[a][i].from]==true) {
      if (to[_trans[a][i].to]==false) {
	changed=true;
	to[_trans[a][i].to]=true;
      }
    }
  } // end loop on i

  return !changed;
}; 

bool nfa::epsilon_closure(bool *from,
			  bool *to) {
  bool done=delta(from,0,to);
  if (!done) {
    //printf("not done yet\n");
    epsilon_closure(to,from);
  }
  return done;
};

void nfa::dot(string &file) {

  FILE *out=NULL;
  out=fopen(file.c_str(),"w");
  if (out==NULL) {
    fprintf(stderr,"cannot write file \"%s\". Aborting dot export.\n",file.c_str());
  } else {
    // header 
    fprintf(out,"/** dot %s -Tpdf > fsa.pdf && acroread fsa.pdf */\n",file.c_str());
    fprintf(out,"digraph \"dfa generated by SPatt 2.0\" {\n");
    fprintf(out,"nodesep=0.3;\n");
    fprintf(out,"rankdir=LR;\n");
    fprintf(out,"center=true;\n");
    fprintf(out,"start0 [shape=plaintext,label=\"\"];\n");
    if (_nstates<=MAX_DOT) {
      // start
      fprintf(out,"%i [peripheries=1];\n",_start);
      fprintf(out,"start0 -> %i [label=\"\"];\n",_start);
      // _final
      fprintf(out,"%i [style=filled,peripheries=2];\n",_final);
      // regular states
      for (unsigned int i=0; i<_nstates; i++) {
	if ((i!=_start)&&(i!=_final))
	  fprintf(out,"%i [peripheries=1];\n",i);
      }
      // epsilon-transitions
      for (unsigned int i=0; i<_ntrans[0]; i++) {
	fprintf(out,"%i -> %i [label=\"\",style=dashed];\n",_trans[0][i].from,_trans[0][i].to);
      }
      // alpha-trans
      // gather all alpha-trans building labels
      map<unsigned int,map<unsigned int,vector<unsigned short> > > gather;
      for (unsigned short a=1; a<=_alphabet_size; a++) {
	for (unsigned int i=0; i<_ntrans[a]; i++) {
	  gather[_trans[a][i].from][_trans[a][i].to].push_back(a-1);
	}
      }
      // write transitions
      for (map<unsigned int,map<unsigned int,vector<unsigned short> > >::iterator it1=gather.begin(); it1!=gather.end(); it1++) {
	for (map<unsigned int,vector<unsigned short> >::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
	  char label[256];
	  unsigned short pos=0;
	  {
	  //if (it2->second.size()==_alphabet_size) {
	  //  label[pos++]='.';
	  //} else {
	    //label[pos++]=_alphabet_label[it2->second[0]];
	    label[pos++]=_Alpha->decode(it2->second[0]);
	    for (unsigned short i=1; i<it2->second.size(); i++) {
	      label[pos++]=',';
	      //label[pos++]=_alphabet_label[it2->second[i]];
	      label[pos++]=_Alpha->decode(it2->second[i]);
	    }
	  }
	  label[pos]='\0';
	  fprintf(out,"%i -> %i [label=\"%s\"];\n",it1->first,it2->first,label);
	}
      }
    } else {
      fprintf(out,"0 [style=filled,shape=circle,color=olivedrab1,peripheries=1];\n");
      fprintf(out,"1 [style=filled,shape=circle,color=tomato1,peripheries=2];\n");
      fprintf(out,"start0 -> 0 [label=\"\"];\n");
      fprintf(out,"0 -> 1 [label=\"'too large(%i); no_display_beyond=%i'\"];\n",_nstates,MAX_DOT);
    }
    // footer
    fprintf(out,"}\n");
  }

  fclose(out);
};

};
