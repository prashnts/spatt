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
#include "dfa.h"

using namespace std;

namespace spatt {

dfa::dfa(nfa &N,bool verbose,bool debug) {

  if (verbose)
    printf(">>> call dfa:dfa(%p)\n",&N);

  _alphabet_size=N._alphabet_size;
  //_alphabet_label=N._alphabet_label;
  //_alphabet_label=N._Alpha->label();
  _Alpha=N._Alpha;

  _renewal=false;
  _m=0;

  /* temp data used during construction */
  vector<vector<unsigned int> > state;


  /* fill _delta with empty maps */
  for (unsigned short a=0; a<_alphabet_size; a++) {
    //map<unsigned long,unsigned long> tmp;
    vector<unsigned long> tmp;
    _delta.push_back(tmp);
  }
  if (debug) {
    for (unsigned short a=0; a<_alphabet_size; a++) {
      printf("size(_delta[%i])=%i\n",a,(int)_delta[a].size());
    }
  }

  { /* subset construction */
    unsigned int n=N._nstates;
    bool *bool_state1=new bool [n];
    bool *bool_state2=new bool [n];
    bool *bool_state3=new bool [n];
    /* initialization */
    for (unsigned int i=0; i<n; i++)
      bool_state1[i]=false;
    bool_state1[N._start]=true;
    N.epsilon_closure(bool_state1,bool_state2);
    _start=0;
    _nstates=1;
    state.push_back(bool2vector(bool_state2,n));
    if (debug) {
      printf("state(0): ");
      print_vector(state[0]);
      printf("\n");
    }
    { // check if final
      bool is_final=false;
      for (vector<unsigned int>::iterator it=state[0].begin(); it!=state[0].end(); it++)
	if (*it==N._final) {
	  is_final=true;
	  break;
	}
      if (is_final) {
	_final.push_back(0);
	_is_final.push_back(true);
      } else {
	_is_final.push_back(false);
      }
      _state.push_back(0);
    } 
    /* main loop */
    unsigned long to;
    for (unsigned long s=0; s<_nstates; s++) {
      if (debug)
	printf("treating state(%i)\n",(int)s);
      // fill bool_state1 with state[s]
      for (unsigned int i=0; i<n; i++)
	bool_state1[i]=false;
      for (vector<unsigned int>::iterator it=state[s].begin(); it!=state[s].end(); it++)
	bool_state1[*it]=true;
      // loop on alphabet
      for (unsigned short a=1; a<=_alphabet_size; a++) {
	if (debug)
	  printf("alpha=%i\n",a);
	// compute delta(state[s],a)
	N.delta(bool_state1,a,bool_state2);
	// take epsilon-closure
	N.epsilon_closure(bool_state2,bool_state3);
	vector<unsigned int> dest=bool2vector(bool_state3,n);
	// check wether dest is an existing  state or not
	bool exists=false;
	unsigned int sdest;
	for (sdest=0; sdest<_nstates; sdest++) {
	  if (state[sdest]==dest) {
	    exists=true;
	    break;
	  }
	}
	if (exists) {
	  to=sdest;
	} else {
	  if (debug) {
	    printf("state(%i): ",(int)_nstates);
	    print_vector(dest);
	    printf("\n");
	  }
	  state.push_back(dest);
	  { // check if final
	    bool is_final=false;
	    
	    for (vector<unsigned int>::iterator it=dest.begin(); it!=dest.end(); it++)
	      if (*it==N._final) {
		is_final=true;
		break;
	      }
	    if (is_final) {
	      _final.push_back(_nstates);
	      _is_final.push_back(true);
	    } else {
	      _is_final.push_back(false);
	    }
	  }
	  _state.push_back(_nstates);
	  to=_nstates;
	  _nstates++;
	}
	_delta[a-1].push_back(to);
	//_delta[a-1][s]=to;
	if (debug) {
	  //printf("%i-(%c)->%i\n",s,_alphabet_label[a-1],to);
	  printf("%i-(%c)->%i\n",(int)s,_Alpha->decode(a-1),(int)to);
	}
      } // end loop on alphabet
    } // end loop on states
    if (debug) {
      printf("_final= ");
      print_lvector(_final);
      printf("\n");
    }	
    delete[] bool_state1;
    delete[] bool_state2;
    delete[] bool_state3;
  } /* end of subset construction */

  if (verbose) {
    printf("after subset construction\n");
    print();
  }

};

dfa::~dfa() {
};

vector<unsigned int> dfa::bool2vector(bool *x,
				      unsigned int size) {
  vector<unsigned int> res;
  for (unsigned int i=0; i<size; i++) {
    if (x[i])
      res.push_back(i);
  }
  return res;
};

void dfa::print_vector(vector<unsigned int> V){

  printf("[ ");
  for (vector<unsigned int>::iterator it=V.begin(); it!=V.end(); it++) 
    printf("%i ",*it);
  printf("]");

};

void dfa::print_lvector(vector<unsigned long> V){

  printf("[ ");
  for (vector<unsigned long>::iterator it=V.begin(); it!=V.end(); it++) 
    printf("%i ",(int)*it);
  printf("]");

};

void dfa::print() {
  printf("nstates = %i\t",(int)_nstates);
  printf("start = %i\t",(int)_start);
  if (_final.size()>1)
    printf("final (%i states) = ",(int)_final.size());
  else
    printf("final (1 state) = ");
  if (_final.size()<=MAX_DISPLAY)
    print_lvector(_final);
  else
    printf(" not printed");
  printf("\n");
  if (_nstates<=MAX_DISPLAY) {
    for (unsigned short a=0; a<_alphabet_size; a++) {
      //printf("delta(%c) : ",_alphabet_label[a]);
      printf("delta(%c) : ",_Alpha->decode(a));
      //      for (map<unsigned long,unsigned long>::iterator it=_delta[a].begin(); it!=_delta[a].end(); it++) 
      //	printf("(%i,%i) ",it->first,it->second);
      for (unsigned long i=0; i<_nstates; i++) 
	printf("(%i,%i) ",(int)i,(int)_delta[a][i]);
      printf("\n");
    }
  } else {
    printf("transitions not printed\n");
  }
};

void dfa::minimize(bool verbose,bool debug) {

  if (verbose)
    printf(">>> call dfa::minimize()\n");

  if (debug)
    printf("compute I:\n");

  // precompute all I[a][s]={p, delta(p,a)=s}
  vector < map<unsigned long,vector<unsigned long> > > I;
  for (unsigned short a=0; a<_alphabet_size; a++) {
    map<unsigned long,vector<unsigned long> > tmp;
    I.push_back(tmp);
    //    for (map<unsigned long,unsigned long>::iterator it=_delta[a].begin(); it!=_delta[a].end(); it++) {
    for (unsigned long i=0; i<_nstates; i++) { 
      //I[a][it->second].push_back(it->first);
      I[a][_delta[a][i]].push_back(i);
    } // end map loop
    if (debug) {
      for (map<unsigned long,vector<unsigned long> >::iterator it=I[a].begin(); it!=I[a].end(); it++) {
	printf("I[%i][%i] = ",(int)a,(int)it->first);
	print_lvector(it->second);
	printf("\n");
      }
    }
  } // end alphabet loop

  // clear the  partition
  _partition.clear();

  if (debug)
    printf("initialization:\n");

  // initialize with F and Q-F
  _partition.push_back(_final);
  _partition.push_back(difference(_state,_final));
  
  if (debug) {
    printf("partition = ");
    for (vector<vector<unsigned long> >::iterator it=_partition.begin(); it!=_partition.end(); it++)
      print_lvector(*it);
    printf("\n");
  }

  if (debug)
    printf("main loop:\n");

  // main loop
  for (unsigned long npart=0; npart<_partition.size(); npart++) {
    // select and remove the last element of working
    if (debug) {
      printf("S = ");
      print_lvector(_partition[npart]);
      printf("\n");
    }
    // loop on alphabet
    for (unsigned short a=0; a<_alphabet_size; a++) {
      if (debug)
	printf("a=%i:\n",a);
      // compute Ia(S)
      vector<unsigned long> aux;
      for (vector<unsigned long>::iterator it=_partition[npart].begin(); it!=_partition[npart].end(); it++) {
	aux=reunion(aux,I[a][*it]);
      }
      if (debug) {
	printf("Ia(S) = ");
	print_lvector(aux);
	printf("\n");
      }
      if (!aux.empty()) {
	// loop on _partition
	unsigned long nmax=_partition.size();
	for (unsigned long n=0; n<nmax; n++) {
	  if (_partition[n].size()>1) { // treat only splittable _partitions
	    vector<unsigned long> aux2=intersection(aux,_partition[n]);
	    if (debug) {
	      printf("R = ");
	      print_lvector(_partition[n]);
	      printf("\n");
	      printf("Ia(S) inter R = ");
	      print_lvector(aux2);
	      printf("\n");
	    }
	    // intersection is non empty and not equal to R
	    if ((!aux2.empty())&&(aux2.size()<_partition[n].size())) {
	      if (debug) {
		printf("split ");
		print_lvector(_partition[n]);
		printf(" in ");
		print_lvector(aux2);
		printf(" and ");
		print_lvector(difference(_partition[n],aux2));
		printf("\n");	    
	      }
	      if ((_partition[n].size()-aux2.size())>aux2.size()) {
		_partition.push_back(aux2);
		_partition[n]=difference(_partition[n],aux2);
	      } else {
		_partition.push_back(difference(_partition[n],aux2));
		_partition[n]=aux2;
	      }
	    } // end if non empty ...
	  } // end if splittable	  
	} // end loop on n
      } // end if (!aux.empty())
      if (debug) {
	printf("partition = ");
	for (vector<vector<unsigned long> >::iterator it=_partition.begin(); it!=_partition.end(); it++)
	  print_lvector(*it);
	printf("\n");	
      }
    } // end loop on alphabet
  } // end loop on npart
  
  if (verbose) {
    printf("minimal partition (%i elements): ",(int)_partition.size());
    if (_partition.size()<=MAX_DISPLAY) {
      for (vector<vector<unsigned long> >::iterator it=_partition.begin(); it!=_partition.end(); it++)
	print_lvector(*it);
      printf("\n");	
    } else {
      printf("not printed\n");
    }
  }  

};

void dfa::rebuild(bool verbose) {

  if (verbose) {
    printf(">>> call dfa::rebuild()\n");
    printf("size reduction: %i -> %i\n",(int)_nstates,(int)_partition.size());
  }

  // get new_nstates
  unsigned long new_nstates=_partition.size();

  if (new_nstates<_nstates) { // work to do
    _nstates=new_nstates;

    //clear _state and _is_final
    _state.clear();
    _is_final.clear();

    // first build a map to get the partition number of each state
    map<unsigned long,unsigned long> get_part;
    for (unsigned long i=0; i<_partition.size(); i++) {
      _state.push_back(i);
      _is_final.push_back(false);
      for (vector<unsigned long>::iterator it=_partition[i].begin(); it!=_partition[i].end(); it++)
	get_part[*it]=i;
    } // end i loop
    
    // update _start
    _start=get_part[_start];

    // update _final (using a map to avoid doubles)
    {
      map<unsigned long,bool> tmp;
      for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++)
	tmp[get_part[*it]]=true;
      _final.clear();
      for (map<unsigned long,bool>::iterator it=tmp.begin(); it!=tmp.end(); it++) {
	_final.push_back(it->first);
	_is_final[it->first]=true;
      }
    }

    // update _delta
    for (unsigned short a=0; a<_alphabet_size; a++) {
      //map<unsigned long,unsigned long> tmp=_delta[a];
      vector<unsigned long> tmp=_delta[a];
      _delta[a].clear();
      for (unsigned long i=0; i<_partition.size(); i++) {
	//_delta[a][i]=get_part[tmp[_partition[i][0]]];
	_delta[a].push_back(get_part[tmp[_partition[i][0]]]);
      }
    }
    
  } // end if work to do

  if (verbose) {
    printf("after rebuilding\n");
    print();
  }

}

vector<unsigned long>
dfa::intersection(vector<unsigned long> &v1,
		  vector<unsigned long> &v2) {

  vector<unsigned long> res;

  if ((v1.empty())||(v2.empty()))
    return res;

  //printf("\n***intersection ");
  //print_lvector(v1);
  //printf(" and ");
  //print_lvector(v2);
  //printf("\n");

  vector<unsigned long>::iterator it1=v1.begin();
  vector<unsigned long>::iterator it2=v2.begin();
  while ((it1!=v1.end())&&(it2!=v2.end())) {
    if (*it1==*it2) {
      res.push_back(*it1);
      it1++;
      it2++;
    } else if (*it1<*it2) {
      it1++;
    } else {
      it2++;
    }
  }

  //printf("***res = ");
  //print_lvector(res);
  //printf("\n");

  return res;
};

vector<unsigned long>
    dfa::difference(vector<unsigned long> &v1,
		    vector<unsigned long> &v2) {

  if (v2.empty())
    return v1;
  
  //printf("\n***difference ");
  //print_lvector(v1);
  //printf(" - ");
  //print_lvector(v2);
  //printf("\n");


  vector<unsigned long> res;
  vector<unsigned long>::iterator it1=v1.begin();
  vector<unsigned long>::iterator it2=v2.begin();
  while ((it1!=v1.end())&&(it2!=v2.end())) {
    if (*it1==*it2) {
      it1++;
      it2++;
    } else if (*it1<*it2) {
      res.push_back(*it1);
      it1++;
    } else {
      it2++;
    }
  }
  if (it2==v2.end()) {
    while (it1!=v1.end()) {
      res.push_back(*it1);
      it1++;
    }
  }

  //printf("***res = ");
  //print_lvector(res);
  //printf("\n");
    
  return res;
};

vector<unsigned long>
dfa::reunion(vector<unsigned long> &v1,
	   vector<unsigned long> &v2) {

  if (v1.empty())
    return v2;
  if (v2.empty())
    return v1;

  //printf("\n***reunion ");
  //print_lvector(v1);
  //printf(" et ");
  //print_lvector(v2);
  //printf("\n");


  vector<unsigned long> res;
  vector<unsigned long>::iterator it1=v1.begin();
  vector<unsigned long>::iterator it2=v2.begin();
  while ((it1!=v1.end())&&(it2!=v2.end())) {
    if (*it1==*it2) {
      res.push_back(*it1);
      it1++;
      it2++;
    } else if (*it1<*it2) {
      res.push_back(*it1);
      it1++;
    } else {
      res.push_back(*it2);
      it2++;
    }
  }
  if (it1==v1.end()) {
    while (it2!=v2.end()) {
      res.push_back(*it2);
      it2++;
    }
  } else {
    while (it1!=v1.end()) {
      res.push_back(*it1);
      it1++;
    }
  }

  //printf("***res = ");
  //print_lvector(res);
  //printf("\n");

  return res;
};

void dfa::dot(string file) {

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
      fprintf(out,"%i [peripheries=1];\n",(int)_start);
      fprintf(out,"start0 -> %i [label=\"\"];\n",(int)_start);
      // _final
      for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++)
	fprintf(out,"%i [style=filled,peripheries=2];\n",(int)*it);
      // regular states
      vector<unsigned long> regular=difference(_state,_final);
      for (vector<unsigned long>::iterator it=regular.begin(); it!=regular.end(); it++) {
	if (*it!=_start)
	  fprintf(out,"%i [peripheries=1];\n",(int)*it);
      }
      // transitions
      // gather all alpha-trans building labels
      map<unsigned long,map<unsigned long,vector<unsigned short> > > gather;
      for (unsigned short a=0; a<_alphabet_size; a++) {
	//	for (map<unsigned long,unsigned long>::iterator it=_delta[a].begin(); it!=_delta[a].end(); it++) {
	for (unsigned long i=0; i<_nstates; i++) { 
	  //gather[it->first][it->second].push_back(a);
	  gather[i][_delta[a][i]].push_back(a);
	}
      }
      // write transitions
      for (map<unsigned long,map<unsigned long,vector<unsigned short> > >::iterator it1=gather.begin(); it1!=gather.end(); it1++) {
	for (map<unsigned long,vector<unsigned short> >::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
	  char label[256];
	  unsigned short pos=0;
	  {
	    //if (it2->second.size()==_alphabet_size) {
	    //label[pos++]='.';
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
	  fprintf(out,"%i -> %i [label=\"%s\"];\n",(int)it1->first,(int)it2->first,label);
	}
      }
    } else {
      fprintf(out,"0 [style=filled,shape=circle,color=olivedrab1,peripheries=1];\n");
      fprintf(out,"1 [style=filled,shape=circle,color=tomato1,peripheries=2];\n");
      fprintf(out,"start0 -> 0 [label=\"\"];\n");
      fprintf(out,"0 -> 1 [label=\"'too large(%i); no_display_beyond=%i'\"];\n",(int)_nstates,(int)MAX_DOT);
    }
    // footer
    fprintf(out,"}\n");
  }

  fclose(out);
};

vector<unsigned long> dfa::locate_occ(string &file,bool verbose) {

  if (verbose)
    printf(">>> call dfa::locate_occ(\"%s\")\n",file.c_str());

  vector<unsigned long> res;
  FILE *in=NULL;
  in=fopen(file.c_str(),"r");
  if (in==NULL) {
    fprintf(stderr,"cannot read file \"%s\". Aborting dot export.\n",file.c_str());
  } else {
//    // build a local coding table
//    char code[256];
//    for (unsigned short i=0; i<256; i++)
//      code[i]=-1;
//    for (unsigned short a=0; a<_alphabet_size; a++)
//      code[(unsigned short)_alphabet_label[a]]=a;
//    code[(unsigned short)'>']=-2;
//    code[(unsigned short)'#']=-2;
//    code[(unsigned short)'%']=-2;
    // build a local map for testing if a state is final
    // map<unsigned long,bool> is_final;
    // for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++) {
    //   //printf("%i\n",*it);
    //   is_final[*it]=true;
    // }
    unsigned long token=_start;
    // buffer read the file
    char buffer[200];
    unsigned long nvalid=0;
    while (fgets(buffer,200,in)) {
      //if (code[(unsigned short)buffer[0]]!=-2) { // not a comment line
      if (_Alpha->code(buffer[0])!=COMMENT) { // not a comment line
	//printf("buffer=%s",buffer);
	for (unsigned int i=0; i<strlen(buffer); i++) {
	  //int c=code[(unsigned short)buffer[i]];
	  int c=_Alpha->code(buffer[i]);
	  //printf("%i->%i\n",i,c);
	  //if (c>=0) {
	  if (c>UNDEF) {
	    nvalid++;
	    token=_delta[c][token];
	    //if (is_final[token])
	    if (_is_final[token])
	      res.push_back(nvalid);
	  }
	}
      } // end if not comment line
    } // end while

    if (verbose) {
      printf("sequence length = %i\n",(int)nvalid);
      printf("number of occurrences = %i\n",(int)res.size());
      printf("at positions : ");
      if (res.size()<=MAX_DISPLAY)
	print_lvector(res);
      else
	printf("not printed");
      printf("\n");
    }

  }
  return res;
};


void dfa::initialize_ambiguity_vectors(bool verbose,bool debug){

  if (verbose)
    printf(">>> call dfa::initialize_ambiguity_vectors()\n");

  // initialize with empty maps
  _m=0;
  _inv_delta.clear();
  _inv_Delta.clear();
  for (unsigned long q=0; q<_nstates; q++) {
    map<string,bool> delta; 
    map<unsigned long,bool> Delta;
    _inv_delta.push_back(delta);
    _inv_Delta.push_back(Delta);
  }

  // loop on transitions filling maps
  string key(1,(char)0);
  unsigned long q;
  for (unsigned long p=0; p<_nstates; p++) {   // loop on states
    for (unsigned short a=0; a<_alphabet_size; a++) { // loop on alphabet
      key[0]=(char)a;
      q=_delta[a][p];
      _inv_delta[q][key]=true;
      _inv_Delta[q][p]=true;
    } // end loop on alphabet
  } // end loop on states

  if (debug)
    print_ambiguity_vectors();

};

// todo todo todo
void dfa::print_ambiguity_vectors() {

  {   // display _inv_delta
    // loop on states
    for (unsigned long q=0; q<_nstates; q++) {
      printf("inv_delta[%i]={ ",(int)q);
      for (map<string,bool>::iterator it=_inv_delta[q].begin(); it!=_inv_delta[q].end(); it++){
	print_string(it->first);
	printf(" ");
      }
      printf("}\n");
    } // end q loop
  } // end display _inv_delta
  {   // display _inv_Delta
    // loop on states
    for (unsigned int q=0; q<_nstates; q++) {
      printf("inv_Delta[%i]={ ",(int)q);
      for (map<unsigned long,bool>::iterator it=_inv_Delta[q].begin(); it!=_inv_Delta[q].end(); it++)
	printf("%i ",(int)it->first);
      printf("}\n");
    } // end q loop
  } // end display _inv_Delta
  
};

string dfa::inv_delta(unsigned short m,unsigned long state) {
  if (m>_m) {
    fprintf(stderr,"Only m<=%i available at this stage in inv_delta. Aborting.\n",(int)_m);
    exit(EXIT_FAILURE);
  }
  // get the string in _inv_delta[state].begin();
  string res;
  if (!_inv_delta[state].empty()) {
    string aux;
    aux=(_inv_delta[state].begin()->first);
    if (aux.size()>=m) {
      res=aux;
      res=res.substr(res.size()-m,m);
    }
  }
  return res;
};


void dfa::remove_local_ambiguity(unsigned short m,unsigned long state,bool verbose,bool debug){

  if (verbose)
    printf(">>> call dfa::remove_local_ambiguity(m=%i,state=%i)\n",(int)m,(int)state);

  // check compatibility with _m
  if (m!=_m+1) {
    fprintf(stderr,"Only m=%i available at this stage in remove_local_ambiguity. Aborting.\n",(int)_m+1);
    exit(EXIT_FAILURE);
  }

  // update _inv_delta[state] if m>1
  if (m>1) {
    map<string,bool> delta;
    // loop on states in _inv_Delta
    for (map<unsigned long,bool>::iterator it=_inv_Delta[state].begin(); it!=_inv_Delta[state].end(); it++) {
      unsigned long p=it->first;
      string key=inv_delta(m-1,p)+inv_delta(1,state);
      delta[key]=true;      
      if (debug) {
	printf("inv_delta(%i,%i)=",(int)m-1,(int)p);
	print_string(inv_delta(m-1,p));
	printf(" + ");
	printf("inv_delta(1,%i)=",(int)p);
	print_string(inv_delta(1,p));
	printf("\n");
      }
    } //end it loop
    // update _inv_delta[state] with delta
    _inv_delta[state]=delta;
  } // end if (m>1)
  
  if (debug) {
    printf("_inv_delta[%i] ={ ",(int)state);
    for (map<string,bool>::iterator it=_inv_delta[state].begin(); it!=_inv_delta[state].end(); it++) {
      print_string(it->first);
      printf(" ");
    }
    printf("}\n");
  }

  // loop on _inv_delta[state]
  if (debug) {
    printf("we have %i ambiguous transitions\n",(int)_inv_delta[state].size());
  }
  while (_inv_delta[state].size()>1) {
    // get a
    string a=_inv_delta[state].begin()->first;
    if (a.size()==m) {
      if (debug) {
	printf("treating transition: ");
	print_string(a);
	printf("\n");
      }
      // add a new state
      unsigned long new_state=_nstates;
      _nstates++;
      _state.push_back(new_state);
      if (debug) {
	printf("adding state %i\n",(int)new_state);
      }
      if (_is_final[state]) {
	if (debug)
	  printf("%i is a new final state\n",(int)new_state);
	_final.push_back(new_state);
	_is_final.push_back(true);
      } else {
	_is_final.push_back(false);
      }
      { // add a term to _inv_delta
	if (debug) {
	  printf("initializing _inv_delta[%i]={ ",(int)new_state);
	  print_string(a);
	  printf(" }\n");
	}
	map<string,bool> delta;
	delta[a]=true;
	_inv_delta.push_back(delta);
      }
      { // add a term to _inv_Delta      
	if (debug)
	  printf("initializing _inv_Delta[%i]={ }\n",(int)new_state);
	map<unsigned long,bool> Delta;
	_inv_Delta.push_back(Delta);
      }
      // set delta(new_state,b)=delta(state,b) for all b
      for (unsigned int b=0; b<_alphabet_size; b++) {
	if (debug)
	  printf("set delta(%i,%i)=%i\n",(int)new_state,(int)b,(int)_delta[b][state]);
	_delta[b].push_back(_delta[b][state]);
	// add new_state to Delta^-1(delta(new_state,b))
	if (debug)
	  printf("adding %i to _inv_Delta[%i]\n",(int)new_state,(int)_delta[b][state]);
	_inv_Delta[_delta[b][new_state]][new_state]=true;
      }
      // loop on Delta^-1(state)
      unsigned short last_of_a=a[a.size()-1];
      for (map<unsigned long,bool>::iterator it=_inv_Delta[state].begin(); it!=_inv_Delta[state].end(); it++) {
	if (_delta[last_of_a][it->first]==state) {
	  if ( (inv_delta(m-1,it->first)+string(1,(char)last_of_a) )==a) {
	    if (debug) printf("change _delta[%i][%i] from %i to %i\n",(int)it->first,(int)last_of_a,(int)_delta[last_of_a][it->first],(int)new_state);
	    _delta[last_of_a][it->first]=new_state;
	    // add the corresponding states to _inv_Delta[new_state]
	    if (debug)
	      printf("add %i to _inv_Delta[%i]\n",(int)it->first,(int)new_state);
	    _inv_Delta[new_state][it->first]=true;
	  } // end if == a
	} // end if == state
      } // end loop on it
      // check Delta^-1(state)
      {
	if (debug)
	  printf("checking _inv_Delta[%i]: ",(int)state);
	// loop on Delta^-1(state)
	vector<map<unsigned long,bool>::iterator> toerase;
	for (map<unsigned long,bool>::iterator it=_inv_Delta[state].begin(); it!=_inv_Delta[state].end(); it++) {
	  bool test=false;
	  // loop on letters
	  for (unsigned short b=0; b<_alphabet_size; b++) {
	    //if (debug) printf("_delta[%i][%i]=%i == %i ?\n",it->first,b,_delta[it->first][b],state);
	    if (_delta[b][it->first]==state) {
	      test=true;
	      break;
	    } // end if
	  } // end loop on letters
	  if (!test) {
	    // remove it from Delta^-1(state)
	    if (debug) printf("remove(%i) ",(int)it->first);
	    toerase.push_back(it);
	  }
	} // end it loop
	// delete for true !
	for (vector<map<unsigned long,bool>::iterator>::iterator it=toerase.begin(); it!=toerase.end(); it++) {
	  _inv_Delta[state].erase(*it);
	}
	if (debug)
	  printf("\n");
      } // en Delta^(-1) check
      
      // remove a from _inv_delta[state]
      if (debug) {
	printf("remove ");
	print_string(a);
	printf(" from _inv_delta[%i]\n",(int)state);
      }
    }
    _inv_delta[state].erase(a);
  } // end while
  
};


void dfa::print_string(string w) {
  for (unsigned short i=0; i<w.size(); i++)
    //printf("%c",_alphabet_label[w[i]]);
    printf("%c",_Alpha->decode(w[i]));
};

void dfa::remove_ambiguity(unsigned short m,bool verbose,bool debug){

  if (verbose)
    printf(">>> call remove_ambiguity(m=%i)\n",m);

  for (unsigned short j=_m+1; j<=m; j++) {
    unsigned long imax=_nstates;
    for (unsigned long i=0; i<imax; i++) {
      remove_local_ambiguity(j,i,debug,debug);
    }
    _m=j;
  }

  if (verbose)
    print();

};

void dfa::renewal(bool verbose){
  if (verbose)
    printf(">>> call dfa::renewal()\n");

  // loop on final states
  for (vector<unsigned long>::iterator it=_final.begin(); it!=_final.end(); it++) {
    // set all delta(f,a)=delta(start,a)
    for (unsigned short a=0; a<_alphabet_size; a++) {
      _delta[a][*it]=_delta[a][_start];
    }
  } // end loop on final states
  _renewal=true;

  if (verbose)
    print();
};

vector<unsigned long> dfa::starts(bool verbose) {

  if (verbose)
    printf(">>> call dfa::starts()\n");

  string path(_m,(char)0);
  vector<unsigned long> res;
  recursive_starts(path,res,_start,_m,0);

  if (verbose) {
    printf("starts = ");
    print_lvector(res);
    printf("\n");
  }

  return res;
};

void dfa::recursive_starts(std::string &path, std::vector<unsigned long> &res, unsigned long pos, unsigned short m, unsigned short i) {
  if (i<m) {
    for (unsigned short a=0; a<_alphabet_size; a++) {
      path[i]=(char)a;
      recursive_starts(path,res,_delta[a][pos],m,i+1);
    }
  }
  if (i==m) {
    res.push_back(pos);
  }
};

unsigned long dfa::code(std::string &word) {
  unsigned long res=0;
  for (unsigned short i=0; i<word.size(); i++) {
    res=_alphabet_size*res+(short)word[i];
  }
  return res;
}

};
