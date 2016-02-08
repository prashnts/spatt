/* $Id: ldstat.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "transition.h"

namespace spatt {

  transition::transition(const spattparameters &params,const alphabet &alpha,const sequence &seq,const count &occ,const markov &M,const input &I): _params(&params), _alpha(&alpha), _seq(&seq), _occ(&occ), _M(&M), _In(&I)
  {
    /* set default values */
    _trans_table=NULL;
    _allocated_size=0;
    _trans_size=0;
    _positions=NULL;
    _probas=NULL;
    _counting_pos_size=0;
    _counting_pos=NULL;
    /* construct the common states */
    // length of the words to consider
    _h=_M->get_order(); 
    if (_h<1)
      _h=1;
    string current(_h,_alpha->code2char(0));
    // loop on all words
    _k=_alpha->get_size();
    long imax=(long)pow((float)_k,_h);
    _n_common_states=imax;
    for (int i=0; i<imax; i++) {
      //cout<<current<<endl;
      _common_states[current]=NULL;	
      int j=1;
      while (j<=_h) {
	int cc;
	cc=_alpha->char2code(current[_h-j]);
	cc++;
	if (cc==_k) {
	  cc=0;
	  current[_h-j]=_alpha->code2char(cc);
	  j=j+1;
	} else {
	  current[_h-j]=_alpha->code2char(cc);
	  j=_h+1;
	}
      } // end while
    } // end for
    /* allocate the corresponding trans_table */
    _common_trans_table=(trans *)malloc(sizeof(trans)*imax*_k);
    if (_common_trans_table==NULL) {
      fprintf(stderr,"not enough memory to construct transition\n");
      exit(EXIT_FAILURE);
    }
    /* tag and link to the states */
    {
      trans *pos=_common_trans_table;
      long number=0;
      map<string,trans*> :: iterator it;
      for (it=_common_states.begin(); it!=_common_states.end(); it++) {
	//cout<<"'"<<(*it).first<<"'"<<endl;
	(*it).second=pos;
	for (int i=0; i<_k; i++) {
	  (*pos).number=number;	  
	}
	pos+=_k;
	number++;
      }
    } // tag and link done
    /* fill and connect trans */
    {
      trans *pos;
      string concat(_h,' ');
      map<string,trans*> :: iterator it;
      for (it=_common_states.begin(); it!=_common_states.end(); it++) {
	pos=(*it).second;
	for (int i=0; i<(_h-1); i++)
	  concat[i]=((*it).first)[i+1];
	//cout<<"*********"<<endl;
	for (int i=0; i<_k; i++) {
	  concat[_h-1]=_alpha->code2char(i);
	  //cout<<(*it).first<<"->"<<concat<<endl;
	  pos[i].to=_common_states[concat];
	  pos[i].proba=_M->get_trans(pos->number,i);
	  //printf("proba(%i,%i)=%f\n",pos->number,pos[i].to->number,pos[i].proba);
	}
      }
    } // fill and connect done

  };

  void transition::init() {
    _states.clear();
    _n_counting_trans=0;
    _counted=0;
  }

  void transition::build(word *w){
    init();
    string s=string(w->get_label());
    add(s);
    tag();
    connect();
    fix_trans();
    counting(s);
  };
  
  void transition::build(word *w1,word *w2){
    init();
    string s1=string(w1->get_label());
    string s2=string(w2->get_label());
    add(s1);
    add(s2);
    tag();
    connect();
    fix_trans();
    counting(s1);
    counting(s2);
  };
  
  void transition::build(pattern *patt){
    init();
    word *w=NULL;
    for (int j=0; j<patt->get_optimized_word_list_size(); j++) {
      w=patt->get_optimized_word_list()[j];
      add(string(w->get_label()));
    }
    tag();
    connect();
    fix_trans();
    for (int j=0; j<patt->get_optimized_word_list_size(); j++) {
      w=patt->get_optimized_word_list()[j];
      counting(string(w->get_label()));
    }
    //print_trans();
  };

  void transition::add(string s) {
    _n_counting_trans++;
    int size=s.size();
    //printf("size=%i, _h=%i\n",size,_h);
    if (size<=_h) {
      fprintf(stderr,"warning ! %s too short will not be added\n",s.c_str());
    } else {
      //cout<<s<<":";
      for (int i=_h+1; i<size; i++) {
	string current(s,0,i);
	//cout<<current<<",";
	_states[current]=NULL;
      }
      //cout<<endl;
    }
  }

  void transition::tag() {
    /* check memory requirements */
    _n_states=_states.size();
    if ((_n_states*_k)>_allocated_size) {
      /* free last alloc */
      if (_trans_table!=NULL)
	free(_trans_table);
      /* and allocate new one */
      _allocated_size=_n_states*_k;
      _trans_table=(trans *)malloc(sizeof(trans)*_allocated_size);
      if (_trans_table==NULL) {
	fprintf(stderr,"not enough memory to construct transition\n");
	exit(EXIT_FAILURE);
      }
    }
    /* tag all */
    trans *pos=_trans_table;
    long number=_n_common_states;
    map<string,trans*> :: iterator it;
    for (it=_states.begin(); it!=_states.end(); it++) {
      //cout<<"'"<<(*it).first<<"'"<<endl;
      (*it).second=pos;
      for (int i=0; i<_k; i++) {
	(*pos).number=number;
	(*pos).proba=0.0;
      }
      pos+=_k;
      number++;
    }    
  }

  void transition::connect() {
    //printf("starting connect\n");
    /* find longuest word size in _states */
    map<string,trans*> :: iterator it;
    trans *pos;
    /* common part first */
    {
    
      string concat(_h+1,' ');
      for (it=_common_states.begin(); it!=_common_states.end(); it++) {
	pos=(*it).second;
	// copy starting word
	for (int i=0; i<_h; i++)
	  concat[i]=((*it).first)[i];
	//cout<<"*********"<<endl;
	for (int i=0; i<_k; i++) {
	  concat[_h]=_alpha->code2char(i);
	  //cout<<(*it).first<<"->"<<concat<<endl;
	  // find the longuest suffix
	  string suffix=longuest_suffix(concat);
	  //cout<<"longuest suffix("<<concat<<")="<<suffix<<endl;
	  if (suffix.size()==_h)
	    pos[i].to=_common_states[suffix];
	  else
	    pos[i].to=_states[suffix];	    
	  pos[i].proba=_M->get_trans(pos->number,i);
	  //printf("proba(%i,%i)=%f\n",pos->number,pos[i].to->number,pos[i].proba);
	}
      }
    } // end common part
    /* specific part then */
    {
      for (it=_states.begin(); it!=_states.end(); it++) {
	pos=(*it).second;
	int word_size=((*it).first).size();
	string concat(word_size+1,' ');      
	// copy starting word
	for (int i=0; i<word_size; i++)
	  concat[i]=((*it).first)[i];
	//cout<<"*********"<<endl;
	for (int i=0; i<_k; i++) {
	  concat[word_size]=_alpha->code2char(i);
	  //cout<<(*it).first<<"->"<<concat<<endl;
	  // find the longuest suffix
	  string suffix=longuest_suffix(concat);
	  //cout<<"longuest suffix("<<concat<<")="<<suffix<<endl;
	  if (suffix.size()==_h)
	    pos[i].to=_common_states[suffix];
	  else
	    pos[i].to=_states[suffix];
	  // get code
	  long code=0;
	  {
	    long factor=1;
	    for (int ii=1; ii<=_h; ii++) {
	      code+=factor*_alpha->char2code( ( (*it).first )[word_size-ii] );
	      factor*=_k;
	    }
	  }
	  //printf("code=%i\n",code);
	  pos[i].proba=_M->get_trans(code,i);
	  //printf("proba(%i,%i)=%f\n",pos->number,pos[i].to->number,pos[i].proba);
	}
      }
    } // end specific part
  }

  string transition::longuest_suffix(string s) {
    int size;
    map<string,trans*> :: iterator it;
    size=s.size();
    for (int i=0; i<size; i++) {
      string current(s,i,size);
      //cout<<"testing "<<current<<endl;
      it=_states.find(current);
      if (it!=_states.end()) {
	//printf("state number %i\n",((*it).second)[0].number);
	return current;
      }
      it=_common_states.find(current);
      if (it!=_common_states.end()) {
	//printf("state number %i\n",((*it).second)[0].number);
	return current;
      }
    }
    return s;
  }

  void transition::print_fixed_trans() {

    if (_positions==NULL || _probas==NULL) {
      fprintf(stderr,"warning ! fix_trans() must be called before print_trans()\n");
    } else {
      printf("state");
      for (int i=0; i<_k; i++)
	printf("\t%c",_alpha->code2char(i));
      printf("\n");
      long *lpos=_positions;
      double *dpos=_probas;
      for (long i=0; i<_order; i++) {
	printf("%i",i);
	for (int j=0; j<_k; j++) 
	  printf("\t%i(%.2f)",lpos[j],dpos[j]);
	printf("\n");
	lpos+=_k;
	dpos+=_k;
      }
    }
  }

  void transition::print_trans() {
    
    trans *pos;
    map<string,trans*> :: iterator it;
    printf("state");
    for (int i=0; i<_k; i++)
      printf("\t%c",_alpha->code2char(i));
    printf("\n");
    for (it=_common_states.begin(); it!=_common_states.end(); it++) {
      pos=(*it).second;
      cout<<(*pos).number<<"("<<(*it).first<<")";
      for (int i=0; i<_k; i++)
	printf("\t%i(%.2f)",pos[i].to->number,pos[i].proba);
      printf("\n");
    }
    for (it=_states.begin(); it!=_states.end(); it++) {
      pos=(*it).second;
      cout<<(*pos).number<<"("<<(*it).first<<")";
      for (int i=0; i<_k; i++)
	printf("\t%i(%.2f)",pos[i].to->number,pos[i].proba);
      printf("\n");
    }
  }

  void transition::fix_trans() {
    /* check memory requirements */
    _order=_n_common_states+_n_states;
    if (_order>_trans_size) {
      /* need more memory */
      if (_positions!=NULL)
	free(_positions);    
      if (_probas!=NULL)
	free(_probas);
      _trans_size=_order;
      _positions=(long *)malloc(sizeof(long)*_trans_size*_k);
      if (_positions==NULL) {
	fprintf(stderr,"not enough memory to construct transition\n");
	exit(EXIT_FAILURE);
      }
      _probas=(double *)malloc(sizeof(double)*_trans_size*_k);
      if (_probas==NULL) {
	fprintf(stderr,"not enough memory to construct transition\n");
	exit(EXIT_FAILURE);
      }      
    }
    /* read map and tables */
    long *lpos=_positions;
    double *dpos=_probas;
    trans *pos;
    map<string,trans*> :: iterator it;
    for (it=_common_states.begin(); it!=_common_states.end(); it++) {
      pos=(*it).second;
      for (int i=0; i<_k; i++) {
	lpos[i]=pos[i].to->number;
	dpos[i]=pos[i].proba;
      }
      lpos+=_k;
      dpos+=_k;
    }
    for (it=_states.begin(); it!=_states.end(); it++) {
      pos=(*it).second;
      for (int i=0; i<_k; i++) {
	lpos[i]=pos[i].to->number;
	dpos[i]=pos[i].proba;
      }
      lpos+=_k;
      dpos+=_k;
    }
    /* _counting_pos stuff */
    {
      /* check memory requirements */
      if (_n_counting_trans>_counting_pos_size) {
	/* need more memory */
	if (_counting_pos!=NULL)
	  free(_counting_pos);
	_counting_pos_size=_n_counting_trans;
	_counting_pos=(long *)malloc(sizeof(long)*_counting_pos_size);
	if (_counting_pos==NULL) {
	  fprintf(stderr,"not enough memory to construct transition\n");
	  exit(EXIT_FAILURE);
	}
      }
    } // end _counting_pos stuff
  }

  void transition::counting(string s) {
    //printf("_counted=%i\t_counting_pos[0]=%i\n",_counted,_counting_pos[0]);
    int size=s.size();
    if (size<=_h) {
      fprintf(stderr,"warning ! %s too short will not be counted\n",s.c_str());      
    } else {
      string prefix(s,0,size-1);
      map<string,trans*> :: iterator it;      
      if ( (size-1)==_h ) {
	/* looking in common part */
	it=_common_states.find(prefix);
	if (it==_common_states.end()) {
	  fprintf(stderr,"warning ! %s must be added before counting call\n",s.c_str());
	} else {
	  trans *pos=(*it).second;
	  _counting_pos[_counted]=pos->number*_k+_alpha->char2code(s[size-1]);
	  _counted++;	  
	}	
      } else {	
	/* looking in specific part */
	it=_states.find(prefix);
	if (it==_states.end()) {
	  fprintf(stderr,"warning ! %s must be added before counting call\n",s.c_str());
	} else {
	  trans *pos=(*it).second;
	  //cout<<prefix<<"->"<<s[size-1]<<endl;
	  _counting_pos[_counted]=pos->number*_k+_alpha->char2code(s[size-1]);
	  _counted++;	  
	}	
      }
    }
  }

  void transition::print_counting() {
    //printf("print_counting\t_n_counting_trans=%i\t_counting_pos_size=%i\n",_n_counting_trans,_counting_pos_size);
    printf("counting:\t");
    for (long i=0; i<_n_counting_trans; i++) {
      //printf("_counting_pos[%i]=%i\n",i,_counting_pos[i]);
      printf("%i(%.2f)\t",_counting_pos[i],_probas[_counting_pos[i]]);
    }
    printf("\n");
  }
  
  // this function multiply all counting proba by t
  void transition::multiply_Q_by(double t) {
    for (long i=0; i<_n_counting_trans; i++) {
      _probas[_counting_pos[i]]*=t;
    }
  }

  // compte res = transition *  x
  void transition::transition_by_vector(double *x,double *res){
    double *dpos=_probas;
    long *lpos=_positions;
    /* loop on lines */
    for (long i=0; i<_order; i++) {
      res[i]=0.0;
      for (int j=0; j<_k; j++) {
	res[i]+=dpos[j]*x[lpos[j]];
      }
      dpos+=_k;
      lpos+=_k;
    }    
  }

  void transition::fill_fortran(double *res){
    // init with zeros
    for (int i=0; i<_order*_order; i++)
      res[i]=0;
    // loop on lines
    {
      double *dpos=_probas;
      long *lpos=_positions;
      for (int i=0; i<_order; i++) {
	for (int j=0; j<_k; j++) {	  	  
	  //res[i][lpos[j]]=dpos[j];
	  res[i+_order*lpos[j]]=dpos[j];
	}
	lpos+=_k;
	dpos+=_k;
      }      
    }
    /* verification */
    //printf("A:\n");
    //for (int i=0; i<_order; i++) {
    //  printf("\t");
    //  for (int j=0; j<_order; j++)
    //    printf("%.2f\t",_AA[i+_order*j]);
    //  printf("\n");
    //}
  }

  // res = P * u
  void transition::P_times(double *u,double *res) {
    // fixme: not optimal here. need another transition structure.
    // res = (P+Q) * u
    transition_by_vector(u,res);
    // res = res - Q * u
    long j=0;
    for (long i=0; i<_n_counting_trans; i++) {
      j=_counting_pos[i]/_k; 
      res[j]-=_probas[_counting_pos[i]]*u[_positions[_counting_pos[i]]];
    }
    //printf("P_times\n");
    //printf("u = ");
    //for (int i=0; i<_order; i++)
    //  printf("%f\t",u[i]);
    //printf("\n");
    //printf("res = ");
    //for (int i=0; i<_order; i++)
    //  printf("%f\t",res[i]);
    //printf("\n");    
  }; 

  // res = log( P * exp(u) )
  void transition::log_P_times(double *u,double *res) {
    // find largest term of u
    double z=mylog(0.0);
    for (int i=0; i<_order; i++) {
      if (u[i]>z)
	z=u[i];
    }
    if (z==mylog(0.0))
      z=0.0;
    // compute res = (P+Q) * exp(u-z)
    {
      double *dpos=_probas;
      long *lpos=_positions;
      /* loop on lines */
      for (long i=0; i<_order; i++) {
	res[i]=0.0;
	for (int j=0; j<_k; j++) {
	  res[i]+=dpos[j]*myexp(u[lpos[j]]-z);
	}
	dpos+=_k;
	lpos+=_k;
      }          
    }
    // res = res - Q * exp(u-z)
    long j=0;
    for (long i=0; i<_n_counting_trans; i++) {
      j=_counting_pos[i]/_k; 
      res[j]-=_probas[_counting_pos[i]]*myexp(u[_positions[_counting_pos[i]]]-z);
    }
    // res = log(res) + z
    for (int i=0; i<_order; i++) {
      res[i]=mylog(res[i])+z;
    }
    //printf("log_P_times: u[0]=%e\tres[0]=%e\n",u[0],res[0]);
  }; 

  // res = P * u + Q * v
  void transition::P_times_plus_Q_times(double *u,double *v,double *res){
    // fixme: not optimal here. need another transition structure.
    P_times(u,res);
    // res = res + Q * v
    long j=0;
    for (long i=0; i<_n_counting_trans; i++) {
      j=_counting_pos[i]/_k; 
      res[j]+=_probas[_counting_pos[i]]*v[_positions[_counting_pos[i]]];
    }
    //printf("P_times_plus_Q_times\n");
    //printf("u = ");
    //for (int i=0; i<_order; i++)
    //  printf("%f\t",u[i]);
    //printf("\n");
    //printf("v = ");
    //for (int i=0; i<_order; i++)
    //  printf("%f\t",v[i]);
    //printf("\n");
    //printf("res = ");
    //for (int i=0; i<_order; i++)
    //  printf("%f\t",res[i]);
    //printf("\n");
  };

  // res = log( P * exp(u) + Q * exp(v) )
  void transition::log_P_times_plus_Q_times(double *u,double *v,double *res){
    //printf("log_P_times_plus_Q_times: start\n");
    // find largest term of u and v
    double z=mylog(0.0);
    for (int i=0; i<_order; i++) {
      if (u[i]>z)
	z=u[i];
      if (v[i]>z)
	z=v[i];
    }
    if (z==mylog(0.0))
      z=0.0;
    //printf("log_P_times_plus_Q_times: z=%e\n",z);
    // compute res = (P+Q) * exp(u-z)
    {
      double *dpos=_probas;
      long *lpos=_positions;
      /* loop on lines */
      for (long i=0; i<_order; i++) {
	res[i]=0.0;
	for (int j=0; j<_k; j++) {
	  res[i]+=dpos[j]*myexp(u[lpos[j]]-z);
	}
	dpos+=_k;
	lpos+=_k;
      }          
    }
    //printf("res[0]=%e\n",res[0]);
    // res = res - Q * exp(u-z) + Q * exp(v-z)
    long j=0;
    for (long i=0; i<_n_counting_trans; i++) {
      j=_counting_pos[i]/_k; 
      res[j]-=_probas[_counting_pos[i]]*(myexp(u[_positions[_counting_pos[i]]]-z)-myexp(v[_positions[_counting_pos[i]]]-z));
    }
    //printf("res[0]=%e\n",res[0]);
    // res = log(res) + z
    for (int i=0; i<_order; i++) {
      res[i]=mylog(res[i])+z;
    }
    //printf("res[0]=%e\n",res[0]);
    //printf("log_P_times_plus_Q_times: u[0]=%e\tv[0]=%e\tres[0]=%e\n",u[0],v[0],res[0]);
    //printf("log_P_times_plus_Q_times: end\n");
  };    


  // res = sum of Q by row (sum of row1, sum of row2, ...)
  void transition::sum_Q(double *res){
    // initialization
    for (long i=0; i<_order; i++)
      res[i]=0.0;
    // parsing Q
    long j=0;
    for (long i=0; i<_n_counting_trans; i++) {
      j=_counting_pos[i]/_k; 
      res[j]+=_probas[_counting_pos[i]];
    }
  }


  transition::~transition() {
    if (_common_trans_table!=NULL)
      free(_common_trans_table);
    if (_trans_table!=NULL)
      free(_trans_table);    
    if (_positions!=NULL)
      free(_positions);    
    if (_probas!=NULL)
      free(_probas);
    if (_counting_pos!=NULL)
      free(_counting_pos);
  };
  
};
