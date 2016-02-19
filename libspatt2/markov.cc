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
#include "markov.h"

using namespace std;

namespace spatt {

  markov::markov(unsigned short alphabet_size,unsigned short m,string markov_file,bool stationary,bool verbose) {

    if (verbose)
      printf(">>> call markov::markov(%i,%i,\"%s\")\n",alphabet_size,m,markov_file.c_str());

    _alphabet_size=alphabet_size;
    _stationary=stationary;
    _m=m;
    unsigned long n=1;
    for (unsigned short i=0; i<=_m; i++)
      n*=_alphabet_size;

    FILE *in=fopen(markov_file.c_str(),"r");
    if (in==NULL) {
      fprintf(stderr,"markov::markov : cannot read markov file \"%s\". Aborting.\n",markov_file.c_str());
      exit(EXIT_FAILURE);
    }

    char buffer[BUFFER_SIZE];
    char *pos;
    unsigned long nread=0;
    while (fgets(buffer,BUFFER_SIZE,in)) {
      if (buffer[0]!=COMMENT_CHAR) {
	pos=strtok(buffer,SEPARATORS);
	if (pos) {
	  _param.push_back(atof(pos));
	  nread++;
	}
	while (pos=strtok(NULL,SEPARATORS)) {
	  _param.push_back(atof(pos));
	  nread++;
	  if (nread>n) {
	    fprintf(stderr,"markov::markov : Wrong number of parameter in the Markov file. Aborting.\n");
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }

    if (nread!=n) {
      fprintf(stderr,"markov::markov : Wrong number of parameter in the Markov file. Aborting.\n");
      exit(EXIT_FAILURE);
    }

    // buffer reading and split (fgets, strtok and atof)
    fclose(in);

    normalize();

    if (_stationary)
      compute_mu0();

    if (verbose)
      print();

  };

  markov::markov(unsigned short alphabet_size,unsigned short m,double* transitions,bool stationary,bool verbose) {
    if (verbose)
      printf(">>> call markov::markov(%i,%i)\n",alphabet_size,m);

    _alphabet_size=alphabet_size;
    _stationary=stationary;
    _m=m;
    unsigned long n=1;
    for (unsigned short i=0; i<=_m; i++)
      n*=_alphabet_size;

    for (unsigned long i=0; i<n; i++)
      _param.push_back(transitions[i]);

    normalize();

    if (_stationary)
      compute_mu0();

    if (verbose)
      print();
  };

  markov::markov(unsigned short alphabet_size,unsigned short m,bool stationary,bool verbose) {
    if (verbose)
      printf(">>> call markov::markov(%i,%i)\n",alphabet_size,m);

    _alphabet_size=alphabet_size;
    _stationary=stationary;
    _m=m;
    unsigned long n=1;
    for (unsigned short i=0; i<=_m; i++)
      n*=_alphabet_size;

    for (unsigned long i=0; i<n; i++)
      _param.push_back(1.0/(double)_alphabet_size);

    normalize();

    if (_stationary)
      compute_mu0();

    if (verbose)
      print();
  };


  markov::markov(unsigned short alphabet_size,bool stationary,bool verbose) {
    if (verbose)
      printf(">>> call markov::markov(%i)\n",alphabet_size);

    _alphabet_size=alphabet_size;
    _stationary=stationary;
    _m=0;

    for (unsigned short a=0; a<_alphabet_size; a++)
      _param.push_back(1.0/(double)_alphabet_size);

    if (_stationary)
      _mu0=_param;

    if (verbose)
      print();
  };

  markov::~markov() {

  };

  void markov::print() {
    printf("alphabet_size=%i\tm=%i\n",_alphabet_size,_m);
    printf("param (size=%i)= [ ",_param.size());
    if (_param.size()<=20)
      for (vector<double>::iterator it=_param.begin(); it!=_param.end(); it++)
	printf("%f ",*it);
    else
      printf("more than 20 terms");
    printf("]\n");
    if (_stationary) {
      printf("mu0 (size=%i)= [ ",_mu0.size());
      if (_mu0.size()<=20)
	for (vector<double>::iterator it=_mu0.begin(); it!=_mu0.end(); it++)
	  printf("%f ",*it);
      else
	printf("more than 20 terms");
      printf("]\n");
    }
  };

  void markov::normalize() {

    unsigned long i=0;
    unsigned long imax=_param.size();
    while (i<imax) {
      double sum=0.0;
      for (unsigned long j=i; j<i+_alphabet_size; j++)
	sum+=_param[j];
      if (sum!=0.0)
	for (unsigned long j=i; j<i+_alphabet_size; j++)
	  _param[j]/=sum;
      i+=_alphabet_size;
    }
  };

  void markov::compute_mu0(bool verbose) {
    /* fixme: replace power computations by something more efficient */

    if (verbose)
      printf(">>> call markov:compute_mu0();\n");

    _mu0.resize((size_t)pow((double)_alphabet_size,(double)_m));
    if (verbose)
      printf("_mu0 size = %i\n",_mu0.size());

    // initialize _mu0 with a random distribution
    {
      double sum=0.0;
      for (vector<double>::iterator it=_mu0.begin(); it!=_mu0.end(); it++) {
	*it=rand()/(double)RAND_MAX;
	sum+=*it;
      }
      for (vector<double>::iterator it=_mu0.begin(); it!=_mu0.end(); it++)
	*it/=sum;
    }

    int iter=0;
    // main loop
    {
      double test=1.0;
      vector<double> aux;
      while ((test>1e-10) && (iter<MAX_ITER)) {
	iter++;
	// save _mu0
	aux=_mu0;
	//printf("aux[0]=%f mu0[0]=%f\n",aux[0],_mu0[0]);
	// reset _mu0 to zero
	for (vector<double>::iterator it=_mu0.begin(); it!=_mu0.end(); it++)
	  *it=0.0;
	//printf("aux[0]=%f mu0[0]=%f\n",aux[0],_mu0[0]);
	// compute _mu0 = aux x Pi
	unsigned long pos1=0;
	unsigned long pos2=0;
	unsigned long pos3=0;
	for (vector<double>::iterator it=_param.begin(); it!=_param.end(); it++) {
	  _mu0[pos1]+=*it*aux[pos3];
	  pos1++;
	  pos2++;
	  if (pos2==_alphabet_size) {
	    pos2=0;
	    pos3++;
	  }
	  if (pos1==_mu0.size()) {
	    pos1=0;
	  }
	}
	// compute test
	{
	  double testbis;
	  test=0.0;
	  for (unsigned long i=0; i<_mu0.size(); i++) {
	    testbis=fabs((_mu0[i]-aux[i])/aux[i]);
	    if (testbis>test)
	      test=testbis;
	  }
	}
	if (verbose)
	  printf("test=%e\n",test);

      }
    } // end main loop
    if (iter==MAX_ITER) {
      fprintf(stderr,"Markov model does not have a stationary distribution. Try again without the --stationary option !\n");
      exit(EXIT_FAILURE);
    }

    //print();
  };

};
