/* $Id: x_WordFam.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

// This file was originally part of the PWCAL package
// author: H. Richard

#include "x_WordFam.h"

//#include "const.h"
//#include "SeqUtils.h" 
//#include "Default.h"

//#include <ctype.h>
//#include <iostream.h>
//#include <fstream.h>
//#include <string.h>
//#include <math.h>

//Constructeur de base, on met tout a NULL

WordFam::WordFam() 
{
  _nalpha = 0 ;
  _nwords = 0 ;
  _size=0;
  _words = NULL ;
  _lwords = NULL ;
  //  _alphabet = NULL ;
  _overlap = NULL ;
}

// constructeur par recopie

WordFam::WordFam(const WordFam & W)
{
  _nalpha=W._nalpha ;
  _nwords=W._nwords ;
  _lwords=new int[_nwords] ;
  for(int i=0; i<_nwords; i++) _lwords[i]=W._lwords[i] ;

  _words = new long[_nwords] ;
  for(int i=0; i<_nwords;i++)
    _words[i] = W._words[i] ;

  if (isOverlap())
  {
    _overlap = new int*[_nalpha*_nalpha] ;
    for(int i=0; i<_nalpha*_nalpha; i++)
	{
	  _overlap[i] = new int[_lwords[i]] ;
	  for(int j=0; j < _lwords[i]; j++)
	     {
		_overlap[i][j] = W._overlap[i][j] ;
	     }
	}
  }
  else
  {
    _overlap=NULL;
  }
}

//Destructeur

WordFam::~WordFam() 
{
  int i;
  delete [] _words ;
  delete [] _lwords ; 
  if (isOverlap())
    {
      for (i=0;i<_nwords*_nwords;i++)
	delete _overlap[i];
      delete []_overlap;
    }
}

//Initialisation ou changement d'un mot
//marche aussi pour une famille de mots separés par des virgules (pas testé)
// not operational atm
int WordFam::SetWord(word *w)
{

  _nalpha=(w->get_alpha())->get_size();
  _nwords = 1 ;
  /* check memory requirements */
  if (_size<_nwords) {
    /* free old alloc if necessary */
    if (_lwords!=NULL)
      delete _lwords;
    if (_words!=NULL)
      delete _words;
    /* and alloc */
    _size=_nwords;
    _lwords=new int[_size];
    _words=new long[_size];
    //_lwords=(int *)malloc(sizeof(int)*_size);
    //_words=(long *)malloc(sizeof(long)*_size);
    //if (_lwords==NULL || _words==NULL) {
    //  fprintf(stderr,"not enough memory in WordFam\n");
    //  exit(EXIT_FAILURE);
    //}
  }

  _lwords[0]=w->get_size();
  _words[0]=w->get_code();

  return(EXIT_SUCCESS) ;
}

int WordFam::SetWord(word *v,word *w)
{

  _nalpha=(w->get_alpha())->get_size();
  _nwords = 2 ;
  /* check memory requirements */
  if (_size<_nwords) {
    /* free old alloc if necessary */
    if (_lwords!=NULL)
      delete _lwords;
    if (_words!=NULL)
      delete _words;
    /* and alloc */
    _size=_nwords;
    _lwords=new int[_size];
    _words=new long[_size];
  }

  _lwords[0]=v->get_size();
  _words[0]=v->get_code();
  _lwords[1]=w->get_size();
  _words[1]=w->get_code();

  return(EXIT_SUCCESS) ;
}

int WordFam::SetWord(pattern *patt)
{

  _nalpha=(patt->get_alpha())->get_size();
  _nwords = patt->get_optimized_word_list_size();
  /* check memory requirements */
  if (_size<_nwords) {
    /* free old alloc if necessary */
    if (_lwords!=NULL)
      delete _lwords;
    if (_words!=NULL)
      delete _words;
    /* and alloc */
    _size=_nwords;
    _lwords=new int[_size];
    _words=new long[_size];
    //_lwords=(int *)malloc(sizeof(int)*_size);
    //_words=(long *)malloc(sizeof(long)*_size);
    //if (_lwords==NULL || _words==NULL) {
    //  fprintf(stderr,"not enough memory in WordFam\n");
    //  exit(EXIT_FAILURE);
    //}
  }

  {
    word *w=NULL;
    for (int i=0; i<_nwords; i++) {
      w=patt->get_optimized_word_list()[i];
      _lwords[i]=w->get_size();
      _words[i]=w->get_code();
    }
  }
  return(EXIT_SUCCESS) ;
}

// ATTENTION: On considere qu'il n'y a pas d'overlap 
// entre "acgt" et "cg"
void WordFam::SetOverlap()
{
  int i,j,k,ind,l1,l2;
  int mot1,mot2;

  _overlap=new int*[_nwords*_nwords];

  for (i=0;i<_nwords;i++)
    for (j=0;j<_nwords;j++)
      {
	mot1=_words[i];
	l1=_lwords[i];
	mot2=_words[j];
	l2=_lwords[j];
	// Si le 2eme mot est plus longs, on enleve
	// les dernieres lettres
	for (k=l1;k<l2;k++)
	  {
	    mot2=mot2/_nalpha;
	    l2--;
	  }
	ind=i*_nwords+j;
	_overlap[ind]=new int[_lwords[i]];
	//construction de l'overlap entre i et j
	for (k=0;k<_lwords[i];k++)
	  {  
	    if  (l1==l2)
	      {
		_overlap[ind][k]=(mot1==mot2)?1:0;
		mot1=mot1%((int)pow((double)_nalpha,(double)(l1-1)));
		l1--;
		mot2=mot2/_nalpha;
		l2--;
	      }
	    else if (l1>l2)
	      {
		_overlap[ind][k]=0;
		mot1=mot1%((int)pow((double)_nalpha,(double)(l1-1)));
		l1--;
	      }
	  }
      }
}


//operateurs
WordFam &WordFam::operator= (const WordFam & W)
{
  if (this!= &W)
    {//Destruction des anciens elements dynamiques
      delete [] _lwords ;
      delete [] _words ;
      _nalpha=W._nalpha ;
      _nwords=W._nwords ;
      _lwords=new int[_nwords] ;
      for(int i=0; i<_nwords; i++) _lwords[i]=W._lwords[i] ;
      //(*_alphabet) = (*(W._alphabet)) ;
      //Affectations
      _words = new long[_nwords] ;
      for(int i=0; i<_nwords;i++)
	{
	  _words[i]=W._words[i] ;
	}
    }
  return(*this) ;
}


void WordFam::Print() const
{
  int i, j, k ;
  
  cout  << _nwords << " mots dans un alphabet a " 
	<< _nalpha << " lettres\n" ;

  cout << "alphabet utilise :\n" ;
  //(*_alphabet) -> Print() ;
  cout << "mots de la famille : " << endl ;
  //w_curr = "toto" ;
  for(i=0; i<_nwords; i++) 
    {
      cout << setw(5) << _words[i] 
	   <<"  longueur : " << _lwords[i] 
	   << endl ;
      //delete w_curr ;
    }

  if (isOverlap())
    {
      cout << "liste des overlaps pour les mots \n" ;
      for(i=0; i<_nwords; i++) 
	{
	  for(j=0; j<_nwords; j++) 
	    {
	      cout << _words[i] << " et " << _words[j] << "   ";
	      for (k=0;k<_lwords[i];k++)
		cout << _overlap[i*_nwords+j][k] << "   ";
	      cout << endl;
	    }
	}
    }

}
