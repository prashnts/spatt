/* $Id: x_WordFam.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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

/////////////*********Gestion des familles de mots***************/////////
/************** le type d'alphabet n'est pas tres fonctionnel**************/
//On ne peut par exemple pas fusionner de maniere simple deux familles de 
//mots... on ne teste pas non plus les differences entre numerotation....


#ifndef X_WORDFAM_H
#define X_WORDFAM_H

#include "word.h"
#include "pattern.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
//#include <fstream>
#include <iomanip>
//#include <string>
//#include <iostream.h>
//#include <fstream.h>
//#include <iomanip>
//#include <glib.h>
#include <math.h>

using namespace std;
using namespace spatt;

class WordFam
{ 
 protected :
  
  int _nalpha;  //redondant...
  int _nwords;  //n° of words in the family
  long *_words; //vector of the words coded in base Alphabet.GetSize()
  int *_lwords; //length of each word
  int **_overlap ; //Bits d'overlap entre les mots. Le tableau
  int _size;
  //entre le mot n° i et j se trouve sur la ligne
  // i*(_nwords-1)+j

 public :

  WordFam() ;
  WordFam(const WordFam &); // constructeur par recopie
  ~WordFam() ;

  //initialisations pour des mots simples.
  
  //Initialise un mot (permet de le changer aussi)
  int SetWord(word *w);
  int SetWord(word *v,word *w);
  int SetWord(pattern *patt);

  //Methodes
  void SetOverlap(); //calcul la matrice d'overlap entre les mots
  inline bool isOverlap() const { return (_overlap!=NULL);} ;
  inline int nw() const { return(_nwords) ;} 
  inline int lw(int i)  const //On démarre le tableau a zero
    {if (i<_nwords) return(_lwords[i]) ;  else return(-1) ;} 
  inline int nalpha() const { return(_nalpha) ;} 

  //decode le char* correspondant au n_eme mot
  //char *GetWord(int n) const ;
  //Operateurs
  WordFam &operator= (const WordFam & W) ; 
  
  //nth word maybe with a define INVALID_WORD
  int operator [] (const int n) const
  {if (n>_nwords) return(-1) ; else return(_words[n]) ; }
  
  //vec of length _lwords[n]
  int* operator() (int n) const 
    {
      if(n>=_nwords)
	{cerr << "subscript of word out of bound\n" ; exit(-1) ;}
      int l = _lwords[n] ;
      int currword ;
      int* res = new int[_lwords[n]] ;
      int i= (int) pow((double)_nalpha,(double)(l-1)) ;
      currword = _words[n] ;
      for (int j=0; j<l ; j++ ) 
	{
	  res[j] = currword / i ;
	  currword %= i ;
	  i/= _nalpha ;
	}
      return(res) ;
    }
  
  // k_th letter of n_th word (prefer the generic [] method)
  int operator() (int n, int k) const
  {
    if(n>=_nwords)
      {cerr << "subscript of word out of bound\n" ; exit(-1) ;}
    else if (k>=_lwords[n])
      {cerr << "subscript of word out of bound\n" ; exit(-1) ;}
    else 
      { 
	int i=_nalpha ;
	int let ;
	for(int j=0; j<k; j++) i*=_nalpha ;
	let = (_words[n] % i ) / (i/_nalpha) ;
	return(let) ;
      }
  }

  
  void Print() const;
  
  //FAIRE UNE METHODE QUI RENVOIE UN TABLEAU DE BOOL
  
  //Does V_i and V_j overlap on u letters ?
  inline bool Epsilon(int i,int j,int u) const ;
  
  //The Same for V_i and W_j
  inline bool Epsilon(int i, const WordFam & W, int j, int u) const ;    
};

//do the u last letters of Vi are the same than the u of Vj ?
bool WordFam::Epsilon(int i,int j,int u) const
{
  int k1,k2,l;
  if(i<_nwords && j<_nwords)
    if(u<=_lwords[i] && u<=_lwords[j])
      {
	k1=k2=1 ;
	for(l=0; l<_lwords[i]-u; l++) k1*=_nalpha ;
	for (l=0; l<u; l++) k2*= _nalpha ;
	return( (_words[i]%k2) == (_words[j]/k1) );
      }
    else 
      {
	cerr << "Error in method Epsilon, length of"
	     << " at least one word out of bound\n";
	exit(-1) ; 
      }
  else 
    {
      cerr << "Error in method Epsilon, n° of word out of bound\n";
      exit(-1) ; 
    }
}

bool WordFam::Epsilon(int i, const WordFam & W, int j, int u) const
{
  int k1,k2,l;
  //test simplement des dim pour l'instant
  if(_nalpha!=W._nalpha)
    cerr << "***********WARNING************\n"
	 << "Dim of alphabets don't seem to agree\n" ;
  if(i<_nwords && j<W._nwords)
    if(u<=_lwords[i] && u<=W._lwords[j])
      {
       	k1=k2=1 ;
	for(l=0; l<_lwords[i]-u; l++) k1*=_nalpha ;
	for (l=0; l<u; l++) k2*= _nalpha ;
	return ( (_words[i]%k2)==(W._words[j]/k1) );
      }
    else 
      {
	cerr << "Error in method Epsilon, length of "
	     << "at least one word out of bound\n";
	exit(-1) ;
      }
  else 
    {
      cerr << "Error in method Epsilon, n° of word out of bound\n" ;
      exit(-1) ;
    }
}


#endif /* WORDFAM_H */
