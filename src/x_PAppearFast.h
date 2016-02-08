/* $Id: x_PAppearFast.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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

//Reimplementation de la classe PAppear dans une version allégée 
//et rapide --> On ne stocke plus que la proba en cours et la proba precedente

#ifndef X_PAPPEARFAST_H
#define X_PAPPEARFAST_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "markov.h"
#include "x_WordFam.h"
#include "x_PSucceed.h"

class PAppearFast
{
  private :
  
  const markov *_M;
  int _rank;
  double *_Mu_W;
  double *_power;
  double _stat;

  int _Nobs ; //Current number of Occurrence calculated
  int _nw ; //number of words
  int _lseq ; //length of the sequence -> initialised to lseqmax
  long _truelseq;
  const WordFam **_Wp ; //family of words
  //powers of the matrix must have be initialized !!
  //two tabs of dim _lseq*_nw :
  double **_ProbCur ; //probabilities for _Nobs obs from 1 to lseqmax 
  double **_ProbPrev ; //same for Nobs - 1
  public :

  //Constructeur
  //Classical initialisation 
  PAppearFast(const markov *M,double *power,double stat,const WordFam * W ,long lseq,long limit);
  
  ~PAppearFast() ;
  //proba of first apparition -> must use this one before the others
  int PFirst(const PSucceed & Q_WW) ;
  //proba of first apparition conditionned by a given word
  int PFirst(const PSucceed & Q_WW, const PSucceed & Q_VW) ;
  
  //Proba of (n+1)_th apparition
  int PNext(const PSucceed & Q_WW);
  //goto the occurrence Nocc by successive calls to PNext
  int GoToNobs(const PSucceed & Q_WW, int Nocc)
    {
      //printf("GoToNobs(%i):\t_Nobs=%i\n",Nocc,_Nobs);
      if (_Nobs >= Nocc) 
	{
	//  cerr << "No need to call GoToNobs here\n" 
	//       << "current _Nobs : " << _Nobs << endl 
	//       << "Nobs asked : " << Nocc << endl ; 
	  return(Nocc - _Nobs) ;
	}
      for (int i= _Nobs; i< Nocc; i++)
	{
	  //printf("GoToNobs(%i):\ti=%i\n",Nocc,i);
	  this -> PNext(Q_WW) ;
	}
      return(_Nobs) ; 
    }

  int Nobs() const { return(_Nobs) ;}
  int lseq() const { return(_lseq) ;}
  
  //Calculus of proba, stat and Nocc for a given ArrayDnaObsSeq
  //les vecteurs doivent avoir ete declares de la bonne taille
  //renvoie 0 si c'est bon et -<nombre d'erreurs de calcul sinon>
  int GetStatsFromSeq(int* Nocc,double* Probas,long nobs);// double* GaussStats) ;
  
  

  //operateurs
  //proba to see a word of the family at least _Nobs times in the sequence
  //of length l ; l must be <= _lseq
  double operator()(int l) const 
    {
      double prob=0 ;
      if (l>_lseq)
	{
	  cerr << "Bad call for calc of proba !!!" << endl ;
	  return(-1) ;
	}

      for (int w=0; w < _nw ; w++)
	for (int i=1; i<=l; i++) 
	  prob += _ProbCur[i][w] ;
      if (_Nobs==0)
	prob = 1-prob ;
      if ((prob <=0) || (prob >= 1)) 
	{
//	  logfile << "*****Problem with words : " ;
//	    for (int i=0; i < (**_Wp).nw()-1; i++) 
//	      logfile << (**_Wp).GetWord(i) <<", " ;
//	  logfile << (**_Wp).GetWord((**_Wp).nw()-1)  << endl ;

	}
      return(prob) ;
    }
  

} ;


#endif //fin P_APPEAR_H
