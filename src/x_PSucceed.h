/* $Id: x_PSucceed.h 440 2005-07-20 12:34:26Z mhoebeke $ */
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

//Calcul de des differents types de probas
//de succession, avec ou sans contrainte sur le
//groupe de mots

/* legal mention copyright H. Richard and G. Nuel */

#ifndef X_PSUCCEED_H
#define X_PSUCCEED_H

#include "markov.h"
#include "x_WordFam.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "math.h"


//Probas of succession from V to W
//from end of word to end of word...
class PSucceed
{
  private :

  const markov *_M;
  double *_power;

  int _lseq ;
  int _rank ; //si rank!=0 on travaille avec Q_V_W
  int _nv ;
  int _nw ;
  //  const Alphabet **_alphabet ;
  const WordFam **_Vp ;
  const WordFam **_Wp ;
  //  const Markov **_Mp ;

  //WARNING, TENSORS OF PROBAS MUST BEGIN AT POS 1
  //(probas[0] should be the complete overlapping of words)
  //dim is _nv *_nw * (_rank != 0) ? _rank : _lseq+1
  double ***_Probas;
  double **TauTab(const WordFam &) const ; 

  public :
  //A-T-ON vraiment besoin de alphabet ???????
  //PSucceed(markov *M,double *power,const WordFam * V,const WordFam * W,long lseq);
  PSucceed(const markov *M,double *power,const WordFam * W,long lseq);
  ~PSucceed() ;

  //**********operateurs****************
  //Extraction de la proba de Vi vers Wj en k etapes
    
  inline double operator() (int i, int j, int x) const 
    {
      if ( (_rank!=0) && (x>_rank) )
	return( _Probas[i][j][_rank]) ;
      else
	return (_Probas[i][j][x]) ;
    }

  //Proba de succession sans contrainte sur les mots de W
  int Q_V_W() ;
  //Proba avec contrainte de premiere 
  //apparition sur les mots de W
  int P_V_W(const PSucceed & Q_VW, const PSucceed & Q_WW) ;
  //on pourrait definir P_W_W qui n'a besoin que d'un argument 
  //(pour les motifs par exemple)

};



#endif
