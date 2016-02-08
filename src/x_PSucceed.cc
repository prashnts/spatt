/* $Id: x_PSucceed.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "x_PSucceed.h"

//constante pour les probas qui n'ont pas encore ete calculees

extern const int PROBA_NOT_DEFINED = -200 ;

/* several call to old Markov class to replace */
/* search "to replace" */


/* not used */
//PSucceed::PSucceed(markov *M,double *power,const WordFam * V,const WordFam * W, long lseq)
//{
//  int i,j,k;
//
//  _M=M;
//  _power=power;
//
//  //initialisation des elements dynamiques
//  _Wp = new const WordFam* ;
//  _Vp = new const WordFam* ;
//
//  *_Wp = W ;
//  *_Vp = V ;
//  
//    
//  _nv = (**_Vp).nw() ;
//  _nw = (**_Wp).nw() ;
//  
//  _lseq = lseq ;
//  _rank = 0 ;
//
//  /* very ugly allocation */
//  _Probas = new double**[_nv] ;
//  for (i=0; i<_nv; i++)
//    {
//      _Probas[i] = new double*[_nw] ;
//      for (j=0; j<_nw; j++)
//	{
//	  _Probas[i][j]=new double[_lseq+1] ; //???????
//	  for(k=0; k<=_lseq;k++)
//	    _Probas[i][j][k] = PROBA_NOT_DEFINED ;
//	}
//    }
//}


PSucceed::PSucceed(const markov *M,double *power,const WordFam * W,long lseq)
{
  int i,j,k;
  int nw;

  _M=M;
  _power=power;
  //printf("_power=%p\n",_power);

  _lseq = lseq ;
  _rank =0 ;
  //initialisation des elements dynamiques
  _Wp = new const WordFam* ;
  _Vp = new const WordFam* ;
  *_Vp = W ;
  *_Wp = W ;

  //WARNING !! TENSOR OF PROBAS BEGINS AT POS 1
  //Must change this type of allocation (prefer MemChunk)
  _nw = (**_Wp).nw() ;
  _nv = _nw ;
  /* very ugly allocation */
  _Probas = new double**[_nw] ;
  for (i=0; i<_nw; i++)
    {
      _Probas[i] = new double*[_nw] ;
      for (j=0; j<_nw; j++)
	{
	  _Probas[i][j]=new double[_lseq+1] ; //???????
	  for(k=0; k<=_lseq;k++)
	    _Probas[i][j][k]= PROBA_NOT_DEFINED ;
	}
    }
  
}


//Destruction
PSucceed::~PSucceed() 
{
  int i,j;
  for (i=0; i<_nv; i++)
    {
      for (j=0; j<_nw; j++)
	{
	  delete [] _Probas[i][j];
	}
      delete [] _Probas[i] ;
    }
  delete _Wp ;
  delete _Vp ;

}
  

//Proba de succession sans contrainte sur les mots de W
int PSucceed::Q_V_W() 
{
  int i,j,k,x ;
  int lw,lv ;
  //fin de V_i (les _order dernieres lettres)
  int V_end ;
  int W_beg ;
  int W_tmp ;
  int left ;
  int right ;
  int ind ;

  /* old stuff */
  //int nalpha = (**_Mp).nalpha() ; // number of letters in the alphabet
  //int order = (**_Mp).order() ; // order of the model
  //int nMu = (**_Mp).nMu() ; // dim of mu (_nalpha^_order)
  //_rank = (**_Mp).rank() ; // how many steps to converge to Mu
  /* end old stuff */
  int nalpha=_M->get_alpha()->get_size();
  int order=_M->get_order();
  int nMu=_M->get_n();
  double **Tau_W = TauTab(**_Wp) ;
  
  /* get _rank, the worse value for the moment */
  _rank = (int)(log(1e-16)/log(_M->get_secondmag()))+1;

  //pour tous les V_i
  for (i=0; i < _nv; i++)
    {
      //cout << "appel de Mp.nMu() pour le mot V_" << i  << endl ;
      V_end = (**_Vp)[i] % nMu ;
      lv = (**_Vp).lw(i) ;
      //pour tous les W_j
      for (j=0; j < _nw; j++)
	{
	  lw = (**_Wp).lw(j) ;
	  W_beg = (**_Wp)[j];
	  //attention, la famille de mots doit 
	  //absolument etre declaree dans l'alphabet phasé
	  for(k=0; k < lw - order ; k++) W_beg /= nalpha ;
	  //cas ou les mots peuvent se recouvrir
	  //Prise en compte du cas ou V et W se recouvrent completement
 	  if ((lw <= lv) && (**_Vp).Epsilon(i,(**_Wp),j,lw))
	    _Probas[i][j][0] = 1.0 ;
	  else 
	    _Probas[i][j][0] = 0.0 ;
	  for(x=1; x<lw; x++)
	    {
	      if ( (x>=lw-lv) && ((**_Vp).Epsilon(i,(**_Wp),j,lw-x)) )
		//si les mots peuvent se recouvrir, on doit faire attention
		//a l'ordre et rajouter les termes necessaires
		{
		  if ((lw-x)<order)
		    {
		      _Probas[i][j][x] = 1.0 ;
		      left = V_end ;
		      //on ne garde que les order - (lw-x) dernieres lettres
		      //de W_beg ;
		      ind = (int) pow((double)nalpha,(double)(order-lw+x)) ;
		      W_tmp = W_beg % ind ;
		      ind /= nalpha ;
		      while(ind != 0)
			{
			  right = W_tmp / ind ;
			  left *= nalpha ;
			  left += right ;
			  /* old stuff */
			  //_Probas[i][j][x] *=(**_Mp)[left] ; // return _Pi[left]
			  _Probas[i][j][x] *=(_M->get_model())[0][left];
			  left %= nMu ;
			  W_tmp %= ind ;
			  ind /= nalpha ;
			}
		      _Probas[i][j][x] *= Tau_W[j][1] ;
		    }
		  else
		    {
		      _Probas[i][j][x]=Tau_W[j][lw-x-order+1] ;
		    }
		}
	      //on impose 0 si W est trop long
	      else 
		_Probas[i][j][x]=0.0 ;
	    }
	  //proba jusqu'a rank
	  for(x=lw; x<_rank ; x++)
	    {
	      //D'apres la def des puissances de la matrice
	      //on dit bien rajouter order pour arriver en W_beg
	      /* old stuff */
	      // (**_Mp)(i,j,step) return _PowPi[step-1][i][j]
	      //_Probas[i][j][x] = (**_Mp)(V_end,W_beg,x-lw+order)* Tau_W[j][1] ;
	      /* careful here, power are stocked in the fortran style: column,row */
	      _Probas[i][j][x]=_power[(x-1)*nMu*nMu+j*nMu+i]*Tau_W[j][1];	      
	    }
	  //quand on arrive en régime stationnaire
	  //(a la longueur de Wj près)
	  _Probas[i][j][_rank] = Tau_W[j][0] ;
	}//fin W_j
      
    }//fin V_i
  for (i=0; i<_nw ; i++)
    delete [] Tau_W[i] ;
  delete [] Tau_W ;

  return(1) ;
  
}


//Proba avec contrainte de premiere 
//apparition sur les mots de W
int PSucceed::P_V_W(const PSucceed & Q_VW, const PSucceed & Q_WW) 
{
  int i,j,l,x,z;
  int lv,lw ;
  double p_sum_rec ;
  /* old stuff */
  //double * Mu_W = (**_Mp)(**_Wp) ; // operator() on a WordFam : stationary proba of it
  /* alloc and fill Mu_W */
  double *Mu_W;
  {
    int nw=(**_Wp).nw();
    Mu_W= new double[nw];
    for (int i=0; i<nw; i++){
      Mu_W[i]=_M->mu((**_Wp)[i],(**_Wp).lw(i));
    }
  }

  //Tests : QWW est-il bien defini ?
  if ((*Q_WW._Wp) != (*Q_WW._Vp))
    {
      cerr << "**********WARNING***************\n"
	   << "La proba de succession sans contrainte n'est pas definie\n"
	   << "sur le meme groupe de mots\nExiting\n"  ;
      exit(-1) ;
    }
  //On teste le W de l'objet contre celui de QWW ici :
  if (_nw == Q_WW._nw)
    {
      for (i=0; i<_nw;i++)
	if ( ((**_Wp)[i])!=((**Q_WW._Wp)[i]))
	  {
	    cerr << "Error in P_V_W of class PSucceed, words " 
		 << i << "don't agree\n" 
		 << "You should define your families with the same WordFam"
		 << " object\n" ;
	    exit(-1) ;
	  }
    }
  else
    {
      cerr << "Error in P_V_W of class PSucceed, dim of motifs "
	   << "don't seem to agree\nExiting" ;
      exit(-1) ;
    }

  //On teste le V de l'objet contre celui de Q_VW ici :
  if (_nv == Q_VW._nv)
    {
      for (i=0; i<_nw;i++)
	if ( ((**_Wp)[i])!=((**Q_VW._Wp)[i]))
	  {
	    cerr << "Error in P_V_W of class PSucceed, words " 
		 << i << "don't agree\n" 
		 << "You should define your families with the same WordFam"
		 << " object\n" ;
	    exit(-1) ;
	  }
    }
  else
    {
      cerr << "Error in P_V_W of class PSucceed, dim of motifs "
	   << "don't seem to agree\nExiting" ;
      exit(-1) ;
    }

  //On pourrait aussi tester les chaines de Markov etc...


  //a chaque position avant rank
  for (x=1; x <= Q_VW._rank ; x++)
    {
      //pour tous les V_i
      for (i=0 ; i< _nv; i++)
	{
	  lv = (**_Vp).lw(i) ;
	  //pour chaque W_j
	  for (j=0 ; j< _nw; j++)
	    {
	      lw = (**_Wp).lw(j) ;
	      
	      _Probas[i][j][x] = Q_VW(i,j,x) ;
	      //terme a soustraire
	      p_sum_rec = 0.0 ;
	      for (z=1 ; z < x ; z++)
		for (l=0 ; l< _nw; l++)
		  //erreur d'arrondi possible ici
		  p_sum_rec += (_Probas[i][l][z]*Q_WW(l,j,x-z)) ;
	      //tester le rapport p_sum_rec terme de dte si besoin
	      
	      _Probas[i][j][x] -= p_sum_rec ;
	    }
	}
    }
  //l'autre partie de la recurrence
  for (x=Q_VW._rank+1 ; x <= _lseq ; x++)
    {
      //pour tous les Vi
      for(i=0; i < _nv; i++)
	{
	  
	  //pour tous les Wj
	  for(j=0; j<_nw; j++)
	    {
	      _Probas[i][j][x] = _Probas[i][j][x-1] ;

	      //recurrence allegee
	      p_sum_rec =0 ;
	      for (z = x-Q_WW._rank + 1; z < x-1 ; z++)
		for (l=0; l< _nw; l++)
		  {
		    //tester le rapport pour les pbmes d'arrondis
		    p_sum_rec+=( (_Probas[i][l][z-1]-_Probas[i][l][z]) \
				*Q_WW(l,j,x-z)) ;
		  }
	      _Probas[i][j][x] += p_sum_rec ;
	      
	      //terme residuel
	      p_sum_rec = 0 ;
	      for (l=0; l < _nw; l++)
		{
		  p_sum_rec += _Probas[i][l][x-Q_WW._rank] ;
		}
	      p_sum_rec *= Mu_W[j] ;
	      
	      _Probas[i][j][x] -= p_sum_rec ;
	    }
	}
    }
  delete Mu_W ;
  return(1) ;  
}

double **PSucceed::TauTab(const WordFam & W) const
{

  int _nalpha=_M->get_alpha()->get_size();
  int _order=_M->get_order();
  int _nMu=_M->get_n();
  int _nPi=_nMu*_nalpha;
  double *_Mu=_M->get_stationary();
  double *_Pi;
  {
    double **tmp=_M->get_model();
    _Pi=tmp[0];
  }

  /* old stuff */
  //int _nalpha;
  //int _order;
  //double *_Mu;
  //int _nMu;
  //int _nPi;
  //double *_Pi;
  /* end old stuff */

  int nw = W.nw(), lw ;
  int word,word_tmp, div;
  int i,j ;
  if ( W.nalpha() != _nalpha ) 
    {
      cerr << "*******Warning*************\n"
	   << "in method TauTab, length of alphabets for Markov chain and\n"
	   << "word family don't seem to agree\n" 
	   << "(" << W.nalpha() << " and " << _nalpha << ")\n"
	   << "(maybe a problem with a phased Markov chain)\n" ;
      exit(-1) ;
    }
  double **TauFamTab = new double*[nw] ;

  for (i=0 ; i<nw ; i++) 
    {
      lw = W.lw(i) ;
      if (lw < _order) 
	{
	  cerr << "Error for stationnary proba of a word :\n"
	       << "length lower than order for the " << i <<"_th word\n" ;
	  exit(-1) ;
	}
      if (lw == _order)
	{
	  //on doit rajouter une valeur fictive pour certains calculs
	  //p exple Q_V_W
	  TauFamTab[i] = new double[2] ;
	  /* old stuff */
	  //TauFamTab[i][0] = (isMu()) ? _Mu[W[i]] : -1 ;
	  TauFamTab[i][0] =_Mu[W[i]];
	  TauFamTab[i][1] = 1.0 ;
	}
      else
	{
	  //ici lw > _order
	  TauFamTab[i] = new double[lw-_order+1] ;
	  word = W[i]  ;
//  	  cout << "word  : " << word 
//  	       << ", length : " << lw <<endl ;
	  //printf("word=%i\tlength=%i\n",word,lw);
	  lw = W.lw(i) ;
	  div = 1 ;
	  //printf("_nPi=%i\n",_nPi);
	  word_tmp = word % _nPi ;
	  //printf("word_tmp=%i\n",word_tmp);
	  word /= _nPi ;
	  //printf("word_tmp=%i\n",word_tmp);
	  TauFamTab[i][lw-_order] = _Pi[word_tmp] ;
	  word_tmp /= _nalpha ;
	  for(j=lw-_order-1; j>0 ; j--)
	    {
	      word_tmp = (word % _nalpha)*_nMu + word_tmp ;
	      TauFamTab[i][j] = _Pi[word_tmp]*TauFamTab[i][j+1] ;
	      word_tmp /= _nalpha ;
	      word /= _nalpha ;
	      //_nMu = _nalpha^order , so we keep only the order last letters
	    }
	  /* old stuff */
	  //TauFamTab[i][0] = (isMu()) ? (_Mu[word_tmp]*TauFamTab[i][1]) : -1 ;
	  TauFamTab[i][0] =(_Mu[word_tmp]*TauFamTab[i][1]) ;
	}
    }
  return(TauFamTab) ;
}

