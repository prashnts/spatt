/* $Id: x_PAppearFast.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

#include "x_PAppearFast.h"

//for Random number generation
//#include <gsl/gsl_rng.h>

//Constructeurs

//RAJOUTER UN TEST SUR LA LONGUEUR DES MOTS ET L'ORDRE DU MARKOV
//Classical initialisation
PAppearFast::PAppearFast(const markov *M,double *power,double stat,const WordFam * W ,long lseq,long limit)
{

  _M=M;
  _power=power;
  _stat=stat;

  //initialisation des elements dynamiques
  _Wp = new const WordFam* ;
  *_Wp =  W ;


  _nw = (**_Wp).nw() ;
  _Nobs = 0 ;  
  _lseq = limit;
  _truelseq= lseq;
  //_lseq = lseq ;
  //printf("in PAppearFast:\n\t_lseq=%i\n",_lseq);
  //Init des 2 tableaux de probas
  _ProbCur = new double*[_lseq+1];
  _ProbPrev = new double*[_lseq+1];
  double* tmpcur = new double[(_lseq+1)*_nw] ;
  double* tmpprev = new double[(_lseq+1)*_nw] ;
  for (int i=0; i <= _lseq ; i++)
    {
      _ProbCur[i] = tmpcur + i*_nw ;
      _ProbPrev[i] = tmpprev + i*_nw ;
      for (int j=0; j<_nw;j++)
	_ProbCur[i][j] =  _ProbPrev[i][j] = 0 ;
    }

}

//destructor

PAppearFast::~PAppearFast()
{
  int j ;

  //il n'y a que de 2 allocations 
  //printf("~PAppearFast() start\n");
  delete [] _ProbCur[0] ;
  delete [] _ProbPrev[0] ;
  delete [] _ProbCur ;
  delete [] _ProbPrev ;
  delete _Wp ;
  delete [] _Mu_W;
  //printf("~PAppearFast() end\n");
}


//Calculus of the first proba
//without first word
int PAppearFast::PFirst(const PSucceed & Q_WW) 
{

  int i,j ;
  int l,k,x,z;//for the recursion
  double p_sum_rec ;
  //temporary variables
  int lw;
  /* old stuff */
  //int rank=(**_Mp).rank() ;
  int rank;
  {
     rank = (int)(log(1e-16)/log(_M->get_secondmag()))+1;
  }
  _rank=rank;
  //  cout << "Calcul de Mu_W" << endl ;
  /* old stuff */
  //double * Mu_W = (**_Mp)(**_Wp) ;
  /* alloc and fill Mu_W */
  double *Mu_W;
  {
    int nw=(**_Wp).nw();
    Mu_W= new double[nw];
    for (int i=0; i<nw; i++){
      Mu_W[i]=_M->mu((**_Wp)[i],(**_Wp).lw(i));
    }
  }
  _Mu_W=Mu_W;

  //cout << "fin du calcul " << endl ;
  //rajouter des tests sur les pointeurs de mots etc.

  if (_Nobs)
    {
      cerr << "Warning, Launching PFirst on au tab already initialized\n" 
	   << "Reseting values...\n" ;
      for(i=0; i<=_lseq; i++)
	for(j=0; j<=_nw; j++)
	  _ProbCur[i][j] = _ProbPrev[i][j] = 0 ;	  
    }

  _Nobs = 1 ;

  //partie 1
  for(x = 1; x <= rank; x++)
    {

      for(i=0; i<_nw; i++)
	{
	  lw = (**_Wp).lw(i) ;
	  _ProbCur[x][i] = (x<lw) ? 0.0 : Mu_W[i] ;
	  p_sum_rec = 0.0 ;
	  
	  for (z=lw; z < x ; z++)
	    for(l=0 ; l<_nw; l++)
	      //problemes d'arrondis possibles ici
	      p_sum_rec += ( _ProbCur[z][l]*Q_WW(l,i,x-z) ) ;
	  
	  _ProbCur[x][i] -= p_sum_rec ;
	}
    }
  
  //approximation
  for(x = rank+1; x <= _lseq; x++)
    {

      for(i=0; i<_nw; i++)
	{
	  lw = (**_Wp).lw(i) ;
	  
	  _ProbCur[x][i] = _ProbCur[x-1][i] ;
	  p_sum_rec = 0.0 ;
	  
	  //terme 1
	  for (z=x-rank+1; z < x ; z++)
	    for(l=0 ; l<_nw; l++)
	      //problemes d'arrondis possibles ici
	      p_sum_rec += ((_ProbCur[z][l] - _ProbCur[z-1][l])\
			    *Q_WW(l,i,x-z)) ;

	  _ProbCur[x][i] -= p_sum_rec ;

	  p_sum_rec = 0.0 ;
	  //terme 2
	  for(l=0 ; l<_nw; l++)
	    //problemes d'arrondis possibles ici
	    p_sum_rec += _ProbCur[x-rank][l] ;
	  p_sum_rec *= Mu_W[i] ;

	  _ProbCur[x][i] -= p_sum_rec ;

	}
    }

  //delete [] Mu_W ;
  return(1) ;
}


//Proba of (n+1)_th apparition, multiple type of calculus
//we don't use the second variable for the moment
int PAppearFast::PNext(const PSucceed & Q_WW)
{
  int i,j ;
  int l,k,x,z;//for the recursion
  double p_sum_rec ;
  //temporary variables
  int lw ;
  /* old stuff */
  //int rank=(**_Mp).rank() ;
  int rank=_rank;
  /* old stuff */
  //double * Mu_W = (**_Mp)(**_Wp) ;
  double * Mu_W=_Mu_W;
  
  //
  _Nobs++ ;
  double** ptrtmp;
  //Substitution des tableaux  PPrev <-> PCur et PCur<- 0
  ptrtmp = _ProbPrev ;
  _ProbPrev = _ProbCur ;
  _ProbCur = ptrtmp ;
  for (i=0 ; i<= _lseq; i++)
    for (j=0; j<_nw; j++)
      _ProbCur[i][j] = 0 ;

  //premiere partie 
  //les _Nobs+lw-1 premieres valeures sont nulles
  for (x=1 ; x <= rank+_Nobs ; x++)
    {
      //j'ai un doute pour les premieres valeurs de W
      for(i=0 ; i <_nw ; i++)
	{
	  lw = (**_Wp).lw(i) ;
	  p_sum_rec = 0.0 ;

	  for (z=1; z<x; z++)
	    for(l=0; l<_nw; l++)
	      p_sum_rec += \
		( (_ProbPrev[z][l] - _ProbCur[z][l])\
		  *Q_WW(l,i,x-z) ) ;
	  
	  _ProbCur[x][i] += p_sum_rec ;

	}
      
    }

  //avec l'approximation
  for (x=rank + _Nobs + 1 ; x <= _lseq ; x++)
    {
      
      for(i=0 ; i <_nw ; i++)
	{
	  _ProbCur[x][i] = _ProbCur[x-1][i] ;

	  p_sum_rec = 0.0 ;
	  
	  for(l=0 ; l< _nw ; l++)
	    p_sum_rec += (_ProbPrev[x-rank][l]\
			  - _ProbCur[x-rank][l]);
	  
	  p_sum_rec *= Mu_W[i] ;
	  
	  _ProbCur[x][i] += p_sum_rec ;

	  p_sum_rec = 0.0 ;

	  for (z=x-rank+1; z < x ; z++)
	    for(l=0 ; l<_nw; l++)
	      p_sum_rec +=( (_ProbPrev[z][l]+_ProbCur[z-1][l]\
			     -_ProbCur[z][l]-_ProbPrev[z-1][l])\
			    *Q_WW(l,i,x-z)) ;	      

	  _ProbCur[x][i] += p_sum_rec ;
	}
      
    }

  //delete Mu_W ;
  //inutile
  return(1) ;
}

//Calcule proba, nombre d'occ et stat pour l'ArrayDnaObsSequence donne
int PAppearFast::GetStatsFromSeq(int* Nocc,double* Probas,long nobs)// double* GaussStats) 
{
  //int nseq = arrseq -> NbSequence() ;
  int nseq = 1;
  int seqcount ;
  int nseqequal ;
  int nerrors = 0 ;
  double prob ;
  //La proba de succession 
  /* old stuff */
  //PSucceed QVW(*_Wp, arrseq->LengthSeqM) ; // length of the sequence
 //printf("_power=%p\n",_power);
  PSucceed QVW(_M,_power,*_Wp, 42) ;
  QVW.Q_V_W() ;
  //pour ranger les sequences dans l'ordre du nombre d'observations
  //printf("nseq=%i\n",nseq);
  size_t* SortIndex=new size_t[nseq] ; // no use if nseq=1
  SortIndex[0]=0;

  
  /* old stuff */
  //int* tmpcount = GetNoccWordFam(**_Wp, arrseq) ; // return a array of size nseq (here 1) with number of occurrences of the wordfam
  int* tmpcount=NULL;
  /* alloc and fill */
  //{
  //  int nw=(**_Wp).nw();
  //  tmpcount=new int[1];
  //  for (int i=0; i<nw; i++){
  //    tmpcount[i]=(_M->get_occ())->get((**_Wp).lw(i),(**_Wp)[i]);
  //  }
  //}
  tmpcount=new int[1];
  tmpcount[0]=(int)nobs;

  /* old stuff */
  //gsl_sort_int_index(SortIndex,tmpcount, 1, nseq) ; // no use

  //Si la longueur ne correspond pas, on realloue ProbCur et ProbPrev
  /* old stuff */
  //if (_lseq != arrseq -> LengthSeqMax())
  if (false)
    {
      for(int i=0; i<=_lseq; i++)
	{
	  delete [] _ProbCur[i] ;
	  delete [] _ProbPrev[i] ;
	}
      delete [] _ProbCur ; delete [] _ProbPrev ;
      /* to replace */
      //_lseq = arrseq -> LengthSeqMax() ;
      _ProbCur = new double*[_lseq+1];
      _ProbPrev = new double*[_lseq+1];
      for (int i=0; i <= _lseq ; i++)
	{
	  /* very ugly way to alloc these array */
	  _ProbCur[i] = new double[_nw];
	  _ProbPrev[i] = new double[_nw];
	  for (int j=0; j<_nw;j++)
	    _ProbCur[i][j] =  _ProbPrev[i][j] = 0 ;
	}      
    }
 
  //Calcul de PFirst pour commencer
  this -> PFirst(QVW) ;

  //printf("PFirst done\n");

  seqcount = 0 ;
  nseqequal = 0  ;
  int nobscurr ;
  int* tmplength ;
  size_t* tmpsort ;
  //Cas de comptages nuls
//  if (tmpcount[SortIndex[0]]==0) {
//    while (tmpcount[SortIndex[nseqequal]]==0 ) {
//      nseqequal++ ;
//      if (nseqequal == nseq)
//	break ;
//    }
//    //si il y a plusieurs sequences de même longueur
//    if (nseqequal>1) {
//      //dans ce cas on "process" les sequences par ordre de taille
//      tmplength = new int[nseqequal] ;
//      tmpsort = new size_t[nseqequal] ;
//      for (int i=0 ; i < nseqequal ; i++) {
//	/* old stuff */
//	//tmplength[i] = arrseq -> LengthSequence(SortIndex[i]) ; // replace by sequence length
//	tmplength[i]=_lseq;
//      }
//      /* old stuff */
//      //gsl_sort_int_index(tmpsort,tmplength, 1, nseqequal) ; // no use
//      prob = 0 ;
//      for (int w=0; w < _nw ; w++)
//	for (int l=1; l<=tmplength[tmpsort[0]]; l++) 
//	  prob += _ProbCur[l][w] ;
//      Probas[SortIndex[tmpsort[0]]] = 1 - (double) prob ;
//      //bizarre, il arrive que certains mots qui n'apparaissent 
//      //pas soient surreprésentés
//      for (int i=1 ; i<nseqequal; i++){
//	for (int l=tmplength[tmpsort[i-1]]+1;	\
//	     l<=tmplength[tmpsort[i]];l++)
//	  for (int w=0; w < _nw ; w++)
//	    prob += _ProbCur[l][w] ;
//	
//	Probas[SortIndex[tmpsort[i]]] = (double) prob ;
//      }
//      seqcount += nseqequal ;
//      delete tmplength ;
//      delete tmpsort ;
//    } else {
//      prob = 0 ;
//      for (int w=0; w < _nw ; w++)
//	/* old stuff */
//	//for (int l=1; l<= arrseq->LengthSequence(SortIndex[0]) ; l++) // replace by seq length
//	for (int l=1; l<= _lseq ; l++) 
//	  prob += _ProbCur[l][w] ;	
//      Probas[SortIndex[0]] = (double) prob ;
//      
//      seqcount++ ;
//    }
//  } // end special case when we have 0 occ

  //On doit faire les Nobs suivant maintenant....
  while (seqcount < nseq)
    {
      nobscurr = tmpcount[SortIndex[seqcount]] ;
      if (_stat<0.0) {
	nobscurr++;
      }
      nseqequal = 0 ;
      //On init les probas 
      this -> GoToNobs(QVW,nobscurr) ;
      //tester la version courte
      while (tmpcount[SortIndex[seqcount+nseqequal]] == nobscurr)
	{
	  nseqequal++ ;
	  if ((seqcount+nseqequal) == nseq)
	    break ;
	}
      if (nseqequal>1)
	{
	  //dans ce cas on "process" les sequences par ordre de taille
	  tmplength = new int[nseqequal] ;
	  tmpsort = new size_t[nseqequal] ;
	  for (int i=0; i<nseqequal; i++) {
	    /* old stuff */
	    //tmplength[i] = arrseq -> LengthSequence(SortIndex[seqcount+i]) ; // seqlenth
	    tmplength[i]=_lseq;
	  }
	  /* old stuff */
	  //gsl_sort_int_index(tmpsort,tmplength, 1, nseqequal) ; // no use
	  prob = 0 ;
	  for (int w=0; w < _nw ; w++)
	    for (int l=1; l<=tmplength[tmpsort[0]]; l++) 
	      prob += _ProbCur[l][w] ;
	  //rajouter le changement de proba ici !!!
	  Probas[SortIndex[seqcount+tmpsort[0]]] = (double) prob ;
	  for (int i=1 ; i<nseqequal; i++)
	    {
	      for (int l=tmplength[tmpsort[i-1]]+1;\
		   l<=tmplength[tmpsort[i]];l++)
		for (int w=0; w < _nw ; w++)
		  prob += _ProbCur[l][w] ;
	      
	      Probas[SortIndex[seqcount+tmpsort[i]]] =  (double) prob ;

	    }
	  seqcount += nseqequal ;
	  delete tmplength ;
	  delete tmpsort ;
	}
      else
	{
	  prob = 0 ;
	  for (int w=0; w < _nw ; w++) {
	    /* old stuff */
	    //for (int l=1; l<= arrseq->LengthSequence(SortIndex[seqcount]) ; l++)  // seq length
	    if (_lseq==_truelseq) {
	      /* regular sum */
	      for (int l=1; l<= _lseq ; l++) 
		prob += _ProbCur[l][w] ;
	    } else {
	      /* tail sum */
	      for (int l=_truelseq+1; l<=_lseq; l++)
		prob += _ProbCur[l][w] ;
	    }
	  }
	  //mettre la proba direct
	  Probas[SortIndex[seqcount]] = (double) prob ;
	  seqcount++ ;
	} 

    }
  double probcompl ;
  double probc ;
  double tmp ;
  double mantisse ;
  long exposant ;
  double stattmp ;
  //rangement des : probas, stats gaussiennes et les nombres d'occurrence
//  if (loglevel >=3)
//    {
//      cerr << "Rangement des probas pour le groupe de mots :" ;
//      for (int i=0; i < _nw; i++) 
//	cerr << (**_Wp).GetWord(i) <<", " ;
//      if (_nw>1)
//	cerr << (**_Wp).GetWord(_nw-1) << endl ;
//    }
  for (int i=0 ; i< nseq; i++)
    {
      Nocc[i] = tmpcount[i] ;
      probcompl = 1 - Probas[i] ;
      probc = Probas[i] ;
      //printf("Probas[0]=%e\t1-Probas[0]=%e\n",probc,probcompl);
      if (_stat>0) {
	Probas[0]=probc;
      } else {
	if (_lseq==_truelseq) {
	  /* use complementary */
	  Probas[0]=-probcompl;
	} else {
	  /* tailsum already computed (just adjust sign) */
	  Probas[0]*=-1.0;
	}
      }
      
//      if (probc<probcompl) {
//	if (prob<0)  nerrors--  ;
//      } else 
//	if (probcompl > 0)
//	  Probas[i] = -probcompl ;
//	else {
//	  Probas[i] = -2000 ;
//	  nerrors-- ;
//	}
      tmp = (log(fabs(Probas[i]))/log(10.0)) ;
      exposant = ((long) tmp)-1 ;
      mantisse = exp((tmp - exposant)*log(10.0)) ;
      /* to replace */
      //inv_fr_nor(mantisse, exposant, &stattmp) ; // no use
      //GaussStats[i] = (Probas[i]>=0) ? stattmp : -stattmp ;
      Nocc[i] = tmpcount[i] ;
//      if (loglevel >=3)
//	{
//	  cerr << Nocc[i] << "\t"  << probc << "\t" << probcompl  
//	       << "\t" << GaussStats[i] << endl ;
//	}
    }
  delete [] SortIndex ;
  delete [] tmpcount ;
  return(nerrors) ;
}


