/* $Id: sspatt.cc 440 2005-07-20 12:34:26Z mhoebeke $ */
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

/*********************************************************************/   
/*  								     */
/*  aspatt main program                                              */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * Ugly fix (part I) to avoid macro redefinition warning messages.
 **/
#undef PACKAGE
#undef VERSION
#include "../tclap-1.0.5/config/config.h"
/**
 * Ugly fix (part II) to avoid macro redefinition warning messages.
 **/
#undef PACKAGE
#undef VERSION
#include <tclap/CmdLine.h>
using namespace TCLAP;

#include "amarkov.h"
#include "timer.h"

int
main(int argc, char *argv[])
{
  // command line parsing
  // third argument is the delimitor
  Arg::setAllowDash(true);
  CmdLine cmdLine("Command description message",' ',VERSION);
  ValueArg<string> MarkovArg("M","Markov","input Markov file (dfa program output)",true,"","filename");
  ValueArg<string> StartingArg("S","Starting-distribution","input starting distribution file",false,"","filename");
  ValueArg<long> lengthArg("l","length","length of the sequence",true,0,"integer"); 
  ValueArg<long> nobsArg("n","nobs","observed number of occurrences of the sequence",true,0,"integer"); 
  cmdLine.add(MarkovArg);
  cmdLine.add(StartingArg);
  cmdLine.add(lengthArg);
  cmdLine.add(nobsArg);
  cmdLine.parse(argc,argv);

  timer T;

  long length=lengthArg.getValue();
  long nobs=nobsArg.getValue();

  // build Markov
  T.start();
  amarkov M(StartingArg.getValue().c_str(),MarkovArg.getValue().c_str());
  //printf("%.2f s to read file and compute stationary distribution\n",T.elapsed_time());

  // get binomial stat
  double nexp=M.nexp(length);
  printf("nexp=%e\n",nexp);
  double bstat=M.bstat(length,nobs,nexp);
  printf("%.2f s to compute binomial statistic=%+f\n",T.elapsed_time(),bstat);

  // get binomial stat
  double cpstat=M.cpstat(length,nobs,nexp);
  printf("%.2f s to compute compound Poisson statistic=%+f\n",T.elapsed_time(),cpstat);

  //printf("length=%i\tnobs=%i\texp=%e\n",length,nobs,nexp);

//  // test part
//  printf("*********\n");
//  printf("test part\n");
//  printf("*********\n");
//
//  // double slow_dcpoi(int n,double lambda,double theta, int alpha,double *p);
//  double v,t;
//  int alpha=9;
//  double p[9]={
//    0.9547189295434349576296995110169518738985061645507812500000000000000000000000000000000000000000000000,
//    0.0429952602806177794358966082199913216754794120788574218750000000000000000000000000000000000000000000,
//    0.0021674800398144671020883578194116125814616680145263671875000000000000000000000000000000000000000000,
//    0.0001121715468990598733303384881843101084086811169981956481933593750000000000000000000000000000000000,
//    0.0000058377001883545242332603127100032480711888638325035572052001953125000000000000000000000000000000,
//    0.0000003041654277310425166620734722944252581555701908655464649200439453125000000000000000000000000000,
//    0.0000000158519985283566015286440544721641154879421264922711998224258422851562500000000000000000000000,
//    0.0000000008261906961546820055735171469887823114675029501086100935935974121093750000000000000000000000,
//    0.0000000000430607093543609368039037841289790736754028444011055398732423782348632812500000000000000000
//  };
//  //double lambda=0.0001;
//  //double lambda=3909.7507;
//  double lambda=3731.7758941451593273086473345756530761718750000000000000000000000000000000000000000000000000000000000000; 
//  //double lambda=10000;
//  double theta=0.0521197219999953512137302880091738188639283180236816406250000000000000000000000000000000000000000000;
//  //for (int n=150; n<=270; n++) {
//  for (int n=100; n<=1200; n++) {
//    printf("%i\t",n);
//
////    T.restart();
////    //v=spatt::complete_log_dcpoi(n,lambda,theta,alpha,p)/log(10.0);
////    //v=log(spatt::barbour_qcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    //v=log(spatt::barbour_qcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    v=log(spatt::fast_dcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    T.stop(); 
////    t=T.elapsed_time();
////    printf("%.10f\t%.2f\t",v,t);
//
//    T.restart();
//    //v=log(spatt::complete_qcpoi(n,lambda,theta,alpha,p))/log(10.0);
//    //v=log(spatt::complete_qcpoi(n,lambda,theta,alpha,p))/log(10.0);
//    v=spatt::fast_log_dcpoi(n,lambda,theta,alpha,p)/log(10.0);
//    T.stop();
//    t=T.elapsed_time();
//    printf("%.10f\t%.2f\t",v,t);
//
////    T.restart();
////    v=log(spatt::fast_pcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    //v=log(spatt::fast_dcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    T.stop();
////    t=T.elapsed_time();
////    printf("%.10f\t%.2f\t",v,t);
//
////    T.restart();
////    v=log(spatt::fast_qcpoi(n,lambda,theta,alpha,p))/log(10.0);
////    T.stop();
////    t=T.elapsed_time();
////    printf("%.10f\t%.2f\t",v,t);
//
//    printf("\n");
//  }

}
