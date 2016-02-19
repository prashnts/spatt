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
#include <cmath>
#include <cstdio>
using namespace std;


/**
 * Ugly fix (part I) to avoid macro redefinition warning messages.
 **/
#undef PACKAGE
#undef VERSION
/**
 * Ugly fix (part II) to avoid macro redefinition warning messages.
 **/
#undef PACKAGE
#undef VERSION
#include <tclap/CmdLine.h>
using namespace TCLAP;

#include <alphabet.h>
#include <sequence.h>
#include <markov.h>
#include <vmarkov.h>
#include <pmc.h>
#include <dfa.h>
#include <nfa.h>
#include <xstat.h>
#include <xwaiting.h>
#include <gstat.h>

using namespace spatt;

int main (int argc,char **argv) {
  
  // global variables
  bool verbose=false;
  bool debug=false;
  string alphabet_label="acgt";
  string pattern_label;
  string nfa_dotfile;
  string dfa_dotfile;
  string sequence_file;
  unsigned long sequence_length;
  short m=-1;
  bool renewal=false;
  string markov_file;
  int rep=UNDER;
  bool repartition=false;
  bool gaussian=false;
  bool sensitivity=false;
  bool exact=true;
  bool have_nobs=false;
  long nobs=-1;
  bool stationary=false;
  bool interrupt=true;
  int offset=0;
  bool presence=false;
  string format_float;

  // setting command line
  CmdLine cmdline("spatt",' ',"2.0");
  SwitchArg verboseArg("V","Verbose","verbose output.",false);
  cmdline.add(verboseArg);
  SwitchArg debugArg("d","debug","debug output.",false);
  cmdline.add(debugArg);
  SwitchArg renewalArg("r","renewal","consider renewal occurrences instead of overlapping ones.",false);
  cmdline.add(renewalArg);
  ValueArg<string> alphabetArg("a","alphabet","define the alphabet to use, ex: \"acgt\" for standard DNA alphabet",false,"acgt","descriptor");
  cmdline.add(alphabetArg);
  ValueArg<string> patternArg("p","pattern","define a pattern, ex: \"ac.(0-2)[agc]acg{a}|gattaca\" means the pattern ac, any char repeated 0 to 2 times, a or g or c, acg, any but a OR the pattern gattaca",true,"","descriptor");
  cmdline.add(patternArg);
  ValueArg<string> nfaArg("","nfa","export the built nfa to a dot file",false,"","dot file");
  cmdline.add(nfaArg);
  ValueArg<string> dfaArg("","dfa","export the built dfa to a dot file",false,"","dot file");
  cmdline.add(dfaArg);
  ValueArg<string> sequenceArg("S","Sequence","provide a sequence in FASTA format",false,"","sequence file");
  cmdline.add(sequenceArg);
  ValueArg<long> seqlenArg("","seqlen","sequence length",false,-1,"integer");
  cmdline.add(seqlenArg);
  ValueArg<short> orderArg("m","","Markov order (m<0 means uniform order 0 model)",false,-1,"integer");
  cmdline.add(orderArg);
  ValueArg<string> markovArg("M","Markov","provide Markov parameters in a file",false,"","Markov file");
  cmdline.add(markovArg);
  SwitchArg overArg("","over","consider over-representation",false);
  cmdline.add(overArg);
  SwitchArg underArg("","under","consider under-representation (default)",false);
  cmdline.add(underArg);
  SwitchArg repartitionArg("","repartition","consider the repartition of pattern occurrences",false);
  cmdline.add(repartitionArg);
  SwitchArg gaussianArg("","gaussian","perform Gaussian computations",false);
  cmdline.add(gaussianArg);
  SwitchArg sensitivityArg("","sensitivity","compute sensitivity to parameter estimation",false);
  cmdline.add(sensitivityArg);
  ValueArg<long> nobsArg("","nobs","number of observed occurrences",false,-1,"integer");
  cmdline.add(nobsArg);
  SwitchArg stationaryArg("","stationary","force spatt to assume that all Markov chains are stationary",false);
  cmdline.add(stationaryArg);
  SwitchArg ignore_interruptArg("","ignore-interrupt","force spatt to ignore sequence interruptions",false);
  cmdline.add(ignore_interruptArg);
  ValueArg<int> offsetArg("","offset","number of position discarded in each sequence when the ignore-interrupt option is activated",false,0,"integer");
  cmdline.add(offsetArg);
  SwitchArg presenceArg("","presence","check for presence/absence of the pattern (rather than on number of occurrences)",false);
  cmdline.add(presenceArg);
  ValueArg<string> formatArg("f","format","define the output format for floats, ex: \"%.4f\"",false,"%e","format");
  cmdline.add(formatArg);

  // parsing command line
  cmdline.parse(argc,argv);
  if (verboseArg.isSet())
    verbose=true;
  if (debugArg.isSet())
    debug=true;
  if (alphabetArg.isSet())
    alphabet_label=alphabetArg.getValue();
  if (patternArg.isSet())
    pattern_label=patternArg.getValue();
  if (nfaArg.isSet())
    nfa_dotfile=nfaArg.getValue();
  if (dfaArg.isSet())
    dfa_dotfile=dfaArg.getValue();
  if (sequenceArg.isSet())
    sequence_file=sequenceArg.getValue();
  if (seqlenArg.isSet())
    sequence_length=seqlenArg.getValue();
  if (orderArg.isSet()) {
    m=orderArg.getValue();
    if (m<0) {
      m=0;
    }
  }
  if (renewalArg.isSet())
    renewal=true;
  if (markovArg.isSet())
    markov_file=markovArg.getValue();
  if (overArg.isSet())
    rep=OVER;
  if (underArg.isSet())
    rep=UNDER;
  if (repartitionArg.isSet())
    repartition=true;
  if (gaussianArg.isSet()) {
    gaussian=true;
    exact=false;
  }
  if (sensitivityArg.isSet())
    sensitivity=true;
  if (nobsArg.isSet())
    nobs=nobsArg.getValue();
  if (offsetArg.isSet())
    offset=offsetArg.getValue();
  if (offset < 0) {
    fprintf(stderr,"negative offset ignored\n");
    offset=0;
  }
  if (stationaryArg.isSet())
    stationary=true;
  if (ignore_interruptArg.isSet())
    interrupt=false;
  if (presenceArg.isSet())
    presence=true;
  format_float=formatArg.getValue();
  
  // build the alphabet
  alphabet my_alphabet(alphabet_label);

  // build the nfa
  //nfa my_nfa(alphabet_label,pattern_label,verbose,debug);
  nfa my_nfa(my_alphabet,pattern_label,verbose,debug);

  // output it
  if (!nfa_dotfile.empty())
    my_nfa.dot(nfa_dotfile);
 
  // build the dfa by subset construction
  dfa my_dfa(my_nfa,verbose,debug);

  // minimize it and rebuild it
  my_dfa.minimize(verbose,debug);
  my_dfa.rebuild(verbose);
  
  // output it
  if (!dfa_dotfile.empty())
    my_dfa.dot(dfa_dotfile);

  if (renewal)
    my_dfa.renewal(verbose);

  if (m>0) {
    my_dfa.initialize_ambiguity_vectors(verbose,debug);
    my_dfa.remove_ambiguity(m,verbose,debug);
    //my_dfa.starts(verbose);
  }

  sequence *seq=NULL;
  if (!sequence_file.empty()) {
    seq = new sequence(my_alphabet,sequence_file);
//    printf("********* testing sequence\n");
//    vector<string> files;
//    files.push_back(sequence_file);
//    files.push_back("/home/gnuel/biodata/ecoli");
//    files.push_back("idiot");
//    files.push_back("/home/gnuel/ecoli");
//    files.push_back(sequence_file);
//    sequence seq(my_alphabet,files);
//    //sequence seq(my_alphabet,sequence_file);
//    int c;
//    while ((c=seq.next())!=SEQEMPTY){
////      if (c>=0)
////	printf("%i",c);
////      else
////	printf("$\n");
//    }
//    printf("nvalidchar=%i\n",seq.nvalidchar());
//    printf("nvalidseq=%i\n",seq.nvalidseq());
//    printf("********* end sequence testing\n");
//    my_dfa.locate_occ(sequence_file,verbose);
  }
  // output it
  if (!dfa_dotfile.empty())
    my_dfa.dot(dfa_dotfile);

  markov *my_markov=NULL;
  if (m>=0) {
    if (!markov_file.empty())
      my_markov = new markov(my_dfa._alphabet_size,my_dfa._m,markov_file,stationary,verbose);
    else if (m==0)
      my_markov = new markov(my_dfa._alphabet_size,stationary,verbose);
    else {
      my_markov = new markov(my_dfa._alphabet_size,my_dfa._m,stationary,verbose);
    }

    if (sensitivity) {
      printf("sensitivity to Markov parameter estimation not fully functional: only mean and covariance of (N_m,N_{m+1}) are computed\n");
      vmarkov my_vmarkov(*my_markov);
      my_vmarkov.compute_mu();
      if (verbose)
	my_vmarkov.print_mu();
      my_vmarkov.compute_Sigma(10000);
      //my_vmarkov.compute_Sigma(580076); // mgenitallium
      //my_vmarkov.compute_Sigma(4639675); // ecoli_K12
      //my_vmarkov.compute_Sigma(4214630); // bsubtilis
      my_vmarkov.dump_Sigma("tmp.vmarkov");
      //if (verbose)
      //my_vmarkov.print_Sigma();
    }

    pmc my_pmc(my_dfa,my_markov->_param,verbose,debug);
    
    my_pmc.sci_export("pmc.sci");
    my_pmc.indexed_export("pmc.markov");

    if (seq || sequence_length) {
        if (!repartition) {
            if (exact) {
                xstat *stat;
                if (seq)
                    stat = new xstat(my_dfa,my_pmc,*seq,rep,nobs,presence,verbose);
                else if (sequence_length)
                    stat = new xstat(my_dfa,my_pmc,sequence_length,rep,nobs,presence,verbose);
                if (stationary)
                    stat->compute(my_markov->_mu0,interrupt,offset,verbose);
                else
                    stat->compute(interrupt,offset,verbose);
                stat->print_distribution(format_float);
                stat->print(pattern_label.c_str(),format_float);
                delete stat;
                //printf("N(%s)=%i\tp-value=%e\n",pattern_label.c_str(),stat.nobs(),stat.pvalue());
                //if (rep==OVER)
                //  printf("P( N(%s) >= %i ) = %e\n",pattern_label.c_str(),stat.nobs(),stat.pvalue());
                //else
                //  printf("P( N(%s) <= %i ) = %e\n",pattern_label.c_str(),stat.nobs(),stat.pvalue());
                //if (rep==OVER)
                //  printf("-log10 P( N(%s) >= %i ) = %f\n",pattern_label.c_str(),stat.nobs(),-log(stat.pvalue())/log(10.0));
                //else
                //  printf("+log10 P( N(%s) <= %i ) = %f\n",pattern_label.c_str(),stat.nobs(),log(stat.pvalue())/log(10.0));
            } // end if (exact)
            if (gaussian) {
                gstat *stat;
                if (seq)
                    stat = new gstat(my_dfa,my_pmc,*seq,rep,nobs,presence,verbose);
                else if (sequence_length)
                    stat = new gstat(my_dfa,my_pmc,sequence_length,rep,nobs,presence,verbose);
                if (stationary)
                    stat->compute(my_markov->_mu0,interrupt,offset,verbose);
                else
                    stat->compute(interrupt,offset,verbose);
                stat->print(pattern_label.c_str(),format_float);
                delete stat;
                //printf("N(%s)=%i\tmean=%f\tsd=%f\tz-score=%f\t",pattern_label.c_str(),stat.nobs(),stat.mean(),stat.sd(),stat.zscore());
                //if (rep==OVER)
                //  printf("P( N(%s) >= %i ) = %.2e\n",pattern_label.c_str(),stat.nobs(),stat.pvalue());
                //else
                //  printf("P( N(%s) <= %i ) = %.2e\n",pattern_label.c_str(),stat.nobs(),stat.pvalue());
                //printf("p-value=%e\n",stat.pvalue());
            } // end if gaussian
        } else {
            xwaiting waiting(my_dfa,my_pmc,*seq,rep,verbose,debug);
            waiting.print();
        }
    }

    if (my_markov)
      delete my_markov;
    if (seq)
      delete seq;
  
  } // end if m>=0

  return EXIT_SUCCESS;
};
