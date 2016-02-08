/* $Id: spattparameterparser.cc 547 2005-11-22 10:30:51Z gnuel $ */
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
/*  Class parameter allows reads the command line and returns the    */
/*  the parameters.                                                  */
/*  								     */
/*  This class uses argtable2                 			     */
/*  								     */
/*  Author: Grégory Nuel					     */
/*  								     */
/*********************************************************************/

#include "spattparameterparser.h"

#include <cstring>
#include <iostream>

using namespace std;
using namespace TCLAP;

namespace spatt {

  spattparameterparser::spattparameterparser(const string &progname,
					     const string &version)
    
  {
    
    Arg::setAllowDash(true);

    _cmdline=new CmdLine(progname,' ',version);
    
    _quietArg=new SwitchArg("q","quiet","quiet output (same effect as --debug-level -1).",false);
  _cmdline->add(_quietArg);


  _verboseArg=new SwitchArg("V","Verbose","verbose output (same effect as --debug-level 1).",false);
  _cmdline->add(_verboseArg);

  _debugArg=new SwitchArg("d","debug","debug output (same effect as --debug-level 2).",false);
  _cmdline->add(_debugArg);


  _debuglevelArg=new ValueArg<int>("","debug-level","set the verbosity level: -1 (quiet), 0 (normal), 1 (verbose), 2 and more (debug)",false,0,"level");
  _cmdline->add(_debuglevelArg);

  _alphabetArg=new ValueArg<string>("a","alphabet","define the alphabet to use, ex: \"acgt:tgca\" for standard DNA alphabet",false,"acgt:tgca","descriptor");
  _cmdline->add(_alphabetArg);

  _sequenceArg=new UnlabeledValueArg<string>("sequence","FASTA sequence","","filename");
  _cmdline->add(_sequenceArg);

  _modelArg=new ValueArg<string>("M","Markov-file","destination file for the Markov model parameters",false,"","filename");
  _cmdline->add(_modelArg);

  _stationaryArg=new ValueArg<string>("S","Stationary-file","destination file for the stationary distribution",false,"","filename");
  _cmdline->add(_stationaryArg);

  _markovArg=new ValueArg<string>("U","Use-markov-file","file containing the Markov model parameters",false,"","filename");
  _cmdline->add(_markovArg);

  _countArg=new ValueArg<string>("C","Count-file","file containing word counts",false,"","filename");
  _cmdline->add(_countArg);

  _patternFileArg=new ValueArg<string>("P","Pattern-file","file containing patterns",false,"","filename");
  _cmdline->add(_patternFileArg);

  _lengthArg=new ValueArg<int>("l","length","length of the longest word counted",false,0,"max");
  _cmdline->add(_lengthArg);

  _patternListArg=new MultiArg<string>("p","pattern","define a pattern (can be used multiple times)",false,"descriptor");
  _cmdline->add(_patternListArg);

  _strandsArg=new SwitchArg("b","both-strands","count patterns on both strands (useless if no complementary alphabet defined).",false);
  _cmdline->add(_strandsArg);

  _orderArg=new ValueArg<int>("m","markov","order of the Markov model",false,-2,"order");
  _cmdline->add(_orderArg);
  
  _allwordsArg=new SwitchArg("","all-words","process all words of the given length",false);
  _cmdline->add(_allwordsArg);

  _normalizeArg=new SwitchArg("n","normalize","statistics are normalized by the length of the sequence",false);
  _cmdline->add(_normalizeArg);

  _maxpvalueArg=new ValueArg<double>("","max-pvalue","only statistics corresponding to a pvalue lower then the threshold are produced",false,1.0,"threshold");
  _cmdline->add(_maxpvalueArg);
  
  _nobsArg=new ValueArg<int>("c","nobs","number of observed pattern occurrences",false,0,"count");
  _cmdline->add(_nobsArg);

}

  void
  spattparameterparser::parse(int argc, char **argv, spattparameters *params) {

  _cmdline->parse(argc,argv);

  if (_quietArg->getValue() == true)
    params->_debug_level=-1;

  if (_verboseArg->getValue() == true)
      params->_debug_level=1;

  if (_debugArg->getValue() == true)
    params->_debug_level=2;

  if (_debuglevelArg->isSet())
    params->_debug_level=_debuglevelArg->getValue();

  vector<string> patterns=_patternListArg->getValue();


  params->_all_words=_allwordsArg->getValue();
  if (params->_all_words == true && patterns.size()>0) 
      cerr << "Warning ! --all-words specified, other patterns ignored" << endl;


  params->_normalize=_normalizeArg->getValue();

  if (_alphabetArg->isSet()) {
    params->_alphabet_label=_alphabetArg->getValue();
    if (params->_debug_level >= 1) {
      cout << "alphabet label =" << params->_alphabet_label << endl;
    }
  }

  params->_sequence_filename=_sequenceArg->getValue();
  params->_use_sequence_file=1;
  if (params->_debug_level >= 1) 
    cout << "sequence filename =" << params->_sequence_filename << endl;


  if (_countArg->isSet()) {
    params->_count_filename=_countArg->getValue();
    params->_use_count_file=true;
    if (params->_debug_level >= 1) {
      cout << "count filename =" << params->_count_filename << endl;    
    }
  }

  if (_markovArg->isSet()) {
    params->_markov_filename=_markovArg->getValue();
    params->_use_markov_file=true;
    if (params->_debug_level >= 1) 
      cout << "Markov filename =" << params->_markov_filename << endl;
  }


  if (_modelArg->isSet()) {
    params->_model_filename=strdup(_modelArg->getValue().c_str());
    params->_use_model_file=true;
    if (params->_debug_level >= 1)
      cout << "Model filename =" << params->_model_filename << endl;
  }

  if (_stationaryArg->isSet()) {
    params->_stationary_filename=strdup(_stationaryArg->getValue().c_str());
    params->_use_stationary_file=1;
    if (params->_debug_level >= 1)
      cout << "Stationary filename =" << params->_stationary_filename << endl;
  }

  if (_patternFileArg->isSet()) {
    params->_pattern_filename=strdup(_patternFileArg->getValue().c_str());
    params->_use_pattern_file=1;
    if (params->_debug_level >= 1)
      cout << "Pattern filename =" << params->_pattern_filename << endl;
  }

  if (_lengthArg->isSet()) {
    params->_length=_lengthArg->getValue();
    if (params->_length > MAX_LENGTH) {
      if (params->_debug_level >= 0)
	cerr << "Warning ! length > " << MAX_LENGTH << " IS set to " << MAX_LENGTH << endl;
      params->_length=MAX_LENGTH;
    }
    if (params->_debug_level >= 1)
      cout << "count length  =" << params->_length << endl;

  }

  for (int i=0;i<patterns.size();i++) {
    if (params->_debug_level >= 1)
      cout << "pattern (" << i << ") = " << patterns[i] << endl;
    params->_pattern_label_list.push_back(patterns[i]);
  }


  if (_strandsArg->isSet()) {
    params->_both_strands=1;
    if (params->_debug_level >= 1)
      cout << "count on both strands enabled" << endl;
  }

  if (_orderArg->isSet()) {
    params->_markov_order=_orderArg->getValue();
    if (params->_markov_order < -1) {
      cerr << "Warning ! Markov order " << params->_markov_order << " not valid and ignored " << endl;
      params->_markov_order=-2;
    } else {
      if (params->_debug_level >= 1)
	cout << "Markov order " << params->_markov_order << " specified" << endl;
    }
  }

  if (params->_markov_order > -2 && params->_length < 1) {
    cerr << "length is set to 1" << endl;
    params->_length=1;
  }

  if (_maxpvalueArg->isSet()) {
    params->_max_pvalue=_maxpvalueArg->getValue();
    if (params->_max_pvalue <= 0.0 ) {
      cout << "pvalue threshold lower than 0.0 specified. No computation to do." << endl;
      exit(EXIT_SUCCESS);
    }
  }
  
  if (_nobsArg->isSet()) {
    params->_nobs=_nobsArg->getValue();
    if (params->_nobs < 0) {
      cerr << "Warning, count must be non-negative in -nobs option" << endl;
      cerr << "count has been set to 0" << endl;
    }
    if (params->_debug_level >= 1)
      cout << "nobs = " << params->_nobs << " specified " << endl;      
  }
  
}

  spattparameterparser::~spattparameterparser() {
    delete _quietArg;
    delete _verboseArg;
    delete _debugArg;
    delete _debuglevelArg;
    delete _patternListArg;
    delete _allwordsArg;
    delete _normalizeArg;
    delete _alphabetArg;
    delete _sequenceArg;
    delete _countArg;
    delete _markovArg;
    delete _modelArg;
    delete _stationaryArg;
    delete _patternFileArg;
    delete _lengthArg;
    delete _strandsArg;
    delete _orderArg;
    delete _cmdline;
    delete _nobsArg;
    delete _maxpvalueArg;
  }

};
