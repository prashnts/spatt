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

#include "alphabet.h"
#include "sequence.h"
#include "markov.h"
#include "pmc.h"
#include "dfa.h"
#include "nfa.h"
#include "xstat.h"
#include "spattc.h"

double spatt_exact(char* alphabet_chars, char* pattern_chars, double* transitions, short markov_order, bool stationary, long nobs, unsigned long sequence_length, char test) {

    // Global variables
    bool verbose=false;
    bool debug=false;
    bool presence=false;
    bool interrupt=true;
    int offset=0;
    int rep;

    if (test=='o')
        rep = OVER;
    else
        rep = UNDER;

    // Cast C -> C++
    std::string alphabet_string(alphabet_chars);
    std::string pattern_string(pattern_chars);

    // Build the alphabet
    spatt::alphabet my_alphabet(alphabet_string);

    // Build the nfa
    spatt::nfa my_nfa(my_alphabet, pattern_string, verbose, debug);

    // Build the dfa by subset construction
    spatt::dfa my_dfa(my_nfa, verbose, debug);

    // Minimize it and rebuild it
    my_dfa.minimize(verbose, debug);
    my_dfa.rebuild(verbose);

    if (markov_order>0) {
        my_dfa.initialize_ambiguity_vectors(verbose, debug);
        my_dfa.remove_ambiguity(markov_order, verbose, debug);
    }

    spatt::markov *my_markov=NULL;
    if (markov_order>=0) {
        if (markov_order==0)
            my_markov = new spatt::markov(my_dfa._alphabet_size, stationary, verbose);
        else
            my_markov = new spatt::markov(my_dfa._alphabet_size, my_dfa._m, transitions, stationary, verbose);

        spatt::pmc my_pmc(my_dfa, my_markov->_param, verbose, debug);

        spatt::xstat stat(my_dfa, my_pmc, sequence_length, rep, nobs, presence, verbose);
        if (stationary)
            stat.compute(my_markov->_mu0, interrupt, offset, verbose);
        else
            stat.compute(interrupt, offset, verbose);
        delete my_markov;
        return stat.pvalue();
    }
    if (my_markov)
        delete my_markov;
}
