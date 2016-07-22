/*################################################################################
# *
# *
# * Copyright Jean-Marc Aury / CEA <jmaury@genoscope.cns.fr>
# *
# * This software is a computer program whose purpose is to do stuff
# *
# * This software is governed by the CeCILL license under French law and
# * abiding by the rules of distribution of free software.  You can  use,
# * modify and/ or redistribute the software under the terms of the CeCILL
# * license as circulated by CEA, CNRS and INRIA at the following URL
# * "http://www.cecill.info".
# *
# * As a counterpart to the access to the source code and  rights to copy,
# * modify and redistribute granted by the license, users are provided only
# * with a limited warranty  and the software's author,  the holder of the
# * economic rights,  and the successive licensors  have only  limited
# * liability.
# *
# * In this respect, the user's attention is drawn to the risks associated
# * with loading,  using,  modifying and/or developing or reproducing the
# * software by the user in light of its specific status of free software,
# * that may mean  that it is complicated to manipulate,  and  that  also
# * therefore means  that it is reserved for developers  and  experienced
# * professionals having in-depth computer knowledge. Users are therefore
# * encouraged to load and test the software's suitability as regards their
# * requirements in conditions enabling the security of their systems
# * and/or data to be ensured and,  more generally, to use and operate it
# * in the same conditions as regards security.
# *
# * The fact that you are presently reading this means that you have had
# * knowledge of the CeCILL license and that you accept its terms.
# *
################################################################################
*/

#include <iostream>
#include <istream>
#include <fstream>
#include <cstdio>

using namespace std;

#define MAX_QUAL_VAL 100

unsigned long int nb_reads=0;
unsigned long long int nb_nt=0;

void usage();
int next(istream*, string&, string&);

int main(int argc, char* argv[])
{
  char *fastaFilename = NULL;
  char c;
  string seq, qual;
  int qual_distrib = 0;
  unsigned long long int qual_val[MAX_QUAL_VAL];
  for(int i=0; i<MAX_QUAL_VAL; i++) { qual_val[i]=0; }

  // Invokes member function `int operator ()(void);'
  while ((c = getopt(argc, argv, "f:hq")) != -1) {
    switch (c) {
    case 'f':
      fastaFilename = optarg;
      break;
    case 'h':
      usage();
      break;
    case 'q':
      qual_distrib=1;
      break;
    default :
      abort();
    }
  }
  istream *_fstrm;
  _fstrm = &cin;
  if(fastaFilename) _fstrm = new ifstream(fastaFilename);

  if (_fstrm==NULL) {
    fprintf(stderr, "failed to open input file '%s'\n", fastaFilename);
    exit(1);
  }

  while( next(_fstrm, seq, qual) != 0 ) {
    nb_reads++;
    nb_nt+=seq.size();
    if(qual_distrib) { for(int j=0; j<qual.size(); j++) { qual_val[qual[j]-64]++; }}
    seq.erase();
    qual.erase();
  }
  fprintf(stdout, "Number of sequences : %15u\n", nb_reads);
  fprintf(stdout, "Number of bases     : %15lld\n", nb_nt);
  if(qual_distrib) {
    cout << endl << "Quality values distribution : " << endl;
    for(int i=0; i<MAX_QUAL_VAL; i++) {
      if(qual_val[i]>0) { fprintf(stdout, "%2d  : %15lld \n", i, qual_val[i]); }
    }
  }
}


void usage() {
  cerr << "usage: fastx_stats [[-h] [-f INFILE] [-q]]" << endl;
  cerr << endl;
  cerr << "   [-h]         = This helpful help screen." << endl;
  cerr << "   [-f INFILE]  = FASTA/Q input file. default is STDIN." << endl;
  cerr << "   [-q]         = Display quality values distribution." << endl;
  exit(0);
}

int next(istream* _fstrm, string& seq, string& qual) {
  char c = _fstrm->get();
  if(_fstrm->eof()) return 0;
  if('>' == c) {
    // fasta format
    _fstrm->ignore(1000, '\n');
    c=' ';
    string seq_buf;
    while(1) {
      c = _fstrm->get();
      _fstrm->unget();
      if('>' == c || _fstrm->eof()) { return 1; }
      *_fstrm>>seq_buf;
      char newline = _fstrm->get();
      seq.append(seq_buf);
      seq_buf.erase();
    }
    cerr << "seq= " << seq << " qual= " << qual << endl;
  } else if('@' == c) {
    // fastq format
    _fstrm->ignore(1000, '\n');
    *_fstrm >> seq;
    _fstrm->ignore(1000, '\n');
    _fstrm->ignore(1000, '\n');
    *_fstrm >> qual;
    _fstrm->ignore(1000, '\n');
    if(seq.size() != qual.size()) {
      cerr<<"fatal error: fq format, sequence length not equal to quality length\n";
      cerr<<seq<< " --> " << qual << endl;
      exit(1);
    }
  } else {
    cerr<<"fatal error: unrecognizable format of input file.\n";
    exit(1);
  }
  return 1;
}
