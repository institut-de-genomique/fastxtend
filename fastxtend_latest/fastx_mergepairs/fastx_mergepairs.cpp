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
/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "sequence_alignment.h"
#include <config.h>
#include "fastx.h"
#include "fastx_args.h"

using namespace std;

const int MAX_READ_SIZE = 2000;

const char* usage=
"usage: fastx_mergepairs [-h] [-l N] [-m N] [-i N] [-s] [-M] [-a INFILE1] [-b INFILE2] [-o OUTFILE]\n" \
"Developped at Genoscope using the " PACKAGE_STRING "\n" \
"\n" \
"   [-h]              = This helpful help screen.\n" \
"   [-l N]            = Fragment size of read 2 used for detecting an overlap, default is 40\n" \
"   [-m N]            = Maximal number of mismatches of the alignment of read 1 and subpart of read2, default is 4\n" \
"   [-i N]            = Minimal identity percent of the alignment of read 1 and subpart of read2, default is 90\n" \
"   [-L N]            = Minimal size of the alignment, default is 15\n" \
"   [-s]              = Silent mode\n" \
"   [-a INFILE1]      = FASTA/Q input file\n" \
"   [-b INFILE2]      = FASTA/Q input file\n" \
"   [-o OUTFILE]      = FASTA/Q output file\n" \
"   [-u OUTFILE1]     = FASTA/Q output file of unpaired reads (e.g. merged reads).\n" \
"   [-p OUTFILE1]     = FASTA/Q output file of paired reads (e.g. non-merged reads).\n" \
"   [-q OUTFILE1a]    = FASTA/Q output file of read1 from paired reads (e.g. non-merged reads).\n" \
"   [-r OUTFILE1b]    = FASTA/Q output file of read2 from paired reads (e.g. non-merged reads).\n" \
"                       Choose between [-o] and [-u],[-p] and [-u],[-q],[-r] arguments.\n" \
"   [-x FILE]         = Print distribution of overlap size between read1 and read2\n" \
"   [-M]              = Only print merged pairs, default is no.\n" \
"\n";

void error(string msg);

int debug = 0;
char *file1= NULL;
char *file2= NULL;
char *output= NULL;
char *distrib_size=NULL;
char *paired_output= NULL;
char *paired_read1_output= NULL;
char *paired_read2_output= NULL;
char *unpaired_output= NULL;
int ali_size = 40;
unsigned int min_len_ali = 15;
int max_mismatches = 4, identity_min = 90;
bool silent = false, only_print_merged=false;

int parse_program_args(int __attribute__((unused)) optind, int c, char* optarg) {
  switch(c) {
  case 'a':
    file1 = optarg;
    break ;
    
  case 'b':
    file2 = optarg;
    break ;
    
  case 'o':
    output = optarg;
    break ;
    
  case 'x':
    distrib_size = optarg;
    break ;
    
  case 'u':
    unpaired_output = optarg;
    break;
    
  case 'p':
    paired_output = optarg;
    break;
    
  case 'q':
    paired_read1_output = optarg;
    break;
    
  case 'r':
    paired_read2_output = optarg;
    break;
    
  case 'd':
    debug = atoi(optarg);
    break;
    
  case 'l':
    if (optarg==NULL) error("[-l] parameter requires an argument value");
    ali_size = atoi(optarg);
    if(ali_size < 9 || ali_size > 100) error("[-l] parameter requires an integer between 10 and 100.");
    break;
    
  case 'i':
    if (optarg==NULL) error("[-i] parameter requires an argument value");
    identity_min = atoi(optarg);
    if(identity_min < 9 || identity_min > 100) error("[-i] parameter requires an integer between 10 and 100.");
    break;
    
  case 'L':
    if (optarg==NULL) error("[Lm] parameter requires an argument value");
    min_len_ali = atoi(optarg);
    break;
    
  case 'm':
    if (optarg==NULL) error("[-m] parameter requires an argument value");
    max_mismatches = atoi(optarg);
    if(max_mismatches < 0 || max_mismatches > 20) error("[-m] parameter requires an integer between 0 and 20.");
    break;
    
  case 'M':
    only_print_merged = true;
    break;
    
  case 's':
    silent = true;
    break;
    
  default:
    errx(1,"Unknown argument (%c)", c );
  }
  return 1;
}

int parse_commandline(int argc, char* argv[]) {  
  fastx_parse_cmdline(argc, argv, "a:b:o:x:u:p:q:r:l:s:i:m:d:L:Mv", parse_program_args);
  return 1;
}

int adapter_cutoff_index ( const SequenceAlignmentResults& alignment_results, int max_mismatches, int identity_min, unsigned int min_ali_size ) {
  int alignment_size = alignment_results.neutral_matches +
    alignment_results.matches + 
    alignment_results.mismatches +
    alignment_results.gaps ;

  if (debug>0) {
    cout << "alignment_size= " << alignment_size << " ; min_alignment_size= " << min_ali_size << endl;
    cout << "gaps= " << alignment_results.gaps << endl;
  }

  if(alignment_size < min_ali_size) return -1;
  if(alignment_results.gaps > 0) return -1;

  int mismatches = alignment_results.mismatches ;
  // add the first bases of the second reads that do not align
  mismatches += alignment_results.target_start ;

   // add the last bases that do not align
  int missing_from_query_end = (alignment_results.query_size - alignment_results.query_end-1);
  int missing_from_target_end = (alignment_results.target_size - alignment_results.target_end-1);
  int missing_from_end = std::min(missing_from_query_end, missing_from_target_end);
  mismatches += missing_from_end;
  
  alignment_size = alignment_results.neutral_matches +
    alignment_results.matches + 
    mismatches +
    alignment_results.gaps ;

  if (debug>0) {
    cout << "mismatches= " << mismatches << " ; max_mismatches= " << max_mismatches << endl;
    cout << "identity= " << (alignment_results.matches * 100 / alignment_size ) << " ; min_identity= " << identity_min << endl;
  }

  if (mismatches > max_mismatches || (alignment_results.matches * 100 / alignment_size ) < identity_min ) return -1;
  return alignment_results.query_start;
}



char reverse_complement_base ( const char input ) {
  switch(input)
    {
    case 'N':
      return 'N';
    case 'n':
      return 'n';
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'a':
      return 't';
    case 't':
      return 'a';
    case 'g':
      return 'c';
    case 'c':
      return 'g';
    default:
      error("Invalid nucleotide value in reverse_complement_base()"); 
    }
  return '0'; //should not get here - just to please the compiler
}



void reverse_complement_fastx(FASTX* pFASTX) {
  int i,j ;
  int length = strlen(pFASTX->nucleotides);
  
  char temp_nuc;
  int  temp_qual;
  
  for (i=0;i<length;i++)
    pFASTX->nucleotides[i] = reverse_complement_base ( pFASTX->nucleotides[i] ) ;
  
  i = 0 ;
  j = length - 1 ;
  while ( i < j ) {
    //Swap the nucleotides
    temp_nuc = pFASTX->nucleotides[i] ;
    pFASTX->nucleotides[i] = pFASTX->nucleotides[j] ;
    pFASTX->nucleotides[j] = temp_nuc;
    
    //Swap the quality scores
    if (pFASTX->read_fastq) {
      temp_qual = pFASTX->quality[i];
      pFASTX->quality[i] = pFASTX->quality[j];
      pFASTX->quality[j] = temp_qual ;
    }
    
    //Advance to next position
    i++;
    j--;
  }
}



int main(int argc, char* argv[]) {

  parse_commandline(argc, argv);

  int correct = 0;
  if(only_print_merged && output && !unpaired_output && !paired_output && !paired_read1_output && !paired_read2_output) correct = 1;
  if(!only_print_merged && output && !unpaired_output && !paired_output && !paired_read1_output && !paired_read2_output) correct = 1;
  if(!only_print_merged && !output && unpaired_output && paired_output && !paired_read1_output && !paired_read2_output) correct = 1;
  if(!only_print_merged && !output && unpaired_output && !paired_output && paired_read1_output && paired_read2_output) correct = 1;
  
  if(!correct) error("need to choose bewteen [-o] OR [u],[-p] OR [-u],[-q],[-r] OR [-M],[-o] arguments");
  if(paired_read1_output && paired_read2_output) paired_output = paired_read1_output;
 
  if(file1 == NULL) error("[-a] parameter requires an argument value");
  if(file2 == NULL) error("[-b] parameter requires an argument value");
 
  FASTX *fastx1 = NULL;
  FASTX *fastx2 = NULL;
  if ((fastx1 = (FASTX*) malloc(sizeof(FASTX))) == NULL) {
	  error("error in fastxmergepairs memory allocation");
  }
  if ((fastx2 = (FASTX*) malloc(sizeof(FASTX))) == NULL) {
	  error("error in fastxmergepairs memory allocation");
  }
  
  HalfLocalSequenceAlignment align;
  fastx_init_reader(fastx1, file1, FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset() );
  fastx_init_reader(fastx2, file2, FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset() );
  if(output) { fastx_init_writer(fastx1, output, OUTPUT_SAME_AS_INPUT, 0); }
  if(unpaired_output) fastx_init_writer(fastx1, unpaired_output, OUTPUT_SAME_AS_INPUT, 0); 
  if(paired_output) fastx_init_writer(fastx1, paired_output, OUTPUT_SAME_AS_INPUT, 0);

  int reads_count = 0, count_input = 0;
  ofstream ofs_u, ofs_p1, ofs_p2, pofsDistrib;
  ofstream *pofsU, *pofsP1, *pofsP2;
  if(output) { 
    ofs_u.open(output);
    pofsU = &ofs_u;
    pofsP1 = &ofs_u;
    pofsP2 = &ofs_u;
  }
  else {
    ofs_u.open(unpaired_output);
    ofs_p1.open(paired_output);
    pofsU = &ofs_u;
    pofsP1 = &ofs_p1;
    pofsP2 = &ofs_p1;
    if(paired_read1_output && paired_read2_output) {
      ofs_p2.open(paired_read2_output);
      pofsP2 = &ofs_p2;
    }
  }
  if(distrib_size) { pofsDistrib.open(distrib_size); }

  unsigned long int nb_tot_pairs = 0;
  unsigned long int nb_merged_pairs = 0;
  unsigned long long int sum_len = 0;
  int min_readsize = MAX_READ_SIZE;
  int max_readsize = 0;
  int size_distrib[MAX_READ_SIZE+1];

  for(int i=0; i<=MAX_READ_SIZE; i++) size_distrib[i]=0;

  while ( fastx_read_next_record(fastx1) ) {
    
    fastx_read_next_record(fastx2);
    reverse_complement_fastx(fastx2);
    nb_tot_pairs++;
    //cout << nb_tot_pairs << endl;

    string s1 = string(fastx1->nucleotides);
    string s2 = string(fastx2->nucleotides);
    string merged = "";
    string merged_qual = "";
    string query = s1;
    string target = s2.substr(0, ali_size);

    align.align( query, target ) ;
    
    if (debug>1) align.print_matrix();
    if (debug>0) align.results().print();
    
    //Find the best match with seed from read2
    int i = adapter_cutoff_index ( align.results(), max_mismatches, identity_min, min_len_ali ) ;

    if ( i!=-1) { // && i>0 ) {
      int overlap = s1.size() - i;

      SequenceAlignmentResults align_res = align.results();
      int deb1=i, deb2=align_res.target_start, fin1=align_res.query_end, fin2=align_res.target_end;
      if(fin2 < s2.size()) {
	int offset = min(s1.size()-fin1, s2.size()-fin2);
	fin2 += offset;
	fin1 += offset;
	//if(fin1 > s1.size()) fin1 = s1.size();
      }
      //cout << "name1= " << fastx1->name << " deb1= " << deb1 << " fin1=" << fin1 << " deb2=" << deb2 << " fin2=" << fin2 << endl;

      // Line 1
      //Print sequence name
      (*pofsU) << fastx1->output_sequence_id_prefix << fastx1->name << endl ;


      // Line 2
      // Retrieve merged sequence
      
      // first part of read 1
      //merged.append(s1.substr(0, deb1));
      
      // common part, by choosing best nucleotide at each position, bad if there is indel
      //int j = 0;
      //for(; deb2+j < fin2 && deb2+j < s2.size(); j++) {
      //if(fastx1->quality[deb1+j] > fastx2->quality[deb2+j]) merged.append(s1.substr(deb1+j, 1));
      //else merged.append(s2.substr(deb2+j, 1));
      //}
      
      // last part of read 2 or read 1 if read 2 is included in read 1
      //if(deb2+j<s2.size()) merged.append(s2.substr(deb2+j, s2.size()-(deb2+j)));
      //if(deb1+j<s1.size()) merged.append(s1.substr(deb1+j, s1.size()-(deb1+j)));
      
      // first part of read 1 and common part, choose the read one (the best in common part ?), a good solution would be to follow the alignment
      merged.append(s1.substr(0, fin1));
      // last part of read 2. Don't keep the last part of read 1 if read 2 is included in read 1.
      if(fin2<s2.size()) merged.append(s2.substr(fin2, s2.size()-fin2));
      
      (*pofsU) << merged << endl;
      size_distrib[merged.size()]++;
      sum_len += merged.size();
      nb_merged_pairs++;
      min_readsize = min((int)merged.size(), min_readsize);
      max_readsize = max((int)merged.size(), max_readsize);


      // Line 3
      // Print quality sequence name
      (*pofsU) <<  "+" << fastx1->name2 << endl;


      // Line 4
      // Retrieve merged quality
      
      // first part of read 1
      //for (j = 0; j < deb1; j++) merged_qual.append(1, char(fastx1->quality[j] + get_fastq_ascii_quality_offset()));
      // common part, by choosing best nucleotide at each position
      //for (j = 0; deb2+j < fin2 && deb2+j < s2.size(); j++) {
      //if(fastx1->quality[deb1+j] > fastx2->quality[deb2+j]) 
	  //merged_qual.append(1, char(fastx1->quality[deb1+j] + get_fastq_ascii_quality_offset()));
      //else 
	  //merged_qual.append(1, char(fastx2->quality[deb2+j] + get_fastq_ascii_quality_offset()));
      //}
      // last part of read 2 or read 1 if read 2 is included in read 1
      //if(deb2+j<s2.size())
      //for(int i=deb2+j ; i<s2.size(); i++) 
      //merged_qual.append(1, char(fastx2->quality[i] + get_fastq_ascii_quality_offset()));
      // if(j+i<s1.size()) 
      //for(int i=deb1+j ; i<s1.size(); i++) 
      // merged_qual.append(1, char(fastx1->quality[i] + get_fastq_ascii_quality_offset()));
      
      // first part of read 1 and common part, choose the read one (the best in common part ?), a good solution would be to follow the alignment
      for (int i = 0; i < fin1; i++) merged_qual.append(1, char(fastx1->quality[i] + get_fastq_ascii_quality_offset()));
      // last part of read 2. Don't keep the last part of read 1 if read 2 is included in read 1.
      for (int i=fin2 ; i<s2.size(); i++) 
    	  merged_qual.append(1, char(fastx2->quality[i] + get_fastq_ascii_quality_offset()));
      
      
      (*pofsU) << merged_qual << endl;

      // DEBUG
      if(debug>0) {
	cout << "********* s1.size= " << s1.size() << " s2.size= " << s2.size() << " merged.size= " << merged.size() << endl;
	cout << "********* s1= " << s1 << endl;
	cout << "********* s2= " << s2 << endl;
	cout << "********* merged= " << merged << endl;
	cout << "********* merged.qual= " << merged_qual << endl;
      }


    } else {
      if(!only_print_merged) {
	(*pofsP1) << fastx1->output_sequence_id_prefix << fastx1->name << endl << s1 << endl << "+" << fastx1->name2 << endl;
	for (int j = 0; j < s1.size(); j++) (*pofsP1) << char(fastx1->quality[j] + get_fastq_ascii_quality_offset());
	(*pofsP1) << endl;
	(*pofsP2) << fastx1->output_sequence_id_prefix << fastx2->name << endl << s2 << endl << "+" << fastx2->name2 << endl;
	for (int j = 0; j < s2.size(); j++) (*pofsP2) << char(fastx2->quality[j] + get_fastq_ascii_quality_offset());
	(*pofsP2) << endl;
      }
    }
  }
  
  if(distrib_size) {
    for(int i=0; i<=MAX_READ_SIZE; i++) {
      if(size_distrib[i]!=0) pofsDistrib << i << " " << size_distrib[i] << endl;
    }
  }

  if(!silent) {
    int med = 0, sum = 0;
    for(int i=0; i<=MAX_READ_SIZE && med==0; i++) { if(sum >= nb_merged_pairs/2) { med = i; } else { sum+=size_distrib[i]; }}
    
    cerr << "Number of initial paired-reads : " << nb_tot_pairs << " (i.e. " << nb_tot_pairs*2 << " reads)" << endl;
    float ratio = ((float)nb_merged_pairs / (float)nb_tot_pairs) * 100.0;
    cerr << "Number of merged paired-reads  : " << nb_merged_pairs << " (" << ratio << "%)" << endl;
    if(nb_merged_pairs > 0) {
      cerr << "Merged reads statistics :" << endl;
      cerr << "  Number     " << nb_merged_pairs << endl;
      cerr << "  Med size   " << med << endl;
      cerr << "  Avg size   " << sum_len/nb_merged_pairs << endl;
      cerr << "  Min size   " << min_readsize << endl;
      cerr << "  Max size   " << max_readsize << endl;
    }
  }

  free(fastx1);
  free(fastx2);
  return 0;
}



void error(string msg) {
  cerr << "[Error] " << msg << endl;
  cerr << "See fastx_mergepairs -h for more details." << endl;
  exit(1);
}



