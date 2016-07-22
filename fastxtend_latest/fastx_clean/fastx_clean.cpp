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

#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>


#include "sequence_alignment.h"

#include <errno.h>
#include <err.h>

#include <config.h>

#include "fastx.h"
#include "fastx_args.h"

using namespace std;

#define MAX_ADAPTER_LEN 100

const char* usage=
"usage: fastx_clean [-h] [-a ADAPTER_FILE] [-D] [-l N] [-n N] [-M N] [-m N] [-p N] [-c] [-C] [-o] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Developped at Genoscope using the " PACKAGE_STRING "\n" \
"\n" \
"   [-h]              = This helpful help screen.\n" \
"   [-a ADAPTER_FILE] = ADAPTER file in fasta format.\n" \

"   [-j]              = Keep the longest sequence before adaptater.\n" \
"   [-l N]            = Discard sequences shorter than N nucleotides. default is 10.\n" \
"   [-q N]            = Quality threshold - nucleotides with lower quality will be trimmed (from both ends of the sequence). \n" \
"                       default value is 2, use 0 to inactivate this trimming.\n" \
"   [-n N]            = Trim sequences after N unknown nucleotides. default is 0 = off. This cleaning is done on the trimmed (adapter + quality) sequence.\n" \
"   [-c]              = Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).\n" \
"   [-C]              = Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).\n" \
"   [-k]              = Report Adapter-Only sequences.\n" \

"   [-v]              = Verbose - report number of sequences.\n" \
"                       If [-o] is specified,  report will be printed to STDOUT.\n" \
"                       If [-o] is not specified (and output goes to STDOUT),\n" \
"                       report will be printed to STDERR.\n" \
"   [-z]              = Compress output with GZIP.\n" \
"   [-D]              = DEBUG output.\n" \

"   [-M N]            = Require minimum adapter alignment length of N.\n" \
"                       If less than N nucleotides aligned with the adapter - don't clip it.\n" \
"   [-m N]            = Maximum number of mismatches allowed, default is 3\n" \
"   [-p N]            = Allow one mismatch every N bases, default is 5 so one mismatch allow every 5 nucleotides\n" \

"   [-r]              = Reverse complement the input fastx file and clip.\n" \
"   [-f]              = Clip the input fastx file in forward, default is true.\n" \

"   [-e]              = Recursive alignment : iterate alignments between sequence and adapters until a match is found.\n"\

"   [-i INFILE]       = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE]      = FASTA/Q output file. default is STDOUT.\n" \
"   [-s STAT_FILE]    = Tabular output file which contains trimming details for each input sequence.\n" \
"\n";


unsigned int min_length = 10;
int discard_unknown_bases = 0;
int discard_non_clipped = 0;
int discard_clipped = 0;
int show_adapter_only = 0;
int junction_adapter = 0;
int debug = 0;
int minimum_adapter_length = 0;
int max_mismatches = 3;
int nb_nt_per_mismatches = 5;
int reverse_complement = 0;
int forward = 0;
int quality_threshold = 2;
int recursive_alignment = 0;
char* adapters_file = NULL;
int searchAdapter = 0;
char* stats_file = NULL;

//Statistics for verbose report
unsigned int count_input=0 ;
unsigned int count_discarded_too_short=0; // see [-l N] option
unsigned int count_discarded_adapter_at_index_zero=0;  //empty sequences (after clipping)
unsigned int count_discarded_no_adapter_found=0; // see [-c] option
unsigned int count_discarded_adapter_found=0; // see [-C] option
unsigned int count_discarded_N=0; // see [-n]
unsigned int count_nb_nt_clipped=0;
unsigned int count_nb_nt_out=0;
unsigned int count_nb_reads_clipped=0;
unsigned int max_length = 0;	// save max length after clipping

FASTX fastx, adapters;
HalfLocalSequenceAlignment align;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg) {
  switch(optc) {
  case 'M':
    if (optarg==NULL) errx(1, "[-M] parameter requires an argument value");
    minimum_adapter_length = atoi(optarg);
    if (minimum_adapter_length<=0) errx(1,"Invalid minimum adapter length (-M %s)", optarg);
    break;

  case 'm':
    if (optarg==NULL) errx(1, "[-m] parameter requires an argument value");
    max_mismatches = atoi(optarg);
    if (max_mismatches<=0) errx(1,"Invalid maximum number of mismatches (-m %s)", optarg);
    break;

  case 'p':
    if (optarg==NULL) errx(1, "[-p] parameter requires an argument value");
    nb_nt_per_mismatches = atoi(optarg);
    if (nb_nt_per_mismatches<=0) errx(1,"Invalid maximum number of nucleotide per mismatch (-m %s)", optarg);
    break;
  case 'j':
    junction_adapter=1;
    break;
  case 'k':
    show_adapter_only=1;
    break;

  case 'r':
    reverse_complement=1;
    break;

  case 'f':
    forward=1;
    break;

  case 'e':
    recursive_alignment=1;
    break;

  case 'D':
    debug++;
    break ;

  case 'c':
    discard_non_clipped = 1;
    break;

  case 'C':
    discard_clipped = 1 ;
    break ;

  case 'a':
	if (optarg==NULL) errx(1,"[-a] parameter requires an argument value");
    adapters_file = optarg;
    struct stat sb;
    if (fopen(adapters_file,"r")!=NULL)
    {
    	stat(adapters_file, &sb);
    	if (sb.st_size>0) searchAdapter=1;
    	else errx(1,"[-a] : input adaptator file empty");
    }
    else errx(1,"[-a] input adaptator file '%s' not exist",optarg);
    break;

  case 's':
    stats_file = optarg;
    break;

  case 'l':
    if (optarg==NULL) errx(1,"[-l] parameter requires an argument value");
    min_length = strtoul(optarg, NULL, 10);
    break;

  case 'q':
    if (optarg==NULL) errx(1,"[-q] parameter requires an argument value");
    quality_threshold = strtoul(optarg, NULL, 10);
    if(quality_threshold < 0 || quality_threshold > 64) errx(1,"Invalid quality threshold (-q %s)", optarg);
    break;

  case 'n':
    if (optarg==NULL) errx(1,"[-n] parameter requires an argument value");
    discard_unknown_bases = strtoul(optarg, NULL, 10);
    if(discard_unknown_bases < 0) errx(1,"Invalid number of unknown bases (-n %s)", optarg);
    break;

  default:
    errx(1,"Unknown argument (%c)", optc ) ;
  }
  return 1;
}



int parse_commandline(int argc, char* argv[]) {
  fastx_parse_cmdline(argc, argv, "M:m:p:kDCcrjfea:q:s:l:n:", parse_program_args);
  return 1;
}



int adapter_cutoff_index ( const SequenceAlignmentResults& alignment_results ) __attribute__ ((const));
int adapter_cutoff_index ( const SequenceAlignmentResults& alignment_results ) {

  int alignment_size = alignment_results.neutral_matches + alignment_results.matches +
    alignment_results.mismatches + alignment_results.gaps ;

  int whole_mismatches = alignment_results.mismatches + alignment_results.neutral_matches;
  whole_mismatches += std::min((alignment_results.query_size - alignment_results.query_end), (alignment_results.target_size - alignment_results.target_end));
  whole_mismatches += std::min((alignment_results.query_start), (alignment_results.target_start));

  //No alignment at all?
  if (alignment_size==0) return -1;
  if (minimum_adapter_length > 0 && alignment_size < minimum_adapter_length) return -1;

  int nb_mismatches_allowed = std::min((alignment_size / nb_nt_per_mismatches), max_mismatches);
  if(whole_mismatches <= nb_mismatches_allowed) return alignment_results.query_start;

  return -1;
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
      errx(1,"Invalid nucleotide value (%c) in reverse_complement_base()", input );
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



void get_longest_unclipped_sequence(std::vector<int>& read, int& start, int& stop) {
  int readlen = read.size() - 1;
  int len=0;

  vector<int>::iterator itS = read.begin()+1, itE;

  while( itE != read.end() ) {
    itS = std::find(itS, read.end(), 1);
    itE = std::find(itS, read.end(), 0);
    itS--;

    if(std::distance(itS, itE)-1 > len) {
      start = itS - read.begin() + 1;
      len = std::distance(itS, itE) - 1;
    }
    itS=itE;
  }

  stop = start + len - 1;
}



/*
int unrecognize_sequence(std::vector<int>& read, string& sequence, string& query) {
  int start=1, stop=1;
  get_longest_unclipped_sequence(read, start, stop);
  query = sequence.substr(start-1, stop-start+2);
  cout << "start= " << start << " stop= " << stop << " query=" << query << endl;
  return start-1;
}
*/


void clip_adapter(std::vector<int>& read, FASTX& fastx, FASTX& adapters, bool isRev, FILE* fhstats) {
  bool recursive = true;
  //std::string sequence = std::string(fastx.nucleotides);
  int step=0;
  std::vector<int> internal_read(read.size(), 1);

  //if(isRev) for(int i = 0; i<read.size(); i++) internal_read[i] = read[read.size()-i-1];
  //else
  for(int i = 0; i<read.size(); i++) internal_read[i] = read[i];

  while ( fastx_read_next_record(&adapters) ) {
    if(isRev) reverse_complement_fastx(&adapters);
    std::string target = std::string(adapters.nucleotides);
    std::string sequence = std::string(fastx.nucleotides);
    recursive = true;

    while(recursive) {
      step++;
      align.align( sequence, target ) ;

      if (debug>1) align.print_matrix();
      if (debug>0) align.results().print();

      //Find the best match with the adapter
      int i = adapter_cutoff_index ( align.results() ) ;
      int readlen = read.size() - 1;
      if(i>=0 && i<readlen) {
    	  int align_end = (align.results()).query_end;
    	  int endpos = std::min(align_end, readlen-1);
    	  //int start = (isRev) ? i : i+1;
    	  int start = i+1;
    	  //int end   = (isRev) ? endpos : endpos+1;
    	  int end = endpos+1;
    	  int strand = (isRev) ? -1 : 1;
    	  if(fhstats != NULL) fprintf(fhstats, "%s\t%s\t%i\t%i\t%i\n", fastx.name, adapters.name, start, end, strand);
    	  if(junction_adapter){
    		  end = readlen;
    	  }
    	  for(int cpt = start; cpt <= end ; cpt++) {
    		  internal_read[cpt]=0;
    		  sequence[cpt-1]='N';
    	  }
      }
      recursive = (recursive_alignment && i>=0);
    }
    //if(isRev) reverse_complement_fastx(&adapters);
  }

  //if(isRev) for(int i = 0; i<read.size(); i++) read[i] = internal_read[read.size()-i-1];
  //else
  for(int i = 0; i<read.size(); i++) read[i] = internal_read[i];

  rewind(adapters.input);
}



void clip_lowqual(std::vector<int>& read, FASTX& fastx, int min_quality_threshold, FILE* fhstats) {
  int readlen = read.size() - 1;
  int i = readlen;
  // from end of the read
  for ( ; i >=1 ; i--)
    if (fastx.quality[i-1] <= min_quality_threshold) read[i]=0;
    else break;
  if(i != readlen && fhstats != NULL) fprintf(fhstats, "%s\tlow_qual\t%i\t%i\t1\n", fastx.name, i+1, readlen);
  if(i==1) return;

  // from start of the read
  for (i=0 ; i < readlen ; i++)
    if (fastx.quality[i] <= min_quality_threshold) read[i+1]=0;
    else break;
  if(i != 0 && fhstats != NULL) fprintf(fhstats, "%s\tlow_qual\t%i\t%i\t1\n", fastx.name, 1, i);
}



void clip_N(std::vector<int>& read, FASTX& fastx, int max_number_of_N, FILE* fhstats) {
  int readlen = read.size() - 1;
  int i = readlen, cptN = 0, lastGoodNT = 0;
  // from end of the read
  for ( ; i >=1 ; i--)
    if (fastx.nucleotides[i-1] == 'N') read[i]=0;
    else break;
  if(i != readlen && fhstats != NULL) fprintf(fhstats, "%s\tN\t%i\t%i\t1\n", fastx.name, i+1, readlen);
  if(i==1) return;

  // from start of the read
  for (i=0 ; i < readlen ; i++)
    if (fastx.nucleotides[i] == 'N') read[i+1]=0;
    else break;
  if(i != 0 && fhstats != NULL) fprintf(fhstats, "%s\tN\t%i\t%i\t1\n", fastx.name, 1, i);

  // allow only N unknown nucleotides
  for (int j=i ; j < readlen ; j++) {
    if(fastx.nucleotides[j] == 'N') cptN++;
    else lastGoodNT = j+1;
    if(cptN >= max_number_of_N) {
      j = lastGoodNT;
      if(fhstats != NULL) fprintf(fhstats, "%s\tN\t%i\t%i\t1\n", fastx.name, j, readlen);
      for(; j < readlen ; j++) read[j+1]=0;
    }
  }
}



void clip_read(std::vector<int>& read, FASTX& fastx, FILE* fhstats) {
  int readlen = read.size() - 1;
  int start=1, len=0, stop=1;

  vector<int>::iterator itS = read.begin()+1, itE;

  while( itE != read.end() ) {
    itS = std::find(itS, read.end(), 1);
    itE = std::find(itS, read.end(), 0);
    itS--;

    if(std::distance(itS, itE)-1 > len) {
      start = itS - read.begin() + 1;
      len = std::distance(itS, itE) - 1;
    }
    itS=itE;
  }

  stop = start + len - 1;
  int i = start;
  for( ; i<=stop ; i++) {
    fastx.nucleotides[i-start] = fastx.nucleotides[i-1];
    fastx.quality[i-start] = fastx.quality[i-1];
  }
  fastx.nucleotides[i-start] = 0;
  fastx.quality[i-start] = 0;
}



int main(int argc, char* argv[]) {
  int i, initial_length;
  int reads_count;

  parse_commandline(argc, argv);
  if(!reverse_complement) forward=1;

  FILE* fh_stats = NULL;
  if(stats_file != NULL) {
    fh_stats = fopen(stats_file, "w");
    if (fh_stats == NULL) err(1, "failed to open input file '%s'", stats_file);
  }

  fastx_init_reader(&fastx, get_input_filename(),
		    FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		    get_fastq_ascii_quality_offset() );

  if(searchAdapter == 1){
  	  fastx_init_reader(&adapters, adapters_file,
  				FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
  				get_fastq_ascii_quality_offset() );
    }

  fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

  while ( fastx_read_next_record(&fastx) ) {
    reads_count = get_reads_count(&fastx);
    count_input+= reads_count;
    initial_length = strlen(fastx.nucleotides);
    std::vector<int> read(initial_length+1, 1);
    read[0]=0;

    
    if(searchAdapter == 1){
    	if(forward) clip_adapter(read, fastx, adapters, 0, fh_stats);
    	if(reverse_complement) clip_adapter(read, fastx, adapters, 1, fh_stats);
    }

    if(quality_threshold > 0) clip_lowqual(read, fastx, quality_threshold, fh_stats);

    clip_read(read, fastx, fh_stats);

    std::vector<int> read2(initial_length+1, 1);
    read2[0]=0;
    if(discard_unknown_bases) clip_N(read2, fastx, discard_unknown_bases, fh_stats);
    clip_read(read2, fastx, fh_stats);


    int clipread_len = strlen(fastx.nucleotides);

    if ( clipread_len == initial_length && discard_non_clipped ) { // adapter not found (i.e. sequence was not clipped) ?
      count_discarded_no_adapter_found += reads_count;
      continue;
    }

    if (clipread_len > max_length) { // save max length
	max_length = clipread_len;
    }

    if(clipread_len != initial_length) {

      if (debug>0) printf("TRIMMED --------------> output = %s\n", fastx.nucleotides);
      count_nb_nt_clipped += (initial_length - clipread_len);
      count_nb_reads_clipped++;

      if (clipread_len == 0) { // empty sequence ? (in which the adapter was found at index 0)
	count_discarded_adapter_at_index_zero += reads_count;
	if (show_adapter_only) fastx_write_record(&fastx);
	continue;
      }

      if (clipread_len < min_length) { // too-short sequence ?
	count_discarded_too_short += reads_count;
	continue;
      }

      if (discard_clipped) { // adapter found, and user requested to keep only non-clipped sequences
	count_discarded_adapter_found += reads_count;
	continue;
      }
    }

    if (!show_adapter_only) {
      //none of the above condition matched, so print this sequence.
      fastx_write_record(&fastx);
      count_nb_nt_out += clipread_len ;
    }
  }

  //Print verbose report
  if ( verbose_flag() ) {
	if(searchAdapter == 1)
		fprintf(get_report_file(), "Clipping Adapter: %s\n", adapters_file );
    fprintf(get_report_file(), "Min. Length: %d\n", min_length) ;

    if (discard_clipped)
      fprintf(get_report_file(), "Clipped reads - discarded.\n"  ) ;
    if (discard_non_clipped)
      fprintf(get_report_file(), "Non-Clipped reads - discarded.\n"  ) ;

    fprintf(get_report_file(), "Input: %u reads.\n", count_input ) ;
    fprintf(get_report_file(), "Output: %u reads (max read length : %u) , %u nt.\n",
	    count_input - count_discarded_too_short - count_discarded_no_adapter_found - count_discarded_adapter_found -
	    count_discarded_N - count_discarded_adapter_at_index_zero, max_length, count_nb_nt_out ) ;
    fprintf(get_report_file(), "clipped: %u reads ( %.2f \% ) , %u nt.\n", count_nb_reads_clipped, ((float)count_nb_reads_clipped / (float)count_input) * 100.0, count_nb_nt_clipped );

    fprintf(get_report_file(), "discarded %u too-short reads.\n", count_discarded_too_short ) ;
    fprintf(get_report_file(), "discarded %u reads of size 0.\n", count_discarded_adapter_at_index_zero );
    if (discard_non_clipped)
      fprintf(get_report_file(), "discarded %u non-clipped reads.\n", count_discarded_no_adapter_found );
    if (discard_clipped)
      fprintf(get_report_file(), "discarded %u clipped reads.\n", count_discarded_adapter_found );
  }
  if(fh_stats != NULL) fclose(fh_stats);
  return 0;
}
