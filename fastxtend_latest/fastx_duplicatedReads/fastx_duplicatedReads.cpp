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
#include <err.h>
#include <getopt.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <stdio.h>

#include "config.h"

#include "fastx.h"

#include "fastx_args.h"

using namespace std;

const char* usage=
"usage: fastx_collapser [-h] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Developped at Genoscope using the " PACKAGE_STRING "\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-s SAMPLE]  = FASTA/Q input file of a sample extract from INFILE\n" \
"   [-i INFILE]  = FASTA/Q input file. Default is STDIN.\n" \
"   [-t SAMPLE2]  = FASTA/Q input file of a sample extract from INFILE2 (Optional)\n" \
"   [-j INFILE2]  = FASTA/Q input file for Read 2 (Optional)\n" \
"   [-c INT]  = Trim sides of reads by a specified percentage (default: 0%)\n" \
"\n";

char* sample1_filename = NULL;
char* sample2_filename = NULL;
char* input1_filename = NULL;
char* input2_filename = NULL;
int trimnucl = 0;
FASTX sample, fastx, sample2, fastx2;

typedef struct occurrence {long long int nbFastx;long long int nbSample;};

#include <tr1/unordered_map>
std::tr1::unordered_map<string, occurrence> collapsed_sequencesA, collapsed_sequencesB, collapsed_sequencesC;


int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg) {
  switch(optc) {
  case 's':
    sample1_filename = optarg;
    break;

  case 't':
    sample2_filename = optarg;
    break;

  case 'j':
	input2_filename = optarg;
	break;

  case 'c':
	trimnucl = atoi(optarg);
	break;

  default:
    errx(1,"Unknown argument (%c)", optc ) ;
  }
  return 1;
}

string id_sequence(string A, string B) {
	if (A.compare(B) < 0){
		  return A+B;
	}else{
		  return B+A;
	}

}

string trim(string A, int P) {
	int x,y;
	if ((P <= 0) || (P >= 50)){
		return A;
	}else {
		x = int(A.length() * (float(P)/100));
		y = A.length() - (2 * x);
		return A.substr(x,y);
	}

}

int parse_commandline(int argc, char* argv[]) {  
  fastx_parse_cmdline(argc, argv, "s:t:j:c:", parse_program_args);
  return 1;
}


int main(int argc, char* argv[]) {
  bool pairs(true);
  std::tr1::unordered_map<string, occurrence>::iterator it;
  parse_commandline(argc, argv);

  if ((string(get_input_filename()) == "-") || (sample1_filename == NULL)){
	    errx(1,"Bad definition of reference or input file") ;
  }
  if (((input2_filename == NULL) && (sample2_filename != NULL)) || ((input2_filename != NULL) && (sample2_filename == NULL))){
	    errx(1,"Bad definition of reference or input file (Read 2)") ;
  }
  if ((input2_filename == NULL) || (sample2_filename == NULL)){
		pairs = false;
  }else{
	    pairs = true;
  }
  if ((trimnucl <= 0) || (trimnucl >= 50)){
	    trimnucl = 0;
  }

  fastx_init_reader(&fastx, get_input_filename(),
		    FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		    get_fastq_ascii_quality_offset() );
  
  fastx_init_reader(&sample, sample1_filename,
		    FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		    get_fastq_ascii_quality_offset() );

  if (pairs == true){
	  fastx_init_reader(&fastx2, input2_filename,
		    FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		    get_fastq_ascii_quality_offset() );

	  fastx_init_reader(&sample2, sample2_filename,
		    FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		    get_fastq_ascii_quality_offset() );
  }

  if (pairs == true){
	  while ( fastx_read_next_record(&sample) && fastx_read_next_record(&sample2) ) {
		  it = collapsed_sequencesA.find(id_sequence(trim(string(sample.nucleotides),trimnucl),trim(string(sample2.nucleotides),trimnucl)));
		  if(it != collapsed_sequencesA.end()) it->second.nbSample++;
		  else
		  {
			  occurrence initOcc; initOcc.nbSample=1; initOcc.nbFastx=0;
			  collapsed_sequencesA[id_sequence(trim(string(sample.nucleotides),trimnucl),trim(string(sample2.nucleotides),trimnucl))]=initOcc;
		  }

		  it = collapsed_sequencesB.find(trim(string(sample.nucleotides),trimnucl));
		  if(it != collapsed_sequencesB.end()) it->second.nbSample++;
		  else
		  {
			  occurrence initOcc; initOcc.nbSample=1; initOcc.nbFastx=0;
			  collapsed_sequencesB[trim(string(sample.nucleotides),trimnucl)]=initOcc;
		  }

		  it = collapsed_sequencesC.find(trim(string(sample2.nucleotides),trimnucl));
		  if(it != collapsed_sequencesC.end()) it->second.nbSample++;
		  else
		  {
			  occurrence initOcc; initOcc.nbSample=1; initOcc.nbFastx=0;
			  collapsed_sequencesC[trim(string(sample2.nucleotides),trimnucl)]=initOcc;
		  }

	  }
	  while ( fastx_read_next_record(&fastx) && fastx_read_next_record(&fastx2) ) {
		  it = collapsed_sequencesA.find(id_sequence(trim(string(fastx.nucleotides),trimnucl),trim(string(fastx2.nucleotides),trimnucl)));
		  if(it != collapsed_sequencesA.end()) it->second.nbFastx++; //get_reads_count(&fastx);
		  it = collapsed_sequencesB.find(trim(string(fastx.nucleotides),trimnucl));
		  if(it != collapsed_sequencesB.end()) it->second.nbFastx++; //get_reads_count(&fastx);
		  it = collapsed_sequencesC.find(trim(string(fastx2.nucleotides),trimnucl));
		  if(it != collapsed_sequencesC.end()) it->second.nbFastx++; //get_reads_count(&fastx);
	  }
  }else{
	  while ( fastx_read_next_record(&sample) ) {
		  it = collapsed_sequencesA.find(trim(string(sample.nucleotides),trimnucl));
		  if(it != collapsed_sequencesA.end()) it->second.nbSample++;
		  else
		  {
			  occurrence initOcc; initOcc.nbSample=1; initOcc.nbFastx=0;
			  collapsed_sequencesA[trim(string(sample.nucleotides),trimnucl)]=initOcc;
		  }
	  }
	  while ( fastx_read_next_record(&fastx) ) {
	    it = collapsed_sequencesA.find(trim(string(fastx.nucleotides),trimnucl));
	    if(it != collapsed_sequencesA.end()) it->second.nbFastx++; //get_reads_count(&fastx);
	  }
  }
  

  int occ_max = 40;
  int occ_max_printed = 10;
  vector<int> nbReads(occ_max+1, 0);
  for(it = collapsed_sequencesA.begin(); it != collapsed_sequencesA.end() ; it++) {
    int indice = (it->second.nbFastx < occ_max) ? it->second.nbFastx : occ_max;
    int nbVu = it->second.nbSample ;
    nbReads[indice]=nbReads[indice]+nbVu;
  }

  int nbDupl = 0 ;
  for(int i=1 ; i<=occ_max ; i++) {
	  nbDupl = nbDupl + nbReads[i]*(i-1)/i;
  }

  double duplRate = ((double)nbDupl / (double)num_input_reads(&sample));

  int nbUniqReads = (int)((double)num_input_reads(&fastx) * (1.0 - duplRate));

  cout << "SideTrimming= " <<  (double)trimnucl << "%" << endl;

  if (pairs == true){
	  cout << "\n-- ESTIMATION OF DUPLICATED READS (PAIRED) \n" << endl;
  }else{
	  cout << "\n--  ESTIMATION OF DUPLICATED READS -- \n" << endl;
  }
  
  cout << "nbReadsInSample= " << num_input_reads(&sample) << " EstimateDuplRate= " << (duplRate * 100.0) << endl;
  cout << "nbReadsInInput= " << num_input_reads(&fastx) << " EstimateUniqReads= " <<  nbUniqReads << "\n" << endl;
  
  for(int i=1 ; i<occ_max_printed ; i++) 
	  cout << "nbOcc= " << i << " rate= " << (((double)nbReads[i] / (double)num_input_reads(&sample)) * 100.0) << endl;
  double rest = 0.0 ;
  for(int i=occ_max_printed ; i<=occ_max ; i++)
      rest = rest + (double)nbReads[i] ;
  cout << "nbOcc> " << occ_max_printed << " rate= " << (((double)rest / (double)num_input_reads(&sample)) * 100.0) << endl;

  if (pairs == true){
	  vector<int> nbReadsB(occ_max+1, 0);
	  vector<int> nbReadsC(occ_max+1, 0);
	  for(it = collapsed_sequencesB.begin() ; it != collapsed_sequencesB.end() ;it++) {
		    int indice = (it->second.nbFastx < occ_max) ? it->second.nbFastx : occ_max;
		    int nbVu = it->second.nbSample ;
		    nbReadsB[indice]=nbReadsB[indice]+nbVu;
	  }
	  for(it = collapsed_sequencesC.begin() ; it != collapsed_sequencesC.end() ;it++) {
		    int indice = (it->second.nbFastx < occ_max) ? it->second.nbFastx : occ_max;
		    int nbVu = it->second.nbSample ;
		    nbReadsC[indice]=nbReadsC[indice]+nbVu;
	  }


	  int nbDuplB = 0 ;
	  int nbDuplC = 0 ;
	  for(int i=1 ; i<=occ_max ; i++) {
		  nbDuplB = nbDuplB + nbReadsB[i]*(i-1)/i;
		  nbDuplC = nbDuplC + nbReadsC[i]*(i-1)/i;
	  }

	  double duplRateB = ((double)nbDuplB / (double)num_input_reads(&sample));
	  int nbUniqReadsB = (int)((double)num_input_reads(&fastx) * (1.0 - duplRateB));

	  double duplRateC = ((double)nbDuplC / (double)num_input_reads(&sample2));
	  int nbUniqReadsC = (int)((double)num_input_reads(&fastx2) * (1.0 - duplRateC));
	  cout << "\n-- ESTIMATION OF DUPLICATED READS (READ 1) " << endl;
	  cout << "\nnbReadsInSample= " << num_input_reads(&sample) << " EstimateDuplRate= " << (duplRateB * 100.0) << endl;
	  cout << "nbReadsInInput= " << num_input_reads(&fastx) << " EstimateUniqReads= " <<  nbUniqReadsB << "\n" << endl;
	  for(int i=1 ; i<occ_max_printed ; i++)
		  cout << "nbOcc= " << i << " rate= " << (((double)nbReadsB[i] / (double)num_input_reads(&sample)) * 100.0) << endl;
	  rest = 0.0 ;
	  for(int i=occ_max_printed ; i<=occ_max ; i++)
	  	  rest = rest + (double)nbReadsB[i] ;
	  cout << "nbOcc> " << occ_max_printed << " rate= " << (((double)rest / (double)num_input_reads(&sample)) * 100.0) << endl;
	  cout << "\n-- ESTIMATION OF DUPLICATED READS (READ 2) " << endl;
	  cout << "\nnbReadsInSample= " << num_input_reads(&sample2) << " EstimateDuplRate= " << (duplRateC * 100.0) << endl;
	  cout << "nbReadsInInput= " << num_input_reads(&fastx2) << " EstimateUniqReads= " <<  nbUniqReadsC << "\n" << endl;
	  for(int i=1 ; i<occ_max_printed ; i++)
		  cout << "nbOcc= " << i << " rate= " << (((double)nbReadsC[i] / (double)num_input_reads(&sample2)) * 100.0) << endl;
	  rest = 0.0 ;
	  for(int i=occ_max_printed ; i<=occ_max ; i++)
		  rest = rest + (double)nbReadsC[i] ;
	  cout << "nbOcc> " << occ_max_printed << " rate= " << (((double)rest / (double)num_input_reads(&sample2)) * 100.0) << endl;
  }
  return 0;
}
