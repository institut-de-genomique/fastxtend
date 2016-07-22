## FASTXTEND

Fastxtend is an extension of FASTX-Toolkit package (http://hannonlab.cshl.edu/fastx_toolkit/index.html). 
Fastxtend is based on the library developed in fastx toolkit, it is written is C and it extend the FASTX-Toolkit with four commands line tools for Short-Reads FASTA/FASTQ files preprocessing:

- fastx_clean allows cleaning (adapters, N, quality) of the reads in fastq files. 
- fastx_duplicatedReads estimates the duplicates rate of reads (single or pair reads). It computes a rapid and accurate estimation of the duplicates rate of an initial read set using a sample of this read set.
- fastx_mergepairs perform the merging of the paired reads and give some statistics (merged size, percent of pairs merged).
- fastx_stats


Fastxtend is distributed open-source under CeCILL 
FREE SOFTWARE LICENSE. Check out http://www.cecill.info/
for more information about the contents of this license.

Fastxtend home on the web is http://www.genoscope.cns.fr/externe/fastxtend/


## COMMAND LINE ARGUMENTS


- All of the tools show usage information with -h 
- The option -Q is not documented in the usages it corresponds to the ASCII offset (generally -Q 33) and performs with all commands.


### fastx_clean


$ fastx_clean -h

usage: fastx_clean [-h] [-a ADAPTER_FILE] [-D] [-l N] [-n N] [-M N] [-m N] [-p N] [-c] [-C] [-o] [-v] [-z] [-i INFILE] [-o OUTFILE]
Developped at Genoscope using the FASTX Toolkit 0.0.13.1
    [-h]              = This helpful help screen.
    [-a ADAPTER_FILE] = ADAPTER file in fasta format.
    [-j]              = Keep the longest sequence before adaptater.
    [-l N]            = Discard sequences shorter than N nucleotides. default is 10.
    [-q N]            = Quality threshold - nucleotides with lower quality will be trimmed (from both ends of the sequence).
                       default value is 2, use 0 to inactivate this trimming.
                       [-n N]            = Trim sequences after N unknown nucleotides. default is 0 = off. This cleaning is done on the trimmed (adapter + quality) sequence.
   [-c]              = Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).
   [-C]              = Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).
   [-k]              = Report Adapter-Only sequences.
   [-v]              = Verbose - report number of sequences.
                       If [-o] is specified,  report will be printed to STDOUT.
                       If [-o] is not specified (and output goes to STDOUT),
                       report will be printed to STDERR.
   [-z]              = Compress output with GZIP.
   [-D]              = DEBUG output.
   [-M N]            = Require minimum adapter alignment length of N.
                       If less than N nucleotides aligned with the adapter - don't clip it.
   [-m N]            = Maximum number of mismatches allowed, default is 3
   [-p N]            = Allow one mismatch every N bases, default is 5 so one mismatch allow every 5 nucleotides
   [-r]              = Reverse complement the input fastx file and clip.
   [-f]              = Clip the input fastx file in forward, default is true.
   [-e]              = Recursive alignment : iterate alignments between sequence and adapters until a match is found.
   [-i INFILE]       = FASTA/Q input file. default is STDIN.
   [-o OUTFILE]      = FASTA/Q output file. default is STDOUT.
   [-s STAT_FILE]    = Tabular output file which contains trimming details for each input sequence.


### fastx_duplicatedReads


$ fastx_duplicatedReads -h
usage: fastx_duplicatedReads [-h] [-v] [-i INFILE] [-o OUTFILE]
Developped at Genoscope using the FASTX Toolkit 0.0.13.1

   [-h]         = This helpful help screen.
   [-s SAMPLE]  = FASTA/Q input file of a sample extract from INFILE
   [-i INFILE]  = FASTA/Q input file. Default is STDIN.
   [-t SAMPLE2]  = FASTA/Q input file of a sample extract from INFILE2 (Optional)
   [-j INFILE2]  = FASTA/Q input file for Read 2 (Optional)
   [-c INT]  = Trim sides of reads by a specified percentage (default: 0%)


### fastx_mergepairs


$ fastx_mergepairs -h
usage: fastx_mergepairs [-h] [-l N] [-m N] [-i N] [-s] [-M] [-a INFILE1] [-b INFILE2] [-o OUTFILE]
Developped at Genoscope using the FASTX Toolkit 0.0.13.1

   [-h]              = This helpful help screen.
   [-l N]            = Fragment size of read 2 used for detecting an overlap, default is 40
   [-m N]            = Maximal number of mismatches of the alignment of read 1 and subpart of read2, default is 4
   [-i N]            = Minimal identity percent of the alignment of read 1 and subpart of read2, default is 90
   [-L N]            = Minimal size of the alignment, default is 15
   [-s]              = Silent mode
   [-a INFILE1]      = FASTA/Q input file
   [-b INFILE2]      = FASTA/Q input file
   [-o OUTFILE]      = FASTA/Q output file
   [-u OUTFILE1]     = FASTA/Q output file of unpaired reads (e.g. merged reads).
   [-p OUTFILE1]     = FASTA/Q output file of paired reads (e.g. non-merged reads).
   [-q OUTFILE1a]    = FASTA/Q output file of read1 from paired reads (e.g. non-merged reads).
   [-r OUTFILE1b]    = FASTA/Q output file of read2 from paired reads (e.g. non-merged reads).
                       Choose between [-o] and [-u],[-p] and [-u],[-q],[-r] arguments.
   [-x FILE]         = Print distribution of overlap size between read1 and read2
   [-M]              = Only print merged pairs, default is no.


### fastx_stats


$ fastx_stats -h
usage: fastx_stats [[-h] [-f INFILE] [-q]]

   [-h]         = This helpful help screen.
   [-f INFILE]  = FASTA/Q input file. default is STDIN.
   [-q]         = Display quality values distribution.


## PRE-REQUISITES

  - A Linux based operating system.
  - Binaries are provided for the following platform : Linux x86_64
  - g++ with gcc 4.1.2 or higher


## INSTALLATION

  1. Clone this GitHub repository    
  2. Compile sources    
  `make;` 
  3. Install binaries  
   `make install`
   This will install the tools into ./bin


## More informations

If you have questions about Fastxtend, you may ask them to sengelen [at] genoscope [.] cns [.] fr and jmaury [at] genoscope [.] cns [.] fr . You may also create an issue to ask questions on github website: https://github.com/institut-de-genomique/fastxtend/issues. 


## ACKNOWLEDGMENTS

Stefan Engelen, Cyril Firmo and Jean-Marc Aury - Fastxtend's authors

This work was financially supported by the Genoscope, 
Institut de Genomique, CEA and Agence Nationale de la 
Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).
