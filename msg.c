#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void PrintMenu(void){
  fprintf(stderr,
  "Usage: smash-contigs [OPTION]... [FILE] [FILE]                           \n"
  "Similarity contigs mapping according to a reference sequence.            \n"
  "A server version, parallelizable, fast, but using high RAM memory.       \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                         give this help,                             \n"
  "  -V                         display version number,                     \n"
  "  -v                         verbose mode (more information),            \n"
  "  -k <k-mer>                 k-mer size [1;20],                          \n"
  "  -m <minimum>               minimum similar block size,                 \n"
  "  -r <repeats>               maximum repeats number,                     \n"
  "  -e <editions>              number of allowed editions,                 \n"
  "  -t <threshold>             maximum gaps allowed between copies,        \n"
  "  -l <link>                  link type between maps [0;4],               \n"
  "  -i                         do NOT map inversions,                      \n"
  "  -n <nThreads>              number of threads,                          \n"
  "  -o <FILE>                  output text filename with positions,        \n"
  "  -x <FILE>                  output image filename with map,             \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>                     contigs filename (FASTA),                   \n"
  "  <FILE>                     reference sequence filename (FASTA).        \n"
  "                                                                         \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                                \n");
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                         \n"
  "                         =====================                           \n"
  "                         | smash-contigs %u.%u |                         \n"
  "                         =====================                           \n"
  "                                                                         \n"
  "smash-contigs: similarity contigs mapping according to a ref sequence.   \n"
  "Copyright (C) 2015-2016 University of Aveiro. This is a Free software.   \n"
  "You may redistribute copies of it under the terms of the GNU - General   \n"
  "Public License v2 <http://www.gnu.org/licenses/gpl.html>. There is NOT   \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by   \n"
  "Diogo Pratas, Armando J. Pinho and Paulo J. S. G. Ferreira.\n\n", VERSION,
  RELEASE);
  }

