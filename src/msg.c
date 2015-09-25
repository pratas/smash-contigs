#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void PrintMenuMap(void){
  fprintf(stderr,
  "Usage: smash-map     [OPTION]... [FILE] [FILE]                           \n"
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
  "  -i                         do NOT map inversions,                      \n"
  "  -n <nThreads>              number of threads,                          \n"
  "  -o <FILE>                  output text filename with positions,        \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>                     contigs filename (FASTA),                   \n"
  "  <FILE>                     reference sequence filename (FASTA).        \n"
  "                                                                         \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                                \n");
  }

void PrintMenuReduce(void){
  fprintf(stderr,
  "Usage: smash-reduce [OPTION]... [FILE]                                   \n"
  "Reduction and clustering of contig maps according to smash-map output.   \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                         give this help,                             \n"
  "  -V                         display version number,                     \n"
  "  -v                         verbose mode (more information),            \n"
  "  -d                         delete input map file,                      \n"
  "  -e <editions>              maximum gaps allowed between maps,          \n"
  "  -o <FILE>                  output filename with reduced maps,          \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>                     contigs filename with positions (.pos),    \n"
  "                                                                         \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                                \n");
  }

void PrintMenuVisual(void){
  fprintf(stderr,
  "Usage: smash-visual [OPTION]... [FILE]                                   \n"
  "Visualisation of contigs map according to smash-(map/reduce) output.     \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                         give this help,                             \n"
  "  -V                         display version number,                     \n"
  "  -v                         verbose mode (more information),            \n"
  "  -l <link>                  link type between maps [0;4],               \n"
  "  -w <width>                 chromosome width,                           \n"
  "  -s <space>                 space between chromosomes,                  \n"
  "  -m <mult>                  color id multiplication factor,             \n"
  "  -b <begin>                 color id beggining,                         \n"
  "  -c <minimum>               minimum block size to consider,             \n"
  "  -i                         do NOT show inversion maps,                 \n"
  "  -r                         do NOT show regular maps,                   \n"
  "  -o <FILE>                  output image filename with map,             \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>                     contigs filename with positions (.pos),    \n"
  "                                                                         \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                                \n");
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                         \n"
  "                         =====================                           \n"
  "                         | SMASH-CONTIGS %u.%u |                         \n"
  "                         =====================                           \n"
  "                                                                         \n"
  "            SMASH-CONTIGS: similarity contigs mapping, reducing          \n"
  "            and visualisation according to a reference sequence          \n"
  "                                                                         \n"
  "                  It includes the following programs:                    \n"
  "                            * smash-map                                  \n" 
  "                            * smash-reduce                               \n" 
  "                            * smash-visual                               \n" 
  "                                                                         \n"
  "Copyright (C) 2015-2016 University of Aveiro. This is a Free software.   \n"
  "You may redistribute copies of it under the terms of the GNU - General   \n"
  "Public License v2 <http://www.gnu.org/licenses/gpl.html>. There is NOT   \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by   \n"
  "Diogo Pratas, Armando J. Pinho and Paulo J. S. G. Ferreira.\n\n", VERSION,
  RELEASE);
  }

