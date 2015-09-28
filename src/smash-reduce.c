#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/uio.h>
#include <sys/mman.h>
#include "mem.h"
#include "lines.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "common.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E D U C E - - - - - - - - - - - - - - - -
void ReduceProjections(char *fn, uint8_t delete){
  FILE *IN = NULL, *OUT = NULL;
  char name[MAX_FILENAME], watermark[MAX_FILENAME];
  sprintf(name, "%s.red", fn);
  int64_t refNBases, conNBases, ri, rf, ci, cf, cx, cy, rx, ry;

  IN  = Fopen(fn,   "r");
  OUT = Fopen(name, "w");

  fprintf(stderr, "Reducing projections ... \n");

  // READ HEADER
  if(fscanf(IN, "%s\t%"PRIi64"\t%"PRIi64"\n", watermark, &conNBases,
  &refNBases) != 3 || watermark[0] != '#' || watermark[1] != 'S' ||
  watermark[2] != 'C' || watermark[3] != 'F'){
    fprintf(stderr, "[x] Error: unknown positions file format!\n");
    exit(1);
    }
  fprintf(OUT, "#SCF\t%"PRIi64"\t%"PRIi64"\n", conNBases, refNBases);

  // READ BODY
  LCACHE *LCache = CreateLCache(100);
  int64_t idx = 0, idxIr = 0, lines = 0;

  while(1){
    char tmp1[MAX_STR] = {'\0'}, tmp2[MAX_STR] = {'\0'};
    if(fscanf(IN, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                  "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n", 
                  tmp1, 
                  &LCache->Lines[LCache->idx].contigs_relative_init_pos, 
                  &LCache->Lines[LCache->idx].contigs_relative_end_pos, 
                  &LCache->Lines[LCache->idx].contigs_absolute_init_pos, 
                  &LCache->Lines[LCache->idx].contigs_absolute_end_pos, 

                  tmp2,
                  &LCache->Lines[LCache->idx].reference_relative_init_pos,
                  &LCache->Lines[LCache->idx].reference_relative_end_pos, 
                  &LCache->Lines[LCache->idx].reference_absolute_init_pos, 
                  &LCache->Lines[LCache->idx].reference_absolute_end_pos) 
                  != 10)
      break; // FOUND UNEXPECTED LINE OR END OF FILE

    if(cf > ci){

      if(lines != 0){ // IT IS NOT THE FIRST REGULAR PATTERN
//      if(posCache[idx][1] - posCache[idx-1][0] <= P->threshold){



          }
        else{


          }

      if(++idx == MAX_POS_CACHE)
        idx = 0;
      }
    else{ // INVERTED
//      posCacheIr[idxIr][0] = cx;
//      posCacheIr[idxIr][1] = cy;
//      posCacheIr[idxIr][2] = rx;
//      posCacheIr[idxIr][3] = ry;

      if(++idxIr == MAX_POS_CACHE)
        idxIr = 0;
      }

    fprintf(OUT, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                 "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n", 
                 tmp1, ci, cf, cx, cy, tmp2, ri, rf, rx, ry);
    ++lines;
    }
  fprintf(stderr, "\rDone!                   \n");

  RemoveLCache(LCache);
  if(!delete)
    unlink(fn);
  fclose(IN);
  fclose(OUT);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t main(int argc, char *argv[]){
  char **p = *&argv;

  uint8_t delete;
  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEF_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenuReduce();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  P->verbose    = ArgsState (DEF_VERBOSE, p, argc, "-v" );
  P->force      = ArgsState (DEF_FORCE,   p, argc, "-F" );
  P->inversion  = ArgsState (DEF_INVE,    p, argc, "-i" );
  delete        = ArgsState (DEF_DELE,    p, argc, "-d" );
  P->minimum    = ArgsNum   (DEF_MINI,    p, argc, "-m", MIN_MINI, MAX_MINI);
  P->threshold  = ArgsNum   (DEF_TSHO,    p, argc, "-t", MIN_TSHO, MAX_TSHO);
  P->nThreads   = ArgsNum   (DEF_THRE,    p, argc, "-n", MIN_THRE, MAX_THRE);
  P->positions  = ArgsFilesReduce        (p, argc, "-o");

  if(P->minimum < P->kmer){
    fprintf(stderr, "  [x] Error: minimum block size must be >= than k-mer!\n");
    exit(1);
    }

  fprintf(stderr, "\n");
  if(P->verbose){
    fprintf(stderr, "==[ CONFIGURATION ]=================\n");
    fprintf(stderr, "Verbose mode ....................... %s\n", P->verbose == 0
    ? "no" : "yes");
    fprintf(stderr, "Force mode ......................... %s\n", P->force == 0 ?
    "no" : "yes");
    fprintf(stderr, "Delete input file .................. %s\n", delete == 1 ? 
    "no" : "yes");
    fprintf(stderr, "Using inversions ................... %s\n", P->inversion ==
     0 ? "no" : "yes");
    fprintf(stderr, "Minimum block size ................. %u\n", P->minimum);
    //fprintf(stderr, "Number of threads .................. %u\n", P->nThreads);
    fprintf(stderr, "Output reduced filename ............ %s\n", P->positions);
    fprintf(stderr, "\n");
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  ReduceProjections(argv[argc-1], delete);
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");
  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
