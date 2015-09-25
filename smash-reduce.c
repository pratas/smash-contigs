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
#include "seq.h"
#include "pos.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "buffer.h"
#include "common.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E D U C E - - - - - - - - - - - - - - - -
void ReduceProjections(char *fn){
  FILE *IN = NULL, *OUT = NULL;
  char name[MAX_FILENAME], watermark[MAX_FILENAME];
  sprintf(name, "%s.red", fn);
  int64_t refNBases, conNBases, ri, rf, ci, cf, cx, cy, rx, ry;

  IN  = Fopen(fn,   "r");
  OUT = Fopen(name, "w");

  fprintf(stderr, "Reducing projections ... \n");

/*
  int64_t posCache[MAX_POS_CACHE][4];
  uint8_t PosUsage[MAX_POS_CACHE];
  int64_t idx = 0;
*/

  // READ HEADER
  if(fscanf(IN, "%s\t%"PRIi64"\t%"PRIi64"\n", watermark, &conNBases,
  &refNBases) != 3 || watermark[0] != '#' || watermark[1] != 'S' ||
  watermark[2] != 'C' || watermark[3] != 'F'){
    fprintf(stderr, "[x] Error: unknown positions file format!\n");
    exit(1);
    }
  fprintf(OUT, "#SCF\t%"PRIi64"\t%"PRIi64"\n", conNBases, refNBases);

  // READ BODY
  while(1){
    char tmp1[MAX_STR] = {'\0'}, tmp2[MAX_STR] = {'\0'};
    if(fscanf(IN, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                  "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n", 
                  tmp1, &ci, &cf, &cx, &cy, tmp2, &ri, &rf, &rx, &ry) != 10)
      break;
/*
    if(cf > ci){
      posCache[idx][0] = ci;
      posCache[idx][1] = cf;
      posCache[idx][2] = ri;
      posCache[idx][3] = rf;

      if(++idx == MAX_POS_CACHE)
        idx = 0;
      }
    else{ // INVERTED

      }

*/ 

    fprintf(OUT, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                 "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n", 
                 tmp1, ci, cf, cx, cy, tmp2, ri, rf, rx, ry);
    }
  fprintf(stderr, "\rDone!                   \n");

  //unlink(fn);
  fclose(IN);
  fclose(OUT);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t main(int argc, char *argv[]){
  char **p = *&argv;

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
  P->minimum    = ArgsNum   (DEF_MINI,    p, argc, "-m", MIN_MINI, MAX_MINI);
  P->threshold  = ArgsNum   (DEF_TSHO,    p, argc, "-t", MIN_TSHO, MAX_TSHO);
  P->nThreads   = ArgsNum   (DEF_THRE,    p, argc, "-n", MIN_THRE, MAX_THRE);
  P->positions  = ArgsFilesReduce        (p, argc, "-o");

  if(P->minimum < P->kmer){
    fprintf(stderr, "  [x] Error: minimum block size must be >= than k-mer!\n");
    exit(1);
    }

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  ReduceProjections(argv[argc-1]);
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");
  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
