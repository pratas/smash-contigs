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
#include "rmodel.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E D U C E - - - - - - - - - - - - - - - -
void ReduceProjections(Threads T){
  FILE *IN = NULL, *OUT = NULL;
  char name[MAX_FILENAME], nameCat[MAX_FILENAME];
  sprintf(name, "%s.t%u", P->positions, T.id+1);
  sprintf(nameCat, "%s.cat", name);
  int64_t ri, rf, ci, cf, cx, cy, rx, ry;

  IN  = Fopen(name, "r");
  OUT = Fopen(nameCat, "w");

/*
  int64_t posCache[MAX_POS_CACHE][4];
  uint8_t PosUsage[MAX_POS_CACHE];
  int64_t idx = 0;
*/

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

    //------------------------------------------------------------------------
*/ 

    fprintf(OUT, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                 "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n", 
                 tmp1, ci, cf, cx, cy, tmp2, ri, rf, rx, ry);
    }

  unlink(name);
  fclose(IN);
  fclose(OUT);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - C   T H R E A D I N G - - - - - - - - - - - - - - -
void *ProjectionsThread(void *Thr){
  Threads *T = (Threads *) Thr;
  ReduceProjections(T[0]);
  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -
void ReduceAction(){
  uint32_t n;
  pthread_t t[P->nThreads];
  Threads  *T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(n = 0 ; n < P->nThreads ; ++n) T[n].id = n; 

  fprintf(stderr, "  [+] Reduce projections ... \n");
  for(n = 0 ; n < P->nThreads ; ++n)
    pthread_create(&(t[n+1]), NULL, ProjectionsThread, (void *) &(T[n]));
  for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
    pthread_join(t[n+1], NULL);
  fprintf(stderr, "\r      Done!                   \n");
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - J O I N E R - - - - - - - - - - - - - - - -
//
// JOINNING THREADS FUNCTION
// IT ALSO ADDS IN THE FIRST LINE A WATERMARK AND THE NUMBER OF BASES FROM REF 
// AND CONTIGS FILES
// 
void ThreadConcatenation(void){
  FILE *OUT = NULL;
  uint32_t n, k;
  uint8_t *buf = (uint8_t *) Malloc(BUFFER_SIZE * sizeof(uint8_t));

  fprintf(stderr, "  [+] Joinning thread files ...\n");

  OUT = Fopen(P->positions, "w");
  fprintf(OUT, "#SCF\t%"PRIi64"\t%"PRIi64"\n", P->Con.nBases, P->Ref.nBases);

  for(n = 0 ; n < P->nThreads ; ++n){
    char tmp[MAX_FILENAME];
    sprintf(tmp, "%s.t%u.cat", P->positions, n+1);
    FILE *IN = Fopen(tmp, "r");
    while((k = fread(buf, 1, BUFFER_SIZE, IN)))
      fwrite(buf, 1, k, OUT);
    fclose(IN);
    unlink(tmp);
    }
  fclose(OUT);

  fprintf(stderr, "      Done!                \n");
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
  P->positions  = ArgsFiles              (p, argc, "-o");
  P->Con.name   = argv[argc-2];
  P->Ref.name   = argv[argc-1];
  P->Con.length = FopenBytesInFile(P->Con.name); 
  P->Ref.length = FopenBytesInFile(P->Ref.name); 

  if(P->minimum < P->kmer){
    fprintf(stderr, "  [x] Error: minimum block size must be >= than k-mer!\n");
    exit(1);
    }

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  ReduceAction();
  ThreadConcatenation();

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
