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
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "buffer.h"
#include "common.h"
#include "rmodel.h"

RCLASS  *Mod;  // MEMORY MODEL SHARED BY THREADING
SEQ     *Seq;  // SEQUENCE SHARED BY THREADING

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I N G - - - - - - - - - - - - - 

void CompressTarget(Threads T){
  FILE        *Reader  = Fopen(P->Con.name, "r");
  uint64_t    nBase = 0, idxPos = 0;
  uint32_t    n, k;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t)), sym;

  FileType(PA, Reader);

  nBase = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1) continue;
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);
      n = 0;
     
      // TODO: NEED SHADOWS?

      StopRM(Mod);
      StartMultipleRMs(Mod, symBuf->buf+symBuf->idx-1);

      printf("%u : %u\n", Mod->nRM, Mod->mRM);

/*
      //RenormWeights(Mod, F);
      ComputeMixture(Mod, F, MX, buf);
      UpdateWeights(C, F, buf, sym);
*/

      UpdateCBuffer(symBuf);
      nBase++;
      }

  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReference(){
  FILE     *Reader = Fopen(P->Ref.name, "r");
  uint32_t n;
  uint64_t nBases = 0, nReads = 0, idx = 0;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  sym, irSym, *readBuf;
  FileType(PA, Reader);
  fclose(Reader);
  struct   stat s;
  size_t   size, k;
  long     fd = open(P->Ref.name, O_RDONLY);

  fstat (fd, & s);
  size = s.st_size;
  readBuf = (uint8_t *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  for(k = 0 ; k < size ; ++k){

    if(ParseSym(PA, (sym = *readBuf++)) == -1) continue;
    sym = DNASymToNum(sym);
    UpdateSeq(Seq, sym);

    if(sym != 4){
      symBuf->buf[symBuf->idx] = sym;
      Mod->P->idx = GetIdxR(symBuf->buf+symBuf->idx-1, Mod);
      InsertKmerPos(Mod, Mod->P->idx, k);
      UpdateCBuffer(symBuf);
      }

    CalcProgress(P->Ref.length, k);
    ++nBases;
    }

  P->Ref.nBases = nBases;
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  close(fd);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void CompressAction(){
  uint32_t n;
  pthread_t t[P->nThreads];
  Threads  *T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(n = 0 ; n < P->nThreads ; ++n) T[n].id = n; 

  fprintf(stderr, "  [+] Building models ...\n");
  Seq = CreateSeq(100000);
  Mod = CreateRC(P->repeats, 1, 0.9, 7, P->kmer, 0.9, P->inversion);
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Loading reference ...\n");
  LoadReference();
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Compressing contigs ... \n");
  for(n = 0 ; n < P->nThreads ; ++n)
    pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[n]));
  for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
    pthread_join(t[n+1], NULL);
  fprintf(stderr, "      Done!\n");

  // TODO: FREE REPEAT MODELS!
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  int32_t     xargc = 0;
  uint32_t    n;
  Threads     *T;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEF_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  P->verbose    = ArgsState  (DEF_VERBOSE, p, argc, "-v" );
  P->force      = ArgsState  (DEF_FORCE,   p, argc, "-F" );
  P->inversion  = ArgsState  (DEF_INVE,    p, argc, "-i" );
  P->kmer       = ArgsNum    (DEF_KMER,    p, argc, "-k", MIN_KMER, MAX_KMER);
  P->minimum    = ArgsNum    (DEF_MINI,    p, argc, "-m", MIN_MINI, MAX_MINI);
  P->repeats    = ArgsNum    (DEF_REPE,    p, argc, "-r", MIN_REPE, MAX_REPE);
  P->window     = ArgsNum    (DEF_WIND,    p, argc, "-w", MIN_WIND, MAX_WIND);
  P->editions   = ArgsNum    (DEF_EDIT,    p, argc, "-e", MIN_EDIT, MAX_EDIT);
  P->nThreads   = ArgsNum    (DEF_THRE,    p, argc, "-n", MIN_THRE, MAX_THRE);
  P->positions  = ArgsFiles  (p, argc, "-o");
  P->Con.name   = argv[argc-2];
  P->Ref.name   = argv[argc-1];
  P->Con.length = FopenBytesInFile(P->Con.name); 
  P->Ref.length = FopenBytesInFile(P->Ref.name); 

  // GET NUMBER OF READS AND NUMBER OF SYMBOLS FOR REFERENCE AND TARGET
  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  CompressAction();
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  fprintf(stderr, "4 positions...\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveSeq(Seq);
  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
