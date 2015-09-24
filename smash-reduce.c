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

HASH     *Hash; // HASH MEMORY SHARED BY THREADING
SEQ      *Seq;  // SEQUENCE SHARED BY THREADING
HEADERS  *Head; // HEADERS SHARED BY THREADING

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I N G - - - - - - - - - - - - - 
void CompressTarget(Threads T){
  FILE      *Reader = Fopen(P->Con.name, "r");
  char      name[MAX_FILENAME];
  sprintf(name, ".t%u", T.id+1);
  FILE      *Writter = Fopen(concatenate(P->positions, name), "w");
  int64_t   nBaseRelative = 0, nBaseAbsolute = 0, idxPos = 0;
  uint32_t  k, r = 0;
  int32_t   action;
  PARSER    *PA = CreateParser();
  CBUF      *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t   *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t)), sym,
            *conName = (uint8_t *) Calloc(MAX_CONTIG_NAME, sizeof(uint8_t));
  RCLASS    *Mod = CreateRClass(P->repeats, P->minimum, P->kmer, P->inversion);

  Mod->nBases = P->Ref.nBases;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      sym = readBuf[idxPos];
      if((action = ParseSym(PA, sym)) < 0){
        switch(action){
          case -1: // IT IS THE BEGGINING OF THE HEADER
            if(PA->nRead > 1 && Mod->nRM > 0)
              ResetAllRMs(Mod, Head, nBaseRelative, nBaseAbsolute, conName, 
              Writter);
            nBaseRelative = 0;
            r = 0;
          break;
          case -2: // IT IS THE '\n' HEADER END
            conName[r] = '\0';
          break;
          case -3: // IF IS A SYMBOL OF THE HEADER
            if(r >= MAX_CONTIG_NAME-1)
              conName[r] = '\0';
            else{ 
              if(sym == ' ' && r == 0) continue;
              conName[r++] = sym;        
              }
          break;
          case -99: // IF IS A SIMPLE FORMAT BREAK
          break;
          default:
            fprintf(stderr, "ERROR: Unknown action!\n");
            exit(1);
          }
        continue; // GO TO NEXT SYMBOL
        }

      if((sym = DNASymToNum(sym)) == 4){
        if(Mod->nRM > 0 && PA->nRead % P->nThreads == T.id) 
          ResetAllRMs(Mod, Head, nBaseRelative, nBaseAbsolute, conName, 
          Writter);
        ++nBaseRelative;
        ++nBaseAbsolute;
        continue;
        }
      
      symBuf->buf[symBuf->idx] = sym;
      GetIdxRM   (symBuf->buf+symBuf->idx, Mod);
      GetIdxRevRM(symBuf->buf+symBuf->idx, Mod);

      if(PA->nRead % P->nThreads == T.id){
        if(nBaseRelative >= Mod->kmer){  // PROTECTING THE BEGGINING OF K-SIZE
          UpdateRMs(Mod, Seq->buf, nBaseRelative, sym);
          StopRMs(Mod, Head, nBaseRelative, nBaseAbsolute, conName, Writter);
          StartMultipleRMs(Mod, Hash, nBaseRelative, nBaseAbsolute);
          }
        }

      UpdateCBuffer(symBuf);
      ++nBaseRelative;
      ++nBaseAbsolute;
      }

  P->Con.nBases = nBaseAbsolute;
  Free(readBuf);
  Free(conName);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  RemoveRClass(Mod);
  fclose(Writter);
  fclose(Reader);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - C   T H R E A D I N G - - - - - - - - - - - - - - -
void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -
void LoadReference(){
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  sym, *readBuf;
  int32_t  action;
  uint64_t r = 0, nBaseRelative = 0;
  struct   stat s;
  size_t   size, k;
  long     fd = open(P->Ref.name, O_RDONLY);
  RCLASS   *Mod = CreateRClass(P->repeats, P->minimum, P->kmer, P->inversion);

  Head->Pos[0].init = 0;
  Mod->nBases = 0;
  fstat (fd, & s);
  size = s.st_size;
  readBuf = (uint8_t *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  for(k = 0 ; k < size ; ++k){
    sym = *readBuf++;
    if((action = ParseSym(PA, sym)) < 0){
      switch(action){
        case -1:
          UpdateHeaders(Head);
          if(Head->iPos != 1){
            Head->Pos[Head->iPos-2].end  = Mod->nBases;
            Head->Pos[Head->iPos-1].init = Mod->nBases + 1;
            }
        break;
        case -2:
          Head->Pos[Head->iPos-1].name[r] = '\0';
          r = 0;
          nBaseRelative = 0;
          Mod->idx = 0;
        break;
        case -3:
          if(r >= MAX_CONTIG_NAME - 1)
            Head->Pos[Head->iPos-1].name[r] = '\0';
          else{
            if(sym == ' ' && r == 0)
              continue;
            Head->Pos[Head->iPos-1].name[r++] = sym;
            }
        break;
        }
      continue; // CASE -99
      }

    sym = DNASymToNum(sym);
    UpdateSeq(Seq, sym);

    if(sym != 4){
      symBuf->buf[symBuf->idx] = sym;
      Mod->idx = GetIdxRM(symBuf->buf+symBuf->idx-1, Mod);
      if(nBaseRelative >= Mod->kmer)
        InsertKmerPos(Hash, Mod->idx, Mod->nBases);
      UpdateCBuffer(symBuf);
      }

    CalcProgress(P->Ref.length, k);
    ++nBaseRelative;
    Mod->nBases++;
    }

  Head->Pos[Head->iPos-1].end = Mod->nBases;
  P->Ref.nBases = Mod->nBases;
  RemoveCBuffer(symBuf);
  RemoveRClass(Mod);
  RemoveParser(PA);
  close(fd);
  }

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
void CompressAction(){
  uint32_t n;
  pthread_t t[P->nThreads];
  Threads  *T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(n = 0 ; n < P->nThreads ; ++n) T[n].id = n; 

  fprintf(stderr, "  [+] Building models ...\n");
  Seq  = CreateSeq(100000);
  Hash = CreateHash();
  Head = CreateHeaders();
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Loading reference ...\n");
  LoadReference();
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Compressing contigs ... \n");
  fprintf(stderr, "      (this may take a while) ");
  for(n = 0 ; n < P->nThreads ; ++n)
    pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[n]));
  for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
    pthread_join(t[n+1], NULL);
  fprintf(stderr, "\r      Done!                   \n");

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
  P->kmer       = ArgsNum   (DEF_KMER,    p, argc, "-k", MIN_KMER, MAX_KMER);
  P->minimum    = ArgsNum   (DEF_MINI,    p, argc, "-m", MIN_MINI, MAX_MINI);
  P->repeats    = ArgsNum   (DEF_REPE,    p, argc, "-r", MIN_REPE, MAX_REPE);
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
  CompressAction();
  ThreadConcatenation();

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveSeq(Seq);
  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
