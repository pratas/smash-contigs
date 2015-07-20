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
  FILE        *Reader  = Fopen(P->Con.name, "r");
  char        name[MAX_FILENAME];
  sprintf(name, ".t%u", T.id+1);
  FILE        *Writter = Fopen(concatenate(P->positions, name), "w");
  uint64_t    nBaseRelative = 0, nBaseAbsolute = 0, nNRelative = 0, idxPos = 0;
  uint32_t    k, r = 0;
  int32_t     action;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t)), sym,
              contigName[MAX_CONTIG_NAME];
  RCLASS      *Mod = CreateRClass(P->repeats, P->editions, P->minimum, P->kmer,
              P->inversion);

  Mod->nBases = P->Ref.nBases;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      sym = readBuf[idxPos];
      if((action = ParseSym(PA, sym)) < 0){
        switch(action){
          case -2:
            contigName[r] = '\0';
            if(Mod->nRM > 0) 
              ResetAllRMs(Mod, Head, nBaseRelative, contigName, Writter);
            nNRelative = 0;
            nBaseRelative = 0;
            r = 0;
          break;
          case -3:
            if(r >= MAX_CONTIG_NAME - 1)
              contigName[r] = '\0';
            else{ 
              if(sym == ' ' && r == 0) 
                continue;
              contigName[r++] = (uint8_t) sym;        
              }
          break;
          }
        continue;
        }
fprintf(stderr, "%c", NumToDNASym(Seq->buf[nBaseAbsolute]));

      if((sym = DNASymToNum(sym)) == 4){
        if(Mod->nRM > 0) // PROTECT Ns IN THE CONTIG SEQUENCE
          ResetAllRMs(Mod, Head, nBaseRelative, contigName, Writter);
        ++nNRelative;
        ++nBaseRelative;
        ++nBaseAbsolute;
        continue;
        }

      
      symBuf->buf[symBuf->idx] = sym;
      GetIdxRM   (symBuf->buf+symBuf->idx-1, Mod);
      GetIdxRevRM(symBuf->buf+symBuf->idx-1, Mod);
//          fprintf(stderr, "%"PRIu64"\n", Mod->idx);

      if(PA->nRead % P->nThreads == T.id){
        if(nBaseRelative > Mod->kmer){  // PROTECTING THE BEGGINING OF K-SIZE
          UpdateRMs(Mod, Seq->buf, nBaseRelative, sym);
          StopRMs(Mod, Head, nBaseRelative, contigName, Writter);
          StartMultipleRMs(Mod, Hash, nBaseRelative);
          }
        }

      UpdateCBuffer(symBuf);
      ++nBaseRelative;
      ++nBaseAbsolute;
      }

  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Writter);
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
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  sym, *readBuf;
  int32_t  action;
  uint64_t r = 0;
  struct   stat s;
  size_t   size, k;
  long     fd = open(P->Ref.name, O_RDONLY);
  RCLASS   *Mod = CreateRClass(P->repeats, P->editions, P->minimum, P->kmer,
           P->inversion);

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
        break;
        case -3:
          if(r >= MAX_CONTIG_NAME - 1)
            Head->Pos[Head->iPos-1].name[r] = '\0';
          else{
            if(sym == ' ' && r == 0)
              continue;
            Head->Pos[Head->iPos-1].name[r++] = (uint8_t) sym;
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

//fprintf(stderr, "%"PRIu64"\n", Mod->idx); // XXX: COMPARE

      //TODO: CONDITION TO LOAD KMER AFTER nBASES & FOR EACH READ RESET IDX
        InsertKmerPos(Hash, Mod->idx, k);
      UpdateCBuffer(symBuf);
      }

//    CalcProgress(P->Ref.length, k);
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
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t main(int argc, char *argv[]){
  char **p = *&argv;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEF_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
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
  P->editions   = ArgsNum   (DEF_EDIT,    p, argc, "-e", MIN_EDIT, MAX_EDIT);
  P->nThreads   = ArgsNum   (DEF_THRE,    p, argc, "-n", MIN_THRE, MAX_THRE);
  P->positions  = ArgsFiles              (p, argc, "-o");
  P->Con.name   = argv[argc-2];
  P->Ref.name   = argv[argc-1];
  P->Con.length = FopenBytesInFile(P->Con.name); 
  P->Ref.length = FopenBytesInFile(P->Ref.name); 
  P->window     = P->kmer;

  if(P->minimum < P->kmer){
    fprintf(stderr, "  [x] Error: minimum block size must be >= than k-mer!\n");
    exit(1);
    }

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  CompressAction();
  // TODO: ReduceProjections();
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
