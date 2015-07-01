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

/*
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I N G - - - - - - - - - - - - - 

void CompressTarget(Threads T){
  FILE        *Reader  = Fopen(P->files[T.id], "r");
  double      *cModelWeight, cModelTotalWeight = 0, bits = 0, instance = 0;
  uint64_t    nBase = 0;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos;
  PModel      **pModel, *MX;
  CModel      **Shadow;
  FloatPModel *PT;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]); 
  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  cModelWeight  = (double   *) Calloc(totModels, sizeof(double));

  for(n = 0 ; n < totModels ; ++n)
    cModelWeight[n] = 1.0 / totModels;

  FileType(PA, Reader);

  nBase = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1) continue;
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        GetPModelIdx(pos, Shadow[cModel]);
        ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->pModelIdx,
        Shadow[cModel]->alphaDen);
        ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
        if(Shadow[cModel]->edits != 0){
          ++n;
          Shadow[cModel]->SUBS.seq->buf[Shadow[cModel]->SUBS.seq->idx] = sym;
          Shadow[cModel]->SUBS.idx = GetPModelIdxCorr(Shadow[cModel]->SUBS.
          seq->buf+Shadow[cModel]->SUBS.seq->idx-1, Shadow[cModel], Shadow
          [cModel]->SUBS.idx);
          ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->SUBS.idx, 
          Shadow[cModel]->SUBS.eDen);
          ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
          }
        ++n;
        }

      MX->sum  = (MX->freqs[0] = 1 + (unsigned) (PT->freqs[0] * MX_PMODEL));
      MX->sum += (MX->freqs[1] = 1 + (unsigned) (PT->freqs[1] * MX_PMODEL));
      MX->sum += (MX->freqs[2] = 1 + (unsigned) (PT->freqs[2] * MX_PMODEL));
      MX->sum += (MX->freqs[3] = 1 + (unsigned) (PT->freqs[3] * MX_PMODEL));
      bits += (instance = PModelSymbolLog(MX, sym));
      nBase++;

      cModelTotalWeight = 0;
      for(n = 0 ; n < totModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        }

      for(n = 0 ; n < totModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(Shadow[cModel]->edits != 0){
          CorrectCModelSUBS(Shadow[cModel], pModel[++n], sym);
          }
        ++n;
        }

      UpdateCBuffer(symBuf);
      }

  Free(cModelWeight);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);

  P->matrix[P->ref][T.id] = nBase == 0 ? 101 : bits / 2 / nBase; // 101 -> nan
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }

*/


//    // COUNT READS
//    if((sym == '@' || sym == '>') && Seq->buf[Seq->idx] == '\n')
//      P->Ref.nReads++;


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

    symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);
    UpdateSeq(Seq, sym);


    Mod->P->idx = GetIdxR(Seq->buf+Seq->idx, Mod);
    InsertKmerPos(Mod, Mod->P->idx, k);
    CalcProgress(P->Ref.length, k);

    ++nBases;
    UpdateCBuffer(symBuf);
    }

  P->Ref.nBases = nBases;
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  close(fd);
  }

/*
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void CompressAction(Threads *T, uint32_t ref){
  uint32_t n, k;
  pthread_t t[P->nThreads];
  P->ref = ref;

  Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Models[n] = CreateCModel(T[ref].model[n].ctx, T[ref].model[n].den, 
    T[ref].model[n].ir, REFERENCE, P->col, T[ref].model[n].edits, 
    T[ref].model[n].eDen);

  fprintf(stderr, "  [+] Loading reference %u ... ", ref+1);
  LoadReference(T[ref]);
  fprintf(stderr, "Done!\n");
  
  fprintf(stderr, "      [+] Compressing %u targets ... ", P->nFiles);
  ref = 0;
  do{
    for(n = 0 ; n < P->nThreads ; ++n)
      pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[ref+n]));
    for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[n+1], NULL);
    }
  while((ref += P->nThreads) < P->nFiles && ref + P->nThreads <= P->nFiles);

  if(ref < P->nFiles){ // EXTRA - OUT OF THE MAIN LOOP
    for(n = ref, k = 0 ; n < P->nFiles ; ++n)
      pthread_create(&(t[++k]), NULL, CompressThread, (void *) &(T[n]));
    for(n = ref, k = 0 ; n < P->nFiles ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[++k], NULL);
    }
  fprintf(stderr, "Done!\n");

  for(n = 0 ; n < P->nModels ; ++n)
    FreeCModel(Models[n]);
  Free(Models);
  }

*/


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
  P->mrepeats   = ArgsNum    (DEF_REPE,    p, argc, "-r", MIN_REPE, MAX_REPE);
  P->window     = ArgsNum    (DEF_WIND,    p, argc, "-w", MIN_WIND, MAX_WIND);
  P->mutations  = ArgsNum    (DEF_MUTA,    p, argc, "-m", MIN_MUTA, MAX_MUTA);
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

  fprintf(stderr, "Building repeats model and sequence ...\n");
  Seq = CreateSeq(1000); 
  Mod = CreateRC(500, 1, 0.9, 7, P->kmer, 0.9, P->inversion);
  fprintf(stderr, "Done!                \n");
  fprintf(stderr, "Loading reference ...\n");
  LoadReference();
  fprintf(stderr, "Done!                \n");


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
