#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED

#include "defs.h"
#include "rmodel.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  char *name;
  U64  nReads;
  U32  nSym;
  U64  length;
  }
FILEP;

typedef struct{
  U8       help;
  U8       verbose;
  U8       force;
  U8       inversion;
  U8       level;
  U32      kmer;
  U32      window;
  U32      mutations;
  U32      nThreads;
  FILEP    Ref;
  FILEP    Con;
  char     *positions;
  }
Parameters;

typedef struct{
  uint32_t id;
  uint32_t tar;
  uint32_t ref;
  uint64_t min;
  //RMODEL   *R;
  }
Threads;

Parameters *P;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
