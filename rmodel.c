#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "rmodel.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATION OF CONTEXT MULTIPLICATOR FOR INDEX FUNCTION USAGE
//
uint64_t CalcMult(uint32_t c){
  uint32_t n;
  uint64_t x[c], p = 1;
  for(n = 0 ; n < c ; ++n){
    x[n] = p;
    p <<= 2;
    }
  return x[c-1];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES THE RCLASS BASIC STRUCTURE 
//
RCLASS *CreateRClass(uint32_t max, uint32_t k, uint8_t ir){
  uint32_t n;

  RCLASS *C = (RCLASS *)  Calloc(1,   sizeof(RCLASS));
  C->RM     = (RMODEL *)  Calloc(max, sizeof(RMODEL));
  C->active = (uint8_t *) Calloc(max, sizeof(uint8_t));
  C->nRM    = 0;
  C->mRM    = max;
  C->rev    = ir;
  C->idx    = 0;
  C->idxRev = 0;
  C->kmer   = k;
  C->mult   = CalcMult(k);
  for(n = 0 ; n < max ; ++n){
    C->RM[n].pos     = 0;
    C->RM[n].nHits   = 0;
    C->RM[n].nTries  = 0;
    C->RM[n].winSize = C->kmer;
    C->RM[n].win     = (uint8_t *) Calloc(C->kmer + 1, sizeof(uint8_t *));
    C->RM[n].rev     = ir;
    }

  return C;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE PERMANENTLY RCLASS 
//
void RemoveRClass(RCLASS *C){
  uint32_t n;

  for(n = 0 ; n < C->mRM ; ++n){
    Free(C->RM[n].win);
    }
  Free(C->RM);
  Free(C->active);
  Free(C);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REVERSE COMPLEMENT INDEX BASED ON PAST SYMBOLS FOR REPEATS
//
uint64_t GetIdxRevR(uint8_t *p, RCLASS *C){
  return (C->idxRev = (C->idxRev>>2)+GetCompNum(*p)*C->mult);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR REPEATS
//
uint64_t GetIdxR(uint8_t *p, RCLASS *C){
  return (C->idx = ((C->idx-*(p-C->kmer)*C->mult)<<2)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START EACH REPEAT MODEL
//
int32_t StartRM(RCLASS *C, HASH *H, uint32_t m, uint64_t i, uint8_t r){
  uint32_t s;
  ENTRY *E;

  if((E = GetHEnt(H, i)) == NULL) 
    return 0;

  if(r == 0){ 
    C->RM[m].pos = E->pos[0];
    }
  else{
    if(E->pos[0] <= C->kmer + 1) 
      return 0;
    C->RM[m].pos = E->pos[0]-C->kmer-1;
    }

  C->RM[m].nHits   = 0;
  C->RM[m].nTries  = 0;
  C->RM[m].rev     = r;
  C->RM[m].winSize = C->kmer;
  C->RM[m].win     = (uint8_t *) Calloc(C->kmer+1, sizeof(uint8_t));

  return 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE REPEAT MODEL
//
void UpdateRM(RMODEL *R, uint8_t *b, uint8_t s){
  //R->lastHit = 1;
  if(R->rev == 0){
    if(b[R->pos++] == s){
      R->nHits++;
      // R->lastHit = 0;
      }
    }
  else{
    if(GetCompNum(b[R->pos--]) == s){
      R->nHits++;
      // R->lastHit = 0;
      }
    }
  R->nTries++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
void StopRM(RCLASS *C, uint64_t iBase, FILE *Writter){
  uint32_t n;

/*
  s = (R->rev == 1) ? GetCompNum(b[R->pos]) : b[R->pos];
  R->probs[s] = (R->nHits+C->P->alpha) / (R->nTries+2*C->P->alpha);
  for(n = 0 ; n < NSYM ; ++n)
    if(n != s) R->probs[n] = (1-R->probs[s])/3;
*/


  for(n = 0 ; n < C->nRM ; ++n){
    }


  for(;;){
    for(n = 0 ; n < C->nRM ; ++n){
    //  if((C->RM[n].acting = C->P->beta * C->RM[n].acting + C->RM[n].lastHit) >
    //  C->P->limit * 1.0 || C->RM[n].pos == 0)
    //    {
    //    if(n != C->nRM-1)
    //      C->RM[n] = C->RM[C->nRM-1];
    //    C->nRM--;
        return;
    //    }
      }
    return;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RESET RMODEL
// 
void ResetRM(RMODEL *R){
  memset(R->win, 0, R->winSize);
  R->pos    = 0;
  R->size   = 0;
  R->nHits  = 0;
  R->nTries = 0;
  R->rev    = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE RM
//
void RemoveRM(RMODEL *R){
  Free(R->win);
  Free(R);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE ALL RM
//
void RemoveAllRM(RCLASS *C){
  uint32_t n;
  for(n = 0 ; n < C->nRM ; ++n){
    Free(&C->RM[n]);
    }
  Free(C->RM);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FORCE STOP REPEAT MODELS DURING END OF READ
//
void ResetAllRM(RCLASS *C, uint64_t iBase, FILE *Writter){
  uint32_t n;
  for(n = 0 ; n < C->nRM ; ++n){
    fprintf(Writter, "%"PRIu64"\t%"PRIu64"\n", iBase, 0);
    ResetRM(&C->RM[n]);
    }
  C->nRM = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START NEW REPEAT MODELS IF THERE IS STILL SPACE
//                         
void StartMultipleRMs(RCLASS *C, HASH *H, uint8_t *b){
  if(C->nRM < C->mRM && StartRM(C, H, C->nRM, GetIdxR(b, C), 0))
    C->nRM++;

  if(C->rev == 1 && C->nRM < C->mRM && StartRM(C, H, C->nRM, GetIdxRevR(b, C), 1))
    C->nRM++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

