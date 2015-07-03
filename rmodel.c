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
// SILVERS HASH
//
static uint64_t XHASH(uint64_t x){
  return (x * 786433 + 196613) % 68719476735;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATION OF CONTEXT MULTIPLICATOR FOR INDEX FUNCTION USAGE
//
uint64_t CalcMult(uint32_t k){
  uint32_t n;
  uint64_t x[k], p = 1;
  for(n = 0 ; n < k ; ++n){
    x[n] = p;
    p <<= 2;
    }
  return x[k-1];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES THE HASH STRUCTURE 
//
HASH *CreateHash(void){
  HASH *H = (HASH *) Calloc(1, sizeof(HASH));
  H->ent  = (ENTRY   **) Calloc(HSIZE, sizeof(ENTRY  *));
  H->size = (uint16_t *) Calloc(HSIZE, sizeof(uint16_t));
  return H;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES THE RCLASS BASIC STRUCTURE 
//
RCLASS *CreateRC(uint32_t max, uint32_t k, uint8_t ir, uint64_t size){
  uint32_t n;

  RCLASS *C = (RCLASS *) Calloc(1,   sizeof(RCLASS));
  C->RM     = (RMODEL *) Calloc(max, sizeof(RMODEL));
  C->mRM    = max;

  for(n = 0 ; n < max ; ++n){
    C->RM[n].pos    = 0;
    C->RM[n].nHits  = 0;
    C->RM[n].nTries = 0;
    C->RM[n].rev    = ir;
    }

  C->rev    = ir;
  C->idx    = 0;
  C->idxRev = 0;
  C->kmer   = k;
  C->mult   = CalcMult(k);
  C->size   = size;
  return C;
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
// GET REPEAT MODEL HASH ENTRY
//
ENTRY *GetHEnt(HASH *H, uint64_t key){
  uint32_t n, h = (uint32_t) (key % HSIZE);
  uint64_t b = (uint64_t) key & 0xfffffff0000;

  for(n = 0 ; n < H->size[h] ; ++n)
    if(((uint64_t) H->ent[h][n].key | b) == key)
      return &H->ent[h][n];

  return NULL;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START EACH REPEAT MODEL
//
int32_t StartRM(RCLASS *C, HASH *H, uint32_t m, uint64_t i, uint8_t r){
  uint32_t s;
  ENTRY *E;

  if((E = GetHEnt(H, i)) == NULL) return 0;
  if(r == 0) C->RM[m].pos = E->pos[0];
  else{
    if(E->pos[0] <= C->kmer+1) return 0;
    C->RM[m].pos = E->pos[0]-C->kmer-1;
    }

  C->RM[m].nHits  = 0;
  C->RM[m].nTries = 0;
  C->RM[m].rev    = r;

  return 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INSERT KMER POSITION INTO HASH TABLE 
//
void InsertKmerPos(HASH *H, uint64_t key, uint32_t pos){
  uint32_t n, h = (uint32_t) key % HSIZE;
  uint64_t b = key & 0xfffffff0000;

  for(n = 0 ; n < H->size[h] ; ++n)
    if(((uint64_t) H->ent[h][n].key | b) == key){
      H->ent[h][n].pos = (PPR *) Realloc(H->ent[h][n].pos, 
      (H->ent[h][n].nPos + 1) * sizeof(PPR));
      H->ent[h][n].pos[H->ent[h][n].nPos++] = pos;           
      return;
      }

  // CREATE A NEW ENTRY
  H->ent[h] = (ENTRY *) Realloc(H->ent[h], (H->size[h]+1) * sizeof(ENTRY));

  // CREATE A NEW POSITION
  H->ent[h][H->size[h]].pos    = (PPR *) Malloc(sizeof(PPR));
  H->ent[h][H->size[h]].pos[0] = pos;
  H->ent[h][H->size[h]].nPos   = 1;
  H->ent[h][H->size[h]].key    = (uint16_t) (key & 0xffff);
  H->size[h]++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE REPEAT MODEL
//
void UpdateRM(RMODEL *R, uint8_t *b, uint8_t s){
  R->lastHit = 1;
  if(R->rev == 0){
    if(b[R->pos++] == s){
      R->nHits++;
      R->lastHit = 0;
      }
    }
  else{
    if(GetCompNum(b[R->pos--]) == s){
      R->nHits++;
      R->lastHit = 0;
      }
    }
  R->nTries++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE REPEAT MODEL PROBABILITIES
//
/*
void ComputeRMProbs(RCLASS *C, RMODEL *R, uint8_t *b){
  uint8_t n, s;
  s = (R->rev == 1) ? GetCompNum(b[R->pos]) : b[R->pos];
  R->probs[s] = (R->nHits+C->P->alpha) / (R->nTries+2*C->P->alpha);
  for(n = 0 ; n < NSYM ; ++n)
    if(n != s) R->probs[n] = (1-R->probs[s])/3;
  }
*/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
void StopRM(RCLASS *C){
  uint32_t n;

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
// START NEW REPEAT MODELS IF THERE IS STILL SPACE
//                         
void StartMultipleRMs(RCLASS *C, HASH *H, uint8_t *b){
  if(C->nRM < C->mRM && StartRM(C, H, C->nRM, GetIdxR(b, C), 0))
    C->nRM++;

  if(C->rev == 1 && C->nRM < C->mRM && StartRM(C, H, C->nRM, GetIdxRevR(b, C), 1))
    C->nRM++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

