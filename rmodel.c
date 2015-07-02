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
RCLASS *CreateRC(uint32_t m, double a, double b, uint32_t l, uint32_t c,
double g, uint8_t i){
  RCLASS *C     = (RCLASS   *) Calloc(1,     sizeof(RCLASS  ));
  C->hash       = (HASH     *) Calloc(1,     sizeof(HASH    ));
  C->P          = (RPARAM   *) Calloc(1,     sizeof(RPARAM  ));
  C->hash->ent  = (ENTRY   **) Calloc(HSIZE, sizeof(ENTRY  *));
  C->hash->size = (uint16_t *) Calloc(HSIZE, sizeof(uint16_t));
  C->RM         = (RMODEL   *) Calloc(m,     sizeof(RMODEL  ));
  C->mRM        = m;
  C->P->rev     = i;
  C->P->alpha   = ((int)(a*65535))/65535.0;
  C->P->beta    = ((int)(b*65535))/65535.0;
  C->P->gamma   = ((int)(g*65535))/65535.0;
  C->P->limit   = l;
  C->P->ctx     = c;
  C->P->mult    = CalcMult(c);
  C->P->idx     = 0;
  C->P->idxRev  = 0;
  return C;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REVERSE COMPLEMENT INDEX BASED ON PAST SYMBOLS FOR REPEATS
//
uint64_t GetIdxRevR(uint8_t *p, RCLASS *C){
  return (C->P->idxRev = (C->P->idxRev>>2)+GetCompNum(*p)*C->P->mult);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR REPEATS
//
uint64_t GetIdxR(uint8_t *p, RCLASS *C){
  return (C->P->idx = ((C->P->idx-*(p-C->P->ctx)*C->P->mult)<<2)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET REPEAT MODEL HASH ENTRY
//
ENTRY *GetHEnt(RCLASS *C, uint64_t key){
  uint32_t n, h = (uint32_t) (key % HSIZE);
  uint64_t b = (uint64_t) key & 0xfffffff0000;

  for(n = 0 ; n < C->hash->size[h] ; ++n)
    if(((uint64_t) C->hash->ent[h][n].key | b) == key)
      return &C->hash->ent[h][n];

  return NULL;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START EACH REPEAT MODEL
//
int32_t StartRM(RCLASS *C, uint32_t m, uint64_t i, uint8_t r){
  uint32_t s;
  ENTRY *E;

  if((E = GetHEnt(C, i)) == NULL) return 0;
  if(r == 0) C->RM[m].pos = E->pos[0];
  else{
    if(E->pos[0] <= C->P->ctx+1) return 0;
    C->RM[m].pos = E->pos[0]-C->P->ctx-1;
    }

  C->RM[m].nHits  = 0;
  C->RM[m].nTries = 0;
  C->RM[m].rev    = r;
  C->RM[m].acting = 0;
  C->RM[m].weight = INIWEIGHT;
  for(s = 0 ; s < NSYM ; ++s)
    C->RM[m].probs[s] = 0;

  return 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INSERT KMER POSITION INTO HASH TABLE 
//
void InsertKmerPos(RCLASS *C, uint64_t key, uint32_t pos){
  uint32_t n, h = (uint32_t) key % HSIZE;
  uint64_t b = key & 0xfffffff0000;

  for(n = 0 ; n < C->hash->size[h] ; ++n)
    if(((uint64_t) C->hash->ent[h][n].key | b) == key){
      C->hash->ent[h][n].pos = (PPR *) Realloc(C->hash->ent[h][n].pos, 
      (C->hash->ent[h][n].nPos + 1) * sizeof(PPR));
      C->hash->ent[h][n].pos[C->hash->ent[h][n].nPos++] = pos;           
      return;
      }

  // CREATE A NEW ENTRY
  C->hash->ent[h] = (ENTRY *) Realloc(C->hash->ent[h], (C->hash->size[h]+1) *
  sizeof(ENTRY));

  // CREATE A NEW POSITION
  C->hash->ent[h][C->hash->size[h]].pos    = (PPR *) Malloc(sizeof(PPR));
  C->hash->ent[h][C->hash->size[h]].pos[0] = pos;
  C->hash->ent[h][C->hash->size[h]].nPos   = 1;
  C->hash->ent[h][C->hash->size[h]].key    = (uint16_t) (key & 0xffff);
  C->hash->size[h]++;
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
void ComputeRMProbs(RCLASS *C, RMODEL *R, uint8_t *b){
  uint8_t n, s;
  s = (R->rev == 1) ? GetCompNum(b[R->pos]) : b[R->pos];
  R->probs[s] = (R->nHits+C->P->alpha) / (R->nTries+2*C->P->alpha);
  for(n = 0 ; n < NSYM ; ++n)
    if(n != s) R->probs[n] = (1-R->probs[s])/3;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
void StopRM(RCLASS *C){
  uint32_t n;
  uint8_t  a;
  do{
    a = 0;
    for(n = 0 ; n < C->nRM ; ++n)
      if((C->RM[n].acting = C->P->beta * C->RM[n].acting + C->RM[n].lastHit) >
      C->P->limit * 1.0 || C->RM[n].pos == 0)
        {
        if(n != C->nRM-1)
          C->RM[n] = C->RM[C->nRM-1];
        C->nRM--;
        a = 1;
        break;
        }
    }
  while(a); //for(;;) & return?
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START NEW REPEAT MODELS IF THERE IS STILL SPACE
//                         
void StartMultipleRMs(RCLASS *C, uint8_t *b){
  if(C->nRM < C->mRM && StartRM(C, C->nRM, GetIdxR(b, C), 0))
    C->nRM++;

  if(C->P->rev == 1 && C->nRM < C->mRM && StartRM(C, C->nRM,
  GetIdxRevR(b, C), 1))
    C->nRM++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

