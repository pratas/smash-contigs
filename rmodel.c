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
    C->RM[n].winSize = k;
    C->RM[n].win     = (uint8_t *) Calloc(k + 1, sizeof(uint8_t));
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
int32_t StartRMs(RCLASS *C, HASH *H, uint64_t idx, uint8_t ir){
  uint32_t n, k;
  ENTRY *E;

  if((E = GetHEnt(H, idx)) == NULL)
    return 0; // NEVER SEEN IN THE HASH TABLE, SO LETS QUIT

  while(C->nRM < C->mRM && n < E->nPos)
    {
    for(k = 0 ; k < C->mRM ; ++k)
      if(C->active[k] == 0)
        break;   // GET NON ACTIVE REPEAT ID
 
    if(ir == 0)
      {
      C->RM[k].pos = E->pos[n];
      }
    else
      {
/*
      if(E->pos[n] <= C->kmer + 1) // IN THE BEGINNIG [THIS MIGHT BE REMOVED]
        {
        ++n;
        continue;
        }
*/

      C->RM[k].pos = E->pos[n] - C->kmer - 1;
      }

    // RESET TO DEFAULTS
    C->RM[k].nHits   = 0;
    C->RM[k].nTries  = 0;
    C->RM[k].rev     = ir;
    memset(C->RM[k].win, 0, C->RM[k].winSize); 

    C->active[k] = 1;  // SET IT ACTIVE
    C->nRM++;          // INCREASE NUMBER OF REPEATS

    ++n;
    }

  return 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE REPEAT MODEL
//
void UpdateRMs(RMODEL *R, uint8_t *b, uint8_t s){
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
void StopRMs(RCLASS *C, uint64_t iBase, FILE *Writter){
  uint32_t n;

  for(n = 0 ; n < C->mRM ; ++n){
    if(C->active[n] == 1)
       break;


    }

  for(;;){
    for(n = 0 ; n < C->nRM ; ++n){
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
// FORCE STOP REPEAT MODELS DURING END OF READ
//
void ResetAllRMs(RCLASS *C, uint64_t iBase, FILE *Writter){
  uint32_t n;

  for(n = 0 ; n < C->mRM ; ++n){
    fprintf(Writter, "%"PRIu64"\t%"PRIu64"\n", iBase, 0);
    C->active[n] = 0;
    }
  C->nRM = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START NEW REPEAT MODELS IF THERE IS STILL SPACE
//                         
void StartMultipleRMs(RCLASS *C, HASH *H, uint8_t *b){

  if(C->nRM < C->mRM)
    StartRMs(C, H, GetIdxR(b, C), 0);

  if(C->rev == 1 && C->nRM < C->mRM)
    StartRMs(C, H, GetIdxRevR(b, C), 1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

