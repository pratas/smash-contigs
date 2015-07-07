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
RCLASS *CreateRClass(uint32_t max, uint32_t editions, uint32_t min, uint32_t k, 
uint8_t ir){
  uint32_t n;

  RCLASS *C   = (RCLASS *)  Calloc(1,   sizeof(RCLASS));
  C->RM       = (RMODEL *)  Calloc(max, sizeof(RMODEL));
  C->active   = (uint8_t *) Calloc(max, sizeof(uint8_t));
  C->nRM      = 0;
  C->mRM      = max;
  C->rev      = ir;
  C->idx      = 0;
  C->idxRev   = 0;
  C->kmer     = k;
  C->mult     = CalcMult(k);
  C->maxFails = editions;
  C->minSize  = min;
  for(n = 0 ; n < max ; ++n){
    C->RM[n].pos     = 0;
    C->RM[n].nFails  = 0;
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
uint64_t GetIdxRevRM(uint8_t *p, RCLASS *C){
  return (C->idxRev = (C->idxRev>>2)+GetCompNum(*p)*C->mult);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR REPEATS
//
uint64_t GetIdxRM(uint8_t *p, RCLASS *C){
  return (C->idx = ((C->idx-*(p-C->kmer)*C->mult)<<2)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET NON ACTIVE RMODEL ID
//
static int32_t GetFirstNonActiveRM(RCLASS *C){
  uint32_t k;
  for(k = 0 ; k < C->mRM ; ++k)
    if(C->active[k] == 0)
      return k;

  fprintf(stderr, "  [x] Error: impossible non active RM!\n");
  exit(1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET ACTIVE RMODEL ID
//
static int32_t GetFirstActiveRM(RCLASS *C){
  uint32_t k;
  for(k = 0 ; k < C->mRM ; ++k)
    if(C->active[k] == 1)
      return k;

  fprintf(stderr, "  [x] Error: impossible active RM!\n");
  exit(1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START EACH REPEAT MODEL
//
int32_t StartRMs(RCLASS *C, HASH *H, uint64_t idx, uint8_t ir){
  uint32_t n, k;
  ENTRY *E;

  if((E = GetHEnt(H, idx)) == NULL)
    return 0; // NEVER SEEN IN THE HASH TABLE, SO LETS QUIT

  while(C->nRM < C->mRM && n < E->nPos){

    k = GetFirstNonActiveRM(C);
 
    if(ir == 0){
      C->RM[k].pos = E->pos[n];
      }
    else{

/*    if(E->pos[n] <= C->kmer + 1){
        ++n;
        continue;
        }     */

      C->RM[k].pos = E->pos[n] - C->kmer - 1;
      }

    // RESET TO DEFAULTS
    C->RM[k].nFails  = 0;
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
    if(b[R->pos++] != s){
      R->nFails++;
      // R->lastHit = 0;
      }
    }
  else{
    if(GetCompNum(b[R->pos--]) != s){
      R->nFails++;
      // R->lastHit = 0;
      }
    }

  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
void StopRMs(RCLASS *C, uint64_t iBase, FILE *Writter){
  uint32_t id, n;

  if(C->nRM > 0){
    for(id = 0 ; id < C->mRM ; ++id){
      if(C->active[id] == 1){

        if(C->RM[id].nFails > C->maxFails){

          if(100 > C->minSize)
            ; //WRITE POS TO FILE

          C->active[id] == 0;
          }

        }
      }
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
    StartRMs(C, H, GetIdxRM(b, C), 0);

  if(C->rev == 1 && C->nRM < C->mRM)
    StartRMs(C, H, GetIdxRevRM(b, C), 1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

