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
  C->n        = 0;
  C->nBases   = 0;
  C->mult     = CalcMult(k);
  C->maxFails = editions;
  C->minSize  = min;
  for(n = 0 ; n < max ; ++n){
    C->RM[n].pos     = 0;
    C->RM[n].init    = 0;
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
  int32_t k;

  for(k = 0 ; k < C->mRM ; ++k){
    if(C->active[k] == 0)
      return k;
    }

  fprintf(stderr, "Impossible state!\n");
  exit(1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START EACH REPEAT MODEL
//
void StartRMs(RCLASS *C, HASH *H, uint64_t iPos, uint64_t idx, uint8_t ir){
  uint32_t n = 0, k = 0;
  ENTRY *E;

  if((E = GetHEnt(H, idx)) == NULL)
    return; // NEVER SEEN IN THE HASH TABLE, SO LETS QUIT

  while(C->nRM < C->mRM && n < E->nPos){
    k = GetFirstNonActiveRM(C);

    if(E->pos[n] >= C->nBases - C->kmer || E->pos[n] < C->kmer){
      ++n;
      continue;
      }

    if(ir == 0)
      C->RM[k].init = C->RM[k].pos = E->pos[n];
    else
      C->RM[k].init = C->RM[k].pos = E->pos[n]-C->kmer-1;

    // RESET TO DEFAULTS
    C->RM[k].nFails  = 0;
    C->RM[k].stop    = 0;
    C->RM[k].write   = 0;
    C->RM[k].rev     = ir;
    C->RM[k].initRel = iPos;
    memset(C->RM[k].win, 0, C->RM[k].winSize); 

    C->active[k] = 1;  // SET IT ACTIVE
    C->nRM++;          // INCREASE NUMBER OF REPEATS
    ++n;
    }

  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// HIT REPEAT
//
static void Hit(RMODEL *R){
  if(R->nFails > 1)
    R->nFails--;
  ShiftBuffer(R->win, R->winSize, 0);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FAIL REPEAT
//
static void Fail(RMODEL *R){
  if(R->nFails < R->winSize)
    R->nFails++;
  else
    fprintf(stderr, "ERROR: nfails!\n");
  ShiftBuffer(R->win, R->winSize, 1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE REPEAT MODEL
//
void UpdateRMs(RCLASS *C, uint8_t *b, uint64_t ePos, uint8_t sym){
  uint32_t n;

  for(n = 0 ; n < C->mRM ; ++n){
    if(C->active[n] == 1){

      C->RM[n].size = labs(ePos-C->RM[n].initRel) + C->kmer;
      if(C->RM[n].win[0] == 1) C->RM[n].nFails--;
      
      if(C->RM[n].rev == 0){ // REGULAR REPEAT
        // HITS & FAILS
        if(b[C->RM[n].pos] != sym) Fail(&C->RM[n]);
        else                       Hit (&C->RM[n]);
        // STOP IF POS <= KMER
        if(C->RM[n].pos < C->nBases-1) C->RM[n].pos++;
        else                           C->RM[n].stop = 1;
        }
      else{ // INVERTED REPEAT
        // PROTECT EXTRA SYMBOLS
        if(b[C->RM[n].pos] == 4){
          Fail(&C->RM[n]); // SEE AFTER: DISCARDING POLITICS
          if(C->RM[n].pos > C->kmer) C->RM[n].pos--;
          else                       C->RM[n].stop = 1;
          continue;
          }
        // HITS & FAILS
        if(GetCompNum(b[C->RM[n].pos]) != sym) Fail(&C->RM[n]);
        else                                   Hit (&C->RM[n]);
        // STOP IF POS <= KMER
        if(C->RM[n].pos > C->kmer) C->RM[n].pos--;
        else                       C->RM[n].stop = 1;
        }
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET INDEX POSITION
//
static uint64_t GetIPoint(HEADERS *Head, uint64_t init){
  uint64_t id;
  for(id = 0 ; id < Head->iPos ; ++id){
    if(Head->Pos[id].init <= init && Head->Pos[id].end >= init){
      return id;   
      }
    }
  return id;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PROTECT NAMES
//
void ProtectVoidName(uint8_t *name, uint8_t type){
  if(name[0] == '\0'){
    name[0] = type == 0 ? 't' : 'r';
    name[1] = type == 0 ? 'a' : 'e';
    name[2] = type == 0 ? 'r' : 'f';
    name[3] = '\0';
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PRINT BLOCK
//
void PrintBlock(RCLASS *C, HEADERS *Head, uint64_t ePos, uint32_t n, uint8_t 
*nm, FILE *W){
  // # ID_TAR INIT_REL_TAR END_REL_TAR ID_REF INIT_ABS_REF END_ABS_REF SIZE
  uint64_t idxPos = 0;

  if(C->RM[n].rev == 0){
    // REGULAR REPEAT
    idxPos = GetIPoint(Head, C->RM[n].init-C->kmer);
    ProtectVoidName(nm, 0);
    ProtectVoidName(Head->Pos[idxPos].name, 1);

    fprintf(W, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t"
    "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",
    nm,                                                  // SAMPLE CONTIG NAME
    C->RM[n].initRel - C->kmer,                          // SAMPLE CONTIG INIT
    ePos,                                                // SAMPLE CONTIG END
    Head->Pos[idxPos].name,                              // TARGET CONTIG NAME
    (C->RM[n].init-C->kmer) - Head->Pos[idxPos].init,    // TARGET CONTIG INIT
    C->RM[n].pos - Head->Pos[idxPos].init,               // TARGET CONTIG END
    C->RM[n].size);
    }
  else{
    // REVERSE REPEAT
    idxPos = GetIPoint(Head, C->RM[n].pos);
    ProtectVoidName(nm, 0);
    ProtectVoidName(Head->Pos[idxPos].name, 1);

    fprintf(W, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t"
    "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",
    nm,                                                  // SAMPLE CONTIG NAME
    C->RM[n].initRel - C->kmer,                          // SAMPLE CONTIG INIT
    ePos,                                                // SAMPLE CONTIG END
    Head->Pos[idxPos].name,                              // TARGET CONTIG NAME

    C->RM[n].init + C->kmer - Head->Pos[idxPos].init,    // TARGET CONTIG INIT
    C->RM[n].pos - Head->Pos[idxPos].init,               // TARGET CONTIG END
    C->RM[n].size);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
static void ResetRM(RCLASS *C, uint32_t id){
  C->RM[id].write = 0;
  C->RM[id].stop  = 0;
  C->active[id]   = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// STOP USELESS REPEAT MODELS
//
void StopRMs(RCLASS *C, HEADERS *Head, uint64_t position, uint8_t *buf, FILE 
*Writter){
  int32_t id, largerRM = -1, largerRMIR = -1;
  uint64_t size = 0, sizeIR = 0;

  if(C->nRM > 0){
    for(id = 0 ; id < C->mRM ; ++id){
      if(C->active[id] == 1){
        if(C->RM[id].nFails > C->maxFails || C->RM[id].stop == 1){

          if(C->RM[id].size >= C->minSize){

            if(C->RM[id].rev == 0 && size < C->RM[id].size){
              size = C->RM[id].size;
              largerRM = id;
              }

            if(C->RM[id].rev == 1 && sizeIR < C->RM[id].size){
              sizeIR = C->RM[id].size;
              largerRMIR = id;
              }
           
            C->RM[id].write = 2;
            continue;
            }

          ResetRM(C, id);
          --C->nRM;
          }
        }
      }

    if(largerRM != -1)
      PrintBlock(C, Head, position, largerRM, buf, Writter);

    if(C->rev == 1 && largerRMIR != -1)
      PrintBlock(C, Head, position, largerRMIR, buf, Writter);

    for(id = 0 ; id < C->mRM ; ++id){
      if(C->RM[id].write == 2){ // IT WAS SMALLER
        if(C->RM[id].rev == 0){
          if(C->RM[largerRM].init > C->RM[id].init || 
             C->RM[largerRM].pos  < C->RM[id].pos){
            PrintBlock(C, Head, position, id, buf, Writter);
            }
          }
        else{
          if(C->RM[largerRMIR].init > C->RM[id].init || 
             C->RM[largerRMIR].pos  < C->RM[id].pos){
            PrintBlock(C, Head, position, id, buf, Writter);
            }
          }
        ResetRM(C, id);
        --C->nRM;
        }
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FORCE STOP REPEAT MODELS DURING END OF READ
//
void ResetAllRMs(RCLASS *C, HEADERS *Header, uint64_t position, uint8_t *buf, 
FILE *Writter){
  uint32_t n;

  for(n = 0 ; n < C->mRM ; ++n) // FORCE STOP RM's
    if(C->active[n] == 1)
      C->RM[n].stop = 1;
  
  StopRMs(C, Header, position, buf, Writter);
  for(n = 0 ; n < C->mRM ; ++n)
    ResetRM(C, n);
  C->nRM = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// START NEW REPEAT MODELS IF THERE IS STILL SPACE
//                         
void StartMultipleRMs(RCLASS *C, HASH *H, uint64_t iPos){
  if(C->nRM < C->mRM) StartRMs(C, H, iPos, C->idx, 0);
  if(C->rev == 1 && C->nRM < C->mRM) StartRMs(C, H, iPos, C->idxRev, 1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

