#ifndef RMODEL_H_INCLUDED
#define RMODEL_H_INCLUDED

#include "defs.h"
#include "buffer.h"
#include "hash.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REPEAT MODELS TO HANDLE LONG SEGMENTS. DATA SUBSTITUTIONS DO NOT AFFECT THE
// PERFORMANCE SO MUCH AS IN CONTEXT MODELS.
//
typedef struct{
  uint32_t pos;      // POSITION OF THE FIRST PREDICTED K-MER 
  uint32_t size;     // CURRENT SIZE OF THE RMODEL
  uint32_t nHits;    // NUMBER OF TIMES THIS MODEL WAS CORRECT
  uint32_t nTries;   // NUMBER OF TIMES THIS MODEL WAS USED
  uint32_t winSize;  // WINDOW SIZE FOR MODEL ACTING
  uint8_t  *win;     // MODEL ACTING ACCORDING TO THE WINDOW SIZE
  uint8_t  rev;      // INVERTED REPEAT MODEL. IF REV=1 THEN IS ON
  }
RMODEL;

typedef struct{
  RMODEL   *RM;      // POINTER FOR EACH OF THE MULTIPLE REPEAT MODELS
  uint8_t  *active;  // THE REPEAT MODEL IS ACTIVE OR NOT
  uint32_t nRM;      // CURRENT NUMBER OF REPEAT MODELS
  uint32_t mRM;      // MAXIMUM NUMBER OF REPEAT MODELS
  uint32_t kmer;     // CONTEXT TEMPLATE SIZE FOR REPEAT MODEL
  uint64_t mult;     // INDEX MULTIPLIER
  uint64_t idx;      // CURRENT CONTEXT INDEX
  uint64_t idxRev;   // CURRENT INVERTED REPEAT INDEX
  uint8_t  rev;      // INVERTED REPETAT MODEL. IF REV='Y' THEN IS TRUE
  }
RCLASS;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t    CalcMult          (uint32_t);
RCLASS      *CreateRClass     (uint32_t, uint32_t, uint8_t);
uint64_t    GetIdxRevR        (uint8_t *, RCLASS *);
uint64_t    GetIdxR           (uint8_t *, RCLASS *);
int32_t     StartRM           (RCLASS *, HASH *, uint32_t, uint64_t, uint8_t);
void        UpdateRM          (RMODEL *, uint8_t *, uint8_t);
void        ResetRM           (RMODEL *);
void        RemoveRM          (RMODEL *);
void        RemoveAllRM       (RCLASS *);
void        StopRM            (RCLASS *, uint64_t, FILE *);
void        ResetAllRM        (RCLASS *, uint64_t, FILE *);
void        StartMultipleRMs  (RCLASS *, HASH *, uint8_t *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
