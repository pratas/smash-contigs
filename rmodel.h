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
  uint64_t initRel;  // FIRST POSITION OF THE REPEAT IN THE CONTIG
  uint64_t size;     // CURRENT SIZE OF THE RMODEL IN THE CONTIG
  uint64_t init;     // FIRST POSITION OF THE FIRST PREDICTED K-MER 
  uint64_t pos;      // CURRENT POSITION OF THE FIRST PREDICTED K-MER 
  int32_t  nFails;   // NUMBER OF TIMES THIS MODEL WAS WRONG
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
  uint32_t minSize;  // MINIMUM BLOCK SIZE
  uint32_t maxFails; // MAXIMUM NUMBER OF FAILS ALLOWED 
  uint64_t mult;     // INDEX MULTIPLIER
  uint64_t idx;      // CURRENT CONTEXT INDEX
  uint64_t idxRev;   // CURRENT INVERTED REPEAT INDEX
  uint8_t  rev;      // INVERTED REPETAT MODEL. IF REV='Y' THEN IS TRUE
  uint64_t n;        // NUMBER OF PROCESSED N's
  }
RCLASS;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t   CalcMult         (uint32_t);
RCLASS     *CreateRClass    (uint32_t, uint32_t, uint32_t, uint32_t, uint8_t);
void       RemoveRClass     (RCLASS *);
uint64_t   GetIdxRevRM      (uint8_t *, RCLASS *);
uint64_t   GetIdxRM         (uint8_t *, RCLASS *);
void       PrintBlock       (RCLASS *, uint64_t, uint32_t, FILE *);
void       StartRMs         (RCLASS *, HASH *, uint64_t, uint64_t, uint8_t);
void       UpdateRMs        (RCLASS *, uint8_t *, uint64_t, uint8_t);
void       StopRMs          (RCLASS *, uint64_t, FILE *);
void       ResetAllRMs      (RCLASS *, uint64_t, FILE *);
void       StartMultipleRMs (RCLASS *, HASH *, uint64_t, uint8_t *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
