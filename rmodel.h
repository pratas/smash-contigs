#ifndef RMODEL_H_INCLUDED
#define RMODEL_H_INCLUDED

#include "defs.h"
#include "buffer.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// HASH TABLE TO STORE REPEATING POSITIONS AND INDEXES ALONG THE SEQUENCE
// DO NOT CHANGE THESE MACRO VALUES UNLESS YOU REALLY KNOW WHAT YOU ARE DOING!
//
#define HSIZE        16777259   // NEXT PRIME AFTER 16777216 (24 BITS)
#define MAX_CTX      20         // ((HASH_SIZE (24 B) + KEY (16 B))>>1) = 20 

typedef uint32_t PPR;  // PRECISION OF THE POSITION POINTER FOR REPEATS

typedef struct{
  uint16_t key;      // THE KEY (INDEX / HASHSIZE) STORED IN THIS ENTRY
  uint16_t nPos;
  PPR      *pos;     // LIST WITH THE REPEATING POSITIONS
  }
ENTRY;

typedef struct{
  uint16_t *size;    // NUMBER OF KEYS FOR EACH ENTRY
  ENTRY    **ent;    // ENTRIES VECTORS POINTERS
  }
HASH;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REPEAT MODELS TO HANDLE LONG SEGMENTS. DATA SUBSTITUTIONS DO NOT AFFECT THE
// PERFORMANCE SO MUCH AS IN CONTEXT MODELS.
//
typedef struct{
  uint32_t pos;      // POSITION OF THE FIRST PREDICTED K-MER 
  uint32_t size;     // CURRENT SIZE OF THE RMODEL
  uint32_t nHits;    // NUMBER OF TIMES THIS MODEL WAS CORRECT
  uint32_t nTries;   // NUMBER OF TIMES THIS MODEL WAS USED
  uint8_t  *act;     // MODEL ACTING ACCORDING TO THE WINDOW SIZE
  uint8_t  rev;      // INVERTED REPEAT MODEL. IF REV=1 THEN IS ON
  }
RMODEL;

typedef struct{
  RMODEL   *RM;      // POINTER FOR EACH OF THE MULTIPLE REPEAT MODELS
  uint32_t nRM;      // CURRENT NUMBER OF REPEAT MODELS
  uint32_t mRM;      // MAXIMUM NUMBER OF REPEAT MODELS
  uint32_t kmer;     // CONTEXT TEMPLATE SIZE FOR REPEAT MODEL
  uint64_t mult;     // INDEX MULTIPLIER
  uint64_t idx;      // CURRENT CONTEXT INDEX
  uint64_t idxRev;   // CURRENT INVERTED REPEAT INDEX
  uint8_t  rev;      // INVERTED REPETAT MODEL. IF REV='Y' THEN IS TRUE
  uint64_t size;     // XXX: MAYBE NOT USED // SIZE OF THE INPUT SEQUENCE
  }
RCLASS;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t    CalcMult          (uint32_t);
HASH        *CreateHash       (void);
RCLASS      *CreateRC         (uint32_t, uint32_t, uint8_t, uint64_t);
uint64_t    GetIdxRevR        (uint8_t *, RCLASS *);
uint64_t    GetIdxR           (uint8_t *, RCLASS *);
ENTRY       *GetHEnt          (HASH   *, uint64_t);
int32_t     StartRM           (RCLASS *, HASH *, uint32_t, uint64_t, uint8_t);
void        InsertKmerPos     (HASH   *, uint64_t, uint32_t);
void        UpdateRM          (RMODEL *, uint8_t *, uint8_t);
void        StopRM            (RCLASS *);
void        StartMultipleRMs  (RCLASS *, HASH *, uint8_t *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
