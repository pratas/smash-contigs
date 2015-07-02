#ifndef RMODEL_H_INCLUDED
#define RMODEL_H_INCLUDED

#include "defs.h"
#include "buffer.h"

#define INIWEIGHT  0.990 //XXX: REMOVE MIXTURE...

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
  uint64_t idx;      // CURRENT CONTEXT INDEX
  uint64_t idxRev;   // CURRENT INVERTED REPEAT INDEX
  uint64_t mult;     // INDEX MULTIPLIER
  uint32_t ctx;      // CONTEXT TEMPLATE SIZE FOR REPEAT MODEL
  uint32_t limit;    // REPEAT PERFORMANCE LIMIT, ASSOCIATED WITH BETA
  double   alpha;    // ALPHA PROBABILITY ESTIMATOR
  double   beta;     // REPEAT PERFORMANCE DECAYING FOR REPEAT MOVE
  double   gamma;    // PERFORMANCE DECAYING PARAMETER
  uint8_t  rev;      // INVERTED REPEAT USAGE [MEMORY/TIME PAINFUL]
  }
RPARAM;

typedef struct{
  uint32_t pos;      // POSITION OF THE FIRST PREDICTED K-MER
  uint32_t nHits;    // NUMBER OF TIMES THIS MODEL WAS CORRECT
  uint32_t nTries;   // NUMBER OF TIMES THIS MODEL WAS USED
  double   probs[4]; // REPEAT MODEL SYMBOL PROBABILITIES
  double   weight;   // WEIGHT OF THE MODEL FOR MIXTURE
  double   acting;   // PERFORMANCE PARAMETER
  double   lastHit;  // PERFORMANCE PARAMETER
  uint8_t  rev;      // INVERTED REPETAT MODEL. IF REV='Y' THEN IS TRUE
  }
RMODEL;

typedef struct{
  HASH     *hash;    // REPEATING KMERS HASH TABLE
  RMODEL   *RM;      // POINTER FOR EACH OF THE MULTIPLE REPEAT MODELS
  RPARAM   *P;       // EXTRA PARAMETERS FOR REPEAT MODELS
  uint32_t nRM;      // CURRENT NUMBER OF REPEAT MODELS
  uint32_t mRM;      // MAXIMUM NUMBER OF REPEAT MODELS
  uint64_t size;     // SIZE OF THE INPUT SEQUENCE
  }
RCLASS;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t        CalcMult          (uint32_t);
RCLASS          *CreateRC         (uint32_t, double, double, uint32_t, 
                                   uint32_t, double, uint8_t);
uint64_t        GetIdxRevR        (uint8_t *, RCLASS *);
uint64_t        GetIdxR           (uint8_t *, RCLASS *);
ENTRY           *GetHEnt          (RCLASS *, uint64_t);
int32_t         StartRM           (RCLASS *, uint32_t, uint64_t, uint8_t);
void            InsertKmerPos     (RCLASS *, uint64_t, uint32_t);
void            UpdateRM          (RMODEL *, uint8_t *, uint8_t);
void            ComputeRMProbs    (RCLASS *, RMODEL *, uint8_t *);
void            StopRM            (RCLASS *);
void            StartMultipleRMs  (RCLASS *, uint8_t *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
