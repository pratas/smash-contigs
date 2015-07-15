#include "pos.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HEADERS *CreateHeaders(void){
  HEADERS *Headers = (HEADERS *) Calloc(1, sizeof(HEADERS));
  Headers->nPos        = 1;
  Headers->Pos         = (POS *) Calloc(1, sizeof(POS));
  return Headers;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateHeaders(HEADERS *Headers){
  Headers->nPos++;
  Headers->Pos = (POS *) Realloc(Headers->Pos, Headers->nPos * sizeof(POS));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveHeaders(HEADERS *Headers){
  uint64_t n;
  for(n = 0 ; n < Headers->nPos ; ++n)
    Free(Headers->Pos);
  Free(Headers);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

