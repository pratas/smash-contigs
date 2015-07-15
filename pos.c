#include "pos.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HEADERS *CreateHeaders(void){
  HEADERS *Headers = (HEADERS *) Calloc(1, sizeof(HEADERS));
  Headers->nPos        = 1;
  Headers->Pos[0].init = 0;
  Headers->Pos[0].end  = 0;
  return Headers;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveHeaders(HEADERS *Headers){
  uint64_t n;

  for(n = 0 ; n < Headers->nPos ; ++n)
    Free(Headers->Pos);
  Free(Headers);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

