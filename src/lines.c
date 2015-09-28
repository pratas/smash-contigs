#include <stdio.h>
#include <stdlib.h>
#include "lines.h"
#include "defs.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LCACHE *CreateLCache(size_t size){
  LCACHE *LC = (LCACHE *) Calloc(1, sizeof(LCACHE));
  LC->max    = size;
  LC->idx    = 0;
  LC->Lines  = (LINE *) Calloc(LC->max, sizeof(LINE));
  return LC;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintLine(LCACHE *LC, FILE *F, uint32_t idx){
  fprintf(F, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                 "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n",
          LC->Lines[idx].contigs_name, 
          LC->Lines[idx].contigs_relative_init_pos,
          LC->Lines[idx].contigs_relative_end_pos,
          LC->Lines[idx].contigs_absolute_init_pos,
          LC->Lines[idx].contigs_absolute_end_pos,
          LC->Lines[idx].reference_name, 
          LC->Lines[idx].reference_relative_init_pos,
          LC->Lines[idx].reference_relative_end_pos, 
          LC->Lines[idx].reference_absolute_init_pos,
          LC->Lines[idx].reference_absolute_end_pos);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetChar2Bar0(char *str){
  uint32_t n;
  for(n = 0 ; n < MAX_STR ; ++n)
    str[n] = '\0';
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateLCacheIdx(LCACHE *LC){
  if(++LC->idx == LC->max)
    LC->idx = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveLCache(LCACHE *LC){
  uint32_t n;

  for(n = 0 ; n < LC->max ; ++n)
    Free(LC->Lines);
  Free(LC);
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
