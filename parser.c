#include <string.h>
#include "parser.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATE PARSER
//
PARSER *CreateParser(void){
  PARSER *PA = (PARSER *) Calloc(1, sizeof(PARSER));
  PA->sym    = 0;
  PA->type   = 0;
  PA->header = 0;
  PA->line   = 0;
  PA->dna    = 0;
  PA->nRead  = 0;
  return PA;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FILE TYPE
//
void FileType(PARSER *PA, FILE *IN){
  rewind(IN);
  PA->sym = fgetc(IN);
    switch(PA->sym){
    case '>': PA->type = 1; PA->nRead = 0; break;
    case '@': PA->type = 2; PA->nRead = 0; break;
    default : PA->type = 0; PA->nRead = 1; break;
    }
  rewind(IN);
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PARSE SYM
//
int32_t ParseSym(PARSER *PA, uint8_t sym){

  switch(PA->type){
    // IS A FASTA FILE
    case 1:
      switch(sym){
        case '>':  PA->header = 1; PA->nRead++;  return -1;
        case '\n': PA->header = 0;               return -1;
        default:   if(PA->header==1) return -1;
        }
    break;

    // IS A FASTQ FILE
    case 2:
      switch(PA->line){
        case 0: if(sym == '\n'){ 
                  PA->line = 1; 
                  PA->nRead++; 
                  PA->dna = 1; 
                  } 
        break;
        case 1: if(sym == '\n'){ 
                  PA->line = 2;
                  PA->dna = 0; 
                  } 
        break;
        case 2: if(sym == '\n'){ 
                  PA->line = 3; 
                  PA->dna = 0; 
                  } 
        break;
        case 3: if(sym == '\n'){ 
                  PA->line = 0; 
                  PA->dna = 0; 
                  } 
        break;
        }
      if(PA->dna == 0 || sym == '\n') return -1;
    break;

    // OTHER (SUCH AS DNA SEQ)
    default: ; // DO NOTHING
    }

  // NUCLEOTIDE PARSE
  if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
    return 4; // MAP EXTRA INTO 'N'

  return sym;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE PARSER
//
void RemoveParser(PARSER *PA){
  Free(PA);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
