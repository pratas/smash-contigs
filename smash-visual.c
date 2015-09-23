#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/uio.h>
#include <sys/mman.h>
#include "mem.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "paint.h"
#include "common.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - P L O T - - - - - - - - - - - - - - - -
void PrintPlot(char *posFile, uint32_t width, uint32_t space, uint32_t mult,
uint32_t start, uint64_t minimum){
  FILE *PLOT = NULL, *POS = NULL;
  char backColor[] = "#ffffff";
  int64_t conNBases = 0, refNBases = 0;
  char watermark[MAX_FILENAME];
  Painter *Paint;

  fprintf(stderr, "  [+] Printing plot ...\n");

  POS  = Fopen(posFile,  "r");
  PLOT = Fopen(P->image, "w");
 
  if(fscanf(POS, "%s\t%"PRIi64"\t%"PRIi64"\n", watermark, &conNBases,
  &refNBases) != 3 || watermark[0] != '#' || watermark[1] != 'S' ||
  watermark[2] != 'C' || watermark[3] != 'F'){
    fprintf(stderr, "  [x] Error: unknown positions file format!\n");
    exit(1);
    }

  if(P->verbose){
    fprintf(stderr, "      Reference number of bases: %"PRIu64"\n", refNBases);
    fprintf(stderr, "      Target number of bases: %"PRIu64"\n", conNBases);
    }

  SetRatio(MAX(refNBases, conNBases) / DEFAULT_SCALE);
  Paint = CreatePainter(GetPoint(refNBases), GetPoint(conNBases), (double) 
  width, (double) space, backColor);

  PrintHead(PLOT, (2 * DEFAULT_CX) + (((Paint->width + Paint->space) * 2) -
  Paint->space), Paint->maxSize + EXTRA);
  Rect(PLOT, (2 * DEFAULT_CX) + (((Paint->width + Paint->space) * 2) -
  Paint->space), Paint->maxSize + EXTRA, 0, 0, backColor);
  RectOval(PLOT, Paint->width, Paint->refSize, Paint->cx, Paint->cy,
  backColor);
  RectOval(PLOT, Paint->width, Paint->tarSize, Paint->cx, Paint->cy,
  backColor);
  Text(PLOT, Paint->cx,                           Paint->cy-15, "REF");
  Text(PLOT, Paint->cx+Paint->width+Paint->space, Paint->cy-15, "CON");

  // IF MINIMUM IS SET DEFAULT, RESET TO BASE MAX PROPORTION
  if(minimum == 0)
    minimum = MAX(refNBases, conNBases) / 100;

  int64_t ri, rf, ci, cf, cx, cy;
  uint64_t regular = 0, inverse = 0, ignored = 0;  
  while(1){
    char tmp1[MAX_STR] = {'\0'}, tmp2[MAX_STR] = {'\0'};
    if(fscanf(POS, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
    "%"PRIi64"\t%"PRIi64"\n", tmp1, &ci, &cf, &cx, &cy, tmp2, &ri, &rf) != 8)
      break;

    if(labs(rf-ri) < minimum || labs(cx-cy) < minimum){
      ++ignored;
      continue;
      }

    if(rf > ri){
      switch(P->link){
        case 1: 
          Line(PLOT, 2, Paint->cx + Paint->width, 
          Paint->cy + GetPoint(ri+((rf-ri)/2.0)), 
          Paint->cx + Paint->space + Paint->width, 
          Paint->cy + GetPoint(cx+((cy-cx)/2.0)), "black");
        break;
        case 2:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri+((rf-ri)/2.0)),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx+((cy-cx)/2.0)), 
          GetRgbColor(start * mult));
        break;
        case 3:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx), "black");
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy), "black");
        break;
        case 4:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx), 
          GetRgbColor(start * mult));
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy), 
          GetRgbColor(start * mult));
        break;
        case 5:
          Polygon(PLOT, 
          Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx),
          GetRgbColor(start * mult), "grey");    
        break;
        default:
        break;
        }        
      
      Rect(PLOT, Paint->width, GetPoint(rf-ri), Paint->cx, Paint->cy +
      GetPoint(ri), GetRgbColor(start * mult));

      Rect(PLOT, Paint->width, GetPoint(cy-cx), Paint->cx + Paint->space + 
      Paint->width, Paint->cy + GetPoint(cx), GetRgbColor(start * mult));

      ++regular; 
      }
    else{ 
      switch(P->link){
        case 1:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf+((ri-rf)/2.0)),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy+((cx-cy)/2.0)), "green");
        break;
        case 2:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf+((ri-rf)/2.0)),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy+((cx-cy)/2.0)), 
          GetRgbColor(start * mult));
        break;
        case 3:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy), "green");
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx), "green");
        break;
        case 4:
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy), 
          GetRgbColor(start * mult));
          Line(PLOT, 2, Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx), 
          GetRgbColor(start * mult));
        break;
        case 5:
          Polygon(PLOT, 
          Paint->cx + Paint->width,
          Paint->cy + GetPoint(rf),
          Paint->cx + Paint->width,
          Paint->cy + GetPoint(ri),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cx),
          Paint->cx + Paint->space + Paint->width,
          Paint->cy + GetPoint(cy),
          GetRgbColor(start * mult), "grey");
        break;
        default:
        break;
        }

      Rect(PLOT, Paint->width, GetPoint(ri-rf), Paint->cx, Paint->cy +
      GetPoint(rf), GetRgbColor(start * mult));

      RectIR(PLOT, Paint->width, GetPoint(cy-cx), Paint->cx + Paint->space + 
      Paint->width, Paint->cy + GetPoint(cx), GetRgbColor(start * mult));

      ++inverse;
      }

    ++start;
    }
  rewind(POS);

  fprintf(stderr, "      Found %"PRIu64" regular regions. \n", regular);
  fprintf(stderr, "      Found %"PRIu64" inverted regions.\n", inverse);
  if(P->verbose)
    fprintf(stderr, "      Ignored %"PRIu64" regions.\n", ignored);

  Chromosome(PLOT, Paint->width, Paint->refSize, Paint->cx, Paint->cy);
  Chromosome(PLOT, Paint->width, Paint->tarSize, Paint->cx + Paint->space +
  Paint->width, Paint->cy);
  PrintFinal(PLOT);
  fclose(POS);

  fprintf(stderr, "      Done!\n");
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t main(int argc, char *argv[]){
  char **p = *&argv;
  uint32_t width, space, mult, start, minimum;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEF_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenuVisual();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersionVisual();
    return EXIT_SUCCESS;
    }

  P->verbose    = ArgsState (DEF_VERBOSE, p, argc, "-v" );
  P->force      = ArgsState (DEF_FORCE,   p, argc, "-F" );
  P->link       = ArgsNum   (DEF_LINK,    p, argc, "-l", MIN_LINK, MAX_LINK);
  width         = ArgsNum   (DEF_WIDT,    p, argc, "-w", MIN_WIDT, MAX_WIDT);
  space         = ArgsNum   (DEF_SPAC,    p, argc, "-s", MIN_SPAC, MAX_SPAC);
  mult          = ArgsNum   (DEF_MULT,    p, argc, "-m", MIN_MULT, MAX_MULT);
  start         = ArgsNum   (DEF_BEGI,    p, argc, "-b", MIN_BEGI, MAX_BEGI);
  minimum       = ArgsNum   (DEF_MINP,    p, argc, "-c", MIN_MINP, MAX_MINP);
  P->image      = ArgsFilesImg           (p, argc, "-x");

  fprintf(stderr, "\n");
  fprintf(stderr, "==[ PROCESSING ]====================\n");
  PrintPlot(argv[argc-1], width, space, mult, start, minimum);
  fprintf(stderr, "\n");

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
