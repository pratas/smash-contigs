#ifndef PAINT_H_INCLUDED
#define PAINT_H_INCLUDED

#include "defs.h"

#define DEFAULT_CX             70.0
#define DEFAULT_CY             75.0
#define DEFAULT_TX             50.0
#define DEFAULT_TY             82.0
#define DEFAULT_WIDTH          25.0
#define LEVEL_SATURATION       160
#define LEVEL_VALUE            160
#define DEFAULT_WIDTH          25.0
#define DEFAULT_SPACE          10.0
#define EXTRA                  150 

uint32_t ratio;

typedef struct
  {
  char    *backColor;
  double  width;
  double  cx; 
  double  cy;
  double  tx;
  double  ty;
  double  refSize;
  double  tarSize;
  double  maxSize;  
  }
Painter;

typedef struct
  {
  uint8_t  r;
  uint8_t  g;
  uint8_t  b;
  } 
RgbColor;

typedef struct
  {
  uint8_t  h;
  uint8_t  s;
  uint8_t  v;
  } 
HsvColor;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Painter   *CreatePainter (double, double, char *);
RgbColor  HsvToRgb       (HsvColor);
HsvColor  RgbToHsv       (RgbColor);
char      *GetRgbColor   (uint8_t);
void      PrintFinal     (FILE *);
void      PrintHead      (FILE *, double, double);
void      RectOval       (FILE *, double, double, double, double, char *);
void      RectOvalIR     (FILE *, double, double, double, double, char *);
void      Rect           (FILE *, double, double, double, double, char *);
void      RectIR         (FILE *, double, double, double, double, char *);
void      Chromosome     (FILE *, double, double, double, double);
void      Text           (FILE *, double, double, char *);
void      Textfloatoat   (FILE *, double, double, double);
double    GetPoint       (uint64_t);
void      SetRatio       (uint32_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
