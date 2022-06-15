
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "primers.h"

#include "parameter.h"
#include "filters.h"
#include "delta.h"
#include "tm.h"

#include "index.h"
#include "hybrid.h"
#include "display.h"

int main(int argc, char *argv[])
{
  FILE *fpri;                       // fichier d'amorces

  int i,j,k,x,y;                    // variables de travail
  char comments[128];
  char c;
  char primer[128];
  char invprimer[128];
  char invcomprimer[128];

  int  sizePrimer, sizeComments, tm;
  int deltaG[10];

  if (argc != 3)
    {
      fprintf (stderr, "usage : %s <primers> <parameters>\n", argv[0]);
      exit (0);
    }

  getParameter(argv[2]);

  if ((fpri=fopen(argv[1],"r"))==NULL) { fprintf (stderr,"cannot open %s\n",argv[1]); exit (0); }

  //while (fscanf(fpri,"%[^\n]%c",comments,&c)!=EOF)
  while ((sizeComments = getLine(fpri,comments))!=0)
    {
      sizePrimer = getLine(fpri,primer);

      k=max(sizePrimer,sizeComments);
      
      printf ("%s",comments);
      i=sizeComments;
      while (i<k) { printf (" "); i++; }

      printf ("\t%cGC",37);
      printf ("\tTm");
      printf ("\t#rep");
      printf ("\tHpin");
      printf ("\tCp_dG");
      printf ("\tCp_dG3'");
      printf ("\tSta_5'");
      printf ("\tSta_3'");
      printf ("\n");

      printf ("%s",primer);
      i=sizePrimer;
      while (i<k) { printf (" "); i++; }

      for (j=0; j<sizePrimer; j++) invprimer[j] = primer[sizePrimer-j-1]; invprimer[sizePrimer]='\0';
      for (j=0; j<sizePrimer; j++) invcomprimer[j] = comp (invprimer[j]); invcomprimer[sizePrimer]='\0';

      printf ("\t%3d",percentGC(primer,sizePrimer));
      tm = (int) oligotm (primer, (double)dnaConc, (double)saltConc,sizePrimer);
      printf ("\t%d",tm);
      printf ("\t%4d",repeat (primer, sizePrimer));
      printf ("\t%4d",hairpin (primer, sizePrimer));
      autocomp2(primer,sizePrimer,invprimer,deltaG);
      printf ("\t%5d\t%7d",deltaG[0],deltaG[1]);
      printf ("\t%6d",deltaGExt5Val (primer));
      printf ("\t%6d",deltaGExt3Val (primer,sizePrimer));
      printf ("\n\n");
    }
  fclose(fpri);
}

