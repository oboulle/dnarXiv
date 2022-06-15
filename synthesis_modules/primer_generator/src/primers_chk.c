
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
  TYPE_BANK  *BK;                   // banque de sequences d'ADN
  TYPE_HYBR  *HY;
  FILE *fpri;
  FILE *fok;                        // fichier r�sultat : amorces sans hybridation 
  FILE *fko;                        // fichier r�sultat : amorces avec hybridation
  FILE *finfo;                      // fichier r�sultat : informations relatives � fko
  FILE *fdic;                       // dictionnaire des amorces
  FILE *fpos;                       // position des hybridations

  int i,j,k,x,y;                    // variables de travail
  char tmpstring[128];
  char filename[1024];
  char cmd[1024];
  char c;
  char lastprimer[128];
  char primer[128];
  char ccseq[128];
  char invprimer[128];
  char invcomprimer[128];

  int  sizePrimer,startSeq,endSeq, nbprimer,delta, idxg;

  int listHybrid[10000];

  int nbHybrid;

  if (argc < 8)
    {
      fprintf (stderr, "usage : %s <genome> <primers> <parameters> <primers_out> <primers_hyb> <primers_dic> <primers_pos>\n", argv[0]);
      exit (0);
    }

  getParameter(argv[3]);

  if ((fpri = fopen(argv[2],"r"))==NULL) { fprintf (stderr,"cannot open %s\n",argv[2]); exit (0); }
  if ((fok  = fopen(argv[4],"w"))==NULL) { fprintf (stderr,"cannot open %s\n",argv[4]); exit (0); }
  if ((fko  = fopen(argv[5],"w"))==NULL) { fprintf (stderr,"cannot open %s\n",argv[5]); exit (0); }
  if ((fdic = fopen(argv[6],"w"))==NULL) { fprintf (stderr,"cannot open %s\n",argv[6]); exit (0); }
  sprintf(filename,"%s.tmp",argv[7]);
  if ((fpos = fopen(filename,"w"))==NULL) { fprintf (stderr,"cannot open %s.tmp\n",argv[7]); exit (0); }
  sprintf(filename,"%s.info",argv[5]);
  if ((finfo =fopen(filename,"w"))==NULL) { fprintf (stderr,"cannot open %s\n",filename); exit (0); }


  if ((BK       = (TYPE_BANK *)     malloc (sizeof(TYPE_BANK)))==NULL) MemError();
  if ((HY       = (TYPE_HYBR *)     malloc (sizeof(TYPE_HYBR)))==NULL) MemError();
  if ((BK->tseq = (int *)           malloc (NB_DIFF_SEED*sizeof(int)))==NULL) MemError();

  k = getSize(argv[1]); // misc.h
  if ((BK->seq  = (char *) malloc(k*sizeof(char)))==NULL) MemError();
  if ((BK->iseq = (int *) malloc(k*sizeof(int)))==NULL) MemError();
  getDnaSeq(BK,argv[1]); // misc.h

  if (BK->nb_seq>1)
    {
      fprintf (stderr,"file %s must contain only one sequence\n",argv[1]);
      exit(0);
    }

  k = indexDnaSeq(BK); // (index.h) indexation de la 1ere banque avec un pas de 1

  startSeq = BK->idx_offset[0];
  endSeq = startSeq + BK->idx_size[0]-1;

  strcpy(lastprimer,"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\0");
  nbprimer = 0;
  while ((k = getLine(fpri,tmpstring))!=0)
    {
      // pour faire joli : nombre de primers trait�s
      if (nbprimer%100==0) fprintf (stderr," %d %c",nbprimer,13);

      idxg = -1;

      // lecture du primer
      sizePrimer = getLine(fpri,primer);

      if (strncmp(primer,lastprimer,sizePrimer)==0) continue;
      strcpy(lastprimer,primer);

      fprintf (fdic,"%4d %s\n",nbprimer,primer);

      for (j=0; j<sizePrimer; j++) invprimer[j] = primer[sizePrimer-j-1]; invprimer[sizePrimer]='\0';
      for (j=0; j<sizePrimer; j++) invcomprimer[j] = comp (invprimer[j]); invcomprimer[sizePrimer]='\0';


      nbHybrid=0;

      // v�rification : 
      //
      //              5' G T T T G C G T G A C G T G T  3'
      // genome          | | | | | | | | | | | | | | |
      //              3' C A A A C G C A C T G C A C A  5'
      //                               | | |
      // amorce                    5'  T G A  3'
      //
      // il faut retrouver la meme s�quence sur le g�nome
      // et calculer l'hybridation par rapport au g�nome compl�ment� (BK->cseq)
      // hybrid ( 5'-ATG-3' , 3'-ACT-5')
      //          primer      comp(g�nome)

      for (i=0; i<sizePrimer-SIZE_SEED+1; i++)
	{
	  x = seed_code(&primer[i]);
	  k = BK->tseq[x];
	  while (k!= -1)
	    {
	      if ((k-i>startSeq)&&(k-i+sizePrimer-1<=endSeq)&&((k-i-startSeq)!=idxg))
		{
		  x=0; j=0;
		  while ((x<nbHybrid)&&(j==0)) { if (listHybrid[x]==k-i) j=1; else x++; }
		  if (x==nbHybrid) 
		    {
		      HY->loc = k-i-startSeq;
		      for (j=0; j<sizePrimer; j++) ccseq[j] = comp(BK->seq[k-i+j]);
		      if (k-i+sizePrimer < endSeq) 
			ccseq[sizePrimer]=comp(BK->seq[k-i+sizePrimer]); 
		      else
			ccseq[sizePrimer]=comp(BK->seq[k-i+sizePrimer-1]); 
		      check_hybrid (primer,sizePrimer,ccseq,comp(BK->seq[k-i-1]),HY); 
		      if ((HY->deltaG_max<maxDeltaGHybrid)||(HY->deltaG_3p<maxDeltaGHybrid3))
			{
			  listHybrid[x]=k-i; nbHybrid++;
			  printHybrid(finfo,nbHybrid,tmpstring,primer,&BK->seq[k-i],sizePrimer,idxg,HY,0);
			}
		    }
		}
	      k = BK->iseq[k];
	    }
	}

      // v�rification : 
      //
      //                                     3'  A G T  5'
      //                                         | | |
      //              5' G T T T G C G T G A C G T C A T T  3'
      // genome          | | | | | | | | | | | | | | | | |
      //              3' C A A A C G C A C T G C A G T A A  5'
      //
      // il faut retrouver l'amorce compl�ment�e et invers�e INV(COMP(TGA) = INV (ACT) = TCA sur le g�nome
      // et calculer l'hybridation par rapport au g�nome invers�
      // hybridation (5'-TGA-3' , 3'-ACT- 5')
      //              primer       inv(g�nome)

      for (i=0; i<sizePrimer-SIZE_SEED+1; i++)
	{
	  x = seed_code(&invcomprimer[i]); 
	  k = BK->tseq[x];
	  while (k!= -1)
	    {
	      if ((k-i>startSeq)&&(k-i+sizePrimer-1<=endSeq)&&((k-i-startSeq)!=idxg))
		{
		  x=0; j=0;
		  while ((x<nbHybrid)&&(j==0)) { if (listHybrid[x]==k-i) j=1; else x++; }
		  if (x==nbHybrid) 
		    { 
		      HY->loc = k-i-startSeq;
		      for (j=0; j<=sizePrimer; j++) ccseq[j] = BK->seq[k-i+sizePrimer-1-j];
		      if (k-i > startSeq) 
			ccseq[sizePrimer]=comp(BK->seq[k-i-1]); 
		      else
			ccseq[sizePrimer]=comp(BK->seq[k-i]); 
		      check_hybrid (primer,sizePrimer,ccseq,BK->seq[k-i+sizePrimer],HY);
		      if ((HY->deltaG_max<maxDeltaGHybrid)||(HY->deltaG_3p<maxDeltaGHybrid3))
			{ 
			  listHybrid[x]=k-i; nbHybrid++;
			  printHybrid(finfo,nbHybrid,tmpstring,primer,&BK->seq[k-i],sizePrimer,idxg,HY,1);
			}
		    }
		}
	      k = BK->iseq[k];
	    }
	}

      if (nbHybrid==0) fprintf (fok,"%s\n%s\n",tmpstring,primer); 
      else 
	{
	  fprintf(fko,"%s\n%s\n",tmpstring,primer); 
	  for (i=0; i<nbHybrid; i++)
	    {
	      k = getNumSeq(BK,listHybrid[i]);
	      k = BK->idx_offset[k];
	      fprintf (fpos,"%d %d\n",listHybrid[i]-k,nbprimer);
	    }
	}
      nbprimer++;
    }
  fclose(fpri);
  fclose(finfo);
  fclose(fok);
  fclose(fko);
  fclose(fdic);
  fclose(fpos);
  sprintf(cmd,"sort -n %s.tmp > %s; rm -f %s.tmp",argv[7],argv[7],argv[7]);
  system(cmd);
}

