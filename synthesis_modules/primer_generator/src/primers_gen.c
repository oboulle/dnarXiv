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
#include "display.h"

char messerror[128] = "usage : primers_gen <fasta file> <parameters> <output file> [-c] [-s]\n\0";

void paramerror(char *s)
{
  fprintf (stderr,"option %s unknown\n",s);
  fprintf (stderr,"%s",messerror);
  exit (0);
}


int main (int argc, char *argv[])
{

  FILE *fseq;		        // s�quences au fromat FASTA
  FILE *fres;                   // r�sultats
  char *seq;                    // s�quence en cours de traitement

  int  i, j, k, x, y, t;
  char c;
  char *s;
  long L;

  char com[1024];               // commentaire associ� a la sequence
  int sizeSeq;                  // taille de la s�quence 
  int idxSeq;                   // index du primer sur la s�quence
  int nbSeq;                    // nombre de s�quences trait�es
  int gencomp;                  // 1: g�nome compl�ment�.
  int sortprimer;               // 1: primer tri� par ordre d'apparition sur le g�nome; 0: par ordre alphab�thique
  int DeltaG[2];                // valeur du deltaG
  int nb_primer;                // nombre de primers par s�quence

  int pCG;                      // pourcentage en GC du primer en cours de traitement
  float tmPrimer;               // valeur du tm du primer en cours de traitement

  char tmp[1024];               // chaine de caract�res � usage temporaire
  char primer[64];              // primer courant
  char Iprimer[64];             // primer invers�

  char **TabPrimers;            // tableau de pointeurs sur les primers

  if (argc < 4) { fprintf (stderr,"%s",messerror); exit (0); }

  sortprimer=0;
  gencomp=0;
  if (argc >= 5)
    {
      if (strcmp(argv[4],"-c") == 0) gencomp=1;
      else
	{
	  if (strcmp(argv[4],"-s") == 0) sortprimer=1;
	  else
	    {
	      paramerror(argv[4]);
	    }
	}
    }
  if (argc >= 6)
    {
      if (strcmp(argv[5],"-c") == 0) gencomp=1;
      else
	{
	  if (strcmp(argv[5],"-s") == 0) sortprimer=1;
	  else
	    {
	      paramerror(argv[5]);
	    }
	}
    }


  getParameter(argv[2]);

  // ouverture du fichier r�sultat

  if ((fres = fopen (argv[3], "w")) == NULL) 
    { fprintf  (stderr, "cannot open %s\nexit...\n", argv[3]); exit (0); }

  // lecture du fichier d'entr�e au format  FASTA

  if ((fseq = fopen (argv[1], "r")) == NULL) 
    { fprintf (stderr, "cannot open %s\nexit...\n", argv[1]); exit (0); }

  nbSeq=0;
  fscanf (fseq,"%c",&c);
  while (fscanf (fseq,"%c",&c)!=EOF)
    {
      // lecture commentaire
      k=0;
      while (c!='\n') 
	{ 
	  if (c>=32) { com[k]=c; k++; }
	  fscanf (fseq,"%c",&c);
	} 
      com[k]=' '; com[k+1]='\0';

      // lecture s�quence : premi�re fois pour obtenir sa taille
      L = ftell(fseq);
      sizeSeq=0;
      while (fscanf (fseq,"%c",&c)!=EOF)
	{
	  if (((c>='a')&&(c<='t')) || ((c>='A')&&(c<='T'))) sizeSeq++;
	  if (c=='>') break;
	}

      TabPrimers = (char **) calloc (sizeSeq+1, sizeof (char *));
      seq = (char *) calloc (sizeSeq+1, sizeof (char));

      // lecture s�quence : seconde fois pour la m�moriser
      fseek(fseq,L,SEEK_SET);
      sizeSeq=0;
      while (fscanf (fseq,"%c",&c)!=EOF)
	{
	  if (((c>='a')&&(c<='t')) || ((c>='A')&&(c<='T'))) seq[sizeSeq++]=c;
	  if (c=='>') break;
	}
      seq[sizeSeq]='\0';

      nbSeq++;

      if (gencomp==1)
	{
	  for (i=0; i<sizeSeq/2; i++) 
	    {
	      j = comp(seq[i]); seq[i] = comp(seq[sizeSeq-i-1]); seq[sizeSeq-i-1] = j; 
	    }
	  if ((sizeSeq%2)==1) seq[sizeSeq/2] = comp(seq[sizeSeq/2]);
	}

      // traitement de la s�quence.

      nb_primer = 0;

      for (idxSeq=0; idxSeq<sizeSeq-sizePrimer+1; idxSeq++)
	{
	  // -----------------------------------------------------------------------------------  filtre A,C,G,T
	  if (checkPrimer(&seq[idxSeq],sizePrimer,primer)==0) continue;
	  
	  // -----------------------------------------------------------------------------------  filtre pourcentage en GC
	  pCG = percentGC(primer,sizePrimer);
	  if  ((pCG < pcGCMin) || (pCG > pcGCMax)) continue;

	  // -----------------------------------------------------------------------------------  filtre tm
	  tmPrimer = oligotm (primer, (double)dnaConc, (double)saltConc, sizePrimer);
	  if ((tmPrimer < oligoTmMin) || (tmPrimer > oligoTmMax)) continue;
	  
	  // -----------------------------------------------------------------------------------  filtre nombre de r�p�titions
	  if (repeat (primer, sizePrimer) > nbRepeat) continue;

	  // -----------------------------------------------------------------------------------  filtre tiges-boucles
	  if (hairpin (primer, sizePrimer) > maxHpDup) continue;

	  // -----------------------------------------------------------------------------------  filtre auto compl�mentarit� 
	  for (k = 0; k < sizePrimer; k++) Iprimer[k] =  primer[sizePrimer-k-1];
	  Iprimer[sizePrimer] = '\0'; 
	  autocomp2(primer,sizePrimer,Iprimer,DeltaG);
	  if ((DeltaG[0]<maxDeltaGAuto)||(DeltaG[1]<maxDeltaGAuto3)) continue; 

	  // -----------------------------------------------------------------------------------  filtre DELTA_G aux extr�mit�s 5'et 3'
	  x = deltaG (primer, sizePrimer);
	  if (x == 0) continue;

	  // -----------------------------------------------------------------------------------  �criture des r�sultats
	  TabPrimers[nb_primer] = &seq[idxSeq];
	  nb_primer++;
	}
      if (sortprimer == 0) 
	{
	  qsort(TabPrimers,nb_primer,sizeof (* TabPrimers),compare_primers);
	  for (i=0; i<nb_primer; i++)
	    {
	      printResults(fres, com, TabPrimers[i], sizePrimer, TabPrimers[i]-seq, gencomp*sizeSeq);
	    }
	}
      else
	{
	  if (gencomp == 0)
	    {
	      for (i=0; i<nb_primer; i++)
		{
		  printResults(fres, com, TabPrimers[i], sizePrimer, TabPrimers[i]-seq, gencomp*sizeSeq);
		}
	    }
	  else
	    {
	      for (i=nb_primer-1; i>=0; i--)
		{
		  printResults(fres, com, TabPrimers[i], sizePrimer, TabPrimers[i]-seq, gencomp*sizeSeq);
		}
	    }
	}


      free (TabPrimers);
      free (seq);
    }

  if (nbSeq>1) printf ("%d sequences processed\n",nbSeq); else printf ("%d sequence processed\n",nbSeq);
  fclose(fseq);
  fclose(fres);
}
