#define sizePrimerMax                 128
#define deltaGAutoInit               1960
#define deltaGHybridInit             1960

// valeurs par defaut pour la taille min/max d'une amorce

int sizePrimer = 25;

// valeurs par defaut pour le pourcentage en GC

int pcGCMin = 40;
int pcGCMax = 60;

// valeurs par defaut pour le calcul du Tm

int oligoTmMin = 57;
int oligoTmMax = 62;
int saltConc   = 500;
int dnaConc    = 50;

// valeurs par defaut pour la detection de tige-boucle

int maxHpDup  = 4;
int maxHpLoop = 4; 

// valeurs par defaut pour le nombre de repetitions

int nbRepeat = 5;

// valeurs par defaut pour l'auto complementarite

int  maxDeltaGAuto   = -10000;
int  maxDeltaGAuto3  = -7000;
int  sizeDeltaGAuto  = 6;
int  sizeDeltaGAuto3 = 8;

// valeurs par defaut pour la stabilite aux extremites

int  sizeExt5 = 5;
int  sizeExt3 = 5;
int  deltaG5  = -4000;
int  deltaG3min = -6000;
int  deltaG3max = -4000;

// valeurs par defaut pour l'hybridation

int  maxDeltaGHybrid3 = -9000;
int  maxDeltaGHybrid  = -16000;
int  sizeDeltaGHybrid3 = 8;

// constantes pour l'indexation

#define SIZE_SEED                        6
#define NB_DIFF_SEED                  4096


// ouverture du fichier de parametres et initialisation des parametres

int getParameter(char *filename)
{
  FILE *fpar;
  char tmp[1024];
  char c;

  if ((fpar = fopen (filename, "r")) == NULL) 
    {
      fprintf (stderr, "cannot open %s\n", filename);
      fprintf (stderr, "exit...\n");
      exit (0);
    }

  while (fscanf (fpar, "%s", tmp) != EOF)
    {
      // lecture des commentaire : ligne qui commence par #
      if (tmp[0]== '#')
	{
	  fscanf(fpar,"%c",&c);
	  while (c!='\n') fscanf(fpar,"%c",&c);
	  continue;
	}

      if (strcmp (tmp, "sizePrimer")            == 0) { fscanf (fpar, "%d", &sizePrimer);            continue; }
      if (strcmp (tmp, "pcGCMin")               == 0) { fscanf (fpar, "%d", &pcGCMin);               continue; }
      if (strcmp (tmp, "pcGCMax")               == 0) { fscanf (fpar, "%d", &pcGCMax);               continue; }
      if (strcmp (tmp, "oligoTmMin")            == 0) { fscanf (fpar, "%d", &oligoTmMin);            continue; }
      if (strcmp (tmp, "oligoTmMax")            == 0) { fscanf (fpar, "%d", &oligoTmMax);            continue; }
      if (strcmp (tmp, "dnaConc")               == 0) { fscanf (fpar, "%d", &dnaConc);               continue; }
      if (strcmp (tmp, "saltConc")              == 0) { fscanf (fpar, "%d", &saltConc);              continue; }
      if (strcmp (tmp, "maxHpDup")              == 0) { fscanf (fpar, "%d", &maxHpDup);              continue; }
      if (strcmp (tmp, "maxHpLoop")             == 0) { fscanf (fpar, "%d", &maxHpLoop);             continue; }
      if (strcmp (tmp, "nbRepeat")              == 0) { fscanf (fpar, "%d", &nbRepeat);              continue; }
      if (strcmp (tmp, "sizeExt5")              == 0) { fscanf (fpar, "%d", &sizeExt5);              continue; }
      if (strcmp (tmp, "sizeExt3")              == 0) { fscanf (fpar, "%d", &sizeExt3);              continue; }
      if (strcmp (tmp, "deltaG5")               == 0) { fscanf (fpar, "%d", &deltaG5);               continue; }
      if (strcmp (tmp, "deltaG3min")            == 0) { fscanf (fpar, "%d", &deltaG3min);            continue; }
      if (strcmp (tmp, "deltaG3max")            == 0) { fscanf (fpar, "%d", &deltaG3max);            continue; }
      if (strcmp (tmp, "maxDeltaGAuto")         == 0) { fscanf (fpar, "%d", &maxDeltaGAuto);         continue; }
      if (strcmp (tmp, "maxDeltaGAuto3")        == 0) { fscanf (fpar, "%d", &maxDeltaGAuto3);        continue; }
      if (strcmp (tmp, "sizeDeltaGAuto")        == 0) { fscanf (fpar, "%d", &sizeDeltaGAuto);        continue; }
      if (strcmp (tmp, "sizeDeltaGAuto3")       == 0) { fscanf (fpar, "%d", &sizeDeltaGAuto3);       continue; }
      if (strcmp (tmp, "maxDeltaGHybrid")       == 0) { fscanf (fpar, "%d", &maxDeltaGHybrid);       continue; }
      if (strcmp (tmp, "maxDeltaGHybrid3")      == 0) { fscanf (fpar, "%d", &maxDeltaGHybrid3);      continue; }
      if (strcmp (tmp, "sizeDeltaGHybrid3")     == 0) { fscanf (fpar, "%d", &sizeDeltaGHybrid3);     continue; }
   
      fprintf (stderr," *** parameter %s unknown\n     exit...\n", tmp);
      exit (0);
    }
  fclose (fpar);
  return (0);
}
