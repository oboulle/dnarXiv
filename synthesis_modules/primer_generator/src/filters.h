
// verification qu'un primer comporte uniquement les 4 lettres ACGT

int checkPrimer(char *seq, int sizeP, char *primer)
{
  int i;

  for (i=0; i<sizeP; i++)
    {
      switch (seq[i])
	{
	case 'a' : primer[i]='A'; break;
	case 'c' : primer[i]='C'; break;
	case 'g' : primer[i]='G'; break;
	case 't' : primer[i]='T'; break;
	case 'A' : primer[i]='A'; break;
	case 'C' : primer[i]='C'; break;
	case 'G' : primer[i]='G'; break;
	case 'T' : primer[i]='T'; break;
	default : return 0;
	}
    }
  primer[i]='\0';
  return 1;
}


// filtre sur le pourcentage en GC
// verifie qu'il est compris entre pcGCMin et pcGCMax
// retourne 1 si la condition est verifiee; 0 sinon

int percentGC (char *primer, int primerLength)
{
  int i, k, x;
  k = 0;			
  for (i = 0; i < primerLength; i++)
    {
      if ((primer[i] == 'G') || (primer[i] == 'C') || (primer[i] == 'c') || (primer[i] == 'g'))
	k++;
    }				
  x = k * 100;
  x = x / primerLength;
  return x;
}

//repeat: indique si le primer contient des repetitions du type (x)^nbRepeat ou (xy)^nbRepeat

int repeat (char *primer, int primerLength)
{
  int i, k, t, l, r, x;
  char pred1, pred2;

  k = 1;
  t = 1;
  l = 0;
  r = 0;
  pred1 = '\n';
  pred2 = '\n';

  for (i = 0; i < primerLength - 1; i++)
    {
      if (primer[i] == primer[i + 1]) 
	k++;
      else if (pred1 == primer[i] && pred2 == primer[i + 1])
	{
	  t++;
	  if ((t%2)==1) l = (t+1)/2 + 1;
	  else l = t / 2 + 1;
	}
      else
	{
	  r = max(r,max(k,l));
	  k = 1; 
	  t = 1;
	}

      pred1 = pred2;
      pred2 = primer[i];
    }
  x = max(r,max(k,l)); 
  return x;
}


// pour tout caractere 'a' etant un nucleotide rend son complementaire

char comp (char a){
  switch (a)
    {
    case 'a' :  return 't';
    case 'c' :  return 'g';
    case 'g' :  return 'c';
    case 't' :  return 'a';
    case 'n' :  return 'n';
    case 'A' :  return 'T';
    case 'C' :  return 'G';
    case 'G' :  return 'C';
    case 'T' :  return 'A';
    case 'N' :  return 'N';
    default  :
      printf ("ERROR function comp (filters.h) -- %d %c\n",a,a);
      exit (0);
    }
}

// indique si l'amorce ne peut pas s'apparier en "epingle a�cheveux"

int hairpin (char *primer, int primerLength)
{
  int i, j, k, x, r;
  int bi, bj1;

  bj1 = 2 * maxHpDup + maxHpLoop;
  bi = primerLength - bj1;
  r=0;

  for (i = 0; i < bi; i++)
    {
      for (j = i + bj1; j < primerLength; j++)
	{
	  k=0;
	  while ((j-k-1-k-i>maxHpLoop)&&(primer[i + k] == comp (primer[j - k - 1]))) k++;
	  r=max(r,k);
	}
    }
  return r;
}



