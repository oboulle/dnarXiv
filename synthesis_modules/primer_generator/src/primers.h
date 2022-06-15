#define SIZE_MAX_PRIMERS 128

typedef struct {
  int size;               // taille totale de la banque : commentaire + nucleotides
  int nb_seq;             // nombre de sequences dans la banque
  int nb_nt;              // nombre de nucleotides
  int *idx_offset;        // tableau d'offset qui indique le debut de chaque sequence (1er caractere)
  int *idx_size;          // tableau qui indique la taille des sequences
  char *seq;              // stockage de la banque
  int  *iseq;             // pointeurs sur les memes graines
  int  *tseq;             // tableau initial de NB_DIF_SEED pointeurs sur *iseq 
} TYPE_BANK;

typedef struct {
  int deltaG_max;
  int deltaG_3p;
  int loc;
} TYPE_HYBR;


int compare_primers (void const *p1, void const *p2)
{
  return strcmp (*(char **)p1,*(char **)p2);
}

int MemError()
{
  fprintf (stderr,"\nMEMORY ERROR: ");
  fprintf (stderr,"nor enougth memory\n\n"); exit (0);
}

int min (int a, int b)
{
  if (a<b) return a; else return b;
}

int max (int a, int b)
{
  if (a>b) return a; else return b;
}

// donne la taille en octet du fichier

int getSize(filename)
     char * filename;
{
  FILE *fseq;
  int k;

  if ((fseq=fopen(filename,"r"))==NULL)
    {
      fprintf (stderr,"cannot open %s\n",filename); exit (0);
    }

  fseek (fseq,0,SEEK_END);
  k = ftell(fseq);
  fclose (fseq);
  return k;
}

char char2char(c)
     char c;
{
  switch (c)
    {
    case 'a' : return 'A';
    case 'c' : return 'C';
    case 'g' : return 'G';
    case 't' : return 'T';
    case 'A' : return 'A';
    case 'C' : return 'C';
    case 'G' : return 'G';
    case 'T' : return 'T';
    default  : return 'N';
    }
}

// lit un fichier de sequences d'ADN au format FASTA
// rempli les champs de la structure TYPE_BANK
//    - int size   : taille totale de la banque : commentaire + nucleotides
//    - int nb_nt  : nombre de nucleotides
//    - int nb_seq : nombre de sequences dans la banque
//    - char *seq  : la banque est stockee dans un tableau de caracteres (commentaire + nt)
//                   >com1 toto>ATGGACCCAGGATTAGC>com2 titi>TTAGGACCAGGATA>
//                   les caracteres < et > dans les commentaires sont remplaces pas [ et ]
//    - int *idx_offset : tableau d'offset qui indique le debut de chaque sequence (1er caractere)
//    - int *idx_size   : tableau qui indique la taille des sequences

int getDnaSeq (bk,filename)
     TYPE_BANK *bk;
     char *filename;
{
  FILE *fseq;
  char c;
  int i,k,comment,nbnt;
  int *tmp1, *tmp2;

  if ((tmp1 = (int *)malloc(1000000*sizeof(int)))==NULL) MemError();
  if ((tmp2 = (int *)malloc(1000000*sizeof(int)))==NULL) MemError();

  fseq=fopen(filename,"r");

  bk->nb_seq=0;
  comment=0; 
  k=0;
  nbnt=0;
  while (fscanf (fseq,"%c",&c)!=EOF)
    {
      if ((c == '>')&&(comment == 0)) { comment = 1; tmp2[bk->nb_seq] = 0; bk->nb_seq++; bk->seq[k++]='>'; }
      else if (comment == 1)
	{
	  if (c == '>') bk->seq[k++]=']';
	  else if (c == '<') bk->seq[k++]='[';
	  else
	    {
	      if (c=='\n') { comment = 0; bk->seq[k++]='>'; tmp1[bk->nb_seq-1]=k; }
	      else bk->seq[k++]=c;
	    }
	}
      else
	{
	  if ((bk->nb_seq>0) && ( ((c>='a')&&(c<='t')) || ((c>='A')&&(c<='T'))))
	    {
	      bk->seq[k++]= char2char(c); 
	      nbnt++; 
	      tmp2[bk->nb_seq-1]++;
	    }
	}
    }
  if ((bk->idx_offset = (int *) malloc((bk->nb_seq+1)*sizeof(int)))==NULL) MemError();
  if ((bk->idx_size   = (int *) malloc((bk->nb_seq+1)*sizeof(int)))==NULL) MemError();

  for (i=0; i<bk->nb_seq; i++) bk->idx_offset[i]=tmp1[i];
  for (i=0; i<bk->nb_seq; i++) bk->idx_size[i]=tmp2[i];

  free(tmp1); free(tmp2);

  bk->seq[k++]='>';
  bk->size=k;
  bk->nb_nt=nbnt;
  fclose(fseq);
  return k;
}

int getLine(FILE *ff, char *line)
{
  int i, endofline;
  char c;

  endofline = 0;
  i=0;
  while (endofline==0)
    {
      if (fscanf(ff,"%c",&c)==EOF) return 0;
      if (c=='\n') endofline=1;
      else if (c>=20) line[i++]=c;
    }
  line [i]='\0';
  return i;
}


