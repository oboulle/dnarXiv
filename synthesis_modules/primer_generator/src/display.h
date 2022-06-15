
void printResults (FILE *fres, char *com, char *primer, int sizePrimer, int idxPrimer, int sizeSeq)
{
  int i;

  if (sizeSeq==0)
    fprintf (fres,">%s  [%d %d]",com, idxPrimer, idxPrimer+sizePrimer-1);
  else
    fprintf (fres,">%s  [%d %d]",com, sizeSeq-idxPrimer-1, sizeSeq-(idxPrimer+sizePrimer));
  fprintf (fres,"\n");
  for (i=0; i<sizePrimer; i++) fprintf (fres,"%c",primer[i]); fprintf (fres,"\n");
}

void printHybrid (FILE *fko, int nbHybrid, char *com, char *primer, char *seq, int sizeprimer, int idxg, TYPE_HYBR *hy, int sens)
{
  int j;

  if (nbHybrid==1) 
    {
      fprintf (fko,"\n%s\n%s\n",com,primer);
      for (j=0; j<sizeprimer; j++) fprintf (fko," "); fprintf (fko,"      ");
      fprintf (fko,"\tseq \t   start      end\tdG max \t dG 3'\n\n");
    }

  if (sens==0)
    {
      fprintf (fko,"5' ");
      for (j=0; j<sizeprimer; j++) fprintf (fko,"%c",primer[j]);
      fprintf (fko," 3' primer\n");      // fprintf (fko," 3' primer\n",primer); 2022-05-11
      fprintf (fko,"   ");
      for (j=0; j<sizeprimer; j++) if (primer[j]==seq[j]) fprintf (fko,"|"); else fprintf (fko," "); 
      fprintf (fko,"\n3' ");
      for (j=0; j<sizeprimer; j++) fprintf (fko,"%c",comp(seq[j]));
      fprintf (fko," 5' ");
    }
  else
    {
      fprintf (fko,"3' ");
      for (j=0; j<sizeprimer; j++) fprintf (fko,"%c",primer[sizeprimer-j-1]);
      fprintf (fko," 5' primer\n"); // fprintf (fko," 5' primer\n",primer); 2022-05-11
      fprintf (fko,"   ");
      for (j=0; j<sizeprimer; j++) if (primer[sizeprimer-j-1]==comp(seq[j])) fprintf (fko,"|"); else fprintf (fko," "); 
      fprintf (fko,"\n5' ");
      for (j=0; j<sizeprimer; j++) fprintf (fko,"%c",seq[j]);
      fprintf (fko," 3' ");
    }

  fprintf (fko,"genome \t%8d %8d",hy->loc,hy->loc+sizeprimer-1);
  fprintf (fko,"\t%6d \t%6d",hy->deltaG_max,hy->deltaG_3p);
  fprintf (fko,"\n\n");
}


