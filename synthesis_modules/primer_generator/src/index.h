int nt_code(char c)
{
  switch (c) 
    {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    default  : return -1;
    }
}

int seed_code(char *s)
{
  int j,x,z;

  x=0;
  for (j=0; j<SIZE_SEED; j++)
    {
      z = nt_code(s[j]); 
      if (z<0) return -1;
      x = x*4 + z; 
    }
  return x;
}

int getNumSeq(bk,ix)
     TYPE_BANK *bk;
     int ix;
{
  int i, imin, imax, ixx;

  if (ix<0) ixx = -ix; else ixx = ix;

  imin = 0;
  imax = bk->nb_seq;
  i = (imax - imin) / 2;

  while ((bk->idx_offset[i] > ixx) || (bk->idx_offset[i]+bk->idx_size[i] < ixx))
    {
      //printf ("*1* %d %d %d -- %d %d / %d\n",imin,imax,i,bk->idx_offset[i],bk->idx_offset[i]+bk->idx_size[i],ixx);
      if (bk->idx_offset[i] < ixx) imin = i; else imax = i;
      i = imin + (imax - imin) / 2;
      //printf ("*2* %d %d %d -- %d %d / %d\n",imin,imax,i,bk->idx_offset[i],bk->idx_offset[i]+bk->idx_size[i],ixx);
    }
  return i;
}



int indexDnaSeq(bk)
     TYPE_BANK *bk;
{
  int i,j,k,x,y,z,ok,ii;
  int start_dna, end_dna, maxi;
  int *cpt;

  // pour determiner la classe qui contient le + grand nb de graines
  if ((cpt = (int *) malloc (NB_DIFF_SEED*sizeof(int)))==NULL) MemError(); 

  for (i=0; i<NB_DIFF_SEED; i++) { bk->tseq[i] = -1; cpt[i]=0; }

  maxi=0;  k=1; 
  for (k=0; k<bk->nb_seq; k++)
    {
      start_dna = bk->idx_offset[k];
      end_dna   = bk->idx_size[k] + start_dna;

      for (i=start_dna; i<end_dna-SIZE_SEED+1; i++) 
	{
	  x = seed_code(&bk->seq[i]);
	  if (x>=0)
	    {
	      bk->iseq[i] = bk->tseq[x];
	      bk->tseq[x] = i;
	      cpt[x]++;
	      if (maxi<cpt[x]) maxi=cpt[x];
	    }
	}
    }

  free (cpt);
  return maxi;
}
