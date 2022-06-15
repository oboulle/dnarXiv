
//dangling5: donne la valeur de l'extremite flottante en 5'
//ref : San2004 : Next Nearest Neighbor model (NNN)

int dangling5 (char ext2, char ext1, char dang2)
{
  if (ext2 == 'A' && ext1 == 'T')
    {
      if (dang2 == 'A') return -510;
      else if (dang2 == 'C') return -420;
      else if (dang2 == 'G') return -620;
      else return -710;
    }
  else if (ext2 == 'C' && ext1 == 'G')
    {
      if (dang2 == 'A') return -960;
      else if (dang2 == 'C') return -520;
      else if (dang2 == 'G') return -720;
      else return -580;
    }
  else if (ext2 == 'G' && ext1 == 'C')
    {
      if (dang2 == 'A') return -580;
      else if (dang2 == 'C') return -340;
      else if (dang2 == 'G') return -560;
      else return -610;
    }
   else if (ext2 == 'T' && ext1 == 'A')
    {
      if (dang2 == 'A') return -500;
      else if (dang2 == 'C')return -20;
      else if (dang2 == 'G')return 480;
      else return -100;
    }
  else
    return 0;
}

//dangling3: donne la valeur de l'extremite flottante en 3'
//ref : San2004 : Next Nearest Neighbor model (NNN)

int dangling3 (char ext2, char ext1, char dang2)
{
  if (ext2 == 'A' && ext1 == 'T')
    {
      if (dang2 == 'A') return -120;
      else if (dang2 == 'C') return 280;
      else if (dang2 == 'G') return -10;
      else return -130;
    }
  else if (ext2 == 'C' && ext1 == 'G')
    {
      if (dang2 == 'A') return -820;
      else if (dang2 == 'C') return -310;
      else if (dang2 == 'G') return -10;
      else return -520;
    }
  else if (ext2 == 'G' && ext1 == 'C')
    {
      if (dang2 == 'A') return -920;
      else if (dang2 == 'C') return -230;
      else if (dang2 == 'G') return -440;
      else return -350;
    }
   else if (ext2 == 'T' && ext1 == 'A')
    {
      if (dang2 == 'A') return -480;
      else if (dang2 == 'C') return -190;
      else if (dang2 == 'G') return -500;
      else return -290;
    }
  else
    return 0;
}

int check_hybrid(char *primer, int sizeprimer, char *seq, char sseq, TYPE_HYBR *hy)
{
  int i,j;

  float Na_mM, saltCorrFloat;
  int   predPri, predSeq, saltCorr, deltaGI;
  int   nbMisMatch, delta3p, deltaT;
  char cx, lcx;


  // initialisation des predecesseurs

  predPri = primer[sizeprimer-1];
  predSeq   = seq[sizeprimer-1];

  // initialisation NbMisMatch

  if (predPri != comp(predSeq)) nbMisMatch = 1; else nbMisMatch = 0;
  
  //correction salt dependancy d'un duplex (San1998)

  Na_mM = (float) saltConc;
  saltCorrFloat = 175.0 * log(Na_mM/1000.0) - 200.0;
  saltCorr = (int) saltCorrFloat;

  //intialisation de deltaGI apres corrections
  
  deltaGI = deltaGHybridInit + saltCorr + dangling3(predSeq, predPri, seq[sizeprimer]);
  deltaT  = deltaGI;
  delta3p = deltaGI;

  // calcul thermodynamique en partant de l'extrémité 3'

  j = sizeprimer -2;


  while (j>=0)
    {
      cx  = comp(seq[j]);
      if (j>0) lcx = comp(seq[j-1]);

      // 0) les 2 nucleotides courants et la paire suivante sont en mismatch mais pas la paire precedente :

      if ((j>0) && (nbMisMatch == 0) && (primer[j] != cx) && (primer[j-1] != lcx))
	{
	  nbMisMatch++;
	  predPri = primer[j];
	  predSeq = seq[j];
	  j--;
	  continue;
	}

      // 1) les 2 nucleotides courants sont en mismatch mais la paire precedente et la suivante sont en match :

      else if ((j>0) && (nbMisMatch == 0) && (primer[j] != cx) && (primer[j-1] == lcx))
	{
	  deltaGI = deltaGI + valComp (predSeq, predPri,seq[j],primer[j]) 
	                    + valComp (primer[j-1],seq[j-1], primer[j],seq[j]);

	  if (deltaT > deltaGI) deltaT = deltaGI;
	  if (sizeprimer - 1 - j <= sizeDeltaGHybrid3) if (delta3p > deltaGI) delta3p = deltaGI;

	  nbMisMatch = 0;
	  predPri = primer[j-1];
	  predSeq = seq[j-1];
	  j = j-2;;
	  continue;
	}

      // 2) les 2 nucleotides courants sont en mismatch ainsi que la paire precedente :

      else if ((nbMisMatch == 1) && (primer[j] != cx))
	{
	  deltaGI= deltaGHybridInit + saltCorr;
	  predPri = primer[j];
	  predSeq = seq[j];
	  nbMisMatch++;
	  j--;
	  continue;
	}

      // 3) les 2 nucleotides courants sont en mismatch ainsi que plus d'une paire precedente :

      else if ((nbMisMatch > 1) && (primer[j] != cx))
	{
	  predPri = primer[j];
	  predSeq =  seq[j];
	  nbMisMatch++;
	  j--;
	  continue;
	}      

      // 4) les 2 nucleotides courants sont en match mais pas la ou les paire precedente :

      else if ((nbMisMatch >= 1) && (primer[j] == cx))
	{
	  predPri = primer[j];
	  predSeq = seq[j];
	  nbMisMatch = 0;
	  j--;
	  continue;
	}

      // 5) les 2 nucleotides courants et la paire precedente sont en match :

      else 
	{
	  deltaGI = deltaGI + valComp (predSeq, predPri, seq[j],primer[j]);
	  nbMisMatch = 0;

	  if (deltaT > deltaGI) deltaT = deltaGI;
	  if (sizeprimer - 1 - j <= sizeDeltaGHybrid3) if (delta3p > deltaGI) delta3p = deltaGI;

	  predPri = primer[j];
	  predSeq = seq[j];
	  j--;
	}
    }

  deltaGI = deltaGI + dangling5(seq[0],primer[0],sseq);
  if (deltaT > deltaGI) deltaT = deltaGI;

  //printf ("--> %d\n",deltaT);

  hy->deltaG_max = deltaT;
  hy->deltaG_3p  = delta3p;

  return deltaT;
}
