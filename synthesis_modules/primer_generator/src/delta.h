// rend pour tout couples de bases voisines une valeur indiquant la stabilite de la liaison
// en fonction de la nature des 2 bases

int valDeltaG (char a, char b)
{
  if (a == 'A')
    {
      if (b == 'A')
	return -1000;
      if (b == 'C')
	return -1440;
      if (b == 'G')
	return -1280;
      if (b == 'T')
	return -880;
    }
  else if (a == 'C')
    {
      if (b == 'A')
	return -1450;
      if (b == 'C')
	return -1840;
      if (b == 'G')
	return -2170;
      if (b == 'T')
	return -1280;
    }
  else if (a == 'G')
    {
      if (b == 'A')
	return -1300;
      if (b == 'C')
	return -2240;
      if (b == 'G')
	return -1840;
      if (b == 'T')
	return -1440;
    }
  else
    {
      if (b == 'A')
	return -580;
      if (b == 'C')
	return -1300;
      if (b == 'G')
	return -1450;
      if (b == 'T')
	return -1000;
    }
  return -1;
}


// calcul du deltaG de l'extremite 5'

int deltaGExt5Val (char *primer)
{
  int i;
  int dgExt5;
  dgExt5 = 0;
  for (i = 0; i < sizeExt5 - 1; i++)
    {
      dgExt5 = dgExt5 + valDeltaG (primer[i], primer[i + 1]);
    }
  return dgExt5;
}

// calcul du deltaG de l'extremite 3'

int deltaGExt3Val (char *primer, int primerLength)
{
  int i;
  int dgExt3;
  dgExt3 = 0;
  for (i = primerLength - sizeExt3; i < primerLength - 1; i++)
    {
      dgExt3 = dgExt3 + valDeltaG (primer[i], primer[i + 1]);
    }
  return dgExt3;
}

// indique si l'amorce a des valeurs de stabilit� pour ses extremit�s 3' et 5'
// compatible avec les valeurs en parametre deltaG5, deltaG3min et deltaG3max

int deltaG (char *primer, int primerLength)
{
  int dgExt5 = deltaGExt5Val (primer);
  int dgExt3 = deltaGExt3Val (primer, primerLength);
  if ((dgExt5 <= deltaG5) && (dgExt3 <= deltaG3max) && (dgExt3 >= deltaG3min)) return 1;
  return 0;
}


//valcomp : table des valeurs thermodynamiques pour un internal single mismatch

int valComp (char pred2, char pred1, char n2, char n1)
{
  if (pred2 == 'G' && pred1 == 'C')
    {
      if (n2 == 'G')
	{
	  if (n1 == 'G')
	    return -1110;
	  else if (n1 == 'C')
	    return -1840;
	  else if (n1 == 'A')
	    return -520;
	  else
	    return 80;
	}
      else if (n2 == 'C')
	{
	  if (n1 == 'G')
	    return -2240;
	  else if (n1 == 'C')
	    return 790;
	  else if (n1 == 'A')
	    return 470;
	  else
	    return 620;
	}
      else if (n2 == 'A')
	{
	  if (n1 == 'G')
	    return -250;
	  else if (n1 == 'C')
	    return 810;
	  else if (n1 == 'A')
	    return 170;
	  else
	    return -1300;
	}
      else
	{
	  if (n1 == 'G')
	    return -590;
	  else if (n1 == 'C')
	    return 980;
	  else if (n1 == 'A')
	    return -1440;
	  else
	    return 450;
	}
    }
  else if (pred2 == 'C' && pred1 == 'G')
    {
      if (n2 == 'G')
	{
	  if (n1 == 'G')
	    return -110;
	  else if (n1 == 'C')
	    return -2170;
	  else if (n1 == 'A')
	    return 110;
	  else
	    return -470;
	}
      else if (n2 == 'C')
	{
	  if (n1 == 'G')
	    return -1840;
	  else if (n1 == 'C')
	    return 700;
	  else if (n1 == 'A')
	    return 790;
	  else
	    return 620;
	}
      else if (n2 == 'A')
	{
	  if (n1 == 'G')
	    return 30;
	  else if (n1 == 'C')
	    return 750;
	  else if (n1 == 'A')
	    return 430;
	  else
	    return -1450;
	}
      else
	{
	  if (n1 == 'G')
	    return -320;
	  else if (n1 == 'C')
	    return 400;
	  else if (n1 == 'A')
	    return -1280;
	  else
	    return -120;
	}
    }
  else if (pred2 == 'A' && pred1 == 'T')
    {
      if (n2 == 'G')
	{
	  if (n1 == 'G')
	    return -130;
	  else if (n1 == 'C')
	    return -1280;
	  else if (n1 == 'A')
	    return 20;
	  else
	    return 710;
	}
      else if (n2 == 'C')
	{
	  if (n1 == 'G')
	    return -1440;
	  else if (n1 == 'C')
	    return 1330;
	  else if (n1 == 'A')
	    return 770;
	  else
	    return 640;
	}
      else if (n2 == 'A')
	{
	  if (n1 == 'G')
	    return 140;
	  else if (n1 == 'C')
	    return 880;
	  else if (n1 == 'A')
	    return 610;
	  else
	    return -1000;
	}
      else
	{
	  if (n1 == 'G')
	    return 70;
	  else if (n1 == 'C')
	    return 730;
	  else if (n1 == 'A')
	    return -880;
	  else
	    return 690;
	}
    }
  else if (pred2 == 'T' && pred1 == 'A')
    {
      if (n2 == 'G')
	{
	  if (n1 == 'G')
	    return 440;
	  else if (n1 == 'C')
	    return -1450;
	  else if (n1 == 'A')
	    return 740;
	  else
	    return 430;
	}
      else if (n2 == 'C')
	{
	  if (n1 == 'G')
	    return -1300;
	  else if (n1 == 'C')
	    return 1050;
	  else if (n1 == 'A')
	    return 1330;
	  else
	    return 970;
	}
      else if (n2 == 'A')
	{
	  if (n1 == 'G')
	    return 420;
	  else if (n1 == 'C')
	    return 920;
	  else if (n1 == 'A')
	    return 690;
	  else
	    return -580;
	}
      else
	{
	  if (n1 == 'G')
	    return 340;
	  else if (n1 == 'C')
	    return 750;
	  else if (n1 == 'A')
	    return -1000;
	  else
	    return 680;
	}
    }
  else
    return 0;
}


int autocomp2 (char *primer, int sizeprimer, char *invprimer, int *DELTAG)
{
  int x, x1, x2;
  int ps, pe, is , ie;
  int p1, i1;
  int deltaG, deltaGT, deltaG3;

  x1 = sizeprimer-sizeDeltaGAuto;
  x2 = sizeDeltaGAuto-1;


  for (x=0; x<(sizeprimer-sizeDeltaGAuto)*2 + 1; x++)
    {
      ps = max(x1,0); 
      pe = min(x1+sizeprimer-1,sizeprimer-1);
      is = max(0,x2-sizeprimer+1);                
      ie = min(x2, sizeprimer-1);

      p1 = pe-1; i1 = ie-1;
      deltaG = deltaGAutoInit;
      deltaGT = deltaG;
      deltaG3 = deltaG;

      while (p1>=ps)
	{
	  if (primer[p1]==comp(invprimer[i1]))
	    {
	      if (primer[p1+1]==comp(invprimer[i1+1]))
		{
		  // la pair courante et la pair pr�cedente matchent 

		  deltaG = deltaG + valComp(invprimer[i1+1],primer[p1+1],invprimer[i1],primer[p1]);
		  deltaGT = min (deltaG,deltaGT);
		  if (p1>=(sizeprimer-sizeDeltaGAuto3)) deltaG3 = min (deltaG,deltaG3); 
		}
	      else
		{
		  if ((p1<=pe-2)&&(primer[p1+2]==comp(invprimer[i1+2])))
		    {
		      // la pair courante match, pas la pr�c�dente, mais celle d'avant oui 
		      
		      deltaG = deltaG + valComp(invprimer[i1+1],primer[p1+1],invprimer[i1],primer[p1]);
		      deltaGT = min (deltaG,deltaGT);
		      if (p1>=(sizeprimer-sizeDeltaGAuto3)) deltaG3 = min (deltaG,deltaG3);
			
		    }
		  else
		    {
		      // la pair courante match, mais pas les 2 pr�c�dentes
		      deltaG = deltaGAutoInit;
		    }
		}
	    }
	  else
	    {
	      if ((p1>ps) && (primer[p1+1]==comp(invprimer[i1+1])) && (primer[p1-1]==comp(invprimer[i1-1])))
		{
		  // la pair courante est en mismatch, mais pas la suivante ni la pr�c�dente

		  deltaG = deltaG + valComp(invprimer[i1+1],primer[p1+1],invprimer[i1],primer[p1]);
		  deltaGT = min (deltaG,deltaGT);
		  if (p1>=(sizeprimer-sizeDeltaGAuto3)) deltaG3 = min (deltaG,deltaG3);
		}
	      else
		{
		  // autre cas : on remet deltaG � sa valeur d'initialisation
		  deltaG = deltaGAutoInit;
		}
	    }
	  p1--; i1--;
	}
      x1--; x2++;
    }
  DELTAG[0]=deltaGT; 
  DELTAG[1]=deltaG3;
  return 1;
}

