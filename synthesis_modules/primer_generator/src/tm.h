/* 
 * Tables of nearest-neighbor thermodynamics for DNA bases from Santalucia 1998.
 */

#define S_A_A 222
#define S_A_C 224
#define S_A_G 210
#define S_A_T 204

#define S_C_A 227
#define S_C_C 199
#define S_C_G 272
#define S_C_T 210
 
#define S_G_A 222
#define S_G_C 244
#define S_G_G 199
#define S_G_T 224

#define S_T_A 213
#define S_T_C 222
#define S_T_G 227
#define S_T_T 222

#define H_A_A  79
#define H_A_C  84
#define H_A_G  78
#define H_A_T  72

#define H_C_A  85
#define H_C_C  80
#define H_C_G 106
#define H_C_T  78

#define H_G_A  82
#define H_G_C  98
#define H_G_G  80
#define H_G_T  84

#define H_T_A  72
#define H_T_C  82
#define H_T_G  85
#define H_T_T  79

#define A_CHAR 'A'
#define G_CHAR 'G'
#define T_CHAR 'T'
#define C_CHAR 'C'

#define CATID5(A,B,C,D,E) A##B##C##D##E
#define CATID2(A,B) A##B
#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dh += CATID5(H,_,LAST,_,THIS); \
     ds += CATID5(S,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

#define STATE(LAST)     \
   CATID2(LAST,_STATE): \
   c = *s; s++;         \
   DO_PAIR(LAST,A)      \
   else DO_PAIR(LAST,T) \
   else DO_PAIR(LAST,G) \
   else DO_PAIR(LAST,C) \
   else if ('\0' == c)  \
             goto DONE; \
   else goto ERROR1 \
   

double 
oligotm(s, DNA_nM, Na_mM, size)
     const char *s;
     double DNA_nM;
     double Na_mM;
     int size;
{
  register int dh = -2, ds = 71; // dh -> init = +0.2 kcal.mol-1 + sym correction = 0.0 pour self dup
                                   // ds -> init = -5.7 e.u. + sym correction = -1.4 pour self dup
  register char c;
  double delta_H, delta_S;

  /* Use a finite-state machine (DFA) to calucluate dh and ds for s. */
  c = *s; s++;
  if (c == 'A') goto A_STATE;
  else if (c == 'G') goto G_STATE;
  else if (c == 'T') goto T_STATE;
  else if (c == 'C') goto C_STATE;
  else goto ERROR2;
  STATE(A);
  STATE(T);
  STATE(G);
  STATE(C);
  
 DONE:  /* dh and ds are now computed for the given sequence. */
  delta_H = dh * -100.0;  
  delta_S = ds * -0.1;     
  
  // justification : san98 for salt correction for an oligonucleotide (delta_H unaffected)
  // and DNA_nM/1.10e9 in place of 4.10e9 because the Tm calculated concerns self-complementary molecules
  
  return (delta_H / ((delta_S+(0.368*(size-1.0)*log(Na_mM/1000.0))) + 1.9872 * log(DNA_nM/4000000000.0))) - 273.15;
  //return delta_H / (delta_S + 1.9872 * log(DNA_nM/4000000000.0)) - 273.15 + 16.6 * log10(Na_mM/1000.0);
  
  
 ERROR1:  
  printf("ERROR #1:\n"); return 0.0;

 ERROR2: 
  printf("ERROR #2: %c \n",c); exit (0); return 0.0;
}
#undef DO_PAIR

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
    dg += CATID5(G,_,LAST,_,THIS);  \
    goto CATID2(THIS,_STATE);	    \
  }

