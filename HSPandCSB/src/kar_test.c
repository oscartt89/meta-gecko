/* Kar_test.c : Pruebas del alg. de karlin
   -------------------------------------[ABR.95]--*/

#include "../otsSoft/global.h"
int fin_db, fin_buffer;
#include "../otsSoft/seqalea.c"
#include "karlin.c"
void impreseq();


main(argc,argv)
int argc;char *argv[];
{

struct Sequence sX, sY;
int n0,n1;
INToFLO PAM[MAXALF][MAXALF];
char Ax[MAXALF],Ay[MAXALF];
int lAx,lAy;
int Tx[MAXALF],Ty[MAXALF];
char nX[MAXLS],nY[MAXLS];
float FrX[MAXALF],FrY[MAXALF], FrA[MAXALF];
struct parametros p;
int i,j;
INToFLO pMax,pMin;
char seqA[MAXLS];
int nA=2000;
char nseqA[MAXLS];
float FFx[MAXALF],FFy[MAXALF];
float *ProbScoX;
  

LeerParametros(argc,argv,&p);
/*--
  printf("\nSeq Horizontal =%s",p.nomfichbd);
  printf("\nsec Vertical   =%s",p.nomfichx);
  printf("\npam  =%s",p.nomfscor);
  printf("\nt_pam=%d",p.TipoPam);
---*/


if(!(n0 = LeeSeqX(p.nomfichx,&sX))) terror("leeSx");
if(!(n1 = LeeSeqX(p.nomfichbd,&sY)))terror("leeY");
printf("\nN0=%d seq=%s",n0,sX.ident); impreseq(sX.datos);
printf("\nN1=%d seq=%s",n1,sY.ident); impreseq(sY.datos);
if (!(LoadScores(p.nomfscor,p.TipoPam,Ax,Ay,&lAx,&lAy,PAM))) terror("pam");;
Traductor(Ax,Tx);
Traductor(Ay,Ty);
if(!SeqToNum(sX.datos,Tx,nX)) terror("al traducir");;
if(!SeqToNum(sY.datos,Ty,nY)) terror("Al trad2");;


MaxMinPam(PAM,lAx,lAy,&pMax,&pMin);
Frequency(nX,n0,Tx,FrX,lAx);
Frequency(nY,n1,Ty,FrY,lAy);
/*---Analisis de las freq positivas y negativas en la otra seq--*/
/*----lo de KARLIN----
boolean	karlin(int low,int high,double *pr,double *lambda,double *K,double *H)
---*/
 
/*--
printf("\nnuevo modulo aleatorios---------\n");
GenSeqRan(seqA, nA, Ax, lAx);
impreseq(seqA);
if(!SeqToNum(seqA,Tx,nseqA)) terror("al traducir");;
Frequency(nseqA,nA,Tx,FrA,lAx);
for (i=0;i<lAx;i++) printf("\n %c = %f ",Ax[i],FrA[i]);
printf("\nnuevo modulo aleatorios [2]---------\n");
GenSeqRanFreq(seqA, nA, Ax, lAx, FrX);
impreseq(seqA);
if(!SeqToNum(seqA,Tx,nseqA)) terror("al traducir");;
Frequency(nseqA,nA,Tx,FrA,lAx);
for (i=0;i<lAx;i++) printf("\n %c = Orig=%f  Resul=%f ",Ax[i],FrX[i],FrA[i]);

--*/

printf("\n-------------------\n");

}
void impreseq(s)
  char *s;
{
  int i;
  printf("\nSeq: ");
  for (i=0; i<50 && s[i]; i++) printf("%c",s[i]);
  printf("...\n");
  return;
}

