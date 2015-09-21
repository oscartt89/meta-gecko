/* test2 de Karlin (Feb.2000)
    Sintaxis: kar2test File.Freq  Score.Mat  Tipo(0=Triangular/1=Cuadrada con 2 alfabetos)

     Freq.file   : per-residue frequency. 
                   Format:
                     Residue    frequency
                   (note: frequency sum must be 1.0)
     Score.mat   : Per-residue score matrix
                   Format:
                     [0] : Alphabeth
                           Score matrix triangular
                     [1] : Alphabeth
                           Alphabeth
                           Score matrix cuadrada
     Tipo de matrix: 0=Triangular
                     1=Cuadrada con 2 alfabetos


      Para obtener detalles de la ejecucion, compilar con 
      la opcion -DDEBUG

   e-mail: ots@ac.uma.es
                                       Feb.96
   -------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./karlin.c"

#define MAXALF   ('Z'-'A'+1)




void terror(char *s) { printf("\n** ERR ** %s **\n",s); exit(-1);}

int **LoadScores(filename, Tipo, AlfabX,Alfab2, lAlfaX,lAlfa2)
    char *filename, *AlfabX, *Alfab2;
    int Tipo, *lAlfaX, *lAlfa2;
    /* Tipo : 0 = triangular con 1 alfabeto
              1 = Rectangular de 2 alfabetos
       ---------------------------------------*/
    {
    int i,j;
    int fx;
    FILE *fich;
    int **P;
    

   for(i=0;i<MAXALF;i++) AlfabX[i]=Alfab2[i]=0x00;

   if ((fich = fopen(filename, "rt"))==NULL)
       terror("can not open Scores file");
   while(1) /* Jun.95-------insert initial comments--*/
   { fscanf(fich,"%s",Alfab2);   
     if (Alfab2[0]!='#') break;
   }
   *lAlfa2=(int)strlen(Alfab2);
   if (Tipo) { fscanf(fich,"%s\n",AlfabX);
               *lAlfaX =(int)strlen(AlfabX); }
       else  { strcpy(AlfabX, Alfab2);
               *lAlfaX = *lAlfa2; }
   if (*lAlfaX >MAXALF) 
      { printf("Alfab. too long in Scores file (%s)",AlfabX );
        return 0; }
   if (*lAlfa2 > MAXALF) 
      { printf("Alfab. too long in Scores file (%s)",Alfab2 );
        return 0; }

  /* memory for Pam matrix ------------------*/

  if ((P=(int**) calloc(*lAlfaX, sizeof(int*)))==NULL) return NULL;
  for (i=0;i<*lAlfaX;i++)
    if ((P[i]=(int*) calloc(*lAlfa2, sizeof(int)))==NULL) return NULL;

   if (Tipo) 
     { for (i=0;i<*lAlfaX; i++)
          for (j=0;j<*lAlfa2 ; j++)
          {
            fscanf(fich,"%d",&fx);
            P[i][j]= fx;
          }
      } 
      else
      { for (i=0;i<*lAlfaX; i++)
          for (j=0;j<=i ; j++)
          {
            fscanf(fich,"%d",&fx);
            P[i][j]=P[j][i]=fx;
          }
      } 
   
/* be sure they are upercase---*/
    for (i=0;i<*lAlfaX;i++) AlfabX[i]=toupper(AlfabX[i]);
    for (i=0;i<*lAlfa2;i++) Alfab2[i]=toupper(Alfab2[i]);

    fclose(fich);
    return P;
}
int main( int argc,char *argv[])
{


/* Open output file */	
	FILE* fout;
	if((fout=fopen(argv[4],"wb"))==NULL){
		printf("***ERROR Opening output file %s\n",argv[4]);
		exit(-1);
	}
	
	
    int **PAM,TipoPam;
    char Ax[MAXALF],Ay[MAXALF], linea[500],le;
    int lA=MAXALF,lAx,lAy;
    int i,j,k=0;
    FILE *fIn;
    double Fr[MAXALF];
    double valtot=0.;
    float valor,suma;
    int xlow=0,xhigh=0, rango;
    double *pr,lambda,K,H;

    if (argc!=5) terror("Use ExpectVal File.freq  Score.mat TipoMat(0=Triang/1=Cuadrada karlin.out");

    TipoPam=atoi(argv[3]);
   /* ----Load Score Matrix-----------------------*/
   if ((PAM=LoadScores(argv[2],TipoPam,Ax,Ay,&lAx,&lAy))==NULL)
       terror("loading score-matrix");;

   for (i=0;i<lAx;i++)
     for (j=0;j<lAy;j++){
      if (i==0 && j==0) { xlow=xhigh=PAM[i][j]; }
      if (PAM[i][j]<xlow) xlow=PAM[i][j];
      if (PAM[i][j]>xhigh) xhigh=PAM[i][j];
     }
   rango=xhigh-xlow+1;
   //printf("lowPAM=%d  highPAM=%d rango=%d\n",xlow,xhigh,rango);


   if ((pr=(double*) calloc(rango, sizeof(double)))==NULL)
      terror("memory for pr");


    /* read per-residue frquencies from file-----*/
   if ((fIn=fopen(argv[1],"rt"))==NULL) terror("opening frequency file");

   for(i=0;i<lA; i++) Fr[i] =0;

   fgets(linea,499,fIn);
   j=0;

   while(!feof(fIn)) {
    i=sscanf(linea,"%c %f \n",&le,&valor);
    if (i!=2 ) { //printf("warning...file format =%s",linea);
                fgets(linea,499,fIn); continue;
               }
    valtot+=valor;
    le=toupper(le);
    j=-1;
    for (i=0;i<lAx;i++) if (le==Ax[i]) {j=i;break;}
    if (j==-1) { printf("letra=%c Alfa=%s\n",le,Ax); 
                 terror("no esta letra en alfabeto");}
    Fr[j]=(double) valor;
    fgets(linea,499,fIn);
   }

   fclose(fIn);
   for(i=0;i<lA; i++) //printf("Freq.%d=%f\n",i,Fr[i]);

   //printf("freq.file reading...done (total=%f)\n",valtot);

   /*--Normalise---*/
   for(i=0;i<lA; i++) Fr[i] /=valtot;
//printf("Normializadas\n");
   for(i=0;i<j; i++) //printf("Freq.%d=%f\n",i,Fr[i]);

//printf("lax=%d lay=%d\n",lAx,lAy);
   /* calcula todas las probabilidades */
   for (i=0;i<lAx;i++)
     for (j=0;j<lAy;j++) {
      k =PAM[i][j]-xlow;
      pr[k]+=(double)Fr[i]*(double)Fr[j];
/*
printf("scoK=%d  Fr(%d)=%f  Fr(%d)=%d  prod=%f\n",k,
*/
     }
   valor=0;
   for (i=0;i<rango ;i++) {
     valor+=pr[i];
     //printf("Prob del valor %3d=%lf\n",(i+xlow),pr[i]);
   } 
   //printf("TOT=%lf\n",valor);
   
   /* score medio positivo (para el scomaximo*/
   valor=0;
   suma=0;
   for (i=0;i<lAx;i++) {
      valor+=PAM[i][i]*Fr[i]*Fr[i];
      suma+=Fr[i]*Fr[i];
   }
   //printf("valor=%f  scomediopositivo=%f \n",valor, valor/suma);


/*----lo de KARLIN----*/
karlin(xlow,xhigh,pr,&lambda,&K,&H);
/* ---*/
 

   //printf("\nLambda=%lf  K=%lf\n",lambda, K);
   
   fprintf(fout,"Lambda=%f  K=%f\n",lambda,K);
   fclose(fout);

   return 0;
}


