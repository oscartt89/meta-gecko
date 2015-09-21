/* fromTxttoMAT.c

Este programa coge un fichero de 3 columnas y lo mete en una matriz de acomulado.
La primera columna es de longitud
la segudna es de similarity
la tercera es el coverage que cada fragmento tiene

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>

int MAXLEN = 1001;

int main(int ac,char** av){

if(ac<4){
	printf("Use: fromTxt2MAT archivotxt sem.MAT cov.MAT \n");
	exit(-1);
}

/** Open files */
FILE* fe,*fsemMAT,*fcovMAT;
// Open features
if((fe=fopen(av[1],"r"))==NULL){
		printf("*****ERROR opening %s file",av[1]);
		exit(-1);
	}
// Open new seq

if((fsemMAT=fopen(av[2],"w"))==NULL){
		printf("*****ERROR opening %s file",av[2]);
		exit(-1);
	}
	
// Open new seq

if((fcovMAT=fopen(av[3],"w"))==NULL){
		printf("*****ERROR opening %s file",av[3]);
		exit(-1);
	}

/** creates a matrix */
// Creamos la matriz
int i,j;
int** matrixSEM;
matrixSEM = (int**)malloc(sizeof(int*)*MAXLEN);
for(i=0;i<MAXLEN;i++){
	matrixSEM[i] = (int*)malloc(sizeof(int)*101);
	for(j=0;j<101;j++)matrixSEM[i][j]=0;
}

int** matrixCOV;
matrixCOV = (int**)malloc(sizeof(int*)*MAXLEN);
for(i=0;i<MAXLEN;i++){
	matrixCOV[i] = (int*)malloc(sizeof(int)*101);
	for(j=0;j<101;j++)matrixCOV[i][j]=0;
}

/** Read file */	
int length,sem,cov;
char linea[1000];
int max=1000;	

fgets(linea,max,fe);
	while(!feof(fe)){
		if(sscanf(linea,"%d\t%d\t%d\n",&length,&sem,&cov)==3){
			//printf("%d\t%d\t%d\n",length,sem,cov);
			if(length<1000){
				matrixSEM[length][sem]++;
				matrixCOV[length][cov]++;
			}
		}
		fgets(linea,max,fe);
		
	}

	
for(i=0;i<MAXLEN;i++){
	for(j=0;j<101;j++){
		fprintf(fsemMAT,"%d\t%d\t%d\n",i,j,matrixSEM[i][j]);
		fprintf(fcovMAT,"%d\t%d\t%d\n",i,j,matrixCOV[i][j]);
	}
}

fclose(fsemMAT);
fclose(fcovMAT);


exit(-1);
}
