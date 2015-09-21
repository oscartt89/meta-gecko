/*
 * getHistogramByLength.c
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	



Parametro de entrada: Inicio, Final.

Muestra por pantalla los fragmentos que estan comprendidos en esa franja

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <dirent.h>
#include <math.h>
#include <inttypes.h>

#include "structs.h"
#include "comparisonFunctions.h"
#include "fragmentv3.h"
#include "fragmentv2.h"

#define MAXLINE   1000

struct Files { 
   char *fname;
   uint64_t nx, ny;;
};


int howManyFiles(char *);
struct Files* LoadListOfFiles(char*fListName, int *nF);

int main(int ac,char** av){

	if(ac<5){
	printf("\n Uso:\n\tgetHistogramByLength frags min max out.txt \n");exit(-1);
	}

	// Get parameters
	int min,max;
	min = atoi(av[2]);
	max = atoi(av[3]);
	//fin = atoi(av[3]);
	int i,j;
	// Read fragments
	struct FragFile* f;
	int nf; // Number of fragments
	uint64_t xtotal,ytotal;

	// Create a vector of sem.
	/*
	int* sim;
	sim = (int*)malloc(sizeof(int)*101);
	for(i=0;i<101;i++)sim[i]=0;
	*/
	
	// Create a matrix of length x sim
	int** matrix;
	matrix = (int**)malloc(sizeof(int*)*max);
	for(i=0;i<max;i++){
		matrix[i] = (int*)malloc(sizeof(int)*101);
		for(j=0;j<101;j++)matrix[i][j]=0;
	}
	
	/************/
	int ij,nFiles;
	struct Files *L;
	
	L = LoadListOfFiles(av[1], &nFiles);
	
	for (ij=0;ij<nFiles; ij++) { // process Frgs files

		printf("%s\n",L[ij].fname);
		nf=0;
		f=readFragments(L[ij].fname,&nf,&xtotal,&ytotal);
		for(i=0;i<nf;i++){
		/*
			if( (f[i].length >= ini) && (f[i].length <= fin) ){
				if((int)f[i].similarity<=100)sim[(int)f[i].similarity]++;
			}
			*/
			if( ((int)f[i].similarity<=100)&&((int)f[i].length<max)&&((int)f[i].length>min) ){
				matrix[(int)f[i].length][(int)f[i].similarity]++;
			}
		}
	}
	/******/
	// Filtro para suavizar
	int N=1;
	int k,l;
	int valor,n;
	
	for(i=N;i<max-N;i++){
		for(j=N;j<101-N;j++){
			valor=0;
			n=0;
			for(k=-N;k<=N;k++){
				for(l=-N;l<=N;l++){
					valor += matrix[i+k][j+l];
					n++;
				}
			}
			matrix[i][j]=valor/n;
		}
	}
	/***********PRINT***********/
	FILE* ftxt;
	if((ftxt=fopen(av[4],"w"))==NULL){
		printf("***ERROR Opening output file %s\n",av[4]);
		exit(-1);
	}
	// Print
	/*
	printf("semejanza\n");
	for(i=0;i<101;i++)printf("%d\n",sim[i]);
	*/
	
	fprintf(ftxt,"length\tsem\tcount\n");
	for(i=0;i<max;i++)for(j=0;j<101;j++)fprintf(ftxt,"%d\t%d\t%d\n",i,j,matrix[i][j]);
	fclose(ftxt);
	
	
	/*************/
	return 0;

}

/*******/
int howManyFiles(char *fname){
	FILE *f;
	char line[MAXLINE];
	int nL=0;

	if ((f=fopen(fname,"rt"))==NULL){
		 printf("Opening frags file LIST");
		 exit(-1);
	}
	   


        fgets(line,MAXLINE,f);
        while(!feof(f)) {
 		if (line[0]!='#' && (int)strlen(line)!=0) nL++;
	        fgets(line,MAXLINE,f);
	}
	fclose(f);
	return nL;
}
// Load to memory a list of files --------------------------------------------------
// datafile format: fileName[tab]nSeq[tab][format][newLINE]

struct Files* LoadListOfFiles(char*fListName, int *nF) {

        FILE *f, *ff;
        struct Files*L=NULL;
        char line[MAXLINE];
        int N=0,nFiles;
	uint64_t xnx,xny;
	

	nFiles = howManyFiles(fListName);
	if ((L=(struct Files*) calloc(nFiles,sizeof(struct Files)))==NULL) {
		printf("memory for list of files");
		exit(-1);
	}
           

        if ((f=fopen(fListName,"rt"))==NULL){ 
			printf("erro opening List of Files");
			exit(-1);
		}
        fgets(line,MAXLINE,f);
        while(!feof(f)) {

 		if (line[(int)strlen(line)-1]=='\n') line[(int)strlen(line)-1]=0x00;
		 
 		if (line[0]!='#' && (int)strlen(line)>2) {
			L[N].fname = strdup(line);         

			if ((ff=fopen(L[N].fname,"rt"))==NULL) {
				fprintf(stderr,"abriendo -->%s<--\n",L[N].fname);
                                   // rterror("Opening frags file from LIST X X XX");
                                continue;
			}

		    readSequenceLength(&xnx, ff);
			readSequenceLength(&xny, ff);

			L[N].nx = xnx;
			L[N].ny = xny;
			/*
			if (N!=0 && L[N].nx !=L[N-1].nx)
			   terror("uncoherent NX values --reference genome is not the same");
			*/
			N++;
			fclose(ff);
		}
		fgets(line,MAXLINE,f);
	}
	fclose(f);
	(*nF) = nFiles;
	return L;
}
