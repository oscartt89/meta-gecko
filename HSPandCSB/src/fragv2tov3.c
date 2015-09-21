/*
 * fragv2tov3.c
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	



Convierte un fichero de fragmentos v2 en v3. Hace falta k y lambda para calcular pvalue

We change fragments structure adding new fields. We also calc p-value 

Input:
	frags file.
	karlin parameters: k and lambda.
Out:
	frags files in a new structure.
	frags files in txt format.

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <dirent.h>
#include <math.h>
#include <inttypes.h>


#include "fragmentv3.h"
#include "fragmentv2.h"




int main(int ac,char** av){

	if(ac<5){
	printf("\n Uso:\n\tfragv2tov3 entrada.fragv2 karlinParameters salida.fragv3 salida.v3.txt \n");exit(-1);
	}

// Get Parameters

float lambda,k; // karlin parameters
/*************************************/
/* Open input file */	
	FILE* fin;
	if((fin=fopen(av[2],"r"))==NULL){
		printf("***ERROR Opening output file %s\n",av[2]);
		exit(-1);
	}
	char linea[1000];
	int max=1000;
	
	fgets(linea,max,fin);
	if(!sscanf(linea,"Lambda=%f  K=%f\n",&lambda,&k)==2){
		printf("****ERROR in getKarlinParameters\n");
	}
/************************************/		


// Read fragments
struct FragFile* f;
Fragmentv3* fv3;
int nf; // Number of fragments
uint64_t xtotal,ytotal;

nf=0;
f=readFragments(av[1],&nf,&xtotal,&ytotal);


fv3=(Fragmentv3*)malloc(sizeof(Fragmentv3)*nf);

int i;
long double pvalue,evalue;


int tik=0; // Turns '1' if there is any fragment statistically relevant.

for(i=0;i<nf;i++){
	// P-value calculation
	pvalue=-expm1l(-((long double)k*4*4*((long double)1+expm1l(-((long double)lambda*(long double)f[i].score))))   );
	evalue=-log(1-pvalue);
// Normal fields
	fv3[i].diag=(unsigned long)f[i].diag;
    fv3[i].xIni=(unsigned long)f[i].xStart;
    fv3[i].yIni=(unsigned long)f[i].yStart;
    fv3[i].xFin=(unsigned long)f[i].xEnd;
    fv3[i].yFin=(unsigned long)f[i].yEnd;
 
	fv3[i].length=(unsigned long)f[i].length;
    fv3[i].ident=(unsigned long)f[i].ident;
    fv3[i].score=(unsigned long)f[i].score;
    fv3[i].similarity=f[i].similarity;
    fv3[i].seqX=(unsigned long)f[i].seqX; //sequence number in the 'X' file
    fv3[i].seqY=(unsigned long)f[i].seqY; //sequence number in the 'Y' file
    fv3[i].block=f[i].block;          //synteny block id
    fv3[i].strand=(char)f[i].strand;        //'f' for the forward strain and 'r' for the reverse

	
// New fields	
	fv3[i].pvalue=pvalue;
// Modelling events
	fv3[i].events=0;
	fv3[i].id=0;
	fv3[i].gap=0;
	fv3[i].translocations=0;
	fv3[i].reversions=0;
	fv3[i].duplications=0;
	fv3[i].events=0;

	
	pvalue=1-exp(-k*exp(-lambda*f[i].score));

/*	// Calc if there is any relevant fragment.
	if(fv3[i].pvalue==0){
		tik=1;
			
	}
*/	
tik=1; // We want all fragments.


	// Le damos la vuelta y nos quitamos de complicaciones.
	// Los fragmentos reversos: yIni > yFin. Ojo con eso

	// Change Y component. Then in fragv2tov3 we revert this modification.

	if(fv3[i].strand=='r'){
		fv3[i].yIni=ytotal-fv3[i].yIni;
		fv3[i].yFin=ytotal-fv3[i].yFin;
	}

	
}


			
if(tik){
    
	FILE* fs;
	if((fs=fopen(av[3],"wb"))==NULL){
		printf("***ERROR Opening output file %s\n",av[4]);
		exit(-1);
	}
	FILE* ftxt;
	if((ftxt=fopen(av[4],"w"))==NULL){
		printf("***ERROR Opening output file %s\n",av[5]);
		exit(-1);
	}
	
	fprintf(ftxt,"xIni\tyIni\txFin\tyFin\tlength\tstrand\tident\n");
	
	//printf("abierto %s\n",av[4]);
	fwrite(&xtotal,sizeof(int),1,fs);
	fwrite(&ytotal,sizeof(int),1,fs);
	for(i=0;i<nf;i++){
		if(fv3[i].pvalue==0){
			fv3[i].block=i;// Deberia ser 0
		}else{
			fv3[i].block=-1; // If the fragment is no relevant, we mark it with -1
		}
			
		fwrite(&fv3[i],sizeof(Fragmentv3),1,fs);
		fprintf(ftxt,"%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%d\n",fv3[i].xIni,fv3[i].yIni,fv3[i].xFin,fv3[i].yFin,fv3[i].length,fv3[i].strand,fv3[i].ident,fv3[i].block);
			
	}
	fclose(fs);
	fclose(ftxt);

}


return 0;

}

