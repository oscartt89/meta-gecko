/*
 * duplicatedMetaFrags.c
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	

 Input:
 ---------
	frags file
	A percentage of overlapping. Integer [0-100]
	
Fragment with block<0 will be ignored.

Process:
---------
If exclusive: fragment.block=1;

If fragment1 overlapp more than #SOL# with fragment2:
	if fragment1.score > Fragment2.score then
		fragment1.block =2
		fragment2.block = -2
	else
		fragment2.block= 2
		fragment1.block=-2

output:
--------
	file.dup.frags
	file.dup.txt
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
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "fragmentv2.h"

/* Programs */
int MIN (long int a, long int b);
int solape (struct FragFile f, struct FragFile g,char c,int sol);
int MAX (long int a, long int b);

void criterio(struct FragFile* f, struct FragFile* g,int i,int j);

/*****************/

int main(int ac,char** av){

	FILE* fs;
	if((fs=fopen(av[3],"wb"))==NULL){
		printf("***ERROR Opening output file%d %s\n",3,av[3]);
		exit(-1);
	}
	FILE* ftxt;
	if((ftxt=fopen(av[4],"wt"))==NULL){
		printf("***ERROR Opening output file%d %s\n",4,av[4]);
		exit(-1);
	}
	
	if(ac<5){
	printf("\n Uso:\ncomparaFrags frags1 SOL-(overlap-percentage)[0-100] output.frags output.frags.txt\n");exit(-1);
	}
	// Input
	int sol;
	sol = atoi(av[2]);
	if( sol <0 || sol>100){
		printf("Second parameter (SOL parameter) must be a number between 0 and 100.\n");
		exit(-1);
	}
	
	

	
	
	
	// Read fragments from frags1 and frags2
	struct FragFile* f;
	
	int nf;
	uint64_t xtotal,ytotal;
	
	f=readFragments(av[1],&nf,&xtotal,&ytotal);

	
	
	
	int i,j;
	
	// Ponemos a cero el campo bloque
	for(i=0;i<nf;i++)if(f[i].block>=0)f[i].block=1;
	
	
	printf("xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\tscore\tseqX\tseqY\toverlapped%%\n");
	/******************/
	for(i=0;i<nf;i++){
		if(f[i].block>0){
			for(j=i+1;j<nf;j++){
				// Si solapan
				if( (f[j].block>0) && ((int)(solape(f[i],f[j],'x',sol))>sol) ){
					//printf("sol: %d\n",solape(f[i],f[j],'x'));
					criterio(f,f,i,j); 
					
				}
			}
		}
	}
	/******************/
	printf(">EXCLUSIVE FRAGMENTS\n");
	printf("xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\tscore\tseqX\tseqY\tblock\n");
	for(i=0;i<nf;i++){
		if(f[i].block>0){
			printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f[i].xStart,(long int)f[i].yStart,(long int)f[i].xEnd,(long int)f[i].yEnd,(long int)f[i].length,f[i].strand,(long int)f[i].ident,(long int)f[i].score,(long int)f[i].seqX,(long int)f[i].seqY,f[i].block);

		}
	}
	
	/********* SAVE FRAGMENTS *********/
	
	
	fprintf(ftxt,"xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\n");
	
	//printf("abierto %s\n",av[4]);
	writeSequenceLength(&xtotal,fs);
	writeSequenceLength(&ytotal,fs);
	for(i=0;i<nf;i++){
		
			
		writeFragment(&f[i], fs);
		fprintf(ftxt,"%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\n",(long int)f[i].xStart,(long int)f[i].yStart,(long int)f[i].xEnd,(long int)f[i].yEnd,(long int)f[i].length,f[i].strand,(long int)f[i].ident,(long int)f[i].block);
			
	}
	fclose(fs);
	fclose(ftxt);
	/********* END SAVE FRAGMENTS *********/
	return 0;	

}

/*
Programas
*/

void criterio (struct FragFile* f, struct FragFile* g,int i, int j){

	if( f[i].score > g[j].score){
		f[i].block=2;
		g[j].block=-2;
		
	}else{
		f[i].block=-2;
		g[j].block=2;
	
	}
}


/*************/
int solape ( struct FragFile f, struct FragFile g,char c,int sol){ // Devuelve porcentaje de solapamiento


double d=0.0;

int plf,plg;

long int gxS,fxS,gxE,fxE;
gxS=(long int)g.xStart;
gxE=(long int)g.xEnd;
fxS=(long int)f.xStart;
fxE=(long int)f.xEnd;



	if(c=='x'){
		d= MIN(gxE,fxE)-MAX(gxS,fxS);
		if(d<0)d=0;
		
	}else{
	/*
		if( g.yStart<=f.yEnd){
			if(g.yEnd<=f.yEnd){
				d=abs((long int)g.yEnd-(long int)f.yEnd);
			}else{
				d=abs((long int)f.yEnd-(long int)g.yStart);
			}
		}else{
			return 0;
		}
	*/
	}

	
	if(d){
		
		plf = (int)(100*(double)d/(double)f.length);
		plg = (int)(100*(double)d/(double)g.length);
		
		d=MAX((long int)plf,(long int)plg);
		if((int)d>sol){
			printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
			printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
			printf("-------\n");
		}
		
		return (int)d ;
	}else{
		return 0;
	}
	
	
	
	
}

/*****/
int MIN(long int a, long int b){if (a>=b)return b;else return a;}
int MAX(long int a, long int b){if (a>=b)return a;else return b;}
