/*
 * genomeAbundance.c
 *
 *  Created on: 17/11/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	Given a fragfile, reads-vs-genomas it takes the best genome in read, taking the best set of non-ooverlapping
 fragments. Then, it calculates the abundance of reads in genome.
 *  
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
#include "fragmentv2.h"


/* Programs */
int MIN (long int a, long int b);
int solape (struct FragFile f, struct FragFile g,char c,int sol);
int MAX (long int a, long int b);

int criterio(struct FragFile* f, struct FragFile* g,int i,int j);
void print(struct FragFile f);
int getMaxIdentGenomes(int* ident_genomes,int ngenomes);
void seleccion(int* ident_genomes,int ngenomes,int* genomes,int* genomes2,int* genomes3);

struct FragFile*  sortFragments(struct FragFile* fG,int nFG);
int simScore(struct FragFile* f,int nFG);
int sum(int* genomes,int ngenomes);

/*****************/
void xterror(char *s) {
	fprintf(stderr,"\nERR **** %s  ***\n",s);
	exit(-1);
}

// Global variables
/* glpbal*/
char OPT;
int Go,Ge;
int CONT=0;
int main(int ac,char** av){
	int U;
	int sol;
	int ngenomes;
	FILE* fs;
	FILE* ftxt;
	struct FragFile* f;
	int nf;
	uint64_t xtotal,ytotal;
	int i,j;
	int readX,genomeY;

	if(ac!=10)
	   	xterror("Usage:\ngenomeAbundance frags.FILE SOL-(overlap-percentage)[0-100] output.frags output.frags.txt UmbralCoverage numberOfGenomes Gopen[100] Gextension[4] OPT('i' ident, 'c' coverage, 's' score, 'g' gapped_score\n");
	
	// Input
	Go=atoi(av[7]);
	Ge=atoi(av[8]);
	sol = atoi(av[2]);
	OPT= av[9][0];
	if( sol <0 || sol>100)
		xterror("2nd param (SOL) must be [0-100]");
	U = atoi(av[5]);
	if( U <0 || U>100)
		xterror("5th param (U) must be [0-100]");
	
	if((fs=fopen(av[3],"wb"))==NULL)
		xterror("***Opening output file 1");
	if((ftxt=fopen(av[4],"wt"))==NULL)
		xterror("***Opening output file 2 (param 4)");
	
	ngenomes = atoi(av[6]);
	// Read MetaGenomefragments from frags1 and frags2
	
	f=readFragments(av[1],&nf,&xtotal,&ytotal);

	
	// set block to zero
	for(i=0;i<nf;i++){if(f[i].block>=0){f[i].block=1;}}
	

	/******************/
	//printf("xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\tscore\tseqX\tseqY\toverlapped\%\n");
	/******************/
	int rX_ant,gY_ant,inicio;
	int k;
	rX_ant=gY_ant=inicio=0;
	for(i=0;i<nf;i++){
		if((int)f[i].block>0){
			readX=(int)f[i].seqX;
			genomeY=(int)f[i].seqY;
			if((rX_ant != readX)||(gY_ant!= genomeY)){
				inicio = i;
				rX_ant=readX;
				gY_ant=genomeY;
			}
			j=inicio;
			while(j<nf && ((int)f[i].block>0)&&(readX==(int)f[j].seqX) && (genomeY==(int)f[j].seqY)){
			
				// If f[j] is not overlapped yet
				if( (f[j].block>0) ){
					// Check if it is overlapped
					if(((int)(solape(f[i],f[j],'t',sol))>0)&&(i!=j)){
						//If it is, we check which one has higher score
						if(criterio(f,f,i,j)){
							// This function returns 1 if f[j] has better score than f[i]. 
							//In this case, we have to re-check all fragment that has been overlapped by f[i]
							for(k=i;k<j;k++){
								if((int)f[k].block==-i){ // In "criterio" we put as a block field a minus X. And X it's the index of the fragment
								// who was overlapped by.
									if(((int)(solape(f[j],f[k],'t',sol))>0)){
										if(criterio(f,f,j,k)){// We know that in this case, f[j].score > f[k]
											//printf("ERROR\n");
											//exit(-1);
										} 
									}else{
										f[k].block=0; // If they are not overlapped, we restore the field block.
									}
								}
							}
							j=i; // We go back. Then with j++ we will be in the rigth place.
						} 
					}
				}
				j++;
			}
		}
	}
	/******************/
	
	
	/*
	printf(">EXCLUSIVE FRAGMENTS\n");
	printf("xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\tscore\tseqX\tseqY\tblock\%\n");
	for(i=0;i<nf;i++){
		if(f[i].block>0){
			printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f[i].xStart,(long int)f[i].yStart,(long int)f[i].xEnd,(long int)f[i].yEnd,(long int)f[i].length,(long int)f[i].strand,(long int)f[i].ident,(long int)f[i].score,(long int)f[i].seqX,(long int)f[i].seqY,f[i].block);

		}
	}
	*/
	
	/**********************/
	// Calc best collection of fragments in one Genome 
	// here we have all non-overlapping reads (fragments???)
	
	//
	struct FragFile* fG; // Fragmentos mapeados en genoma
	
	int* ident_genomes;
	int* score_genomes;
	int* cov_genomes;
	int* sim_score_genomes;
	int* genomes;
	int* genomes2;
	int* genomes3;
	genomes=(int*)malloc(sizeof(int)*ngenomes);
	genomes2=(int*)malloc(sizeof(int)*ngenomes);
	genomes3=(int*)malloc(sizeof(int)*ngenomes);
	ident_genomes=(int*)malloc(sizeof(int)*ngenomes);
	cov_genomes=(int*)malloc(sizeof(int)*ngenomes);
	score_genomes=(int*)malloc(sizeof(int)*ngenomes);
	sim_score_genomes=(int*)malloc(sizeof(int)*ngenomes);
	
	
	
	for(i=0;i<ngenomes;i++)genomes[i]=genomes2[i]=genomes3[i]=ident_genomes[i]=0;
	
	double suma1,suma2,suma3;
	int readX_ant,genomeY_ant;
	int maxRead=0;
	int l;
	float cov_read;
	int indice;
	int maxX;
	float Totcov;
	int fin;
	inicio=fin=0;
	maxX=0;
	cov_read=0.0;
	readX_ant=genomeY_ant=0;
	int score=0;
        // OTS vars--------
        int nFG; // number of non-overlapping frags 
        int totIDENT; // total number of identities in the matching
        int totLEN;   // total length of the fragments 
		int totCOV;
		int totSCORE;
       // printf("R\tG\tCover\tLenR\ttotSco\ttotID\ttotLEN\t%%ID\tnF\tfrags\n");
	for(i=0;i<nf;i++){
	maxX=MAX(maxX,f[i].xEnd);
		if(f[i].block>0){
			
			readX=(int)f[i].seqX;
			genomeY=(int)f[i].seqY;
			maxRead=MAX(maxRead,readX);
		// Si no esta solapado
		
                        totIDENT=0;
                        totLEN  =0;
						totCOV=0;
						totSCORE=0;
			if( (readX_ant!=readX) || (genomeY_ant != genomeY) ){
				// end the previous
				
				Totcov = cov_read/(double)maxX;
				
				if(abs(cov_read-maxX)==1)Totcov=1.0;
				fin=i;
				/*** OUT  fragments ***/
				
				nFG=0;
				for(i=inicio;i<fin;i++)
					if(f[i].block>0) {
						nFG++;
						totIDENT+=f[i].ident;
						totLEN  +=f[i].length;
						totCOV += (int)100*Totcov;
						totSCORE += f[i].score;
					}	
				//debug	 printf("ident:%d\t%d\t%5.2f\t%d\n",totIDENT,totLEN,100.*(float)totIDENT/(float)totLEN, nFG);
				// create an array with fragments
				fG = (struct FragFile*)malloc(sizeof(struct FragFile)*nFG);
				nFG=0;
				for(l=inicio;l<fin;l++){
					if(f[l].block>0){
						memcpy(&fG[nFG++],&f[l],sizeof(struct FragFile));
						//print(fG[nFG-1]);
					}
				}
				fG=sortFragments(fG,nFG);
				//debug for(l=0;l<nFG;l++)print(fG[l]);
				//debug printf("*****\n");
				
				
				//if(CONT>4)exit(-1);
				/** Calculo de valores **/
				/*********/
				ident_genomes[genomeY_ant]=totIDENT;
				cov_genomes[genomeY_ant]=totCOV;
				score_genomes[genomeY_ant]=totSCORE;
				sim_score_genomes[genomeY_ant]=simScore(fG,nFG);
				/*********/
				// Borramos fragmentos
				free(fG);
				fG=NULL;
					
				
				// new round, new read
				if(readX_ant != readX){
				CONT++;
				
					switch (OPT){
					
						case 'i': seleccion(ident_genomes,ngenomes,genomes,genomes2,genomes3); break;
						case 'c': seleccion(cov_genomes,ngenomes,genomes,genomes2,genomes3); break;
						case 's': seleccion(score_genomes,ngenomes,genomes,genomes2,genomes3); break;
						case 'g': seleccion(sim_score_genomes,ngenomes,genomes,genomes2,genomes3); break;
						default: seleccion(ident_genomes,ngenomes,genomes,genomes2,genomes3);
					}
					
					
					
					
					
					
					maxX=f[i].xEnd;
				
				}
				
				inicio=i;
				readX_ant= readX;
				genomeY_ant =genomeY;
				
				score= (int)f[i].score;
				
				cov_read = (int)f[i].length;
				//printf("**************\n");
				//print(f[i]);
			}else{
				//continue 
				
				score += (int)f[i].score;
				cov_read += (int)f[i].length;
				//print(f[i]);
			}
		}
	
	}
	
	/* print genomes and abundance */
	printf("Genome\tAbundance %s\n",av[1]);
	printf("Genome\topc1\t%%opc1\topc2\t%%opc2\topc3\t%%opc3\n");
	suma1=(double)sum(genomes,ngenomes);
	suma2=(double)sum(genomes2,ngenomes);
	suma3=(double)sum(genomes3,ngenomes);
	for(i=0;i<ngenomes;i++)printf("%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n",i,genomes[i],100.0*genomes[i]/suma1,genomes2[i],100.0*genomes2[i]/suma2,genomes3[i],100.0*genomes3[i]/suma3);
	printf("reads asignados : %d\t reads totales: %d\n",sum(genomes,ngenomes),maxRead);
		
	
	/********* SAVE FRAGMENTS *********/
	
	
	fprintf(ftxt,"xStart\tyStart\txEnd\tyEnd\tlength\tstrand\tident\n");
	
	//printf("abierto %s\n",av[4]);
	writeSequenceLength(&xtotal,fs);
	writeSequenceLength(&ytotal,fs);
	for(i=0;i<nf;i++){
		
			
		writeFragment(&f[i], fs);
		fprintf(ftxt,"%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\n",(long int)f[i].xStart,(long int)f[i].yStart,(long int)f[i].xEnd,(long int)f[i].yEnd,f[i].length,f[i].strand,(long int)f[i].ident,(long int)f[i].block);
			
	}
	fclose(fs);
	fclose(ftxt);
	/********* END SAVE FRAGMENTS *********/
	return 0;
}

/*
Programas... functions???
*/

void print(struct FragFile f){
	printf("%d\t%d\t%d\t%d\t%c\t%d\t%.2f\t%d\t%d\t%d\t%d\n",(int)f.xStart,(int)f.yStart,(int)f.xEnd,(int)f.yEnd,f.strand,(int)f.length,f.similarity,(int)f.ident,(int)f.seqX,(int)f.seqY,(int)f.block);

}
int criterio (struct FragFile* f, struct FragFile* g,int i, int j){

	if( (int)f[i].score >= (int)g[j].score){
		g[j].block=-1*i;
		f[i].block=2;
		return 0;
		
		
		
	}else{
		f[i].block=-1*j;
		g[j].block=2;
		return 1; // El candidato es peor que el nuevo
	}
	
}


/*************/
// returns the percentage of overlapping
// Devuelve porcentaje de solapamiento
int solape ( struct FragFile f, struct FragFile g,char c,int sol){ 


double d=0.0;
d=sol;
d=0.0;

int plf,plg;

long int gxS,fxS,gxE,fxE;
gxS=(long int)g.xStart;
gxE=(long int)g.xEnd;
fxS=(long int)f.xStart;
fxE=(long int)f.xEnd;

long int gyS,fyS,gyE,fyE;
gyS=MIN(g.yStart,g.yEnd);
gyE=MAX(g.yStart,g.yEnd);
fyS=MIN(f.yStart,f.yEnd);
fyE=MAX(f.yStart,f.yEnd);



	if(c=='x'){
		d= MIN(gxE,fxE)-MAX(gxS,fxS);
		if(d<0)d=0;
		
	}else if(c=='y'){
		d= MIN(gyE,fyE)-MAX(gyS,fyS);
		if(d<0)d=0;
	}else{
		MAX(MIN(gxE,fxE)-MAX(gxS,fxS),MIN(gyE,fyE)-MAX(gyS,fyS));
		if(d<0)d=0;
	}

	
	
	if(d){
		
		plf = (int)(100*(double)d/(double)f.length);
		plg = (int)(100*(double)d/(double)g.length);
		
		d=MAX((long int)plf,(long int)plg);
//		printf("-- SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");		
		return (int)d ;
	}else{
	
//		printf("-- NO SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");
	
		return 0;
	}
	
	
	
	
}
/******/
int getMaxIdentGenomes(int* ident_genomes,int ngenomes){

	int value=0;
	int indice=-1;
	int i;
	i=0;
	for(i=0;i<ngenomes;i++){
		
		if(value<MAX(value,ident_genomes[i])){
			value=MAX(value,ident_genomes[i]);
			indice=i;
		}
			
	}
	ident_genomes[indice]=0;
	return indice;
}
/*****/
int MIN(long int a, long int b){if (a>=b)return b;else return a;}
int MAX(long int a, long int b){if (a>=b)return a;else return b;}
/*****************/

void seleccion(int* ident_genomes,int ngenomes,int* genomes,int* genomes2,int* genomes3){

int j;
int indice;
	//printf("%d\t Nuevo read %d\n",CONT++,readX);
	indice=getMaxIdentGenomes(ident_genomes,ngenomes);
	if(indice>=0){// Si el indice es >=0 es porque hay algun ganador
		genomes[indice]++;
	//	printf("Ganador: %d\n",indice);
	}
	
	indice=getMaxIdentGenomes(ident_genomes,ngenomes);
	if(indice>=0){// Si el indice es >=0 es porque hay algun ganador
		genomes2[indice]++;
	//	printf("Ganador: %d\n",indice);
	}
	
	indice=getMaxIdentGenomes(ident_genomes,ngenomes);
	if(indice>=0){// Si el indice es >=0 es porque hay algun ganador
		genomes3[indice]++;
	//	printf("Ganador: %d\n",indice);
	}
	// Set 0 ident_genomes
	for(j=0;j<ngenomes;j++)ident_genomes[j]=0;
}

/**************/
struct FragFile* sortFragments(struct FragFile* fG,int nFG){

int i,j;
int ymax=pow(2,8*sizeof(int));
int ymin;
int yvalue;
int indice;

struct FragFile* f;

f = (struct FragFile*)malloc(sizeof(struct FragFile)*nFG);
	
	for(i=0;i<nFG;i++){
		ymin=MIN(fG[i].yStart,fG[i].yEnd);
		for(j=0;j<nFG;j++){
			yvalue=MIN(fG[j].yStart,fG[j].yEnd);
			if(yvalue<ymin){
				indice=j;
			}
		}
		memcpy(&f[i],&fG[indice],sizeof(struct FragFile));
		fG[indice].yStart=fG[indice].yEnd=ymax;
	
	}
	
free(fG);
fG=NULL;
return f;
}

/************/
int simScore(struct FragFile* f,int nFG){

int i;
int score_acom=0;
int ymax;
int ymin;
int dif;

	ymax=MAX(f[0].yStart,f[0].yEnd);
	for(i=1;i<nFG;i++){
		ymin=MIN(f[i].yStart,f[i].yEnd);
		dif= ymin-ymax;
		score_acom += MAX(f[i-1].score-Go,f[i-1].score-Ge*dif);
		
	}
	score_acom += f[i-1].score;
	
	return score_acom;

}
/******/
int sum(int* genomes,int ngenomes){
int i;
int suma=0;

for(i=0;i<ngenomes;i++){
	suma += genomes[i];

}
return suma;

}
