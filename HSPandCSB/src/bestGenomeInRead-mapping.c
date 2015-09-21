/*
 * bestGenomeInRead.c
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	
 *  history -----------------31 Oct OTS
 *  use xterror function
 *  translate to english
 *  include a column with the number of non-overlapping frags in the OUT
 *  include totNumberOfIdentities
 *  include totnonoverlapping Leng
 *  computes % of identidad
 *  include a header in the OUT file
 *  OUT---> Read: Genome: Coverage: LengthRead: TotalScore: 
 *               totIDENT: totLEN :%Ident :Nu.Frag : frags:
 *          line 213 -- -strand printted as (long int) ;;; chaged to char
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
/*****************/
void xterror(char *s) {
	fprintf(stderr,"\nERR **** %s  ***\n",s);
	exit(-1);
}

int main(int ac,char** av){
	int U;
	int sol;
	FILE* fs;
	FILE* ftxt;
	struct FragFile* f;
	int nf;
	uint64_t xtotal,ytotal;
	int i,j;
	int readX,genomeY;

	if(ac<6)
	   	xterror("Usage:\nbestGenomeInRead frags.FILE SOL-(overlap-percentage)[0-100] output.frags output.frags.txt UmbralCoverage Header(0=NO, 1=YES; OPTIONAL)\n");
	
	// Input
	int HEADER;
	if(ac==7){HEADER=atoi(av[6]);}else{HEADER=0;}
	sol = atoi(av[2]);
	if( sol <0 || sol>100)
		xterror("2nd param (SOL) must be [0-100]");
	U = atoi(av[5]);
	if( U <0 || U>100)
		xterror("5th param (U) must be [0-100]");
	
	if((fs=fopen(av[3],"wb"))==NULL)
		xterror("***Opening output file 1");
	if((ftxt=fopen(av[4],"wt"))==NULL)
		xterror("***Opening output file 2 (param 4)");
	
	
	// Read MetaGenomefragments from frags1 and frags2
	
	f=readFragments(av[1],&nf,&xtotal,&ytotal);

	
	// set block to zero
	for(i=0;i<nf;i++){if(f[i].block>=0){f[i].block=1;}}
	
	/*
	for(i=1;i<nf;i++){
		printf("%d\t%d\n",(int)f[i].seqX,(int)f[i].seqY);
		if(f[i].seqX<f[i-1].seqX){
			printf("Error en %d,%d",i,i+1);
			exit(-1);
		}else if(f[i].seqX==f[i-1].seqX){
			if(f[i].seqY<f[i-1].seqY){
				printf("Error en %d,%d",i,i+1);
				exit(-1);
			}
			
		}
		
	}
	exit(-1);
	*/
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
					if(((int)(solape(f[i],f[j],'x',sol))>0)&&(i!=j)){
						//If it is, we check which one has higher score
						if(criterio(f,f,i,j)){
							// This function returns 1 if f[j] has better score than f[i]. 
							//In this case, we have to re-check all fragment that has been overlapped by f[i]
							for(k=i;k<j;k++){
								if((int)f[k].block==-i){ // In "criterio" we put as a block field a minus X. And X it's the index of the fragment
								// who was overlapped by.
									if(((int)(solape(f[j],f[k],'x',sol))>0)){
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
	
	int readX_ant,genomeY_ant;
	float cov_read;
	int maxX;
	float Totcov;
	int fin;
	int contador=0;
	inicio=fin=0;
	maxX=0;
	cov_read=0.0;
	readX_ant=genomeY_ant=0;
	int score=0;
        // OTS vars--------
        int nFG=0; // number of non-overlapping frags 
        int totIDENT; // total number of identities in the matching
        int totLEN;   // total length of the fragments 
       if(HEADER) printf("R\tGenomes\n");
	  
	for(i=0;i<nf;i++){
	maxX=MAX(maxX,f[i].xEnd);
		if(f[i].block>0){
			readX=(int)f[i].seqX;
			genomeY=(int)f[i].seqY;
		// Si no esta solapado
		
                        totIDENT=0;
                        totLEN  =0;
			if( (readX_ant!=readX) || (genomeY_ant != genomeY) ){
				// end the previous
				Totcov = cov_read/(double)maxX;
				if( (Totcov) ){ // *100> U
					if(abs(cov_read-maxX)==1)Totcov=1.0;
					fin=i;
					/*** OUT  fragments ***/
					// Read: Genome: Coverage: LengthRead: TotalScore: totIDENT totLEN %Ident Nu.Frag  frags:
					//printf("R:%d\tG: %d\tCov: %.2f\t LRead:%d\tScore: %d\tfrags: ",readX_ant,genomeY_ant,Totcov,maxX,score);
					if(contador){
						contador=0;
						printf("%d\t",readX_ant);
					}
					printf("%d\t",genomeY_ant);
                                        // how many frags to be reported
                                        nFG++;
					for(i=inicio;i<fin;i++)
						if(f[i].block>0) {
							//nFG++;
							totIDENT+=f[i].ident;
							totLEN  +=f[i].length;
						}	
					 //printf("%d\t%d\t%5.2f\t%d\t",totIDENT,totLEN,100.*(float)totIDENT/(float)totLEN, nFG);
					/*
					for(i=inicio;i<fin;i++){
						if(f[i].block>0)
							 //printf("%d\t",i);
					}
					 */
					/*********/
					
				}
				// new round
				if(readX_ant != readX){
					maxX=f[i].xEnd;
					contador=1;
					printf("NGenomes: %d\n",nFG);
					nFG=0;
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

/*****/
int MIN(long int a, long int b){if (a>=b)return b;else return a;}
int MAX(long int a, long int b){if (a>=b)return a;else return b;}
