/*
getCoverage4HSP.c

For each HSP it calculates its coverage for x and y sequence (using genes.c)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>

#include "gene.h"
#include "fragmentv2.h"

int solape ( struct FragFile f, GeneFeatures g);
int MIN(long int a, long int b);
int MAX(long int a, long int b);

int main(int ac,char** av){

if(ac<4){
	printf("Use: getCoverage4LocusTaginHSP input.frags seq1.genes seq2.genes \n");
	exit(-1);
}


// Read fragments
struct FragFile* f;
int nf; // Number of fragments
uint64_t xtotal,ytotal;

nf=0;
f=readFragments(av[1],&nf,&xtotal,&ytotal);





/******/
int i;
int tmp;


// change components ( reverse fragments currently has yStart<yEnd
for(i=0;i<nf;i++){
	if(f[i].strand=='r'){
		tmp= f[i].yStart;
		f[i].yStart = f[i].yEnd;
		f[i].yEnd = tmp;
	}

}


/**************/
// Open gene files
int ng1,ng2;
GeneFeatures *G1,*G2;




G1= readGenes(av[2],&ng1);
G2= readGenes(av[3],&ng2);

//printf("nf1:%d\tnf2:%d\n",ng1,ng2);

// Create secundary vector.
float* COVLocusX;
float* COVLocusY;

COVLocusX = (float*)malloc(sizeof(float)*ng1);
COVLocusY = (float*)malloc(sizeof(float)*ng2);
// Set COVlocus
for(i=0;i<ng1;i++){COVLocusX[i]=0.0;}
for(i=0;i<ng2;i++){COVLocusY[i]=0.0;}


//nf=nf-1;

//printf("xStart\txEnd\tyStart\tyEnd\tstrand\tSeqX\tSeqY\tLocusX\tLocusY\tLength\n");

float locusX,locusY;

for(i=1;i<nf;i++){
//printf("%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%f\t%f\t%ld\t%d\n",f[i].xStart,f[i].xEnd,f[i].yStart,f[i].yEnd,f[i].strand,f[i].seqX,f[i].seqY,locusX,locusY,f[i].length,f[i].block);

	if(f[i].block>0){

		f[i].seqX=getLocusLong(G1, (int)f[i].xStart, (int)f[i].xEnd,ng1,f[i].strand);
		//printf("total: Ini:%ld\t End:%ld \t%ld\n",f[i].xStart,f[i].xEnd,f[i].seqX);
		f[i].seqY=getLocusLong(G2, (int)f[i].yStart, (int)f[i].yEnd,ng2,f[i].strand);
		//printf("total: %ld\n",f[i].seqY);

	//getchar();


		if(f[i].length!=0){
			locusX=(float)f[i].seqX/f[i].length;
			locusY=(float)f[i].seqY/f[i].length;
		}else{
			locusX=(float)f[i].seqX/(f[i].xEnd-f[i].xStart);
			locusY=(float)f[i].seqY/(f[i].yEnd-f[i].yStart);
		}

		
		if(locusX>1){
			locusX=1;
		}
		if(locusY>1){
			locusY=1;
		}

		COVLocusX[i]=locusX;
		COVLocusY[i]=locusY;
	}
//printf("%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%f\t%f\t%ld\t%d\n",f[i].xStart,f[i].xEnd,f[i].yStart,f[i].yEnd,f[i].strand,f[i].seqX,f[i].seqY,locusX,locusY,f[i].length,f[i].block);
		
}

	/***********************************/



	/**********************/
	int j;
	int covG1,covG2;
	covG1=covG2=0;
	int totalG1,totalG2;
	int covTemp=0;
	totalG1=totalG2=0;
	float overlapG1,overlapG2;
	overlapG1=0.0;
	overlapG2=0.0;
	
	//fprintf(fs,"Type,id,locus_tag,gene,strand,start,end,length,cov\n");
		
		
	for(i=0;i<ng1;i++){
	totalG1+=G1[i].length;// length of cds
	covTemp=0;
		for(j=0;j<G1[i].length;j++){
			covTemp=+*G1[i].cov;
			covG1=covG1+covTemp; 
		}
		overlapG1 = (100*covTemp/(double)G1[i].length);
		COVLocusX[i]=overlapG1;
		// print (all) or (covTemp) locus
		//if(1)fprintf(fs,"GX,%s,%s,%s,%c,%d,%d,%d,%f\n",G1[i].id,G1[i].locus,G1[i].gene,G1[i].strand,G1[i].start,G1[i].end,G1[i].length,overlapG1);
		
	}
	for(i=0;i<ng2;i++){
	totalG2+=G2[i].length;
	covTemp=0;
		for(j=0;j<G2[i].length;j++){
			covTemp=covTemp+*G2[i].cov;
			covG2=covG2+covTemp;
			
		}
		overlapG2 = (100*covTemp/(double)G2[i].length);
		COVLocusY[i]=overlapG2;
		// print (all) or (covTemp) locus
		//if(1)fprintf(fs,"GY,%s,%s,%s,%c,%d,%d,%d,%f\n",G2[i].id,G2[i].locus,G2[i].gene,G2[i].strand,G2[i].start,G2[i].end,G2[i].length,overlapG2);
		
	}
	

	
	
	/*********************************/
	//printf("covG1: %d\t total G1: %d\t %f\t Cov2: %d\t total G2 %d\t %f\n",covG1,totalG1,100*(double)covG1/totalG1,covG2,totalG2,100*(double)covG2/totalG2);

	//writeFragmentsv3 (f,av[4],nf,xtotal, ytotal);
	
	
	/************** Aqui sacamos tripleta de fragmentos y su solape ***********************/
	int xStart,xEnd,yStart,yEnd;
	int xOrigen,yOrigen;
	xStart=xEnd=yStart=yEnd=0;
	xOrigen=0;
	yOrigen=0;
	j=0;
	int k;
	int solapeX,solapeY,maxSol;
	solapeX=solapeY=maxSol=0;
	

	//printf("length\tsimilarity\tcov\n");
		
	for(i=0;i<nf;i++){
		// Genes.1
		for(k=0;k<ng1;k++)solapeX += solape(f[i],G1[k]);
		// Genes.2
		for(k=0;k<ng2;k++)solapeY += solape(f[i],G2[k]);
		maxSol = MIN (100*((double)MAX(solapeX,solapeY)/f[i].length),100);
		
		printf("%d\t%d\t%d\n",(int)f[i].length,(int)f[i].similarity,maxSol);
		
		solapeX=solapeY=0;
	}
	
	

exit(-1);
}

/***********************/
/*************/
int solape ( struct FragFile f, GeneFeatures g){ // Devuelve porcentaje de solapamiento

double d;



	long int gxS,fxS,gxE,fxE;
	gxS=(long int)g.start;
	gxE=(long int)g.end;
	fxS=(long int)f.xStart;
	fxE=(long int)f.xEnd;
	if(f.strand=='r'){
		fxS=MIN((long int)f.yStart,(long int)f.yEnd);
		fxE=MAX((long int)f.yEnd,(long int)f.yStart);
	}


long int glength;

glength=abs(g.start-g.end);

	d= MIN(gxE,fxE)-MAX(gxS,fxS);
	if(d<0)d=0;
	

	return d;
			
}

/*****/
int MIN(long int a, long int b){
	if (a>=b){return b;}else{return a;}}
int MAX(long int a, long int b){if (a>=b){return a;}else{ return b;}}

