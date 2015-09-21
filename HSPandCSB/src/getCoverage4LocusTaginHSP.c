/*
getCoverage4LocusTaginHSP.c

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

if(ac<5){
	printf("Use: getCoverage4LocusTaginHSP input.frags seq1.genes seq2.genes salida \n");
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
	if(f[i].strand=='r' ){
		tmp= f[i].yStart;
		f[i].yStart = f[i].yEnd;
		f[i].yEnd = tmp;
	}

}


/**************/
// Open gene files
int ng1,ng2;
GeneFeatures *G1,*G2;

// Open csv output filelength
FILE* fs;


if((fs=fopen(av[4],"w"))==NULL){
	printf("***ERROR Opening output file");
	exit(-1);
}


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
//printf("%ld,%ld,%ld,%ld,%c,%ld,%ld,%f,%f,%ld,%d\n",f[i].xStart,f[i].xEnd,f[i].yStart,f[i].yEnd,f[i].strand,f[i].seqX,f[i].seqY,locusX,locusY,f[i].length,f[i].block);

	if(f[i].block>0){

		f[i].seqX=getLocusLong(G1, (int)f[i].xStart, (int)f[i].xEnd,ng1,f[i].strand);
		//printf("total: Ini:%ld, End:%ld ,%ld\n",f[i].xStart,f[i].xEnd,f[i].seqX);
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

	//	COVLocusX[i]=locusX;
	//	COVLocusY[i]=locusY;
	}
//printf("%ld,%ld,%ld,%ld,%c,%ld,%ld,%f,%f,%ld,%d\n",f[i].xStart,f[i].xEnd,f[i].yStart,f[i].yEnd,f[i].strand,f[i].seqX,f[i].seqY,locusX,locusY,f[i].length,f[i].block);
		
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
	//printf("covG1: %d, total G1: %d, %f, Cov2: %d, total G2 %d, %f\n",covG1,totalG1,100*(double)covG1/totalG1,covG2,totalG2,100*(double)covG2/totalG2);

	//writeFragmentsv3 (f,av[4],nf,xtotal, ytotal);
	
	/************** Aqui sacamos tripleta de fragmentos y genes q solapan ***********************/
	int xStart,xEnd,yStart,yEnd;
	int xOrigen,yOrigen;
	xStart=xEnd=yStart=yEnd=0;
	xOrigen=0;
	yOrigen=0;
	j=0;
	int k;
	int ant=0;
	int ant_i=0;
	int block=0;
	
	// Valores CSB
	float score,similarity,ident;
	int length=0;
	int open=0;
	score=ident=similarity=0.0;
	
	fprintf(fs,"Total frags      : %d\n", nf);
	fprintf(fs,"x========================================================\n");
	fprintf(fs,"Type,xStart,yStart,xEnd,yEnd,strand,block,length,score,ident,similarity,porc_ident,SeqX,SeqY,gene,locus_tag,product\n");
		
	for(i=0;i<nf;i++){
		// Solo los buenos
		block=(int)f[i].block;
		
		if(block>0){
			open=1;
			if(ant_i==0){
				xOrigen=(int)f[ant_i].xStart;
				yOrigen=(int)f[ant_i].yStart;
			}
			
			if( block>ant && i){// block>ant && i
			
				// We have a new CBS
				// We write all genes between xOrigen and xEnd (last fragment)
				// Genes.1
				
				// Pintamos CSB
				length=MAX(abs(xOrigen-f[ant_i].xEnd),abs(yOrigen-f[ant_i].yEnd))+1;
				
				fprintf(fs,"CSB,%d,%d,%d,%d,%c,%d,%d,%.2f,%.2f,%.2f,%.2f,%d,-,-\n",
					xOrigen,yOrigen,(int)f[ant_i].xEnd,(int)f[ant_i].yEnd,f[ant_i].strand,(int)f[ant_i].block,length,  100.0*score/ (length*4.0)  ,(double)ident,100.0*(ident)/(double)length,similarity,(int)score);

				// Now we change xOrigen and yOrigen because we have a new CSB
				
				ident = 0;
				similarity = 0;
				score = 0;
				
				xOrigen=(int)f[i].xStart;
				yOrigen=(int)f[i].yStart;
				ant=block;
				
				j++;
				open=0;
			}
			ident += (float)f[i].ident;
			similarity += (float)f[i].similarity;
			score += (float)f[i].score;
				
				
			xStart=(int)f[i].xStart;
			xEnd=(int)f[i].xEnd;
			yStart=(int)f[i].yStart;
			yEnd=(int)f[i].yEnd;
			if(f[i].strand=='r'){
				yEnd=(int)f[i].yStart;
				yStart=(int)f[i].yEnd;
			}
			ant_i=i;
			// Pintamos frag
			fprintf(fs,"Frag,%d,%d,%d,%d,%c,%d,%d,%.2f,%.2f,%.2f,%.2f,%d,-,-\n",
			xStart,yStart,xEnd,yEnd,f[i].strand,
			(int)f[i].block,(int)f[i].length,
			100.0*(int)f[i].score/ ((int)f[i].length*4.0)  ,
			(double)f[i].ident,
			100.0*((int)f[i].ident)/(double)f[i].length,
			f[i].similarity,
			(int)f[i].score);
			

		//	printf("FRAGMENTO\t%d\t%d\t%d\t%d\t%d\t%d\t%c\n",i,(int)f[i].block,(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand);
		
		}else{
			xStart=(int)f[i].xStart;
			xEnd=(int)f[i].xEnd;
			yStart=(int)f[i].yStart;
			yEnd=(int)f[i].yEnd;
			if(f[i].strand=='r'){
				yEnd=(int)f[i].yStart;
				yStart=(int)f[i].yEnd;
			}
			fprintf(fs,"Frag,%d,%d,%d,%d,%c,%d,%d,%.2f,%.2f,%.2f,%.2f,%d,-,-\n",
			xStart,yStart,xEnd,yEnd,f[i].strand,
			(int)f[i].block,(int)f[i].length,
			100.0*(int)f[i].score/ ((int)f[i].length*4.0)  ,
			(double)f[i].ident,
			100.0*((int)f[i].ident)/(double)f[i].length,
			f[i].similarity,
			(int)f[i].score);
		}
		if( f[i].block!=-1){
			// Pintamos anotaciones
			for(k=0;k<ng1;k++){
			//	printf("TODOS\t %d\t%d\t%d\t%d\n",(int)G1[k].start,0,(int)G1[k].end,0);
				if( solape(f[i],G1[k]) ){
				//	printf("ENTRA\t %d\t%d\t%d\t%d\n",(int)G1[k].start,0,(int)G1[k].end,0);
					fprintf(fs,"GX,%d,%d,%d,%d,%c,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s\n",(int)G1[k].start,0,(int)G1[k].end,0,G1[k].strand,(int)f[i].block,(int)G1[k].end-(int)G1[k].start,0,0,0,0,0,0,G1[k].gene,G1[k].locus,G1[k].product);
				}
				
			}
			// Genes.2
			for(k=0;k<ng2;k++){
			//	printf("TODOS\t %d\t%d\t%d\t%d\n",0,(int)G2[k].start,0,(int)G2[k].end);
				if( solape(f[i],G2[k]) ){
				//	printf("ENTRA\t %d\t%d\t%d\t%d\n",0,(int)G2[k].start,0,(int)G2[k].end);
					fprintf(fs,"GY,%d,%d,%d,%d,%c,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s\n",0,(int)G2[k].start,0,(int)G2[k].end,G2[k].strand,(int)f[i].block,(int)G2[k].end-(int)G2[k].start,0,0,0,0,0,0,G2[k].gene,G2[k].locus,G2[k].product);
			
				}
				
			}
		}
	}
	
	/**************** Last CSB **************/
	// The last one
	// The last one
	if(open){
		
		
		fprintf(fs,"CSB,%d,%d,%d,%d,%c,%d,%d,%.2f,%.2f,%.2f,%.2f,%d,-,-\n",
					xOrigen,yOrigen,(int)f[ant_i].xEnd,(int)f[ant_i].yEnd,f[ant_i].strand,(int)f[ant_i].block,length,  100.0*score/ (length*4.0)  ,(double)ident,100.0*(ident)/(double)length,similarity,(int)score);

		//printf("----------------------\n");
	}		
	// We have a new CBS
	// We write all genes between xOrigen and xEnd (last fragment)
	// Genes.1
	for(k=0;k<ng1;k++){
		if( solape(f[i],G1[k]) ){
		//if((int)G1[k].start<xEnd && (int)G1[k].end>xOrigen){
			fprintf(fs,"GX,%d,%d,%d,%d,%c,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s\n",(int)G1[k].start,0,(int)G1[k].end,0,G1[k].strand,(int)f[ant_i].block,(int)G1[k].end-(int)G1[k].start,0,0,0,0,0,0,G1[k].gene,G1[k].locus);
		}
		
	}
	// Genes.2
	for(k=0;k<ng1;k++){

		if( solape(f[i],G2[k]) ){
		//if((int)G2[k].start<yEnd && (int)G2[k].end>yOrigen){
			fprintf(fs,"GY,%d,%d,%d,%d,%c,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s\n",(int)G2[k].start,0,(int)G2[k].end,0,G2[k].strand,(int)f[ant_i].block,(int)G2[k].end-(int)G2[k].start,0,0,0,0,0,0,G2[k].gene,G2[k].locus);
		}
		
	}

fclose(fs);

// Cat del inf con este nuevo.
char linea[1000];
sprintf(linea,"cat %s.INF %s > %s.csv\n",av[1],av[4],av[4]);
system(linea);
exit(-1);
}

/***********************/
/*************/
int solape ( struct FragFile f, GeneFeatures g){ // Devuelve porcentaje de solapamiento

double d;
int plf,plg;


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
	

	if(d){
		
		plf = (int)(100*(double)d/(double)f.length);
		plg = (int)(100*(double)d/(double)g.length);
		
		d=MAX((long int)plf,(long int)plg);
			
		return (int)d ;
	}else{
		return 0;
	}
			
}

/*****/
int MIN(long int a, long int b){
	if (a>=b){return b;}else{return a;}}
int MAX(long int a, long int b){if (a>=b){return a;}else{ return b;}}

