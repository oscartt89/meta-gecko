/*

getInfoCSB.c

This program takes a fragment file with CSB marked and extract CSB composition

Use: getInfoCSB file.frags fragment_composition(0 no, 1 yes) \n

example:
./getInfo S3chr2R-S2chr2R-L200-S40-K8.csb.frags 1 

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"
#include "fragmentv2.h"

#include "lista.h"
#define  MAX_LEVELS  900000


int main(int ac,char** av){

	
	if(ac<5){
		printf("Use: getInfoCSB file.frags fragment_composition(0 no, 1 yes) outfile.txt SeqX SeqY \n");
		exit(-1);
	}
	
	// Read fragments
	struct FragFile* f;
	int nf; // Number of fragments
	uint64_t xtotal,ytotal;
	nf=0;
	f=readFragments(av[1],&nf,&xtotal,&ytotal);
	
	/******************************************/
	int i;
	int block,ant;
	
	
	int xIni,xFin,yIni,yFin;
	int length,score,nblocks;
	char strand;
	
	length=0;
	nblocks=0;
	ant=0;
	// Calculate number of blocks
	for(i=0;i<nf;i++){
		
		if(f[i].block>0){
			block=(int)f[i].block;
			if(ant!=block){
				nblocks++;
			}
			ant=block;
		}
	}
	
	// Print header. Frags file info
	printf("Total CSB      : %d\n", nblocks);
	printf("x========================================================\n");
	
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
	//for(i=0;i<nf;i++)printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,0,0,(int)f[i].seqX,(int)f[i].seqY);
	//exit(-1);			

	
	/** Borrar lo de antes**/
	length = 0;
	xIni=0;
	yIni=0;
	strand =0;
	xFin=0;
	yFin=0;
	score =  0;

	double similarity,likeness;
	likeness=0;
	int ident;
	int open=0;
	int i_ant=0;
	ident=0;
	
	// Name it.
	/*
	for(i=0;i<nf;i++){
		//if(f[i].block>0)f[i].block=i;
		printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,similarity,likeness,(int)f[i].seqX,(int)f[i].seqY);
	}	
*/	
			
	for(i=0;i<nf;i++){
		block=(int)f[i].block;
		
		if(block >=0){ // It's part of a CSB
			
			
			if(ant==block ){// Last valid fragment and current fragment share CSB
				similarity=(((double)f[i].score)/((double)f[i].length*4.0));
				likeness=(((double)f[i].ident)/((double)f[i].length));
				
				if(atoi(av[2])){// Print fragments
				printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,f[i].similarity,likeness,(int)f[i].seqX,(int)f[i].seqY);
				}
				// Update last fragment
				strand = f[i].strand;
				xFin=(int)f[i].xEnd;
				yFin=(int)f[i].yEnd;
				length = length+ (int)f[i].length;
				score = score+ (int)f[i].length * (int)f[i].score;
				ident = (int)f[i].ident;
				open=1;
				i_ant=i;
			}else{
				if(ant){
					// We finish CSB
					similarity=(((double)score)/(length*4.0));
					if(length)printf("CSB,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,,\n",xIni,yIni,xFin,yFin,strand,ant,length,(int)(score/(double)length),(int)ident,similarity,(int)ident/(double)length);
					if(atoi(av[2])){// Print fragments
						//printf("----------------------\n");
					}
					
					// New CSB
					length = (int)f[i].length;
					xIni=(int)f[i].xStart;
					yIni=(int)f[i].yStart;
					strand = f[i].strand;
					xFin=(int)f[i].xEnd;
					yFin=(int)f[i].yEnd;
					score =  (int)f[i].length * (int)f[i].score;
					ident += (int)f[i].ident;
					i_ant=i;// nuevo
					if(atoi(av[2])){// Print fragments
						printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,similarity,likeness,(int)f[i].seqX,(int)f[i].seqY);
					}
					open=1;
				}
			}
			
			ant=block;
			
			
		}else{
		
			similarity=(((double)f[i].score)/((double)f[i].length*4.0));
			likeness=(((double)f[i].ident)/((double)f[i].length));
			
			printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,similarity,likeness,(int)f[i].seqX,(int)f[i].seqY);
				
		}
	}
	
	// The last one
	if(open){
		strand = f[i_ant].strand;
		xFin=(int)f[i_ant].xEnd;
		yFin=(int)f[i_ant].yEnd;
		similarity=(((double)score)/(length*4.0));
		
		if(length)printf("CSB,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,,\n",xIni,yIni,xFin,yFin,strand,ant,length,(int)(score/(double)length),(int)ident,similarity,(int)ident/(double)length);
		//printf("----------------------\n");
	}			
			
				
	

	
	
	
	return 0;
}
