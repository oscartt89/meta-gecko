/*
 * gene.c
 *
 *  Created on: 21/02/14
 *      Author: jarjonamedina
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "gene.h"



GeneFeatures* readGenes(char* s,int* nf){

	FILE* fe;

	GeneFeatures* fs;
	int n;
	int j;

	if((fe=fopen(s,"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}


	n=0;
	

	long int longFile;
	fseek(fe,0,SEEK_END);
	longFile=ftell(fe);
	n=(int)longFile/sizeof(GeneFeatures);
	//printf("\n ReadFragments Complete\nnum: %d\n",n);
	rewind(fe);

	fs=(GeneFeatures*)malloc(sizeof(GeneFeatures)*n);
	if(fs==NULL){
		printf("****ERROR: Out of memory\n");
		exit(-1);
	}
	//*nf=n-1;
	*nf=n;
	n=0;
	int nlen;
	nlen=0;
	while(!feof(fe)){
		fread(&fs[n],sizeof(GeneFeatures),1,fe);
	
		fs[n].length=fs[n].end-fs[n].start+1;
		fs[n].cov=(int*)malloc(sizeof(int)*fs[n].length);	
		if(fs[n].cov==NULL){
			printf("**** Out of memory in cov gene vector[%d]. Length %d-%d. See gene.c line 72\n",n,fs[n].start,fs[n].length);exit(-1);		
		}

	
		for(j=0;j<fs[n].length;j++){
			fs[n].cov[j]=0;
			nlen++;		
		}
		
		//printf("%s\t%s\t%c\t%d\t%d\t%d\n",fs[n].gene,fs[n].locus,fs[n].strand,fs[n].start,fs[n].end,fs[n].length);
		n++;
	}
	//printf("nlen: %d\n",nlen);


	fclose(fe);
	return fs;
}



/*********************************************
**************************************************
*********************************************/

void setCOV(GeneFeatures* G, int ini,int end,int n){

	int i;

	for(i=ini;i<=end;i++){
		// el indice va de 0 a length. no de ini a fin
		
		G[n].cov[i-G[n].start]=1;
		
		
	}
}

/*********************************************
**************************************************
*********************************************/


int getLocusLong(GeneFeatures* g, int ini, int end,int max,char c){


int longtotal=0;
int i=0;




	//printf("comienzo while: %d > %d  max: %d\n",ini,g[i].end,max);
	while (ini > (int)g[i].end && i<max )i++; // Find start locus
	
		//printf("g[%d].start: %d \t g[%d].end: %d  \n",i,g[i].start,i,g[i].end);
		//printf("ini: %d \t end: %d  \n",ini,end);

		// Check if there is something in the last locus
		//if(ini < g[i-1].end) longtotal+=g[i-1].end-ini;


		/*

			   <=============>
		 <-->

		*/
		if(ini<g[i].start && end < g[i].start){ 
			if(sameStrand(g[i].strand,c)){
				
				longtotal+=0;
				return longtotal;
			}
		i++;
		}

		/*

			  <=============>
		   <----->

		*/
		if(ini<g[i].start && end < g[i].end && end >= g[i].start){ 
			if(sameStrand(g[i].strand,c)){
				longtotal+=end-g[i].start;
				setCOV(g,g[i].start,end,i);
				
				return longtotal;
			}
		i++;
		}
		/*

			  <=============>
				  <----->

		*/
		if(ini>=g[i].start && end <= g[i].end){ 
			if(sameStrand(g[i].strand,c)){
				longtotal+=end-ini;
				setCOV(g,ini,end,i);
		//		printf("longtotal: %d\n",longtotal);
				return longtotal;
			}
		i++;
		}
		/*

			  <=============>
		   <------------------>

		*/
		if(ini<=g[i].start && end >= g[i].end){ 
			if(sameStrand(g[i].strand,c)){
				longtotal+=g[i].end-g[i].start;
				setCOV(g,g[i].start,g[i].end,i);
			}
		i++;
		}
		/*

			  <=============>
						 <-----

		*/
		if(ini<=g[i].end && ini>g[i].start && end > g[i].end){ 
			if(sameStrand(g[i].strand,c)){
				longtotal+=g[i].end-ini;
				setCOV(g,ini,g[i].end,i);
			}
		i++;
		}

		/*

				<=============>   <=============>
			------------------------------>

		*/

		while(end >= g[i].end && i<max){ 
			if(sameStrand(g[i].strand,c)){
				longtotal+=g[i].end - g[i].start;
				setCOV(g,g[i].start,g[i].end,i);
			}
			i++;
		}

		/*

				  <=============>
			----------->

		*/
		// Check the last
		if(end > g[i].start && i<max) {
			if(sameStrand(g[i].strand,c)){
				longtotal+=end-g[i].start;//g[i].strand==c
				setCOV(g,g[i].start,end,i);
			}
		}
	//printf("Devuelve valor: %d / %d\n",longtotal,end-ini);
	return longtotal;

}

int sameStrand(char g,char s){

if(g=='n')g='f';

if(g==s){

	return 1;
}else{
	return 0;
}
}
