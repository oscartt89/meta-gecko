/*
 * parseGeneBankTab.c
 *
 *  Created on: 20/11/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	
		
Given a file with format:

Location         Strand  Length  PID     Gene    Synonym Code    COG     Product
1	1401	+	466	148377269	dnaA	MAG_0010	-	-	chromosomal replication initiation protein

it stores into a binary file with this struct:

GeneFeatures
	id
	gene
	locus_tag
	strand
	Start
	End
	


 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <time.h>

#include "gene.h"

/*
typedef struct {
	char id[50];
	char locus[50];
	char gene[50];
	char strand;
	int start;
	int end;
	// NewField
	int length;	
	int* cov;
	char product[1000];
}GeneFeatures;
*/



/*****/
int main(int ac,char** av){

if(ac<3){
	printf("Use: parseGeneFeatureTab input.gbk.features output.gene output.csv \n");
	exit(-1);
}



FILE* fe,*fs,*fcsv;



// Open features
if((fe=fopen(av[1],"r"))==NULL){
		printf("*****ERROR opening %s file",av[1]);
		exit(-1);
	}
// Open new seq

if((fs=fopen(av[2],"wb"))==NULL){
		printf("*****ERROR opening %s file",av[2]);
		exit(-1);
	}
	
// Open new seq

if((fcsv=fopen(av[3],"w"))==NULL){
		printf("*****ERROR opening %s file",av[3]);
		exit(-1);
	}

/******************************/


GeneFeatures gene;

char linea[1000];
int max=1000;

int number;
char basurilla[50];
	
//CM000363.1_gene_1	Dsim\GD13543	Dsim_GD13543	2226..3544	
fprintf(fcsv,"id,gene,locus,strand,start,end\n");
	
fgets(linea,max,fe);
	while(!feof(fe)){
		
		// For genes
		// Normal line
		//1	1401	f	466	148377269	dnaA	MAG_0010	-	-	chromosomal replication initiation protein
		if(sscanf(linea,"%d\t%d\t%c\t%d\t%d\t%s\t%s\t%s\t%s\t%[^\n]\n",
		&gene.start,&gene.end,&gene.strand,&number,&number,gene.gene,gene.locus,basurilla,basurilla,gene.product)==10){
			
			
			fwrite(&gene,sizeof(GeneFeatures),1,fs);
			fprintf(fcsv,"%d,%d,%c,%s,%s,%s\n",gene.start,gene.end,gene.strand,gene.gene,gene.locus,gene.product);
			// debug printf("%d,%d,%c,%s,%s,%s\n",gene.start,gene.end,gene.strand,gene.gene,gene.locus,gene.product);
		}
		
		fgets(linea,max,fe);
		
		
	}
	


/*


GeneFeatures* g;
int nf;
g=readGenes(av[2],&nf);


int i;
for(i=0;i<nf;i++){
 //printf("%s\t%s\t%c\t%d\t%d\n",g[i].gene,g[i].locus,g[i].strand,g[i].start,g[i].end);
}

//printf("long 2845,4797 : %d\n",getLocusLong(g, 2759,4797));
//printf("long 2845,8547 : %d\n",getLocusLong(g, 2845,8547));
//printf("long 2845,8547 : %d\n",getLocusLong(g, 2758,8552));
*/


	
fclose(fe);
fclose(fs);
fclose(fcsv);

return 0;


}
