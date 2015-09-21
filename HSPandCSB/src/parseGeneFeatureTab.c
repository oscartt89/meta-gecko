/*
 * parsegeneFeatureTab.c
 *
 *  Created on: 06/10/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	
		
Given a file with format:

CM000363.1_gene_1	Dsim\GD13543	Dsim_GD13543	2226..3544
CM000363.1_gene_8	Dsim\GD13539	Dsim_GD13539	complement(119623..122214)

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
	int cov;

}GeneFeatures;
*/



/*****/
int main(int ac,char** av){

if(ac<3){
	printf("Use: parseGeneFeatureTab input.fasta.features output.gene output.csv \n");
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
	
//CM000363.1_gene_1	Dsim\GD13543	Dsim_GD13543	2226..3544	
fprintf(fcsv,"id,gene,locus,strand,start,end\n");
	
fgets(linea,max,fe);
	while(!feof(fe)){
		sprintf(gene.product,"-");
		// For genes
		// Normal line
		if(sscanf(linea,"%s\t%s\t%s\t%d..%d\n",gene.id,gene.gene,gene.locus,&gene.start,&gene.end)==5){
			
			gene.strand='f';
			fwrite(&gene,sizeof(GeneFeatures),1,fs);
			fprintf(fcsv,"%s,%s,%s,%c,%d,%d\n",gene.id,gene.gene,gene.locus,gene.strand,gene.start,gene.end);

		}
		// Reverse
		if(sscanf(linea,"%s\t%s\t%s\tcomplement(%d..%d)\n",gene.id,gene.gene,gene.locus,&gene.start,&gene.end)==5){
			gene.strand='r';
			fwrite(&gene,sizeof(GeneFeatures),1,fs);
			fprintf(fcsv,"%s,%s,%s,%c,%d,%d\n",gene.id,gene.gene,gene.locus,gene.strand,gene.start,gene.end);

		}
		// For locus
		// Normal line
		if(sscanf(linea,"%s\t%s\t%d..%d\n",gene.id,gene.locus,&gene.start,&gene.end)==4){
			sprintf(gene.gene,"-");
			gene.strand='f';
			fwrite(&gene,sizeof(GeneFeatures),1,fs);
			fprintf(fcsv,"%s,%s,%s,%c,%d,%d\n",gene.id,gene.gene,gene.locus,gene.strand,gene.start,gene.end);

		}
		// Reverse
		if(sscanf(linea,"%s\t%s\tcomplement(%d..%d)\n",gene.id,gene.locus,&gene.start,&gene.end)==4){
			sprintf(gene.gene,"-");
			gene.strand='r';
			fwrite(&gene,sizeof(GeneFeatures),1,fs);
			fprintf(fcsv,"%s,%s,%s,%c,%d,%d\n",gene.id,gene.gene,gene.locus,gene.strand,gene.start,gene.end);

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
