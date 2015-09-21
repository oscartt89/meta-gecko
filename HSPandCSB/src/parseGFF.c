/*
 * parseGFF.c
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

void parse(char* ptr,GeneFeatures* gene);
void parseLine2 (GeneFeatures* gene, char* linea2);
/*****/
int main(int ac,char** av){

if(ac<3){
	printf("Use: parseGFF input.gbk.features output.gene output.csv \n");
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


GeneFeatures gene,gene_ant;

gene_ant.start=0;
gene_ant.end=0;

char linea[2000];
char linea2[2000];
int max=1000;
int salto=1;

char basurilla[50];
	
//CM000363.1_gene_1	Dsim\GD13543	Dsim_GD13543	2226..3544	
fprintf(fcsv,"id,gene,locus,strand,start,end\n");
	
fgets(linea,max,fe);
	while(!feof(fe)){
		
		// For genes
		// Normal line
		// NC_007332.1	RefSeq	region	1	920079	.	+	.	ID=id0;Dbxref=taxon:262722;Is_circular=true;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=7448

		if ( (sscanf(linea,"%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%[^\n]\n",
		basurilla,basurilla,basurilla,&gene.start,&gene.end,&gene.strand,linea2)==7) ||
			(sscanf(linea,"%s\t%s\t%s\t%d\t%d\t.\t%c\t0\t%[^\n]\n",
		basurilla,basurilla,basurilla,&gene.start,&gene.end,&gene.strand,linea2)==7) ){
		
			if(gene.strand=='+'){
				gene.strand='f';
			}else{
				gene.strand='r';
			}
			
			sprintf(gene.gene,"-");
			sprintf(gene.locus,"-");
			sprintf(gene.product,"-");
			
			parseLine2(&gene,linea2);
			
			if(salto==0){
				if(gene_ant.start==gene.start && gene_ant.end==gene.end){
					sprintf(gene_ant.product,"%s",gene.product);
					salto=1;
				}
				fwrite(&gene_ant,sizeof(GeneFeatures),1,fs);
				fprintf(fcsv,"%d,%d,%c,%s,%s,%s\n",gene_ant.start,gene_ant.end,gene_ant.strand,gene_ant.gene,gene_ant.locus,gene_ant.product);
				
			}else{
				salto=0;
				
			}	
			memcpy(&gene_ant,&gene,sizeof(GeneFeatures));
				//			 printf("%d,%d,%c,%s,%s,%s\n",gene.start,gene.end,gene.strand,gene.gene,gene.locus,gene.product);
		}
		
		fgets(linea,max,fe);
		
		
	}
	




	
fclose(fe);
fclose(fs);
fclose(fcsv);

return 0;


}

/************************/
void parseLine2 (GeneFeatures* gene, char* linea){


char* fin=";=";

char *ptr;


ptr = (char*)malloc(sizeof(char)*2000);

// printf( "%s\n\n", linea );
 
 ptr = strtok( linea, fin );    // Primera llamada => Primer token
// printf( "%s\n", ptr );

	parse(ptr,gene);
	while( (ptr = strtok( NULL, fin )) != NULL ){    // Posteriores llamadas
		parse(ptr,gene);
		
//		printf( "%s\n", ptr );	 
	}

	
}
/******/
void parse(char* ptr,GeneFeatures* gene){

char geneS[4] = "gene";
char locus_tag[9] = "locus_tag";
char product[7] = "product";
char* fin=";=";
 

 
	if(!strcmp( geneS, ptr )){
		ptr = strtok( NULL, fin );
		
		sprintf(gene->gene,"%s",ptr);
		//printf("geneS: %s\n",gene->gene);
	}
   if(!strcmp( locus_tag, ptr )){
   ptr = strtok( NULL, fin );
	sprintf(gene->locus,"%s",ptr);
   }
   if(!strcmp( product, ptr )){
   ptr = strtok( NULL, fin );
	sprintf(gene->product,"%s",ptr);
   }
   


}
