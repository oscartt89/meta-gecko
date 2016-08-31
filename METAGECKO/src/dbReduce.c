/*********

File		dbReduce.c
Author		EPW <estebanpw@uma.es>
Description	Cuts out a list of genomes from a database to reduce the comparison

USAGE		<genomes_database>		The .fasta database containing the genomes
			<sequence_hits_histogram>	The binary file containing the sequences histogram of hits produced by dictMgMem

**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "metag_common.h"

//goes [a, b)

#define LINE_BUFFER 5000

int main(int argc, char ** av){
	if(argc != 3) terror("USE: dbReduce <genome_database> <sequence_hits_histogram");
	FILE * f, * histdb, * dbout;
	

	//Open files
	f = fopen64(av[1], "rt");
	if(f == NULL) terror("Error opening genome DB file");

	histdb = fopen64(av[2], "rb");
	if(histdb == NULL) terror("Could not open histogram of sequences");

	char buffer[LINE_BUFFER];
	int bufferIdx;



	strcat(name, av[1]);
	strcat(name, "_cut.fasta");
	fprintf(stdout, "[INFO] Opening output file :%s\n", name);
	dbout = fopen64(name, "wt");
	
	//Get number of sequences
    fseeko64(histdb, 0L, SEEK_END);
    uint64_t totalSize = ftello64(histdb);
    fseeko64(histdb, 0L, SEEK_SET);


    //Allocate sequence histogram
    uint64_t nSeqs = totalSize/sizeof(uint64_t);
    fprintf(stdout, "[INFO] Loading %"PRIu64" sequences\n", nSeqs);
	uint64_t * seqHist = (uint64_t *) malloc(nSeqs * sizeof(uint64_t));
	if(seqHist == NULL) terror("Could not allocate memory for sequences histogram");

	//Load sequence histogram into memory
	uint64_t seqsRead = fread(seqHist, sizeof(uint64_t), nSeqs, histdb);
	if(seqsRead != nSeqs) terror("Different amount of read sequences");


	//Read database
	char c;
	uint64_t position;
	int64_t genomeCounter=-1;
	c = getc(f);
	while(!feof(f)){
		
		if(c=='>'){
			genomeCounter++;
			position = ftell(f);
			//If we want the genome CONDITION TODO
			if(seqHist[genomeCounter] > 0){
				//USE A BUFFER
			}
		}
	}
	


	fclose(f);
	fclose(histdb);
	fclose(dbout);

	free(seqHist);
	
}


