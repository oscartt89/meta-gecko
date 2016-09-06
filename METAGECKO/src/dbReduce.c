/*********

File		dbReduce.c
Author		EPW <estebanpw@uma.es>
Description	Cuts out a list of genomes from a database to reduce the comparison

USAGE		<genomes_database>		The .fasta database containing the genomes
			<sequence_hits_histogram>	The binary file containing the sequences histogram of hits produced by dictMgMem
			<working_mode> <n>	Can be [1,2,3...] See below.

WORKING MODES
1	If it has at least n hits [DEFAULT OPTION; Default value n=100]
2	Select only those genomes with at least n% of hits from the total hits achieved by the most-hitted genome [Default n=2]
3	Select only the sequences with hits number n times above the average hits [Default 1]
4	Select only the sequences who satisfy the following expression:
				Be 'Al' the average of the log10 of the hits. Then, a sequence 'Si' will be selected if 
				the Log2 of its hits 'Hi' is higher or equal than 'Al'. Therefore:
					f(Si) 	|	log2(Hi)  >= Al(Hi) 		1
							|	otherwise 					0
				Note that n is not used.
5 	Use a normal distribution of mean being the average and standard deviation computed from the sample. Use n as minimum
	pvalue to filter genomes


**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "metag_common.h"
#include "common.h"

// Returns the probability of x, given the distribution described by mu and sigma.
long double pdf(long double x, long double mu, long double sigma){
	return expl( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
}

// Computes the (approximated) cumulative distribution function in the interval [-inf,x] of a gaussian distribution
long double cdf(long double x, long double mu, long double sigma){
	return 0.5 * (1 + erfl((x - mu) / (sigma * sqrtl(2.))));
}


#define LINE_BUFFER 5000

int main(int argc, char ** av){
	if(argc < 4) terror("USE: dbReduce <genome_database> <sequence_hits_histogram> <working_mode> <n>");
	FILE * f, * histdb, * dbout;
	

	//Open files
	f = fopen64(av[1], "rt");
	if(f == NULL) terror("Error opening genome DB file");

	histdb = fopen64(av[2], "rb");
	if(histdb == NULL) terror("Could not open histogram of sequences");

	// Working mode selection
	int workingMode = atoi(av[3]);
	uint64_t workingMode_value;
	long double workingMode_value_float;

	if(workingMode == 5){ //N is a probability in this case
		workingMode_value_float = atof(av[4]);
		workingMode_value = 1; //To skip warnings...
	}else{ //All other cases use unsigned integers
		workingMode_value = asciiToUint64(av[4]);
		workingMode_value_float = 1.0; //To skip warnings...
	}

	//Reading and writing variables
	char buffer[LINE_BUFFER];
	char name_dbout[LINE_BUFFER];
	name_dbout[0] = '\0';



	strcat(name_dbout, av[1]);
	strcat(name_dbout, ".presec");
	fprintf(stdout, "[INFO] Opening output file :%s\n", name_dbout);
	dbout = fopen64(name_dbout, "wt");
	
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

	//Create mask to tell which genomes are to be used
	unsigned char * mask = (unsigned char *) malloc(seqsRead*sizeof(unsigned char));
	if(mask == NULL) terror("Could not allocate mask");

	//Fill mask depending on working mode
	uint64_t i;
	switch(workingMode){

		//2	Select only those genomes with at least n% of hits from the total hits achieved by the most-hitted genome
		case 2: {
			uint64_t maxFound = 0;
			if(argc == 4) workingMode_value = 2;
			for(i=0;i<seqsRead;i++){
				if(seqHist[i] > maxFound) maxFound = seqHist[i];
			}
			for(i=0;i<seqsRead;i++){
				mask[i] = (((long double)100*seqHist[i]/maxFound >= (long double) workingMode_value)) ? (1) : (0);
			}
		}
		break;

		//3	Select only the sequences with hits number above the average hits
		case 3: {
			uint64_t hitSum = 0;
			uint64_t average;
			if(argc == 4) workingMode_value = 1;
			workingMode_value = ((workingMode_value != 1)) ? (workingMode_value) : (1);
			for(i=0;i<seqsRead;i++){
				hitSum += seqHist[i];
			}
			average = hitSum / seqsRead;
			for(i=0;i<seqsRead;i++){
				mask[i] = ((seqHist[i] >=  average*workingMode_value)) ? (1) : (0);
			}
		}
		break;

		//4	Select only if log2 of the hits is higher than log10 of average hits
		case 4: {
			long double log2average;
			uint64_t hitSum = 0;
			uint64_t average;
			uint64_t * computedValues = (uint64_t *) malloc(seqsRead*sizeof(uint64_t));
			for(i=0;i<seqsRead;i++){
				computedValues[i] = 0;
				if(seqHist[i] > 0){
					computedValues[i] = (uint64_t) log2l((long double)seqHist[i]);
				}
				hitSum += computedValues[i];
			}
			average = hitSum / seqsRead;
			log2average = log2l((long double) average);

			for(i=0;i<seqsRead;i++){
				if(seqHist[i] > 0) mask[i] = (( computedValues[i] >=  log2average)) ? (1) : (0); else mask[i] = 0;
			}
			free(computedValues);

		}
		break;

		case 5: {
			long double stdev = 0;
			long double average = 0;
			for(i=0;i<seqsRead;i++){
				average += (long double) seqHist[i];
			}
			average /= seqsRead;
			for(i=0;i<seqsRead;i++){
				stdev += ((long double)seqHist[i]-average)*((long double)seqHist[i]-average); 
			}
			stdev *= (1/(long double)seqsRead);
			stdev = sqrtl(stdev);

			long double pvalue;
			uint64_t minHits = 0xFFFFFFFFFFFFFFFF; //Max in uint64_t

			for(i=0;i<seqsRead;i++){


				//Bryc, W., 2002. “A Uniform Approximation to the Right Normal Tail Integral,” Applied Mathematics and Computation, Vol. 127, 365-374
				//pvalue = (zscore + 3.333) / (  sqrtl(2*M_PI*zscore*zscore) + 7.32 * zscore + 6.666 ) * expl( - (zscore*zscore) / 2 );

				//Zogheib, Bashar, and Myron Hlynka. Approximations of the Standard Normal Distribution. University of Windsor, Department of Mathematics and Statistics, 2009.
				//pvalue = (0.5 - 0.398942*zscore - 0.066490*zscore*zscore*zscore + 0.09974*zscore*zscore*zscore*zscore*zscore);


				

				if(minHits < seqHist[i]){
					mask[i] = 1; //If a smaller amount of hits yielded a cdf bigger than the filter in a previous occasion, we can use this to skip the calculation again
					continue;
				}else{
					//A Poisson distribution can be approximated using a normal if lambda is sufficiently large
					pvalue = cdf(seqHist[i], average, sqrtl(average));
				}

				if(pvalue >= workingMode_value_float){
					mask[i] = 1;
					if(minHits > seqHist[i]) minHits = seqHist[i]; //Look-up table
				}else{
					mask[i] = 0;
				}
			}

		}

		break;

		//1	If it has at least n hits [DEFAULT]
		default: {
			if(argc == 4) workingMode_value = 100;
			for(i=0;i<seqsRead;i++){
				mask[i] = ((seqHist[i] >= workingMode_value)) ? (1) : (0);
			}
		}
	}


	
	//Only for display of values
	/*
	for(i=0;i<seqsRead;i++){
		fprintf(stdout, "%"PRIu64" :->: %"PRIu64" [MASK:%d]\n", i, seqHist[i], (int) mask[i]);
	}
	*/
	


	//Read database
	buffer[0] = 'N'; //Not a '>'
	int64_t genomeCounter=-1;
	while(!feof(f)){
		
		if(buffer[0]=='>'){ //If there is a sequence
			genomeCounter++;

			//If the genome satisfies the previously calculated condition
			if(mask[genomeCounter] == 1){ 

				fprintf(dbout, "%s", buffer); //Skip first '>'
				fgets(buffer, LINE_BUFFER, f);
				while(buffer[0] != '>' && !feof(f)){
					fprintf(dbout, "%s", buffer);
					fgets(buffer, LINE_BUFFER, f);
				}
			}else{
				fgets(buffer, LINE_BUFFER, f); //Skip the '>' if we dont want it				
			}
		}

		if(buffer[0] != '>') fgets(buffer, LINE_BUFFER, f);
	}
	


	fclose(f);
	fclose(histdb);
	fclose(dbout);

	free(seqHist);
	free(mask);
	return 0;
}


