/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes a process to read and take statistical information
 *    from a metagenome dictionary file.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include <stdio.h> // IO functions
#include <stdint.h> // unit64_t ...
#include <inttypes.h> 
#include <stdbool.h> // Boolean varaibles

int main(int ac, char** av){
	// Variables
	uint16_t wl;
	uint32_t aux32, minReps, maxReps;
	uint64_t numWords = 0, aux64;
	float meanReps = 0;
	FILE *dict;

	// Check
	if(ac!=2){
		fprintf(stderr, "Bad call error.\nUSE: metagdictstats d2hWDict\n");
		return -1;
	}

	// Open metagenome file
	if((dict = fopen(av[1],"rt"))==NULL){
		fprintf(stderr, "Error opening metagenome dictionary file.\n");
		return -1;
	}
	fread(&wl,sizeof(uint16_t),1,dict);

	char kmer[wl/4];

	// Start to read
		// Read first kmer
		fread(&kmer[0],sizeof(unsigned char),wl/4,dict);
		fread(&aux64,sizeof(uint64_t),1,dict);
		fread(&aux32,sizeof(uint32_t),1,dict);

		// Init stats
		minReps = aux32;
		maxReps = aux32;

	while(!feof(dict)){
		meanReps = (meanReps*numWords + aux32)/(numWords+1);
		numWords++;

		if(aux32 < minReps) minReps = aux32;
		if(aux32 > maxReps) maxReps = aux32;

		fread(&kmer[0],sizeof(unsigned char),wl/4,dict);
		fread(&aux64,sizeof(uint64_t),1,dict);
		fread(&aux32,sizeof(uint32_t),1,dict);
	}


	// Free&Close unnecessary varaibles
	fclose(dict);

	// Print info
	fprintf(stdout, "Metagenome dict stats of %s\n", av[1]);
	fprintf(stdout, "\t Words: %"PRIu64"\tMeanReps: %f\n", numWords,meanReps);
	fprintf(stdout, "\t MinReps: %"PRIu32"\tMaxReps: %"PRIu32"\n", minReps,maxReps);

	// End
	return 0;
}