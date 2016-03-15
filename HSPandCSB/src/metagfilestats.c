/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes a process to read and take statistical information
 *    from a metagenome file.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include <stdio.h> // IO functions
#include <stdint.h> // unit64_t ...
#include <inttypes.h> 
#include <stdbool.h> // Boolean varaibles

int main(int ac, char** av){
	// Variables
	uint32_t numReads = 0;
	uint64_t minLength=0, maxLength=0, currLength = 0;
	float meanLength = 0;
	char c;
	FILE *metag;
	bool firstSeq = true;

	// Check
	if(ac!=2){
		fprintf(stderr, "Bad call error.\nUSE: metagstats metagenome.file\n");
		return -1;
	}

	// Open metagenome file
	if((metag = fopen(av[1],"rt"))==NULL){
		fprintf(stderr, "Error opening metagenome file.\n");
		return -1;
	}

	// Start to read
	c = fgetc(metag);
	while(!feof(metag)){
		// Check if it's a special line
		if(!isupper(toupper(c))){ // Comment, empty or quality (+) line
			if(c=='>'){ // Comment line
				c = fgetc(metag);
				while(c != '\n') // Avoid comment line
					c = fgetc(metag);
				if(currLength > 0) firstSeq = false;
				if(minLength == 0 || currLength < minLength) minLength = currLength;
				if(currLength > maxLength) maxLength = currLength;
				if(!firstSeq)
					meanLength = (meanLength*numReads + currLength)/(numReads+1);
				
				numReads++;
				currLength = 0;
			}
			c=fgetc(metag); // First char of next sequence
			continue;
		}
		// Count nucleotides
		switch (c) {
			case 'A':  
			case 'C': 
			case 'G': 
			case 'T': 
				currLength++;
				break;
		}

		c = fgetc(metag);
	}

	// Free&Close unnecessary varaibles
	fclose(metag);

	// Print info
	fprintf(stdout, "Metagenome stats of %s\n", av[1]);
	fprintf(stdout, "\t Reads: %"PRIu32"\tMeanLength: %f\n", numReads,meanLength);
	fprintf(stdout, "\t MinLength: %"PRIu64"\tMaxLength: %"PRIu64"\n", minLength,maxLength);

	// End
	return 0;
}