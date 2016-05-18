/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes a process to read and take statistical information
 *    from a metagenome file.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "metag_common.h"


int exists(char *file){
	if(access(file,F_OK) != (-1)) return 1;
	else return 0;
}

int main(int ac, char** av){
	// Variables
	uint32_t numReads = 0;
	uint64_t minLength=0, maxLength=0, currLength = 0;
	float meanLength = 0, size;
	char c;
	FILE *metag, *out;
	bool firstSeq = true;

	// Check
	if(ac!=3){
		fprintf(stderr, "Bad call error.\nUSE: metagstats metagenome.file outputFile\n");
		return -1;
	}

	// Open metagenome file
	if((metag = fopen(av[1],"rt"))==NULL){
		fprintf(stderr, "Error opening metagenome file.\n");
		return -1;
	}

	// Open output stats file
	if(exists(av[2])){
		if((out = fopen(av[2],"at"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[4]);
			return -1;
		}
	}else{
	 	if((out = fopen(av[2],"wt"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[4]);
			return -1;
		}
		fprintf(out, "Metagenome\tSize\tReads\tMaxL\tMinL\tMeanL\n");
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

	// Take size
	size = ftell(metag)/1000000; // MB

	// Free&Close unnecessary varaibles
	fclose(metag);

	// Print info
	fprintf(out, "%s\t%f\t%"PRIu32"\t%"PRIu64"\t%"PRIu64"\t%f\n", av[1],size,numReads,maxLength,minLength,meanLength);

	fclose(out);

	// End
	return 0;
}