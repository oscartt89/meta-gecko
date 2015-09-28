/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 28-Sept-2015
 * @description This file contains necessary functions for 
 *     handle metagenome files.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"

int main(int ac, char **av){
	FILE *meta, *out; //Metagenome file
	int numReads;	
	int readIndex = atoi(av[3]);

	// Check bad call error
	if(ac!=4)
		fprintf(stderr, "Bad call error. Use: %s metagenomeFile outputFile readIndex\n", av[0]);

	// Open input stream
	if((meta = fopen(av[1],"r"))==NULL) {
		fprintf(stderr,"Error opening metagenome file. [%s]\n",av[1]);
		return -1;
	}

	// Count reads
	numReads = countReads(meta);

	// Reset file read pointer
	fseek(meta,0,SEEK_SET);

	// Check index out of bounds
	if(readIndex <= 0 | readIndex > numReads){
		fprintf(stderr, "Error: index of read given (%d) is out of bounds.\nNumber of reads=(%d)\n",
				readIndex,numReads );
		return -1;
	}

	

	// Take specific read 
	char *read;

		// Reserve memmory 
		if((read = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
			fprintf(stderr, "Error: couldn't allocate disk space for READ sequence.[%d]\n",sizeof(char)*MAXREADLENGTH);
			return -1;
		}

	// Seek specific read
	int pointerPosition = seekRead(meta,readIndex,read);

	// Open output stream
	if((out = fopen(av[2],"wt"))==NULL) {
		fprintf(stderr,"Error opening output file. [%s]\n",av[2]);
		return -1;
	}

	// Create pseudo-genome file
	fprintf(out, "%s\n", read);

	// Close streams
	fclose(meta);
	fclose(out);

	return 0;
}