/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 28-Sept-2015
 * @function read2PGenome (read to pseudo-genome)
 * @description This file contains the implementation necessary for generate
 *      fasta genome format files using specified reads from a metagenome
 *      file given.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"

int main(int ac, char **av){
	FILE *meta, *out; //Metagenome file
	int numReads, readIndex;	
	int useAllReads = 0;

	if(strcmp(av[3],"-all"))
		useAllReads = 1;
	else
		readIndex = atoi(av[3]);

	// Check bad call error
	if(ac!=4)
		fprintf(stderr, "Bad call error. Use: %s metagenomeFile outputFileName readIndex\nOR: metagenomeFile outputFileName -all", av[0]);

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
				useAllReads? numReads : readIndex,
				numReads );
		return -1;
	}

	// Take specific read 
	char *read;

		// Reserve memmory 
		if((read = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
			fprintf(stderr, "Error: couldn't allocate disk space for READ sequence.[%d]\n",sizeof(char)*MAXREADLENGTH);
			return -1;
		}

	
	int initValue = useAllReads? 0 : readIndex;
	int limit = useAllReads? numReads : readIndex+1;
	char *outputF;

	if((outputF = (char*) malloc(MAXFILENAMELENGTH))==NULL){
		fprintf(stderr, "Error: couldn't allocate space on disk for outputF char*. [%d]\n", MAXFILENAMELENGTH);
		return -1;
	}

	for(int i = initValue; i < limit; ++i){
		// Reset value
		strcpy(outputF,av[2]);

		// Seek specific read
		seekRead(meta,readIndex,read);

		// Open output stream
		if((out = fopen(useAllReads? strcat(strcat(outputF,itoa(i)),".fasta") : strcat(outputF,".fasta"),"wt"))==NULL) {
			fprintf(stderr,"Error opening output file. [%s.fasta]\n",av[2]);
			return -1;
		}

		// Create pseudo-genome file
		fprintf(out, "%s\n", read);

		// Close stream
		fclose(out);
	}

	// Close streams
	fclose(meta);

	return 0;
}