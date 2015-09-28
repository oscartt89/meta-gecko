/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 24-Sept-2015
 * @description This file contains necessary functions for 
 *     handle metagenome files.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"
#include "stdlib.h" // Por alguna razon no lo coge del header

#define MAXREADLENGTH 100000 // Por alguna razon no lo coge del header

/* This function returns the number of reads on a metagenome
 * fasta file. Extension isn't checked but format must be:
 *    ">READ_ID \n READ_SEQUENCE"
 * Also other information about read could be expressed starting
 * line with ">" but will be obviated.
 * If it isn't the format (-1) will be returned.
 * If a read appears without an ID header it will be obviated.
 * @param metagenome is a file pointer to metagenome file.
 * @return the number of reads on a metagenome file or a
 *    negative number if couldn't count it.
 */
int countReads(FILE *metagenome){
	int numReads = 0;
	char *currentLine;
	int readingRead = 0;

	// Reserve memory on disk
	if((currentLine = (char*) malloc(MAXREADLENGTH))==NULL){
			fprintf(stderr, "Error: problem reserving space from reads\n");
			return -1;
	}

	// Read
	while(!feof(metagenome)){
		fscanf(metagenome,"%s\n",currentLine); // Read line
		if(readingRead){
			if(strchr(currentLine,'>') == NULL){ // Not ">" found == Read_Seq
				numReads++;
				readingRead = 0;
			}else if(strcmp(currentLine,"\n")==0) //Empty line while searching a read
				return -1; //Ilegal format
		}else if(strchr(currentLine,'>') != NULL){ // Header found == Read_ID
			readingRead = 1;
		}
	}

	if(readingRead==1) return -1; //File ended while searching a read.

	return numReads;
}