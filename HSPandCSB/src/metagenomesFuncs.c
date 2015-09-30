/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 24-Sept-2015
 * @description This file contains necessary functions for 
 *     handle metagenome files.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"


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
	if((currentLine = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
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


/* This function returns the position of the pointer when an specific
 * read is found. Extension isn't checked but format must be:
 *    ">READ_ID \n READ_SEQUENCE"
 * Also other information about read could be expressed starting
 * line with ">" but will be obviated.
 * If it isn't the format (-1) will be returned.
 * If a read appears without an ID header it will be obviated.
 * @param metagenome is a file pointer to metagenome file.
 * @param readIndex is the index position on the metagenome file given.
 * @param read is a string where read sequence will be stored. If it's 
 *     null, storage will be obviated.
 * @return the position of file pointer just before take the read from
 *     file. This position is an absolute position.
 */
int seekRead(FILE *metagenome,int readIndex,char *read){
	int saveReadSeq = read==NULL? 0 : 1;
	int currentRead = 0, posBeforeRead = -1,readingRead = 0;
	char *currentLine;

	// Reserve memory on disk
	if((currentLine = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
			fprintf(stderr, "Error: problem reserving space from reads\n");
			return -1;
	}

	// Read
	while(!feof(metagenome) && currentRead < readIndex){
		posBeforeRead = ftell(metagenome); // Save position before read
		fscanf(metagenome,"%s\n",currentLine); // Read line
		if(readingRead){
			if(strchr(currentLine,'>') == NULL){ // Not ">" found == Read_Seq
				currentRead++;
				readingRead = 0;
			}else if(strcmp(currentLine,"\n")==0) //Empty line while searching a read
				return -1; //Ilegal format
		}else if(strchr(currentLine,'>') != NULL){ // Header found == Read_ID
			readingRead = 1;
		}
	}

	if(currentRead != readIndex) // Couldn't find read
		return -1;
	else if(saveReadSeq) // Found!
		strcpy(read,currentLine);

	free(currentLine);

	return posBeforeRead;
}// END FUNCTION



/* This function is used to take a read sequence from a metagenome
 * file given . Extension isn't checked but format must be:
 *    ">READ_ID \n READ_SEQUENCE"
 * Also other information about read could be expressed starting
 * line with ">" but will be obviated.
 * If it isn't the format (-1) will be returned.
 * If a read appears without an ID header it will be obviated.
 * @param metagenome is a file pointer to metagenome file.
 * @param read is a string where read sequence will be stored. If it's 
 *     null, storage will be obviated.
 * @return a positive number if the rocess ended without errors and 
 *     a negative number if any error happend. Posible errors are:
 *          - File pointer given is in the end of the file.
 *          - Format error.
 */
int takeRead(FILE *metag,char *read){
	if(feof(metag)){
		return -1;
	}

	// Prepare necessary variables
	int readingRead = 0, searching = 1; // Boolean checker
	char *currentLine;

	if((currentLine = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
		fprintf(stderr, "Error: problem reserving space from reads. [takeRead]\n");
		return -1;
	}

	// Read
	while(!feof(metagenome) && searching){
		fscanf(metagenome,"%s\n",currentLine); // Read line
		if(readingRead){
			if(strcmp(currentLine,"\n")==0){ //Empty line while searching a read
				return -1; //Ilegal format
			}else if(strchr(currentLine,'>') == NULL){ // Not ">" found == Read_Seq
				strcpy(read,currentLine);
				searching = 0;
				//readingRead = 0;
			}
		}else if(strchr(currentLine,'>') != NULL){ // Header found == Read_ID
			readingRead = 1;
		}
	}

	return 1;
}