/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes the workflow of GECKO for create
 *    metagenome dictionaries.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "frags.h"

int main(int ac, char** av){
	// Check arguments
	if(ac!=6){
		fprintf(stderr, "Bad call error.\nUSE: frags metagP metagW genoP genoW out\n");
		return -1;
	}

	// Variables
	FILE *mW,*mP,*gW,*gP; // Dictionaries
	FILE *hIndx,*hts; // Intermediate files
	FILE *fr; // Fragments file
	Hit *buffer;
	uint64_t hitsInBuffer = 0;
	uint16_t gWL = 32,mWL;
	uint16_t BytesGenoWord = 8, BytesMetagWord, MinBytes, MaxBytes;
	buffersWritten = 0;
	char *fname;
	//int GoM; // If >0 --> gWL > mWL; <0 --> gWL < mWL; =0 --> gWL == mWL
	

	// Allocate necessary memory
	// Memory for buffer
	if((buffer = (Hit*) malloc(sizeof(Hit)*MAX_BUFF))==NULL){
		fprintf(stderr, "Error allocating memory for hits buffer.\n");
		return -1;
	}

	// Memory for file names handler
	if((fname = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for file names handler.\n");
		return -1;
	}

	// Open current necessary files
	// Open metagenome positions file
	if((mP = fopen(av[1],"rb"))==NULL){
		fprintf(stderr, "Error opening metagenome positions dictionaries.\n");
		return -1;
	}

	// Open metagenome words file
	if((mW = fopen(av[2],"rb"))==NULL){
		fprintf(stderr, "Error opening metagenome words dictionaries.\n");
		return -1;
	}
	// Read words header = WordLength
		// Check
		if(fread(&mWL,sizeof(uint16_t),1,mW)!=1){ 
			fprintf(stderr, "Error, couldn't find word length.\n");
			return -1;
		}else if(mWL % 4 != 0){
			fprintf(stderr, "Error, word length of metagenome dictionary isn't a 4 multiple.\n");
			return -1;
		}else{
			BytesMetagWord = mWL/4;
			if(BytesMetagWord != BytesGenoWord){
				fprintf(stderr, "Error: metagenome and genome dictionaries have differents word lengths.\n");
				return -1;
			}
		}

			// Select minimum and maximum WL
//			if(gWL == mWL){
//				GoM = 0;
//				MinBytes = BytesGenoWord;
//				MaxBytes  = BytesGenoWord;
//			}else{
//				GoM = gWL > mWL? 1 : -1;
//				MinBytes = GoM > 0? BytesMetagWord : BytesGenoWord;
//				MaxBytes = GoM < 0? BytesMetagWord : BytesGenoWord;
//			}


	// Open genome postions file
	if((gP = fopen(av[3],"rb"))==NULL){
		fprintf(stderr, "Error opening genome positions dictionaries.\n");
		return -1;
	}

	// Open genome words file
	if((gW = fopen(av[4],"rb"))==NULL){
		fprintf(stderr, "Error opening genome words dictionaries.\n");
		return -1;
	}

	//Open intermediate files
	strcpy(fname,av[5]); // Copy outDic name
	if((hIndx = fopen(strcat(fname,".hindx"),"wb"))==NULL){
		fprintf(stderr, "Error opening buffer index file.\n");
		return -1;
	}

	// Open hits repo
	strcpy(fname,av[5]);
	if((hts = fopen(strcat(fname,".hts"),"wb"))==NULL){
		fprintf(stderr, "Error opening hits repository.\n");
		return -1;
	}

	// Search hits
		// Prepare necessary variables
		WordEntry we[2]; // [0]-> Metagenome [1]-> Genome
		// Take memory
		if((we[0].seq = (unsigned char *) malloc(sizeof(unsigned char)*BytesMetagWord))==NULL){
			fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
			return -1;
		}else we[0].WB = BytesMetagWord;

		if((we[1].seq = (unsigned char *) malloc(sizeof(unsigned char)*BytesGenoWord))==NULL){
			fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
			return -1;
		}else we[1].WB = BytesGenoWord;

	// Read first entrances
	if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;
	readHashEntry(&we[1],gW);

	// Search
	int cmp;
	while(!feof(mW) && !feof(gW)){
		if((cmp = wordcmp(we[0].seq,we[1].seq,BytesGenoWord))==0) // Hit
			generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,&hitsInBuffer);

		// Load next word
		if(cmp >= 0) // New genome word is necessary
			readHashEntry(&we[1],gW);
		if(cmp <= 0) // New metagenome word is necessary
			if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;
	}

	// Write buffered hits
	if(hitsInBuffer > 0){
		writeHitsBuff(buffer,hIndx,hts,hitsInBuffer);
		buffersWritten++;
	}

	// Free auxiliar buffers
	free(we[0].seq);
	free(we[1].seq);
	free(buffer);

	// Close files
	fclose(mW); fclose(gW);
	fclose(mP); fclose(gP);
	fclose(hIndx);
	fclose(hts);

	// Open necessary files



	// Free space

	// Everything finished OK
	return 0;
}