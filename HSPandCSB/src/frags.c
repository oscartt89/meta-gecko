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
	Hit *buff;
	uint64_t hitsInBuffer = 0;
	uint16_t gWL = 32,mWL;
	uint16_t BytesGenoWord = 8, BytesMetagWord, MinBytes, MaxBytes;
	int GoM; // If >0 --> gWL > mWL; <0 --> gWL < mWL; =0 --> gWL == mWL
	

	// Allocate necessary memory
	// Memory for buffer
	if((buff = (hit*) malloc(sizeof(hit)*MAX_BUFF))==NULL){
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
	fread(&mWL,sizeof(uint16_t),1,mW);

			// Check
			if(mWL == NULL){ 
				fprintf(stderr, "Error, couldn't find word length.\n");
				return -1;
			}else if(mWL % 4 != 0){
				fprintf(stderr, "Error, word length of metagenome dictionary isn't a 4 multiple.\n");
				return -1;
			}else
				BytesMetagWord = mWL/4;

			// Select minimum and maximum WL
			if(gWL == mWL){
				GoM = 0;
				MinBytes = BytesGenoWord;
				MaxBytes  = BytesGenoWord;
			}else{
				GoM = gWL > mWL? 1 : -1;
				MinBytes = GoM > 0? BytesMetagWord : BytesGenoWord;
				MaxBytes = GoM < 0? BytesMetagWord : BytesGenoWord;
			}


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
		wordEntry we[2]; // [0]-> Metagenome [1]-> Genome
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
	if(readHashEntry(&we[1],gW)<0) return -1;

	// Search
	int cmp
	bool loadMetag,loadGeno;
	while(!feof(mW) && !feof(gW)){
		loadGeno = false;
		loadMetag = false;

		// Compare sequences searching match
		if(GoM == 0){ // Same length
			if((cmp = wordcmp(we[0].seq,we[1].seq,BytesGenoWord))==0){ // Hit
				generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,hitsInBuffer,0,0);
				loadGeno = true;
				loadMetag = true;
			}
		}else if(GoM > 0){ // Metag word is larger
			uint16_t i;
			uint16_t diff = BytesMetagWord - BytesGenoWord;
			for(i=0; i<=diff; ++i)
				if((cmp = wordcmp(&we[0].seq[i],we[1].seq,BytesGenoWord))==0){ // Hit
					generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,hitsInBuffer,i,0);
					
				}
		}else{ // Geno word is larger
			uint16_t i;
			uint16_t diff = BytesGenoWord - BytesMetagWord;
			for(i=0; i<=diff; ++i)
				if((cmp = wordcmp(&we[0].seq,we[1].seq[i],BytesMetagWord))==0) // Hit
						generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,hitsInBuffer,0,i);
		}

		// Load next word
		if()
	}

	// Close files

	// Free space

	// Everything finished OK
	return 0;
}