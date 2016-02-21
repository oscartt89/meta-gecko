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
	if(ac!=7){
		fprintf(stderr, "Bad call error.\nUSE: frags metagP metagW genoP genoW out minS\n");
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
	S_Threshold = (uint64_t) atoi(av[6]);
	char *fname;
	bool removeIntermediataFiles = true;	

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
	if(hitsInBuffer > 0)
		writeHitsBuff(buffer,hIndx,hts,hitsInBuffer);

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
	// Open intermediate files
	// Index file
	strcpy(fname,av[5]); // Copy outDic name
	if((hIndx = fopen(strcat(fname,".hindx"),"rb"))==NULL){
		fprintf(stderr, "Error opening buffer index file (read).\n");
		return -1;
	}

	// Open hits repo
	strcpy(fname,av[5]);
	if((hts = fopen(strcat(fname,".hts"),"rb"))==NULL){
		fprintf(stderr, "Error opening hits repository (read).\n");
		return -1;
	}

	// Open final fragments file
	strcpy(fname,av[5]); // Copy outDic name
	if((fr = fopen(strcat(fname,".frags"),"wb"))==NULL){
		fprintf(stderr, "Error opening fragments final file.\n");
		return -1;
	}

	// Prepare necessary variables
	node *hitsList = NULL;
	int64_t hitsUnread[buffersWritten];
	uint64_t positions[buffersWritten];
	uint64_t lastLoaded, activeBuffers = buffersWritten;
	FragFile frag;


	// Read buffers info
	uint64_t i = 0, aux64;
	do{
		fread(&positions[i],sizeof(uint64_t),1,hIndx);
		fread(&aux64,sizeof(uint64_t),1,hIndx);
		hitsUnread[i] = (int64_t) aux64;
		++i;
	}while(i < buffersWritten);

	// Load first hits
	node *currNode;
	for(i=0 ;i<buffersWritten; ++i){
		currNode = (node*) malloc(sizeof(node));
		currNode->next = hitsList;
		fseek(hts,positions[i],SEEK_SET);
		loadHit(&currNode->hits,hts);
		currNode->buff = i;
		// Update info
		positions[i] = (uint64_t) ftell(hts);
		hitsUnread[i]--;
		lastLoaded = i;
	}

	// Assign head
	hitsList = currNode;

	// Sort hits
	sortList(&hitsList);	

	// First fragment = hit[0]
	frag.diag = hitsList->hits.diag;
	frag.xStart = hitsList->hits.posX;
	frag.yStart = hitsList->hits.posY;
	frag.xEnd = hitsList->hits.posX + hitsList->hits.length;
	frag.yEnd = hitsList->hits.posY + hitsList->hits.length;
	frag.length = hitsList->hits.length;
	frag.ident = hitsList->hits.length;
	frag.score = frag.ident;
	frag.similarity = 100;
	frag.seqX = hitsList->hits.seqX;
	frag.seqY = hitsList->hits.seqY;
	frag.block = 0;
	frag.strand = 'f';

	// Load new hit
	if(hitsUnread[hitsList->buff] > 0){
		if(hitsList->buff != lastLoaded){
			fseek(hts,positions[hitsList->buff],SEEK_SET);
			lastLoaded = hitsList->buff;
		}
		loadHit(&hitsList->hits,hts);
		positions[hitsList->buff] = (uint64_t) ftell(hts);
		hitsUnread[hitsList->buff]--;
		checkOrder(&hitsList,false);
	}else{
		checkOrder(&hitsList,true);
		activeBuffers--;
	}

	// Search new fragments
	float newSimilarity;
	int64_t dist;

	while(activeBuffers > 0){
		if(hitsList->hits.seqX == frag.seqX && 
				hitsList->hits.seqY == frag.seqY &&
				hitsList->hits.diag == frag.diag){ // Possible fragment
			// Check if are collapsable
			dist = hitsList->hits.posX - frag.xStart + frag.length;
			if(dist >= 0){
				newSimilarity = (100*hitsList->hits.length + frag.length * frag.similarity)/(hitsList->hits.length + frag.length + dist);
				if(newSimilarity >= S_Threshold){ // Collapse fagments
					frag.length = hitsList->hits.length + hitsList->hits.posX;
						frag.xEnd = frag.xStart + frag.length;
						frag.yEnd = frag.yStart + frag.length;
					frag.ident += hitsList->hits.length;
					frag.score += hitsList->hits.length - dist; // Equal +1; Difference -1
					frag.similarity = newSimilarity;
				}else{ // Else write fragment and load next frag
					writeFragment(frag,fr);
					// Upload new fragment
					frag.xStart = hitsList->hits.posX;
					frag.yStart = hitsList->hits.posY;
					frag.xEnd = hitsList->hits.posX + hitsList->hits.length;
					frag.yEnd = hitsList->hits.posY + hitsList->hits.length;
					frag.length = hitsList->hits.length;
					frag.ident = hitsList->hits.length;
					frag.score = frag.ident;
					frag.similarity = 100;
				}
			}// Else it's collapsable
		}else{ // New fragment
			// Write fragment
			writeFragment(frag,fr);
			// Upload new frag
			frag.diag = hitsList->hits.diag;
			frag.xStart = hitsList->hits.posX;
			frag.yStart = hitsList->hits.posY;
			frag.xEnd = hitsList->hits.posX + hitsList->hits.length;
			frag.yEnd = hitsList->hits.posY + hitsList->hits.length;
			frag.length = hitsList->hits.length;
			frag.ident = hitsList->hits.length;
			frag.score = frag.ident;
			frag.similarity = 100;
			frag.seqX = hitsList->hits.seqX;
			frag.seqY = hitsList->hits.seqY;
		}
		// Load new hit
		if(hitsUnread[hitsList->buff] > 0){
			if(hitsList->buff != lastLoaded){
				fseek(hts,positions[hitsList->buff],SEEK_SET);
				lastLoaded = hitsList->buff;
			}
			loadHit(&hitsList->hits,hts);
			positions[hitsList->buff] = (uint64_t) ftell(hts);
			hitsUnread[hitsList->buff]--;
			checkOrder(&hitsList,false);
		}else{
			checkOrder(&hitsList,true);
			activeBuffers--;
		}
	}

	// Close files
	fclose(hIndx);
	fclose(hts);
	fclose(fr);

	// Remove intermediate files
	if(removeIntermediataFiles){
		strcpy(fname,av[5]);
		remove(strcat(fname,".hts"));
		strcpy(fname,av[5]);
		remove(strcat(fname,".hindx"));
	}

	// Everything finished OK
	return 0;
}