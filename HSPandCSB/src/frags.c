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
	if(hitsInBuffer > 0){
		if(buffersWritten > 0)
			writeHitsBuff(buffer,hIndx,hts,hitsInBuffer);
		else{ // Only one buffer
			// Sort buffer
			quicksort_H(buffer,0,hitsInBuffer-1);
			
			// Close unnecesary files
			fclose(mW); fclose(gW);
			fclose(mP); fclose(gP);
			fclose(hIndx);
			fclose(hts);

			// Free unnecesary variables 
			free(we[0].seq);
			free(we[1].seq);

			// Declare necessary variables
			FragFile frag;
			uint64_t index;
			float newSimilarity;
			int64_t dist;
			
			// Init first fragment
			frag.diag = buffer[0].diag;
			frag.xStart = buffer[0].posX;
			frag.yStart = buffer[0].posY;
			frag.xEnd = buffer[0].posX + buffer[0].length;
			frag.yEnd = buffer[0].posY + buffer[0].length;
			frag.length = buffer[0].length;
			frag.ident = buffer[0].length;
			frag.score = frag.ident;
			frag.similarity = 100;
			frag.seqX = buffer[0].seqX;
			frag.seqY = buffer[0].seqY;
			frag.block = 0;
			frag.strand = 'f';

			// Open final files
			// Open final fragments file
			strcpy(fname,av[5]); // Copy outDic name
			if((fr = fopen(strcat(fname,".frags"),"wb"))==NULL){
				fprintf(stderr, "Error opening fragments final file.\n");
				return -1;
			}

			// Generate fragments
			for(index=1; index < hitsInBuffer; ++index){
				if( buffer[index].diag == frag.diag &&
						buffer[index].seqX == frag.seqX && 
						buffer[index].seqY == frag.seqY){ // Possible fragment
					// Check if are collapsable
					dist = (int64_t)buffer[index].posX - (int64_t)frag.xEnd;
					if(dist >= 0){
						newSimilarity = (100*buffer[index].length + frag.length * frag.similarity)/(buffer[index].length + frag.length + dist);
						if(newSimilarity >= S_Threshold){ // Collapse fagments
							frag.length = buffer[index].length + buffer[index].posX - frag.xStart;
								frag.xEnd = frag.xStart + frag.length;
								frag.yEnd = frag.yStart + frag.length;
							frag.ident += buffer[index].length;
							frag.score += buffer[index].length - dist; // Equal +1; Difference -1
							frag.similarity = newSimilarity;
						}else{ // Else write fragment
							writeFragment(frag,fr);
							// Upload new fragment
							frag.xStart = buffer[index].posX;
							frag.yStart = buffer[index].posY;
							frag.xEnd = buffer[index].posX + buffer[index].length;
							frag.yEnd = buffer[index].posY + buffer[index].length;
							frag.length = buffer[index].length;
							frag.ident = buffer[index].length;
							frag.score = frag.ident;
							frag.similarity = 100;
						}
					}else{ // Else it's collapsable
						frag.xEnd = buffer[index].posX + buffer[index].length - frag.xStart;
						frag.yEnd = buffer[index].posY + buffer[index].length - frag.yStart;
						uint64_t oldLength = frag.length;						
						frag.length += buffer[index].length + dist;
						frag.ident += buffer[index].length + dist;
						frag.score += buffer[index].length + dist;
						frag.similarity = frag.ident == frag.length? 100: (frag.ident*100 / frag.length);
					}
				}else{ // New fragment
					// Write fragment
					writeFragment(frag,fr);
					// Upload new frag
					frag.diag = buffer[index].diag;
					frag.xStart = buffer[index].posX;
					frag.yStart = buffer[index].posY;
					frag.xEnd = buffer[index].posX + buffer[index].length;
					frag.yEnd = buffer[index].posY + buffer[index].length;
					frag.length = buffer[index].length;
					frag.ident = buffer[index].length;
					frag.score = frag.ident;
					frag.similarity = 100;
					frag.seqX = buffer[index].seqX;
					frag.seqY = buffer[index].seqY;
				}
			}

			// Last fragment
			writeFragment(frag,fr);

			// Close output file
			fclose(fr);

			// Free unnecesary memory
			free(buffer);

			// Remove intermediate files
			if(removeIntermediataFiles){
				strcpy(fname,av[5]);
				remove(strcat(fname,".hts"));
				strcpy(fname,av[5]);
				remove(strcat(fname,".hindx"));
			}

			free(fname);

			// End program
			return 0;
		}
	}else if(hitsInBuffer == 0 && buffersWritten == 0){
		// Free auxiliar buffers
		free(we[0].seq);
		free(we[1].seq);
		free(buffer);

		// Close files
		fclose(mW); fclose(gW);
		fclose(mP); fclose(gP);
		fclose(hIndx);
		fclose(hts);

		// Remove intermediate files
		if(removeIntermediataFiles){
			strcpy(fname,av[5]);
			remove(strcat(fname,".hts"));
			strcpy(fname,av[5]);
			remove(strcat(fname,".hindx"));
		}

		free(fname);

		// End program
		return 0;
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
	uint64_t hitsUnread[buffersWritten];
	uint64_t positions[buffersWritten];
	uint64_t lastLoaded, activeBuffers = buffersWritten;
	Hit *HitsBlock;
	FragFile frag;

	// Take memory for hits
	if((HitsBlock = (Hit*) malloc(sizeof(Hit)*activeBuffers*READ_BUFF_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for hits block.\n");
		return -1;
	}


	// Read buffers info
	uint64_t i = 0;
	do{
		fread(&positions[i],sizeof(uint64_t),1,hIndx);
		fread(&hitsUnread[i],sizeof(uint64_t),1,hIndx);
		++i;
	}while(i < activeBuffers);

	// Load first hits
	node *currNode;
	uint64_t read, blockIndex = 0;
	for(i=0 ;i<activeBuffers; ++i, blockIndex += READ_BUFF_LENGTH){
		currNode = (node*) malloc(sizeof(node));
		currNode->next = hitsList;
		currNode->hits = &HitsBlock[blockIndex];
		currNode->buff = i;
		fseek(hts,positions[i],SEEK_SET);
		read = loadHit(&currNode->hits,hts,hitsUnread[i]);
		currNode->index = 0;
		currNode->hits_loaded = read;
		// Update info
		positions[i] = (uint64_t) ftell(hts);
		hitsUnread[i]-=read;
		lastLoaded = i;
		hitsList = currNode;
	}

	// Assign head
	hitsList = currNode;

	// Sort hits
	sortList(&hitsList);	

	// First fragment = hit[0]
	frag.diag = hitsList->hits[hitsList->index].diag;
	frag.xStart = hitsList->hits[hitsList->index].posX;
	frag.yStart = hitsList->hits[hitsList->index].posY;
	frag.xEnd = hitsList->hits[hitsList->index].posX + hitsList->hits[hitsList->index].length;
	frag.yEnd = hitsList->hits[hitsList->index].posY + hitsList->hits[hitsList->index].length;
	frag.length = hitsList->hits[hitsList->index].length;
	frag.ident = hitsList->hits[hitsList->index].length;
	frag.score = frag.ident;
	frag.similarity = 100;
	frag.seqX = hitsList->hits[hitsList->index].seqX;
	frag.seqY = hitsList->hits[hitsList->index].seqY;
	frag.block = 0;
	frag.strand = 'f';

	hitsList->index +=1;

	// Load new hit
	if(hitsList->index >= hitsList->hits_loaded){
		if(hitsUnread[hitsList->buff] > 0){
			if(hitsList->buff != lastLoaded){
				fseek(hts,positions[hitsList->buff],SEEK_SET);
				lastLoaded = hitsList->buff;
			}
			read = loadHit(&hitsList->hits,hts,hitsUnread[hitsList->buff]);
			hitsList->index = 0;
			hitsList->hits_loaded = read;
			positions[hitsList->buff] = (uint64_t) ftell(hts);
			hitsUnread[hitsList->buff]-=read;
			checkOrder(&hitsList,false);
		}else{
			checkOrder(&hitsList,true);
			activeBuffers--;
		}
	}else checkOrder(&hitsList,false);

	// Search new fragments
	float newSimilarity;
	int64_t dist;

	while(activeBuffers > 0){
		if(hitsList->hits[hitsList->index].diag == frag.diag &&
			hitsList->hits[hitsList->index].seqX == frag.seqX && 
				hitsList->hits[hitsList->index].seqY == frag.seqY){ // Possible fragment
			// Check if are collapsable
			dist = hitsList->hits[hitsList->index].posX - (frag.xStart + frag.length);
			if(dist >= 0){
				newSimilarity = (100*hitsList->hits[hitsList->index].length + frag.length * frag.similarity)/(hitsList->hits[hitsList->index].length + frag.length + dist);
				if(newSimilarity >= S_Threshold){ // Collapse fagments
					frag.length = hitsList->hits[hitsList->index].length + hitsList->hits[hitsList->index].posX;
						frag.xEnd = frag.xStart + frag.length;
						frag.yEnd = frag.yStart + frag.length;
					frag.ident += hitsList->hits[hitsList->index].length;
					frag.score += hitsList->hits[hitsList->index].length - dist; // Equal +1; Difference -1
					frag.similarity = newSimilarity;
				}else{ // Else write fragment and load next frag
					writeFragment(frag,fr);
					// Upload new fragment
					frag.xStart = hitsList->hits[hitsList->index].posX;
					frag.yStart = hitsList->hits[hitsList->index].posY;
					frag.xEnd = hitsList->hits[hitsList->index].posX + hitsList->hits[hitsList->index].length;
					frag.yEnd = hitsList->hits[hitsList->index].posY + hitsList->hits[hitsList->index].length;
					frag.length = hitsList->hits[hitsList->index].length;
					frag.ident = hitsList->hits[hitsList->index].length;
					frag.score = frag.ident;
					frag.similarity = 100;
				}
			}else{// Else it's collapsable
				frag.xEnd = hitsList->hits[hitsList->index].posX + hitsList->hits[hitsList->index].length - frag.xStart;
				frag.yEnd = hitsList->hits[hitsList->index].posY + hitsList->hits[hitsList->index].length - frag.yStart;
				uint64_t oldLength = frag.length;
				frag.length += hitsList->hits[hitsList->index].length - dist;
				frag.ident += hitsList->hits[hitsList->index].length - dist;
				frag.score += hitsList->hits[hitsList->index].length - dist;
				frag.similarity = frag.ident*100 / frag.length;
			}
		}else{ // New fragment
			// Write fragment
			writeFragment(frag,fr);
			// Upload new frag
			frag.diag = hitsList->hits[hitsList->index].diag;
			frag.xStart = hitsList->hits[hitsList->index].posX;
			frag.yStart = hitsList->hits[hitsList->index].posY;
			frag.xEnd = hitsList->hits[hitsList->index].posX + hitsList->hits[hitsList->index].length;
			frag.yEnd = hitsList->hits[hitsList->index].posY + hitsList->hits[hitsList->index].length;
			frag.length = hitsList->hits[hitsList->index].length;
			frag.ident = hitsList->hits[hitsList->index].length;
			frag.score = frag.ident;
			frag.similarity = 100;
			frag.seqX = hitsList->hits[hitsList->index].seqX;
			frag.seqY = hitsList->hits[hitsList->index].seqY;
		}
		hitsList->index +=1;
		// Load new hit
		if(hitsList->index >= hitsList->hits_loaded){
			if(hitsUnread[hitsList->buff] > 0){
				if(hitsList->buff != lastLoaded){
					fseek(hts,positions[hitsList->buff],SEEK_SET);
					lastLoaded = hitsList->buff;
				}
				read = loadHit(&hitsList->hits,hts,hitsUnread[hitsList->buff]);
				hitsList->index = 0;
				hitsList->hits_loaded = read;
				positions[hitsList->buff] = (uint64_t) ftell(hts);
				hitsUnread[hitsList->buff]-=read;
				checkOrder(&hitsList,false);
			}else{
				checkOrder(&hitsList,true);
				activeBuffers--;
			}
		}else checkOrder(&hitsList,false);
	}

	// Close files
	free(HitsBlock);
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

	free(fname);

	// Everything finished OK
	return 0;
}