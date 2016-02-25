/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "frags.h"

/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
int wordcmp(unsigned char *w1, unsigned char *w2, int n){
	int i;
	for(i=0;i<n;i++)
		if(w1[i] < w2[i]) return -1;
		else if(w1[i] > w2[i]) return +1;

	return 0;
}


/* Function used to compare two Hit variables. It function sort with this
 * criterion:
 *    1 - Sequence X
 *    2 - Sequence Y
 *    3 - Diagonal
 *    4 - Position on X
 *    5 - Length 
 *  @param h1 word to be compared.
 *  @param h2 word to be compared
 *  @return zero if both are equal, a positive number if w1 is greater or a
 *     negative number if w2 is greater.
 */
int HComparer(Hit h1, Hit h2){
	if(h1.seqX > h2.seqX) return 1;
	else if(h1.seqX < h2.seqX) return -1;

	if(h1.seqY > h2.seqY) return 1;
	else if(h1.seqY < h2.seqY) return -1;

	if(h1.diag > h2.diag) return 1;
	else if(h1.diag < h2.diag) return -1;

	if(h1.posX > h2.posX) return 1;
	else if(h1.posX < h2.posX) return -1;

	if(h1.length > h2.length) return 1;
	else if(h1.length < h2.length) return -1;

	return 0;
}


/* This function is used to load a hashentry from genome dictionary and
 * store it in a wordEntry.
 *  @param we word entry where info will be stored.
 *  @param wD genome dictionary file.
 *  @return zero if everything finished well or a negative number in other cases.
 */
void readHashEntry(WordEntry *we,FILE *wD){
	hashentry h;
	// Read hashentry
	if(fread(&h,sizeof(hashentry),1,wD)!=1){
		if(!feof(wD))
			fprintf(stderr, "Couldn't read hashentry.\n");
		return; // End process
	}

	// Store info
	we->metag = false;
	// we->WB <- It is written out of the function 
	we->pos = h.pos;
	we->reps = (uint32_t) h.num;
	int i;
	for(i=0;i<8;++i)
		we->seq[i] = h.w.b[i];
}


/* This function is used to load a hashentry from metagenome dictionary and
 * store it in a wordEntry.
 *  @param we word entry where info will be stored.
 *  @param wD metagenome dictionary file.
 *  @return zero if everything finished well or a negative number in other cases.
 */
int readWordEntrance(WordEntry *we,FILE *wD,uint16_t SeqBytes){
	we->metag = true;
	// Read sequence
	uint16_t i;
	for(i=0;i<SeqBytes;++i)
			if(fread(&we->seq[i],sizeof(unsigned char),1,wD)!=1){
				if(feof(wD)){
					return 1;
				}
				fprintf(stderr, "readWordEntrance:: Error loading sequence\n");
				return -1;
			}
	// Read position
	if(fread(&we->pos,sizeof(uint64_t),1,wD)!=1){
		fprintf(stderr, "readWordEntrance:: Error loading position.\n");
		return -1;
	}
	// Read repetitions
	if(fread(&we->reps,sizeof(uint32_t),1,wD)!=1){
		fprintf(stderr, "readWordEntrance:: Error loading repetitions.\n");
		return -1;
	}
	return 0;
}


/* This function is used to generate hits from two WordEntry coincidents.
 *  @param buff buffer where store hits.
 *  @param X word entry coincident with Y.
 *  @param Y word entry coincident with X.
 *  @param XPFile X sequence locations dictionary file.
 *  @param YPFile Y sequence locations dictionary file.
 *  @param outIndx intermediate file necessary if buffer get filled.
 *  @param outBuff intermediate file necessary if buffer get filled.
 *  @param hitsInBuff words stored in buffer.
 *  @retun a non-negative number if the process finished without errors or a
 *     a negative number in other cases.
 */
int generateHits(Hit* buff,WordEntry X,WordEntry Y,FILE* XPFile,FILE* YPFile,FILE* outIndx, FILE* outBuff, uint64_t* hitsInBuff){
	// Positionate on locations files
	if(fseek(XPFile,X.pos,SEEK_SET)!=0){
		fprintf(stderr, "generateHits:: Error positioning on X file.\n");
		return -1;
	}
	if(fseek(YPFile,Y.pos,SEEK_SET)!=0){
		fprintf(stderr, "generateHits:: Error positioning on Y file.\n");
		return -1;
	}

	// Prepare necessary variables
	uint64_t hitLength = (uint64_t)(X.WB*4);
	LocationEntry X_Arr[X.reps];
	LocationEntry Y_Arr[Y.reps];

	// Load entrances
	loadLocationEntrance(&X_Arr[0],XPFile,X.reps,X.metag);
	loadLocationEntrance(&Y_Arr[0],YPFile,Y.reps,Y.metag);

	// Check buffer space
	if(*hitsInBuff == MAX_BUFF){
		writeHitsBuff(buff,outIndx,outBuff,*hitsInBuff);
		*hitsInBuff=0;
	}

	// Generate all hits
	uint32_t i,j;
	for(i=0; i<X.reps; ++i)
		for(j=0; j<Y.reps; ++j){
			storeHit(&buff[*hitsInBuff],X_Arr[i],Y_Arr[j],hitLength);
			*hitsInBuff+=1;
			// Check buffer space
			if(*hitsInBuff == MAX_BUFF){
				writeHitsBuff(buff,outIndx,outBuff,*hitsInBuff);
				*hitsInBuff=0;
			}			
		}
	return 0;
}


/* This function is used to load a set of locations entry from a location dictionary file.
 *  @param arr array where locations will be stored.
 *  @param PFile from load the locations.
 *  @param reps number of locations to be loaded.
 *  @param metagenome if it's true PFile must be a new format file, else will be read as 
 *     an old format genome dictionary file.
 */
void loadLocationEntrance(LocationEntry* arr, FILE* PFile, uint32_t reps, bool metagenome){
	uint32_t i;
	if(metagenome){
		for(i=0; i<reps;++i){
			fread(&arr[i].seq,sizeof(uint32_t),1,PFile);
			fread(&arr[i].pos,sizeof(uint64_t),1,PFile);
		}
	}else{
		location aux;
		for(i=0; i<reps;++i){
			fread(&aux,sizeof(location),1,PFile);
			arr[i].seq = (uint32_t) aux.seq;
			arr[i].pos = aux.pos;
		}
	}
}


/* This function is used to load necessary info in a Hit variable.
 *  @param hit where info will be stored.
 *  @param X location of hit.
 *  @param Y location of hit.
 *  @param HitLength matched sequence legnth.
 */
inline void storeHit(Hit* hit,LocationEntry X,LocationEntry Y,uint64_t HitLength){
	hit->diag = X.pos - Y.pos;
	hit->posX = X.pos;
	hit->seqX = X.seq;
	hit->posY = Y.pos;
	hit->seqY = Y.seq;
	hit->length = HitLength;
}


/* This function is used to write a buffer in the intermediate files. The order
 * in each dictioanry is:
 *   - Index: each entrance: Pos<uint64_t> HitsInBuff<uint64_t>
 *   - Hits: each entrance: X<uint32_t> Y<uint32_t> Diag<uint64_t> PosX<uint64_t> PosY<uint64_t> Length<uint64_t>
 *  @param buff buffer t be written.
 *  @param index intermediate file.
 *  @param hits intermediate file.
 *  @param hitsInBuff number of words stored on buffer.
 */
void writeHitsBuff(Hit* buff,FILE* index,FILE* hits,uint64_t hitsInBuff){
	// Sort buffer
	quicksort_H(buff,0,hitsInBuff-1);

	// Write info on index file
	uint64_t pos = (uint64_t) ftell(hits);
	uint64_t numHits = 0;
	Hit lastHit;
	fwrite(&pos,sizeof(uint64_t),1,index);
	
	// Write first hit
	fwrite(&buff[0].seqX,sizeof(uint32_t),1,hits);
	fwrite(&buff[0].seqY,sizeof(uint32_t),1,hits);
	fwrite(&buff[0].diag,sizeof(int64_t),1,hits);
	fwrite(&buff[0].posX,sizeof(uint64_t),1,hits);
	fwrite(&buff[0].posY,sizeof(uint64_t),1,hits);
	fwrite(&buff[0].length,sizeof(uint64_t),1,hits);
	numHits++;
	lastHit = buff[0];
		
	// Write hits in hits file
	for(pos=1; pos<hitsInBuff; ++pos){
		if(buff[pos].seqX == lastHit.seqX && buff[pos].seqY == lastHit.seqY && buff[pos].diag == lastHit.diag && buff[pos].posX < lastHit.posX + lastHit.length){
			continue; // Collapsable
		}
		fwrite(&buff[pos].seqX,sizeof(uint32_t),1,hits);
		fwrite(&buff[pos].seqY,sizeof(uint32_t),1,hits);
		fwrite(&buff[pos].diag,sizeof(int64_t),1,hits);
		fwrite(&buff[pos].posX,sizeof(uint64_t),1,hits);
		fwrite(&buff[pos].posY,sizeof(uint64_t),1,hits);
		fwrite(&buff[pos].length,sizeof(uint64_t),1,hits);
		lastHit = buff[pos];
		numHits++;
	}
	// Write final number of hits
	fwrite(&numHits,sizeof(uint64_t),1,index);
	buffersWritten++;
}


/* Function used to compare two Hit variables. It function sort with this
 * criterion:
 *    1 - Sequence X
 *    2 - Sequence Y
 *    3 - Diagonal
 *    4 - Position on X
 *    5 - Length 
 *  @param h1 word to be compared.
 *  @param h2 word to be compared
 *  @return zero if w2 are greater or equal and a positive number if
 *     w1 is greater.
 */
int GT(Hit h1, Hit h2){
	if(h1.seqX > h2.seqX) return 1;
	else if(h1.seqX < h2.seqX) return 0;

	if(h1.seqY > h2.seqY) return 1;
	else if(h1.seqY < h2.seqY) return 0;

	if(h1.diag > h2.diag) return 1;
	else if(h1.diag < h2.diag) return 0;

	if(h1.posX > h2.posX) return 1;
	else if(h1.posX < h2.posX) return 0;

	if(h1.length > h2.length) return 1;
	return 0;
}


/* This function is necessary for quicksort functionality.
 *  @param arr array to be sorted.
 *  @param left inde of the sub-array.
 *  @param right index of the sub-array.
 */
int partition(Hit* arr, int left, int right){
   int i = left;
   int j = right + 1;
   Hit t;

   // Pivot variable
   int pivot = (left+right)/2;

   if(GT(arr[pivot],arr[right]))
		 SWAP_H(&arr[pivot],&arr[right],t);

   if(GT(arr[pivot],arr[left]))
		 SWAP_H(&arr[pivot],&arr[left],t);

   if(GT(arr[left],arr[right]))
		 SWAP_H(&arr[left],&arr[right],t);

	while(1){
		do{
			++i;
		}while(!GT(arr[i],arr[left]) && i <= right);

		do{
			--j;
		}while(GT(arr[j],arr[left]) && j >= left);

		if(i >= j) break;

		SWAP_H(&arr[i],&arr[j],t);
	}

	SWAP_H(&arr[left],&arr[j],t);

	return j;
}


/* This function is used to sort a Hit array.
 *  @param arr array to be sorted.
 *  @param left index where start to sort.
 *  @param right index where end sorting action.
 *
 */
void quicksort_H(Hit* arr, int left,int right){
	int j;

	if(left < right){
		// divide and conquer
		j = partition(arr,left,right);
		quicksort_H(arr,left,j-1);
		quicksort_H(arr,j+1,right);
   }
}


/* This function is used to load a hit from a hits intermediate file.
 *  @param hit varaible where loaded hit will be stored.
 *  @param hFile pointer to hits intermediate file.
 *  @param unread rest of hits on intermediate file.
 *  @return Number of hits read from intermediate file.
 */
uint64_t loadHit(Hit **hit,FILE* hFile, int64_t unread){
	uint64_t i,j;
	for(j=0; j<READ_BUFF_LENGTH && unread > 0; ++j){
		fread(&(*hit)[j].seqX,sizeof(uint32_t),1,hFile);
		fread(&(*hit)[j].seqY,sizeof(uint32_t),1,hFile);
		fread(&(*hit)[j].diag,sizeof(int64_t),1,hFile);
		fread(&(*hit)[j].posX,sizeof(uint64_t),1,hFile);
		fread(&(*hit)[j].posY,sizeof(uint64_t),1,hFile);
		fread(&(*hit)[j].length,sizeof(uint64_t),1,hFile);
		unread--;
	}
	return j;
}


/* This function is used to write a fragment in fragment file. The order 
 * of a fragment entrance is:
 *    X<uint32_t> Y<uint32_t> diag<int64_t> xStart<uint64_t> yStart<uint64_t> xEnd<uint64_t> 
 *        yEnd<uint64_t> length<uint64_t> ident<uint64_t> score<uint64_t> similarity<float> 
 *        block<int64_t> strand<char>
 *  @param frag fragment to be written.
 *  @param fr fragment file where fragment will be written.
 */
inline void writeFragment(FragFile frag,FILE *fr){
	fwrite(&frag.seqX,sizeof(uint32_t),1,fr);
	fwrite(&frag.seqY,sizeof(uint32_t),1,fr);
	fwrite(&frag.diag,sizeof(int64_t),1,fr);
	fwrite(&frag.xStart,sizeof(uint64_t),1,fr);
	fwrite(&frag.yStart,sizeof(uint64_t),1,fr);
	fwrite(&frag.xEnd,sizeof(uint64_t),1,fr);
	fwrite(&frag.yEnd,sizeof(uint64_t),1,fr);
	fwrite(&frag.length,sizeof(uint64_t),1,fr);
	fwrite(&frag.ident,sizeof(uint64_t),1,fr);
	fwrite(&frag.score,sizeof(uint64_t),1,fr);
	fwrite(&frag.similarity,sizeof(float),1,fr);
	fwrite(&frag.block,sizeof(int64_t),1,fr);
	fwrite(&frag.strand,sizeof(char),1,fr);
}


/* This function is used to swap two hits variables
 *  @param h1 hit to be swapped.
 *  @param h2 hit to be swapped.
 *  @param t auxiliar hit.
 */
void SWAP_H(Hit* h1, Hit* h2, Hit t){
	copyHit(&t,*h1);
	copyHit(h1,*h2);
	copyHit(h2,t);
}


/* This function is used to copy a hit variable in another hit variable.
 *  @param toCopy where hit will be copied.
 *  @param copy hit to be copied.
 */
void copyHit(Hit* toCopy, Hit copy){
	toCopy->diag = copy.diag;
	toCopy->posX = copy.posX;
	toCopy->posY = copy.posY;
	toCopy->seqX = copy.seqX;
	toCopy->seqY = copy.seqY;
	toCopy->length = copy.length;
}


/* This method push node B after A (A->C ==PUSH==> A->B->C)
 *  @param A node after B will be pushed.
 *  @param B node to be pushed.
 */
void push(node **A,node **B){
	(*B)->next = (*A)->next;
	(*A)->next = *B;
}


/* Move node after B to after A position and make linked list consistent.
 *  @param A reference node.
 *  @param B node after it will be moved.
 */
void move(node **A,node **B){
	node *temp = (*B)->next->next;
	push(A,&(*B)->next);
	(*B)->next = temp;
}


/* This emthod sort a linked list
 *  @param first node of the linked list.
 */
void sortList(node **first){
	if((*first)->next == NULL) return; // Linked list with only one element

	node *current = *first;
	node *aux;
	bool sorted = false;
	// Do until end
	while(!sorted){
		if(current->next == NULL) sorted = true;
		else if(GT(current->next->hits[current->index],current->hits[current->next->index])==0){ // Next is smaller
			// Search position
			if(GT(current->next->hits[current->next->index],(*first)->hits[(*first)->index])==0){ // New first node
				aux = current->next->next;
				current->next->next = *first;
				*first = current->next;
				current->next = aux;
			}else{ // Search position
				aux = *first;			
				while(GT(aux->next->hits[aux->next->index],current->next->hits[current->next->index])==1)
					aux = aux->next;
				move(&aux,&current);
				// Chekc if it's the last node
				if(current->next == NULL) sorted = true;
			}
		}else{ // Go next
			current = current->next;
			if(current->next == NULL){ // End of the list
				// Search position
				if(GT(current->next->hits[current->next->index],(*first)->hits[(*first)->index])==0){ // New first node
					aux = current->next->next;
					current->next->next = *first;
					*first = current->next;
					current->next = aux;
				}else{
					aux = *first;			
					while(aux->next == NULL && GT(current->next->hits[current->next->index],aux->next->hits[aux->next->index])==1)
						aux = aux->next;
					move(&aux,&current);
				}
				// List sorted
				sorted = true;
			}
		}
	}
}


/* This function is used to check the correc order of the first node of a linked list.
 * If it's incorrect, this function sort it.
 *  @param list linked list to be checked.
 *  @param discardFirst a boolean value that indicate if first node should be deleted.
 */
void checkOrder(node** list,bool discardFirst){
	node *aux;
	if(discardFirst){
		aux = *list;
		*list = (*list)->next;
		free(aux);
	}else if((*list)->next != NULL){ // Check new position
		// Search new position
		if(GT((*list)->hits[(*list)->index],(*list)->next->hits[(*list)->next->index])==1){
			node *curr = (*list)->next;
			while(1){
				if(curr->next == NULL) break; // End of list
				else if(GT((*list)->hits[(*list)->index],curr->next->hits[curr->next->index])==0) break; // position found
				else curr = curr->next;
			}
			aux = (*list)->next;
			(*list)->next = curr->next;
			curr->next = *list;
			*list = aux;
		}
	}
}