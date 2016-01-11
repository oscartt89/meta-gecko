/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "dict.h"

/* This function is used to shift bits in a unsigned char array
 *	@param w: word structure where char array to be shifted is stored.
 */
void shift_word(word* w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}


/* This method is used to store a wentry instance into an array of wentry. 
 * @param wArr is the array of wentry** where the new wentry will be stored.
 * @param word is the new word to be stored.
 * @param length is the new length of the array.
 * @param maxStored is a variable used to know if memory allocation is necessary or not.
 * @return 0 if the process ends correctly. If something got wrong return a
 * 		negative value.
 */
int storeWord(wentry *wArr,wentry word,uint64_t length,uint64_t *maxStored){
	// Allocate memory if it's necessary
	if((length > *maxStored) == 0){ // Limit of words
		if((wArr[length-1].w.b = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
			fprintf(stderr, "storeWord:: Error allocating space for word\n");
			return -1;
		}
		*maxStored++;
	}
	// Copy the new word
	memcpy(&wArr[length-1],&word,sizeof(wentry));
	return 0;
}


/* This function is used to write a given set of words in two intermediate files.
 *  @param buff set of words to be written.
 *  @param index file of intermediate files.
 *  @param positions file of intermediate files.
 *  @param numWords number of instances on buff set.
 *  @return a negative number if any error happens or a non negative number if
 *     the process finished correctly.
 */
int writeBuffer(wentry* buff,FILE* index,FILE* positions,uint64_t numWords){
	// Sort buffer
	quickSort_W(buff,0,numWords);
	// Search where start to write
	bool alreadyExists;
	long int writeIn = 0;
	uint64_t aux;
	uint16_t auxNum;
	int i,j;

	for(i=0; i<numWords; i=j){

		// Search word
		writeIn = searchInIndex(buff[i].w,index,writeIn,&alreadyExists);
		// Positioante on file
		fseek(index,writeIn,SEEK_SET);
		// Read FLAG and word if it's necessary
		if(alreadyExists){
			// Read flag
			if(fread(&aux,sizeof(uint64_t),1,index)!=1){
				fprintf(stderr, "writeBuffer:: Error searching FLAG.\n");
				return -1;
			}
			// Search next flag or end of file
			bool control = true;
			while(!feof(index) && control){
				// Read number of repetitions
				if(fread(&auxNum,sizeof(uint16_t),1,index)!=1){
					fprintf(stderr, "writeBuffer:: Error searching rep.num.\n");
					return -1;
				}
				// Read position on positions file
				if(fread(&aux,sizeof(uint64_t),1,index)!=1){
					fprintf(stderr, "writeBuffer:: Error searching position value.\n");
					return -1;	
				}
				// Take current position
				writeIn = ftell(index);
				// Search FLAG
				if(fread(&aux,sizeof(uint64_t),1,index)!=1){ // Could be the end of the file
					if(!feof(index)){
						fprintf(stderr, "writeBuffer:: couldn't found next flag before end of file.\n");
						return -1;
					}else control = false; //end of file
				}else if(aux == (uint64_t)FLAG) // FLAG found
					control = false;				
			}
			// Positionate on write point
			fseek(index,writeIn,SEEK_SET);

		}
		// Start to write
		j = i+1;
		// Search equal words
		while(wordcmp(buff[i].w.b,buff[j].w.b,BYTES_IN_WORD)==0 && j < numWords-1){
			fprintf(stderr, "%d\n",j);
			++j;
		}
		
		fprintf(stderr, "OUT\n");
		// Write words
		writeWords(buff,index,positions,i,j,alreadyExists);
	}

	// Buffer correctly written on intermediate files
	return 0;
}


/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n){
	int i;
	for(i=0;i<n;i++)
		if(w1[i] < w2[i]) return -1;
		else if(w1[i] > w2[i]) return +1;

	return 0;
}


/* This function is used to swap/interchange two wentry instances.
 *  @param w1: wentry that will be swapped.
 *  @param w2: wentry that will be swapped.
 */
inline void SWAP_W(wentry *w1,wentry *w2){
	wentry t;
	if((t.w.b = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
		fprintf(stderr, "SWAP_W:: Error allocating memory.\n");
	}else{
		memcpy(&t,w1,sizeof(wentry));
		memcpy(w1,w2,sizeof(wentry));
		memcpy(w2,&t,sizeof(wentry));
		free(t.w.b);
	}
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequences (alphabetically).
 * 		2 - Compare position on sequence.
 *      3 - Compare Read index.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
	int wComp;
	if((wComp = wordcmp(&w1->w.b[0],&w2->w.b[0],BYTES_IN_WORD)) != 0) return wComp;

	if(w1->pos > w2->pos) return 1;
	else if(w1->pos < w2->pos) return -1;

	if(w1->seq > w2->seq) return 1;
	else if(w1->seq < w2->seq) return -1;

	return 0;
}


/* Function used to sort an array of wentrys.
 *  @param words: array of wentry that will be sorted.
 *  @param start: index where start to sort.
 *  @param length: length of the array to sort (starting on "start" index).
 */
void quickSort_W(wentry *words, uint64_t start, uint64_t length){
	//Check exceptions
  if (length < 2) return;
  else if(length == 2){
  	if(wordComparator(&words[start],&words[start+1])>0)
  		SWAP_W(&words[start],&words[start+1]);
  	return;
  }
  if(start < 0) return;
	uint64_t right, left, pivot, end = start + length - 1;
	int changed = 0;

  pivot = rand() % (end + 1 - start) + start;
  right = start;
  left = end;

  while(right != pivot || left != pivot) {
  	// While right is lower than pivot, continue
      while(wordComparator(&words[right],&words[pivot])<0 && right <= pivot) right++;
      // While left y greater than pivot, continue
      while(wordComparator(&words[pivot],&words[left])<0 && left >= pivot) left--;
      // All array checked -> end
      if(right>=left)
          break;
      // Swap selected values
      SWAP_W(&words[right],&words[left]);
      changed = 1;
      if(right < pivot) right++;
      if(left > pivot) left--;
  }

  if(changed && right<left){ // If anything changed, don't continue, it's sorted
  	quickSort_W(words,start,right);
  	quickSort_W(words,right,length-right);
  }

  return;
}


/* Function used to find the correct place to write the word given.
 *  @param toSearc word to be searched.
 *  @param index file where index are written.
 *  @param start place to start to search on file.
 *  @param alreadyExists boolean value that indicates if the word are already indexed or not.
 *  @return a negative value if some error happens or the position to write the word. If 
 *    alreadyExists are true the position returned is just before the flag of this word, in
 *    other cases the position returned is the place to directly write.
 */
long int searchInIndex(word toSearch,FILE *index,long int start, bool *alreadyExists){
	// Postionate on start point
	if(fseek(index,start,SEEK_SET)!=0){
		fprintf(stderr, "searchInIndex:: Error positioning at START value.\n");
		return -1;
	}
	// Search the correct position
	unsigned char word[BYTES_IN_WORD];
	long int pos = start;
	uint64_t aux; // To read flags and positions
	uint16_t aux2; // To read number of repetitions
	int comp;
	*alreadyExists = false;

	// Start reading flag
	if(fread(&aux,sizeof(uint64_t),1,index)!=1){ // Could be the end of the file
		if(!feof(index)){
			fprintf(stderr, "searchInIndex:: couldn't found flag on START.\n");
			return -1;
		}else return start;
		
	}else if(aux != (uint64_t)FLAG){
		fprintf(stderr, "searchInIndex:: FLAG isn't after START point.\n");
		return -1;
	}

	// Now seek the right place to write the word.
	while(!feof(index)){
		// Read word
		fread(&word,sizeof(unsigned char),BYTES_IN_WORD,index);
		// Are the same word?
		if((comp = memcmp(&word,toSearch.b,BYTES_IN_WORD)) == 0){ // Word found
			*alreadyExists = true;
			return pos;
		}else if(comp < 0) // Word is before than that. This is the correct place
			return pos;
		 // Else -> continue searching
		// Search the next FLAG or end of file
		bool control = true;
		while(!feof(index) && control){
			// Read number of repetitions
			if(fread(&aux2,sizeof(uint16_t),1,index)!=1){
				fprintf(stderr, "searchInIndex:: Error searching rep.num.\n");
				return -1;
			}
			// Read position on positions file
			if(fread(&aux,sizeof(uint64_t),1,index)!=1){
				fprintf(stderr, "searchInIndex:: Error searching position value.\n");
				return -1;	
			}
			// Take current position
			pos = ftell(index);
			// Search FLAG
			if(fread(&aux,sizeof(uint64_t),1,index)!=1){ // Could be the end of the file
				if(!feof(index)){
					fprintf(stderr, "searchInIndex:: couldn't found next flag before end of file.\n");
					return -1;
				}else return pos;
			}else if(aux == (uint64_t)FLAG) // Read next word
				control = false;				
		}
	}

	// File ended and word doesn't found
	return pos;
} 


/* This function is used to write a given set of words in the intermediate files.
 *  @param words set of words.
 *  @param index intermediate file where words are idnexed.
 *  @param positions interediate file where positions are encoded.
 *  @param start index of first word to be written.
 *  @param end index of first word to NOT be written.
 *  @param alreadyExists a boolean value that indicates if the word is already indexed on index file.
 */
inline void writeWords(wentry *words, FILE *index, FILE *positions, int start, int end, bool alreadyExists){
	// Write entrance header if it's necessary
	if(!alreadyExists){
		// Write flag
		uint64_t flag = (uint64_t) FLAG;
		fwrite(&flag,sizeof(uint64_t),1,index);
		// Write word
		fwrite(words[start].w.b,sizeof(unsigned char),BYTES_IN_WORD,index);
	}
	// Update index
	// Write num of repetitions
	uint16_t reps = (uint16_t) end-start;
	fwrite(&reps,sizeof(uint16_t),1,index);
	// Write position on positions file 
	uint64_t pos = (uint64_t)ftell(positions);
	fwrite(&pos,sizeof(uint64_t),1,index);
	// Write words on postions file
	int i;
	for(i=start; i<end; ++i){
		// Write ReadIndex
		fwrite(&words[i].seq,sizeof(uint32_t),1,positions);
		// Write SeqPos
		fwrite(&words[i].pos,sizeof(uint64_t),1,positions);
	}
}


/* This function is used to deallocate a buffer from disk
 */
void freeBuffer(wentry* buff, uint64_t maxStored){
	int k;
	for(k=0; k<maxStored; ++k)
		free(buff[k].w.b);
	free(buff);
}