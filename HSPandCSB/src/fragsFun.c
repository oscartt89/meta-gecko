/* @file fragsFun.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description this file contains the functions necesary for
 * 		correct operation of frags.c code. 
 */

#include "frags.h"

/* This function takes all genome dictionaries from a directory given
 * @param genomeSetPath path to genome dictionaries folder.
 * @param genomes is a dictionaryG array pointer where detected dictionaries
 *    will be stored. Memory on this array will be deallocated before start
 *    to store the new values.
 * @return the number of dictionaries detected or a negative number if
 *    something got wrong.
 */
int readGenomeSet(char* genomeSetPath,dictionaryG** genomes){
	// Variables
	DIR *gFolder;
	struct dirent *ent;
	char *wD, *pD;
  char file[MAX_NAME_L];
	int numGenomes = 0, currentMax = MAX_GENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  free(*genomes);
  if((*genomes = (dictionaryG*) malloc(sizeof(dictionaryG)*MAX_GENOME_SET))==NULL){
    fprintf(stderr, "Error allocating memory for genome dictionary set.\n");
    return -1;
  }

	// Open genomes folder
	if((gFolder = opendir(genomeSetPath))==NULL){
		fprintf(stderr, "Error opening genomes folder.\n");
		return -1;
  }

	// Take genome dictionary files
	while((ent = readdir(gFolder))!=NULL){
    // Check for memory
    if(numGenomes>=currentMax){
      if((*genomes = realloc(*genomes,sizeof(dictionaryG)*(currentMax + MAX_GENOME_SET)))==NULL){
        fprintf(stderr, "Error reallicatin memory for genome dictionary set.\n");
        return -1;
      }
      currentMax += MAX_GENOME_SET;
    }
		// Files are sorted alphabetically
    // Check it's not a metagenome dictionary
    if(strstr(ent->d_name,".metag.d2h")!=NULL) continue;
		// Should appear first d2hP than d2hW
		if(strstr(ent->d_name,".d2hP")!=NULL){ // New dictionary
			// Save name
			memcpy(genomes[numGenomes]->name,ent->d_name,strlen(ent->d_name)-5);
			// Save location dictionary
      strcpy(&file[0],genomeSetPath);
      strcat(&file[0],ent->d_name);
			strcpy(genomes[numGenomes]->P,&file[0]);
			//Next file should be d2hW dictionary
			if((ent = readdir(gFolder))==NULL){
				fprintf(stderr, "Error: incomplete genome pair dictionary. End of file list.\n");
				return -1;
			}
			if(strstr(ent->d_name,".d2hW")!=NULL && strstr(ent->d_name,".metag.d2h")==NULL){
				// Save word dictionary
        strcpy(&file[0],genomeSetPath);
        strcat(&file[0],ent->d_name);
				strcpy(genomes[numGenomes]->W,&file[0]);
				numGenomes++;
			}else if(strstr(ent->d_name,".metag.d2h")!=NULL){
        fprintf(stderr, "Error: it's a metagenome dictionary.\n");
        return -1;
      }else{
				fprintf(stderr, "Error: incomplete genome pair dictionary.\n");
				return -1;
			}
		}
	}

  // Close dir
  closedir(gFolder);

  return numGenomes;
}


/* This function takes all metagenome dictionaries from a directory given.
 * @param metagSetPath path to metagenome dictionaries folder.
 * @param metagenomes is a dictionaryM array pointer where detected dictionaries
 *    will be stored. Memory on this array will be deallocated before start
 *    to store the new values.
 * @return the number of dictionaries detected or a negative number if
 *    something got wrong.
 */
int readMetagenomeSet(char* metagSetPath,dictionaryM** metagenomes){
  // Variables
  DIR *mFolder;
  struct dirent *ent;
  char *wD, *pD, *rD;
  char file[MAX_NAME_L];
  int numMetags = 0, currentMax = MAX_METAGENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  free(*metagenomes);
  if((*metagenomes = (dictionaryM*) malloc(sizeof(dictionaryM)*MAX_METAGENOME_SET))==NULL){
    fprintf(stderr, "Error allocating memory for metagenome dictionary set.\n");
    return -1;
  }

  // Open metagenomes folder
  if((mFolder = opendir(metagSetPath))==NULL){
    fprintf(stderr, "Error opening metagenomes folder.\n");
    return -1;
  }

  // Take metagenome dictionary files
  while((ent = readdir(mFolder))!=NULL){
    // Check for memory
    if(numMetags>=currentMax){
      if((*metagenomes = realloc(*metagenomes,sizeof(dictionaryM)*(currentMax + MAX_METAGENOME_SET)))==NULL){
        fprintf(stderr, "Error reallocatin memory for metagenome dictionary set.\n");
        return -1;
      }
      currentMax += MAX_METAGENOME_SET;
    }
    // Files are sorted alphabetically
    // Should appear first d2hP, then d2hR and d2hW
    if(strstr(ent->d_name,".metag.d2hP")!=NULL){ // New dictionary
      // Save name
      memcpy(&metagenomes[numMetags]->name[0],ent->d_name,strlen(ent->d_name)-11);
      // Save location dictionary
      strcpy(&file[0],metagSetPath);
      strcat(&file[0],ent->d_name);
      strcpy(&metagenomes[numMetags]->P[0],&file[0]);
      //Next file should be d2hR dictionary
      if((ent = readdir(mFolder))==NULL){
        fprintf(stderr, "Error: incomplete metagenome triplet dictionary. End of file list.\n");
        return -1;
      }
      if(strstr(ent->d_name,".metag.d2hR")!=NULL){
        // Save read dictionary
        strcpy(&file[0],metagSetPath);
        strcat(&file[0],ent->d_name);
        strcpy(&metagenomes[numMetags]->R[0],&file[0]);
        // Now should appear words dictionary
        if((ent = readdir(mFolder))==NULL){
          fprintf(stderr, "Error: incomplete metagenome triplet dictionary. End of file list.\n");
          return -1;
        }
        if(strstr(ent->d_name,".metag.d2hW")!=NULL){
          // Save words dictionary
          strcpy(&file[0],metagSetPath);
          strcat(&file[0],ent->d_name);
          strcpy(&metagenomes[numMetags]->W[0],&file[0]);
          numMetags++; 
        }else{
          fprintf(stderr, "Error: incomplete metagenome triple dictionary. Word dictionary not found.\n");
          return -1;  
        }
      }else{
        fprintf(stderr, "Error: incomplete metagenome triple dictionary. Read dictionary not found.\n");
        return -1;
      }
    }
  }

  // Close dir
  closedir(mFolder);

  return numMetags;
}


/* This function takes a read from read dictionary given and
 * load all the words of this read in a wentry array given too.
 *  @param dR: reads dictionary.
 *  @param dW: words dictionary.
 *  @param dP: positions dictionary.
 *  @param kmers: array of wentry where words will be allocated.
 *  @param WL: word length.
 *  @return: number wentry instances allocated on kmers.
 * WANING: kmers memory are allocated inside of this function.
 */
int loadRead(FILE *dR,FILE *dW,FILE *dP,wentry** kmers,int WL){
  // Variables
  READ r;
  hashentry he;
  int numKmers = 0;

  // Load read
  fread(&r,sizeof(READ),1,dR); 

  // KMERS space
  if((*kmers = (wentry*) malloc(sizeof(wentry)*MAX_WORDS))==NULL){
    fprintf(stderr , "Error allocating space for metagenome KMERS.\n");
    return -1;
  }

  // Positionate on word dictionary
  fseek(dW,r.pos,SEEK_SET);

  // Load kmers in wentry array
  int i,j;
  uint64_t pos;

  fprintf(stderr, "R.num:%d\n",r.num);
  for(i=0;i<r.num;++i){
    if(feof(dW)){
      fprintf(stderr, "Error reading hashentry. Premature end of file.\n");
      return -1;
    }
    // Read hashentry
    fread(&he,sizeof(hashentry),1,dW);

    // Positionate on positions dictionary
    fseek(dP,he.pos,SEEK_SET);

    fprintf(stderr, "H.num:%d\n",he.num);
    // Take locations
    for(j=0;j<he.num;++j){
      if(feof(dW)){
        fprintf(stderr, "Error reading position. Premature end of file.\n");
        return -1;
      }
      fprintf(stderr, "TEST\n");
      // Store sequence
      memcpy(&kmers[numKmers]->w.b[0],&he.w.b[0],WL); // ## Programa se para en ejecución al cabo de unas iteraciones 
      //Store seq index
      kmers[numKmers]->seq = r.readIndex;
      // Take pos
      fread(&pos,sizeof(uint64_t),1,dP);
      //Store pos
      kmers[numKmers]->pos = pos;
      // Update num kmers
      numKmers++;
    }
  }

  return numKmers;
}


/* This function all words from a genome dictionary and load
 * it in a wentry array given.
 *  @param genome: genome dictionary structure.
 *  @param kmers: array of wentry where words will be allocated.
 *  @param WL: word length.
 *  @return: number wentry instances allocated on kmers.
 * WANING: kmers memory are allocated inside of this function.
 */
int loadGenome(dictionaryG genome,wentry** kmers,int WL){
  // Variables
  hashentry he;
  location lo;
  int numKmers = 0, currentMax = MAX_WORDS;
  FILE *dW, *dP;

  // KMERS space
  if((*kmers = (wentry*) malloc(sizeof(wentry)*MAX_WORDS))==NULL){
    fprintf(stderr , "Error allocating space for metagenome KMERS.\n");
    return -1;
  }
  // Open dictionaries
  if((dW = fopen(&genome.W[0],"rb"))==NULL){
    fprintf(stderr, "Error opening genome word dictionary. [%s]\n",&genome.W[0]);
    return -1;
  }
  if((dP = fopen(&genome.P[0],"rb"))==NULL){
    fprintf(stderr, "Error opening genome posotion dictionary. [%s]\n",&genome.P[0]);
    return -1;
  }

  // Read first word
  fread(&he,sizeof(hashentry),1,dW);

  // Read words
  int i;
  while(!feof(dW)){
    // Num of kmers exceeded
    if(numKmers >= currentMax){
      if((*kmers = realloc(*kmers,sizeof(wentry)*(currentMax + MAX_WORDS)))==NULL){
        fprintf(stderr, "Error reallocating memory for genome wentry set.\n");
        return -1;
      }
      currentMax += MAX_WORDS;
    }

    for(i=0;i<he.num;++i){
      // Store sequence
      memcpy(&kmers[numKmers]->w.b[0],&he.w.b[0],WL);
      // Read location
      fread(&lo,sizeof(location),1,dP);
      // Store sequence index
      kmers[numKmers]->seq = lo.seq;
      // Store position
      kmers[numKmers]->pos = lo.pos;
      // Update num of kmers
      numKmers++;
    }
    fread(&he,sizeof(hashentry),1,dW);    
  }

  return numKmers;
}


/*
 */
int hits(wentry* w1,wentry* w2,hit** hits,uint64_t numW1, uint64_t numW2,int WL){
  // Variables
  int numHits = 0;

  // Memory for hits
  if((*hits = (hit*) malloc(sizeof(hit)*numW1*numW2))){
    fprintf(stderr, "Error allocating memory for hits array.\n");
    return -1;
  }

  // Compare all reads
  int i,j;
  for(i=0;i<numW1;++i)
    for(j=0;j<numW2;++j)
      if(wordcmp(&w1[i].w.b[0],&w2[j].w.b[0],WL)==0){
        // Store position
        hits[numHits]->start1 = w1[i].pos;
        hits[numHits]->start2 = w2[j].pos;
        // Store length
        hits[numHits]->length = (WL * 4)/ BITS_NUCLEOTIDE; 
        // Store sequences indexes
        hits[numHits]->seq1 = w1[i].seq;
        hits[numHits]->seq2 = w2[j].seq;
        // Update number of hits
        numHits++;
      }

  return numHits;
}


/*
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n) {
  int i;
  for (i=0;i<n;i++) {
    if (w1[i]<w2[i]) return -1;
    if (w1[i]>w2[i]) return +1;
  }
  return 0;
}


/*  
 */
int quickSort(hit* hits, int left, int right){
  int j;

  if(left < right){
    // divide and conquer
    if((j = partition(hits, left, right))<0) return -1;
    quickSort(hits, left, j-1);
    quickSort(hits, j+1, right);
  }
  return 0;
}


/*
 */
int partition(hit* hits, int left, int right){
  int i = left;
  int j = right+1;
  hit *t;

  if((t = (hit*) malloc(sizeof(hit)))==NULL){
    fprintf(stderr, "Error allocating memory for auxiliar variable.\n");
    return -1;
  }

  // left sera el pivote
  // y contendra la mediana de left, right y (left+right)/2
  int mid = (int) (left+right)/2;

  if(hitComparator(&hits[mid],&hits[right]))
    SWAP(&hits[mid],&hits[right],t);

  if(hitComparator(&hits[mid],&hits[left]))
    SWAP(&hits[mid],&hits[left],t);

  if(hitComparator(&hits[left],&hits[right]))
    SWAP(&hits[left],&hits[right],t);

  while(1){
    do{
      ++i;
    }while(!hitComparator(&hits[i],&hits[left]) && i <= right);

    do{
      --j;
    }while(hitComparator(&hits[j],&hits[left]) && j >= left);

    if( i >= j ) break;

    SWAP(&hits[i],&hits[j],t);
  }

  SWAP(&hits[left],&hits[j],t);

  free(t);

  return j;
}


/* This function is used to compare two hit instances. The criterion
 * used is:
 *    1 - Compare sequence 1 index.
 *    2 - Compare sequence 2 index.
 *    3 - Compare sequence 1 start.
 *    4 - Compare sequence 2 start.
 *    2 - Compare length.
 * @param h1 hit to be compared.
 * @param h2 hit to be compared.
 * @return a positive number if h1 is greater than h2, a negative number
 *    if h2 is greater than h1 and zero if both are equal.
 */
int hitComparator(hit* h1,hit* h2){
  if(h1->seq1 > h2->seq1) return 1;
  else if(h1->seq1 < h2->seq1) return -1;

  if(h1->seq2 > h2->seq2) return 1;
  else if(h1->seq2 < h2->seq2) return -1;

  if(h1->start1 > h2->start1) return 1;
  else if(h1->start1 < h2->start1) return -1;

  if(h1->start2 > h2->start2) return 1;
  else if(h1->start2 < h2->start2) return -1;

  if(h1->length > h2->length) return 1;
  else if(h1->length < h2->length) return -1;

  return 0;
}


/*
 */
inline void SWAP(hit *h1,hit *h2,hit *t){
  memcpy(t,h1,sizeof(hit));
  memcpy(h1,h2,sizeof(hit));
  memcpy(h2,t,sizeof(hit));
}


/*
 */
int calculateFragments(hit* hits, FragFile** frags, int numHits, int SThreshold, int minLength){
  // Variables
  int numFragments = 0;
  
  // Take space for frags
  if((*frags = (FragFile*) malloc(sizeof(FragFile)*numHits))==NULL){
    fprintf(stderr, "Error allocating space for fragments array.\n");
    return -1;
  }

  // Calculate fragments
    // First instance
    frags[numFragments]->xStart = hits[0].start1;
    frags[numFragments]->yStart = hits[0].start2; 
    frags[numFragments]->seqX = hits[0].seq1;
    frags[numFragments]->seqY = hits[0].seq2;
    frags[numFragments]->length = hits[0].length;
    frags[numFragments]->similarity = 100;
    frags[numFragments]->diag = hits[0].start1 - hits[0].start2;
    frags[numFragments]->ident = hits[0].length;


  int i;
  uint64_t fragEnd,dist,oldLength,newLength;
  uint8_t oldS, newS;
  for(i=1; i<numHits; i++){
    // If hits are equal discard one
    if(hitComparator(&hits[i],&hits[i-1])==0) continue;

    // Sequence change
    if(hits[i].seq1 != frags[numFragments]->seqX | hits[i].seq2 != frags[numFragments]->seqY){
      if(frags[numFragments]->length >= minLength)
        numFragments++;
      // Init new fragment
      frags[numFragments]->xStart = hits[i].start1;
      frags[numFragments]->yStart = hits[i].start2; 
      frags[numFragments]->seqX = hits[i].seq1;
      frags[numFragments]->seqY = hits[i].seq2;
      frags[numFragments]->length = hits[i].length;
      frags[numFragments]->similarity = 100;
      frags[numFragments]->diag = hits[i].start1 - hits[i].start2;
      frags[numFragments]->ident = hits[0].length;
      // Go to next hit
      continue;
    }

    // IF same sequences -> same diagonal
    fragEnd = frags[numFragments]->xStart + frags[numFragments]->length;
    newLength = hits[i].start1 + hits[i].length - frags[numFragments]->xStart;
    oldLength = frags[numFragments]->length;
    oldS = frags[numFragments]->similarity;
    dist = fragEnd - hits[i].start1;
      dist = dist < 0 ? 0 : dist;
    newS = (newLength - (dist + oldLength - oldLength*oldS))/newLength * 100;

    if(newS >= SThreshold){ // Correct fragment, update
      frags[numFragments]->length = newLength;
      frags[numFragments]->similarity = newS;
      frags[numFragments]->ident += newLength - oldLength - dist;
    }else{ // Not enough quality, create new framgent
      if(frags[numFragments]->length >= minLength)
        numFragments++;
      //else overwrite actual fragment (invalid length)
      frags[numFragments]->xStart = hits[i].start1;
      frags[numFragments]->yStart = hits[i].start2;
      frags[numFragments]->diag = hits[i].start1 - hits[i].start2;
      frags[numFragments]->length = hits[i].length;
      frags[numFragments]->seqX = hits[i].seq1;
      frags[numFragments]->seqY = hits[i].seq2;
      frags[numFragments]->similarity = 100;
      frags[numFragments]->ident = hits[i].length;
    }
  }

  return numFragments;
} 