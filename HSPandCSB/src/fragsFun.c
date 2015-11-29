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
  char file[MAX_FILE_LENGTH];
	int numGenomes = 0, currentMax = MAX_GENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  free(*genomes);
  if((*genomes = (dictionaryG*) malloc(sizeof(dictionaryG)*MAX_GENOME_SET))==NULL){
    fprintf(stderr, "readGenomeSet:: Error allocating memory for genome dictionary set.\n");
    return -1;
  }

	// Open genomes folder
	if((gFolder = opendir(genomeSetPath))==NULL){
		fprintf(stderr, "readGenomeSet:: Error opening genomes folder.\n");
		return -1;
  }

	// Take genome dictionary files
	while((ent = readdir(gFolder))!=NULL){
    // Check for memory
    if(numGenomes>=currentMax){
      if((*genomes = realloc(*genomes,sizeof(dictionaryG)*(currentMax + MAX_GENOME_SET)))==NULL){
        fprintf(stderr, "readGenomeSet:: Error reallocating memory for genome dictionary set.\n");
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
				fprintf(stderr, "readGenomeSet:: Error: incomplete genome pair dictionary. End of file list.\n");
				return -1;
			}
			if(strstr(ent->d_name,".d2hW")!=NULL && strstr(ent->d_name,".metag.d2h")==NULL){
				// Save word dictionary
        strcpy(&file[0],genomeSetPath);
        strcat(&file[0],ent->d_name);
				strcpy(genomes[numGenomes]->W,&file[0]);
				numGenomes++;
			}else if(strstr(ent->d_name,".metag.d2h")!=NULL){
        fprintf(stderr, "readGenomeSet:: Error: it's a metagenome dictionary.\n");
        return -1;
      }else{
				fprintf(stderr, "readGenomeSet:: Error: incomplete genome pair dictionary.\n");
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
  char file[MAX_FILE_LENGTH];
  int numMetags = 0, currentMax = MAX_METAGENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  if((*metagenomes = (dictionaryM*) malloc(sizeof(dictionaryM)*MAX_METAGENOME_SET))==NULL){
    fprintf(stderr, "readMetagenomeSet:: Error allocating memory for metagenome dictionary set.\n");
    return -1;
  }


  // Open metagenomes folder
  if((mFolder = opendir(metagSetPath))==NULL){
    fprintf(stderr, "readMetagenomeSet:: Error opening metagenomes folder.\n");
    return -1;
  }

  // Take metagenome dictionary files
  while((ent = readdir(mFolder))!=NULL){
    // Check for memory
    if(numMetags>=currentMax){
      if((*metagenomes = realloc(*metagenomes,sizeof(dictionaryM)*(currentMax + MAX_METAGENOME_SET)))==NULL){
        fprintf(stderr, "readMetagenomeSet:: Error reallocating memory for metagenome dictionary set.\n");
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
        fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triplet dictionary. End of file list.\n");
        return -1;
      }
      if(strstr(ent->d_name,".metag.d2hR")!=NULL){
        // Save read dictionary
        strcpy(&file[0],metagSetPath);
        strcat(&file[0],ent->d_name);
        strcpy(&metagenomes[numMetags]->R[0],&file[0]);
        // Now should appear words dictionary
        if((ent = readdir(mFolder))==NULL){
          fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triplet dictionary. End of file list.\n");
          return -1;
        }
        if(strstr(ent->d_name,".metag.d2hW")!=NULL){
          // Save words dictionary
          strcpy(&file[0],metagSetPath);
          strcat(&file[0],ent->d_name);
          strcpy(&metagenomes[numMetags]->W[0],&file[0]);
          numMetags++; 
        }else{
          fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triple dictionary. Word dictionary not found.\n");
          return -1;  
        }
      }else{
        fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triple dictionary. Read dictionary not found.\n");
        return -1;
      }
    }
  }

  // Close dir
  closedir(mFolder);

  return numMetags;
}


/* This function takes a read from read dictionary given and
 * load all the words of this read in a HE array given too.
 *  @param dR: reads dictionary.
 *  @param dW: words dictionary.
 *  @param dP: positions dictionary.
 *  @param kmers: array of HE where words will be allocated.
 *  @param WL: word length.
 *  @return: number HE instances allocated on kmers.
 * WANING: kmers memory are allocated inside of this function.
 */
uint64_t loadRead(FILE *dR,FILE *dW,FILE *dP,HE** kmers,int WL){
  // Variables
  Read r;
  hashentry he;
  uint64_t numKmers = 0;

  // Load read
  fread(&r,sizeof(read),1,dR); 

  // KMERS space
  if((*kmers = malloc(sizeof(HE)*MAX_WORDS))==NULL){
    fprintf(stderr , "loadRead:: Error allocating space for metagenome KMERS.\n");
    return -1;
  }

  // Positionate on word dictionary
  fseek(dW,r.pos,SEEK_SET);

  // Load kmers in HE array
  int i,j;
  uint64_t loc;

  for(i=0;i<r.num;++i){
    if(feof(dW)){
      fprintf(stderr, "loadRead:: Error reading hashentry. Premature end of file.\n");
      return -1;
    }
    // Read hashentry
    fread(&he,sizeof(hashentry),1,dW);

    // Positionate on positions dictionary
    fseek(dP,he.pos,SEEK_SET);

    // Store sequence
    memcpy(&(*kmers)[numKmers].w.b[0],&he.w.b[0],(WL*BITS_NUCLEOTIDE)/8);
    // Store seq index
    (*kmers)[numKmers].seq = r.readIndex;
    // Store num instances
    (*kmers)[numKmers].num = he.num;
    // Locations array space
    if(((*kmers)[numKmers].locations = (uint64_t*) malloc(sizeof(uint64_t)*he.num))==NULL){
      fprintf(stderr , "loadRead:: Error allocating space for metagenome KMER locations.(%i)\n",he.num);
      return -1;
    }

    // Take locations
    int read;
    if(feof(dW)){
      fprintf(stderr, "loadRead:: Error reading position. Premature end of file.\n");
      return -1;
    }
    // Take pos
    if((read=fread((*kmers)[numKmers].locations,sizeof(uint64_t),he.num,dP)) != he.num){
      fprintf(stderr, "loadRead:: Error reading locations (%i/%i)\n", read, he.num);
    }
    // Update num kmers
    numKmers++;
  }

  return numKmers;
}


/* This function read all words from a genome dictionary and load
 * it in a wHE array given.
 *  @param genome: genome dictionary structure.
 *  @param kmers: array of HE where hashentry read will be allocated.
 *  @param locations: array of array of integers where locations will be allocated.
 *  @param WL: word length.
 *  @return: number HE instances allocated on kmers.
 * WANING: kmers memory are allocated inside of this function.
 */
uint64_t loadGenome(dictionaryG genome,HE** kmers,int WL){
  // Variables
  hashentryOld he;
  location lo;
  uint64_t numKmers = 0, currentMax = MAX_WORDS;
  FILE *dW, *dP;

  // KMERS space
  if((*kmers = (HE*) malloc(sizeof(HE)*MAX_WORDS))==NULL){
    fprintf(stderr , "loadGenome:: Error allocating space for metagenome KMERS.\n");
    return -1;
  }

  // Open dictionaries
  if((dW = fopen(&genome.W[0],"rb"))==NULL){
    fprintf(stderr, "loadGenome:: Error opening genome word dictionary. [%s]\n",&genome.W[0]);
    return -1;
  }
  if((dP = fopen(&genome.P[0],"rb"))==NULL){
    fprintf(stderr, "loadGenome:: Error opening genome posotion dictionary. [%s]\n",&genome.P[0]);
    return -1;
  }

  // Read first word
  fread(&he,sizeof(hashentry),1,dW);

  // Read words
  int i;
  while(!feof(dW)){
    // Num of kmers exceeded
    if(numKmers >= currentMax){
      if((*kmers = realloc(*kmers,sizeof(hashentry)*(currentMax + MAX_WORDS)))==NULL){
        fprintf(stderr, "loadGenome:: Error reallocating memory for genome hashentry set.\n");
        return -1;
      }
      currentMax += MAX_WORDS;
    }

    // Store sequence
    memcpy(&(*kmers)[numKmers].w.b[0],&he.w.b[0],(WL*BITS_NUCLEOTIDE)/8);
    (*kmers)[numKmers].seq = numKmers;
    (*kmers)[numKmers].num = he.num;

    // Locations array space
    if(((*kmers)[numKmers].locations = (uint64_t*) malloc(sizeof(uint64_t)*he.num))==NULL){
      fprintf(stderr , "loadGenome:: Error allocating space for metagenome KMER locations.\n");
      return -1;
    }

    // Take locations
    for(i=0;i<he.num;++i){
      // Read location
      fread(&lo,sizeof(location),1,dP);
      
      if(i == 0){
        // Store sequence index
        (*kmers)[numKmers].seq = lo.seq;
      }
      // Store position
      (*kmers)[numKmers].locations[i] = lo.pos;
    }
    // Update num of kmers
    numKmers++;
    fread(&he,sizeof(hashentry),1,dW);    
  }

  fclose(dW);
  fclose(dP);

  return numKmers;
}


/* This function calculate all matchs between two HE arrays given.
 *  @param w1: HE array to be compared.
 *  @param w2: HE array to be compared.
 *  @param numW1: Number of instances on w1.
 *  @param numW2: number of instances on w2.
 *  @param WL: words length.
 *  @return: number of hits/matches founded.
 */
uint64_t hits(HE* w1,HE* w2,hit** hits,uint64_t numW1, uint64_t numW2,int WL){
  // Variables
  uint64_t numHits = 0;
  int aux = 10;
  uint64_t currentSize = MAX_HITS;

  // Memory for hits
  if((*hits = (hit*) malloc(sizeof(hit)*currentSize))==NULL){
    fprintf(stderr, "hits:: Error allocating memory for hits array.\n");
    return -1;
  }

  // Compare all reads
  int i,j,comp,lastIndex=0;
  for(i=0;i<numW1;++i){
    comp = 0; // Reset value
    for(j=lastIndex;j<numW2 & comp > 0;++j){
      if(numHits >= currentSize){ // Realloc memory if it's necessary 
        if((*hits = (hit*) realloc(*hits,sizeof(hit)*(currentSize+MAX_HITS)))==NULL){
          fprintf(stderr, "hits:: Error reallocating memory for hits array.\n");
          return -1;
        }else
          currentSize += MAX_HITS; // Update "current size"
      }

      comp = wordcmp(&w1[i].w.b[0],&w2[j].w.b[0],(WL*BITS_NUCLEOTIDE)/8);

      if(comp==0){ // Same word (hit)
        int k,h;
        for(k=0; k < w1[i].num;++i)
          for(h=0; h < w2[j].num;++i){
            // Store position
            (*hits)[numHits].start1 = w1[i].locations[k];
            (*hits)[numHits].start2 = w2[j].locations[h];
            // Store length
            (*hits)[numHits].length = WL; 
            // Store sequences indexes
            (*hits)[numHits].seq1 = w1[i].locations[k];
            (*hits)[numHits].seq2 = w2[j].locations[h];
            // Update number of hits
            numHits++;
          }
        lastIndex = j; // Next iterance starts here
      }
    }
  }

  return numHits;
}


/* Function used to group hits that are in the same diagonal and don't have gaps
 * between them.
 *  @param hits: array of hits that will be studied.
 *  @param numHits: number of hits stored on hits array.
 *  @return: new number of hits (grouped hits).
 */
uint64_t groupHits(hit* hits,uint64_t numHits){
  // Variables
  uint64_t newNumHits = 0;

  // Hits must be ordered
  // Group hits
  int i, collapsable;
  uint64_t newEnd;

  for(i=1;i<numHits;++i){    

    if(hits[i].seq1 == hits[newNumHits].seq1 & // Same seqX
       hits[i].seq2 == hits[newNumHits].seq2 & // Same seqY => Same diagonal
       (hits[newNumHits].start1 + hits[newNumHits].length) >= hits[i].start1){ // End of current hit is after next hit start => Collapsable
        // Calc new length
        newEnd = hits[i].start1 + hits[i].length;
        hits[newNumHits].length = newEnd - hits[newNumHits].start1;   
    }else{ // No collapsable
      // Change current hit
      newNumHits++;
      hits[newNumHits].start1 = hits[i].start1;
      hits[newNumHits].start2 = hits[i].start2;
      hits[newNumHits].length = hits[i].length;
      hits[newNumHits].seq1 = hits[i].seq1;
      hits[newNumHits].seq2 = hits[i].seq2;
    }
  }// ROF

  return newNumHits+1;
}


/* This function is used to calculate fragments that satisfy some thresholds given.
 *  @param hits. array of hits that will be studied.
 *  @param numHits: number of instances on hits array.
 *  @param SThreshold: value of similarity threshold.
 *  @param minLength: value of length threshold.
 *  @param out: file where fragments will be written.
 */
int calculateFragments(hit* hits, uint64_t numHits, int SThreshold, int minLength, FILE *out){
  // Variables
  int numFragments = 0;
  FragFile frag;
  
  // Calculate fragments
    // First instance
    storeFragFile(&frag,&hits[0],100);
    frag.block = 0; // Don't change for now
    frag.strand = 'f'; // Reverse not implemented yet


  int i;
  uint64_t fragEnd,dist,oldLength,newLength;
  uint8_t oldS, newS;

  // Exception
  if(numHits == 1 & SThreshold <= 100 & minLength <= frag.length){
    // Write fragment
      fwrite(&frag,sizeof(FragFile),1,out);
      numFragments++;
  }


  for(i=1; i<numHits; i++){
    // If hits are equal discard one
    if(hitComparator(&hits[i],&hits[i-1])==0) continue;

    // Sequence change
    if(hits[i].seq1 != frag.seqX | hits[i].seq2 != frag.seqY){ // Different sequences
      if(frag.length >= minLength){ // Good enough -> write
        // Update values
        frag.xEnd = frag.xStart + frag.length;
        frag.yEnd = frag.yStart + frag.length;
        frag.score = (2*frag.ident - frag.length)*4;
        // Write fragment
        fwrite(&frag,sizeof(FragFile),1,out);
        numFragments++;
      }      
      // Init new fragment
      storeFragFile(&frag,&hits[i],100);
      // Dont alter block and strand for now
      continue;
    }

    // IF same sequences
    if((hits[i].start1 - hits[i].start2) == frag.diag){ // Same diag
      fragEnd = frag.xStart + frag.length;
      newLength = hits[i].start1 + hits[i].length - frag.xStart;
      oldLength = frag.length;
      oldS = frag.similarity;
      dist = fragEnd - hits[i].start1;
      newS = (newLength - (dist + oldLength - oldLength*oldS))/newLength * 100;

      if(newS >= SThreshold){ // Correct fragment, update
        frag.length = newLength;
        frag.similarity = newS;
        frag.ident += newLength - oldLength - dist;
      }else{ // Not enough quality (of new fragment), create new framgent
        if(frag.length >= minLength){ // 
          // Update values
          frag.xEnd = frag.xStart + frag.length;
          frag.yEnd = frag.yStart + frag.length;
          frag.score = (2*frag.ident - frag.length)*4;
          // Write fragment
          fwrite(&frag,sizeof(FragFile),1,out);
          numFragments++;
        }
        //else overwrite actual fragment (invalid length)
        storeFragFile(&frag,&hits[i],100);
      }
    }else{ // Diff diag
      if(frag.length >= minLength){ // 
          // Update values
          frag.xEnd = frag.xStart + frag.length;
          frag.yEnd = frag.yStart + frag.length;
          frag.score = (2*frag.ident - frag.length)*4;
          // Write fragment
          fwrite(&frag,sizeof(FragFile),1,out);
          numFragments++;
        }
        //else overwrite actual fragment (invalid length)
        storeFragFile(&frag,&hits[i],100);
    }
  }

  return numFragments;
} 


/* This function is used to store necessary info in a FragFile instance.
 *  @param frag: FragFile instance where info will be stored.
 *  @param h: hit used to calculate necessary info.
 *  @param s: similarity that will be stored on FragFile instance.
 */
inline void storeFragFile(FragFile* frag,hit* h,float s){
  frag->xStart = h->start1;
  frag->yStart = h->start2;
  frag->xEnd = h->start1 + h->length;
  frag->yEnd = h->start2 + h->length; 
  frag->seqX = h->seq1;
  frag->seqY = h->seq2;
  frag->length = h->length;
  frag->similarity = s;
  frag->diag = h->start1 - h->start2;
  frag->ident = h->length;
  frag->score = SCORE * h->length;
}