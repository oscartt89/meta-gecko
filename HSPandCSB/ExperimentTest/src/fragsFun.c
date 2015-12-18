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
  int numGenomes = 0, currentMax = MAX_GENOME_SET;
  char **files;

  int numFiles = listFiles(genomeSetPath,&files);

  // Sort files alphabetically
  quickSort_S(files,0,numFiles);

  // Allocate memory for genome dirs
  free(*genomes);
  if((*genomes = (dictionaryG*) malloc(sizeof(dictionaryG)*MAX_GENOME_SET))==NULL){
    fprintf(stderr, "readGenomeSet:: Error allocating memory for genome dictionary set.\n");
    return -1;
  }

  // Take genome dictionary files
  int i;
  for(i=0; i<numFiles; ++i){
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
    if(strstr(files[i],".metag.d2h")!=NULL) continue;
    // Should appear first d2hP than d2hW
    if(endsWith(files[i],".d2hP")){ // New dictionary
      // Save name
      memcpy(genomes[numGenomes]->name,files[i],strlen(files[i])-5);
      // Save location dictionary
      strcpy(genomes[numGenomes]->P,files[i]);
      ++i;
      //Next file should be d2hW dictionary
      if(i >= numFiles){
        fprintf(stderr, "readGenomeSet:: Error: incomplete genome pair dictionary. End of file list.\n");
        return -1;
       }
      
      if(endsWith(files[i],".d2hW") && strstr(files[i],".metag.d2h")==NULL){
        // Save word dictionary
        strcpy(genomes[numGenomes]->W,files[i]);
        numGenomes++;
      }else if(strstr(files[i],".metag.d2h")!=NULL){
        fprintf(stderr, "readGenomeSet:: Error: it's a metagenome dictionary.\n");
        return -1;
      }else{
        fprintf(stderr, "readGenomeSet:: Error: incomplete genome pair dictionary.\n");
        return -1;
      }
    }
  }

  free_Files(files,numFiles);

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
  int numMetags = 0, currentMax = MAX_METAGENOME_SET;
  char **files;

  int numFiles = listFiles(metagSetPath,&files);

  // Sort files alphabetically
  quickSort_S(files,0,numFiles);

  // Allocate memory for genome dirs
  if((*metagenomes = (dictionaryM*) malloc(sizeof(dictionaryM)*MAX_METAGENOME_SET))==NULL){
    fprintf(stderr, "readMetagenomeSet:: Error allocating memory for metagenome dictionary set.\n");
    return -1;
  }

  // Take metagenome dictionary files
  int i;
  for(i=0; i<numFiles; ++i){
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
    if(endsWith(files[i],".metag.d2hP")){ // New dictionary
      // Save name
      memcpy(&metagenomes[numMetags]->name[0],files[i],strlen(files[i])-11);
      // Save location dictionary
      strcpy(&metagenomes[numMetags]->P[0],files[i]);
      //Next file should be d2hR dictionary
      ++i;
      if(i >= numFiles){
        fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triplet dictionary. End of file list searching read dictionary.\n");
        return -1;
      }

      if(endsWith(files[i],".metag.d2hR")){
        // Save read dictionary
        strcpy(&metagenomes[numMetags]->R[0],files[i]);
        // Now should appear words dictionary
        ++i;
        if(i >= numFiles){
          fprintf(stderr, "readMetagenomeSet:: Error: incomplete metagenome triplet dictionary. End of file list searching word dictionary .\n");
          return -1;
        }

        if(endsWith(files[i],".metag.d2hW")){
          // Save words dictionary
          strcpy(&metagenomes[numMetags]->W[0],files[i]);
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
  free_Files(files,numFiles);

  return numMetags;
}


/* This function takes a string and returns if the string 
 * starts with a dot character.
 *  @param string: string to be studied.
 *  @return: zero if the string doesn't starts with a dot or
 *    a positive number in other cases.
 */
int startsWithDot(char* str){
  size_t len = strlen(str);
  return len >= 1 ? strncmp(str,".",1) == 0 : 0;
}


/* This function takes a string and return if the strign
 * ends with a suffix given.
 *  @param str: string to be studied.
 *  @param suffix: string to be matched on str end.
 *  @return: zero if the string doesn't ends with the specified
 *    suffix or a positive number in other case.
 */
int endsWith(char *str, char *suffix){
    size_t l_str = strlen(str),
           l_suf = strlen(suffix);
    
    // Suffix bigger than string to compare
    if(l_suf > l_str)
      return 0;

    return strncmp(str + l_str - l_suf, suffix, l_suf) == 0;
}


/* This function is used to list all files in a given directory.
 *  @param path: path of the directory.
 *  @param files: array of strings where directories will be stored.
 *  @return: the length of the final files array.
 */
int listFiles(char *path, char*** files){
  // Variables
  DIR *mFolder;
  struct dirent *ent;
  char f[MAX_FILE_LENGTH];
  int numFiles = 0;
  int currentMax = MAX_FILES;

  // Open metagenomes folder
  if((mFolder = opendir(path))==NULL){
    fprintf(stderr, "listFiles:: Error opening path folder.\n");
    return -1;
  }

  // Take memory for files
  if((*files = (char**) malloc(sizeof(char)*currentMax))==NULL){
    fprintf(stderr, "listFiles:: Error allocating memory for files set.\n");
    return -1;
  }

  // Take files from directory
  while((ent = readdir(mFolder))!=NULL){
    // Realloc if it's necessary
    if(numFiles >= currentMax){
      if(*files = realloc(*files,sizeof(char)*(currentMax + MAX_FILES))==NULL){
        fprintf(stderr, "listFiles:: Error reallocating memory for files set.\n");
        free(*files);
        return -1;
      }
      currentMax += MAX_FILES;
    }

    // Prepare file name
    strcpy(&f,path);
    strcat(&f,ent->d_name);

    // It it's a directory, avoid it
    if(opendir(&f)!=NULL) continue;

    // Check that it's not a hidden file or a system directory
    if(startsWithDot(ent->d_name)) continue;

    // It's a file, save it
    //Allocate necessary memory
    if(((*files)[numFiles] = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
      fprintf(stderr, "listFiles:: Error allocating space for a file path.\n");
      return -1;
    }
    strcpy((*files)[numFiles],&f);
    numFiles++;
  }

  // Close and finish
  closedir(mFolder);
  return numFiles;
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
  fread(&r,sizeof(Read),1,dR); 

  // KMERS space
  if((*kmers = malloc(sizeof(HE)*MAX_WORDS))==NULL){
    fprintf(stderr , "loadRead:: Error allocating space for metagenome KMERS.\n");
    return -1;
  }

  // Positionate on word dictionary
  fseek(dW,r.pos,SEEK_SET);

  // Load kmers in HE array
  int i;
  for(i=0;i<r.num;++i){
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
    int r;
    // Take pos
    if((r=fread(&(*kmers)[numKmers].locations[0],sizeof(uint64_t),he.num,dP)) != he.num){
      fprintf(stderr, "loadRead:: Error reading locations (%i/%i)\n", r, he.num);
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
  fread(&he,sizeof(hashentryOld),1,dW);

  // Read words
  int i;
  while(!feof(dW)){
    // Num of kmers exceeded
    if(numKmers >= currentMax){
      if((*kmers = realloc(*kmers,sizeof(HE)*(currentMax + MAX_WORDS)))==NULL){
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
      fprintf(stderr , "loadGenome:: Error allocating space for genome KMER locations.\n");
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
  uint64_t currentSize = MAX_HITS;

  // Memory for hits
  if((*hits = (hit*) malloc(sizeof(hit)*currentSize))==NULL){
    fprintf(stderr, "hits:: Error allocating memory for hits array.\n");
    return -1;
  }

  // Compare all reads
  int i,j=0,comp=1;
  for(i=0;i<numW1;++i){
    comp = 1; // Reset value
    for(;j<numW2 && comp > 0;++j){
      comp = wordcmp(&w1[i].w.b[0],&w2[j].w.b[0],(WL*BITS_NUCLEOTIDE)/8);
      if(comp==0){ // Same word (hit)
        int k,h;
        for(k=0; k < w1[i].num;++k){
          for(h=0; h < w2[j].num;++h){
            if(numHits >= currentSize){ // Realloc memory if it's necessary 
              if((*hits = (hit*) realloc(*hits,sizeof(hit)*(currentSize+MAX_HITS)))==NULL){
                fprintf(stderr, "hits:: Error reallocating memory for hits array.\n");
                return -1;
              }else{
                currentSize += MAX_HITS; // Update "current size"
              }
            }

            // Store position
            (*hits)[numHits].start1 = w1[i].locations[k];
            (*hits)[numHits].start2 = w2[j].locations[h];
            // Store length
            (*hits)[numHits].length = WL; 
            // Store sequences indexes
            (*hits)[numHits].seq1 = w1[i].seq;
            (*hits)[numHits].seq2 = w2[j].seq;
            // Update number of hits
            numHits++;
          }
        }
        break;
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
  int i;
  uint64_t newEnd,oldEnd;

  for(i=1;i<numHits;++i){    
    if(hits[i].seq1 == hits[newNumHits].seq1 && // Same seqX
       hits[i].seq2 == hits[newNumHits].seq2 && // Same seqY => Same diagonal
       (hits[newNumHits].start1 + hits[newNumHits].length - 1) >= hits[i].start1){ // End of current hit is after next hit start => Collapsable
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\n\tCollapse -> OldSt=%" PRIu64 " OldL:%" PRIu64 " HitS=%" PRIu64 " HitL=%" PRIu64 "\n", hits[newNumHits].start1, hits[newNumHits].length, hits[i].start1, hits[i].length);
/////////////////////////////////////////////////////////////////// 
        // Calc new length
        oldEnd = hits[newNumHits].start1 + hits[newNumHits].length - 1;
        newEnd = hits[i].start1 + hits[i].length - 1;
        if(oldEnd < newEnd) // Only if new end is greater
          hits[newNumHits].length = newEnd - hits[newNumHits].start1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tCollapse -> NStart=%" PRIu64 " NLength:%" PRIu64 "\n", hits[newNumHits].start1, hits[newNumHits].length);
///////////////////////////////////////////////////////////////////  
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

//fprintf(stdout, "\t");
  
  // Calculate fragments
    // First instance
    storeFragFile(&frag,&hits[0],100);
    frag.block = 0; // Don't change for now
    frag.strand = 'f'; // Reverse not implemented yet


  int i;
  uint64_t fragEnd,dist,oldLength,newLength;
  uint8_t oldS, newS;

  // Exception
  if(numHits == 1 && minLength <= frag.length){
    // Write fragment
      fwrite(&frag,sizeof(FragFile),1,out);
      return 1;
  }

  for(i=1; i<numHits; i++){
//fprintf(stdout, "!");
    // If hits are equal discard one
    if(hitComparator(&hits[i],&hits[i-1])==0) continue;

    // Sequence change
    if(hits[i].seq1 != frag.seqX || hits[i].seq2 != frag.seqY){ // Different sequences
      if(frag.length >= minLength){ // Good enough -> write
        // Update values
        frag.xEnd = frag.xStart + frag.length - 1;
        frag.yEnd = frag.yStart + frag.length - 1;
        frag.score = (2*frag.ident - frag.length)*4;
        // Write fragment
        fwrite(&frag,sizeof(FragFile),1,out);
        numFragments++;
//fprintf(stdout, "+");
      }      
      // Init new fragment
      storeFragFile(&frag,&hits[i],100);
      // Dont alter block and strand for now
      continue;
    }
//fprintf(stdout, "!");
    // IF same sequences
    if((hits[i].start1 - hits[i].start2) == frag.diag){ // Same diag
      fragEnd = frag.xStart + frag.length - 1;
      newLength = hits[i].start1 + hits[i].length - frag.xStart - 1;
      oldLength = frag.length;
      oldS = frag.similarity;
      dist = fragEnd - hits[i].start1;
//fprintf(stdout, "-");
//fprintf(stdout, "NL:%" PRIu64 "  HS:%" PRIu64 "  HL:%" PRIu64 "  FS:%" PRIu64 "\n"
//,newLength,hits[i].start1,hits[i].length,frag.xStart);
      newS = (newLength - (dist + oldLength - oldLength*oldS))/newLength * 100;
//fprintf(stdout, "-");
      if(newS >= SThreshold){ // Correct fragment, update
        frag.length = newLength;
        frag.similarity = newS;
        frag.ident += newLength - oldLength - dist;
//fprintf(stdout, ";");
      }else{ // Not enough quality (of new fragment), create new framgent
        if(frag.length >= minLength){ // 
          // Update values
          frag.xEnd = frag.xStart + frag.length - 1;
          frag.yEnd = frag.yStart + frag.length - 1;
          frag.score = (2*frag.ident - frag.length)*4;
          // Write fragment
          fwrite(&frag,sizeof(FragFile),1,out);
          numFragments++;
//fprintf(stdout, "+");
        }
        //else overwrite actual fragment (invalid length)
        storeFragFile(&frag,&hits[i],100);
      }
    }else{ // Diff diag
//fprintf(stdout, ".");
      if(frag.length >= minLength){ // 
          // Update values
          frag.xEnd = frag.xStart + frag.length - 1;
          frag.yEnd = frag.yStart + frag.length - 1;
          frag.score = (2*frag.ident - frag.length)*4;
          // Write fragment
          fwrite(&frag,sizeof(FragFile),1,out);
          numFragments++;
//fprintf(stdout, "+");
        }
//fprintf(stdout, "_");
        //else overwrite actual fragment (invalid length)
        storeFragFile(&frag,&hits[i],100);
    }
  }

  // Store last fragment
  if(frag.length >= minLength){ // Good enough -> write
    // Update values
    frag.xEnd = frag.xStart + frag.length - 1;
    frag.yEnd = frag.yStart + frag.length - 1;
    frag.score = (2*frag.ident - frag.length)*4;
    // Write fragment
    fwrite(&frag,sizeof(FragFile),1,out);
    numFragments++;
//fprintf(stdout, "_");
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
  frag->xEnd = h->start1 + h->length - 1;
  frag->yEnd = h->start2 + h->length - 1; 
  frag->seqX = h->seq1;
  frag->seqY = h->seq2;
  frag->length = h->length;
  frag->similarity = s;
  frag->diag = h->start1 - h->start2;
  frag->ident = h->length;
  frag->score = SCORE * h->length;
}


/* This functions is used to deallocate memory space of an HE array.
 *  @param heArray: HE array that will be deallocated.
 *  @param length: length of heArray
 */
inline void free_HE(HE* heArray,int length){
  int i;
  for(i=0 ; i<length; ++i)
    free(heArray[i].locations);

  free(heArray);
}


/* This functions is used to deallocate memory space of an files array.
 *  @param fArray: char** array that will be deallocated.
 *  @param length: length of fArray.
 */
inline void free_Files(char** fArray,int length){
  int i;
  for(i=0 ; i<length; ++i)
    free(fArray[i]);
  free(fArray);
}


/* This function is used to swap/interchange two strings.
 *  @param s1: string that will be swapped.
 *  @param s2: string that will be swapped.
 */
inline void SWAP_S(char *s1,char *s2){
  char t[MAX_FILE_LENGTH];
  strcpy(&t,s1);
  strcpy(s1,s2);
  strcpy(s2,&t);
}


/* Function used to sort an array of strings.
 *  @param strs: array of string that will be sorted.
 *  @param start: index where start to sort.
 *  @param length: length of the array to sort (starting on "start" index).
 */
void quickSort_S(char** strs, uint64_t start, uint64_t length){
  //Check exceptions
  if (length < 2) return;
  else if(length == 2){
    if(strcmp(strs[start],strs[start+1])>0)
      SWAP_S(strs[start],strs[start+1]);
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
      while(strcmp(strs[right],strs[pivot])<0 && right <= pivot) right++;
      // While left y greater than pivot, continue
      while(strcmp(strs[pivot],strs[left])<0 && left >= pivot) left--;
      // All array checked -> end
      if(right>=left)
          break;
      // Swap selected values
      SWAP_S(strs[right],strs[left]);
      changed = 1;
      if(right < pivot) right++;
      if(left > pivot) left--;
  }
  if(changed && right<=left){ // If anything changed, don't continue, it's sorted
    quickSort_S(strs,start,right);
    quickSort_S(strs,right,length-right);
    if(length == 3)
      quickSort_S(strs,start,length);
  }

  return;
}