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
int readGenomeSet(char* genomeSetPath,dictionaryG* genomes){
	// Variables
	DIR *gFolder;
	struct dirent *ent;
	char *wD, *pD;
	int numGenomes = 0, currentMax = MAX_GENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  free(genomes);
  if((genomes = (dictionaryG*) malloc(sizeof(dictionaryG)*MAX_GENOME_SET))==NULL){
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
      if((genomes = realloc(genomes,sizeof(dictionaryG)*(currentMax + MAX_GENOME_SET)))==NULL){
        fprintf(stderr, "Error reallicatin memory for genome dictionary set.\n");
        return -1;
      }
      currentMax += MAX_GENOME_SET;
    }
		// Files are sorted alphabetically
		// Should appear first d2hP than d2hW
		if(strstr(ent->d_name,".d2hP")!=NULL){ // New dictionary
			// Save name
			memcpy(genomes[numGenomes].name,ent->d_name,strlen(ent->d_name)-5);
			// Save location dictionary
			memcpy(genomes[numGenomes].P,ent->d_name,sizeof(ent->d_name));
			//Next file should be d2hW dictionary
			if((ent = readdir(gFolder))==NULL){
				fprintf(stderr, "Error: incomplete genome pair dictionary. End of file list.\n");
				return -1;
			}
			if(strstr(ent->d_name,".d2hW")!=NULL){
				// Save word dictionary
				memcpy(genomes[numGenomes].W,ent->d_name,sizeof(ent->d_name));
				numGenomes++;
			}else if(strstr(ent->d_name,".metag.d2hR")!=NULL){
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
int readMetagenomeSet(char* metagSetPath,dictionaryM* metagenomes){
  // Variables
  DIR *mFolder;
  struct dirent *ent;
  char *wD, *pD, *rD;
  int numMetags = 0, currentMax = MAX_METAGENOME_SET;
  FILE *WD, *PD;

  // Allocate memory for genome dirs
  free(metagenomes);
  if((metagenomes = (dictionaryM*) malloc(sizeof(dictionaryM)*MAX_METAGENOME_SET))==NULL){
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
      if((metagenomes = realloc(metagenomes,sizeof(dictionaryM)*(currentMax + MAX_METAGENOME_SET)))==NULL){
        fprintf(stderr, "Error reallicatin memory for metagenome dictionary set.\n");
        return -1;
      }
      currentMax += MAX_METAGENOME_SET;
    }
    // Files are sorted alphabetically
    // Should appear first d2hP, then d2hR and d2hW
    if(strstr(ent->d_name,".metag.d2hP")!=NULL){ // New dictionary
      // Save name
      memcpy(metagenomes[numMetags].name,ent->d_name,strlen(ent->d_name)-11);
      // Save location dictionary
      memcpy(metagenomes[numMetags].P,ent->d_name,sizeof(ent->d_name));
      //Next file should be d2hR dictionary
      if((ent = readdir(mFolder))==NULL){
        fprintf(stderr, "Error: incomplete metagenome triplet dictionary. End of file list.\n");
        return -1;
      }
      if(strstr(ent->d_name,".metag.d2hR")!=NULL){
        // Save read dictionary
        memcpy(metagenomes[numMetags].R,ent->d_name,sizeof(ent->d_name));
        // Now should appear words dictionary
        if((ent = readdir(mFolder))==NULL){
          fprintf(stderr, "Error: incomplete metagenome triplet dictionary. End of file list.\n");
          return -1;
        }
        if(strstr(ent->d_name,".metag.d2hW")!=NULL){
          // Save words dictionary
          memcpy(metagenomes[numMetags].W,ent->d_name,sizeof(ent->d_name));
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