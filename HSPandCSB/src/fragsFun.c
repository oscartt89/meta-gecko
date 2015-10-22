/* @file fragsFun.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description this file contains the functions necesary for
 * 		correct operation of frags.c code. 
 */

#include "frags.h"


int readGenomeSet(char* genomeSetPath,wentry** genomeDicts,int* lengths){
	// Variables
	dictionaryG genomes[MAX_GENOME_SET];
	DIR *gFolder;
	struct dirent *ent;
	char *wD, *pD;
	int numGenomes = 0;
  FILE *WD, *PD;

	// Open genomes folder
	if((gFolder = opendir(genomeSetPath))==NULL){
		fprintf(stderr, "Error opening genomes folder.\n");
		return -1;
  }

  	// Take genome dictionary files
  	while((ent = readdir(gFolder))!=NULL){
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

  	// For each genome dictionary, creat wentry array
  	// Alocate dictionaries array
  	if((genomeDicts = (wentry**) malloc(sizeof(wentry*)*numGenomes))==NULL){
  		fprintf(stderr, "Error allocating array of KMERS\n");
  		return -1;
  	}
    // Alocate dictionaries length array
    if((lengths = (int*) malloc(sizeof(int)*numGenomes))==NULL){
      fprintf(stderr, "Error allocating array of lengths.\n");
      return -1;
    }

    // Create wentry array for each dictionary
    int i;
    for(i=0; i<numGenomes; ++i){

    }

  	return numGenomes;
}