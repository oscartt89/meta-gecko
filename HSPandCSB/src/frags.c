/* @file frags.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description: this file contains de code to compare and create
 * 		common frags between a genome set and a metagenome using 
 * 		their dictionaries.
 */

#include "frags.h"

int main(int ac, char** av){
	if(ac!=6){
		fprintf(stderr,"USE: %s genomeSetFolder metagenomeFolder similarityThreshold minLength out\n",av[0]);
		return -1;
	}

	// Variables
	dictionaryG* dicGSet;
	dictionaryM* dicMSet;
	hit *hitsA;
	frag *frags;
	int numGenomes,numMetags,numHits,numFrags;
	char *outputName;

	// Allocate space for outputName variable
	if((outputName = malloc(sizeof(char)*MAX_NAME_L))==NULL){
		fprintf(stderr, "Error allocating space for outputFile name.\n");
		return -1;
	}

	strcpy(outputName,av[5]);

	// Load genome set dictionaries
	if((numGenomes=readGenomeSet(av[1],&dicGSet))<0) return -1;
	else if(numGenomes==0){
		fprintf(stderr, "No genome dictionaries detected.\n");
		return -1;
	}
	// Load metagenome set dictionaries
	if((numMetags=readMetagenomeSet(av[2],&dicMSet))<0) return -1;
	else if(numMetags==0){
		fprintf(stderr, "No metagenome dictionaries detected.\n");
		free(dicGSet);
		return -1;
	}

	// Compare each read with each genome
		//Necessary variables
		int i=0,j,numWM,numWG;
		FILE *dR, *dW, *dP, *fOut;
		wentry *metag, *geno;

	while(i<numMetags){
		// Open dictionaries
		if((dR = fopen(&dicMSet[i].R[0],"rb"))==NULL){
			fprintf(stderr, "Error opening read dictionary. [%s]\n", &dicMSet[i].R);
			return -1;
		}
		if((dW = fopen(&dicMSet[i].W[0],"rb"))==NULL){
			fprintf(stderr, "Error opening words dictionary. [%s]\n", &dicMSet[i].W);
			return -1;
		}
		if((dP = fopen(&dicMSet[i].P[0],"rb"))==NULL){
			fprintf(stderr, "Error opening locations dictionary. [%s]\n", &dicMSet[i].P);
			return -1;
		}

		if(feof(dR) | feof(dW) | feof(dP)){
			fprintf(stderr, "Any of the metagenome dictionaries are empty.\n");
			return -1;
		}
		
		// Compare each read with each genome
		while(!feof(dR)){
			// Load read
			if((numWM = loadRead(dR,dW,dP,&metag))<0) return -1;
			// Compare with each genome
			for(j=0; j<numGenomes; ++j){
				// Load genome
				if((numWG = loadGenome(dicGSet[j],&geno))<0) return -1;
				// Calc hits
					// For now only 100% are allowed on hits
				if((numHits = hits(metag,geno,&hitsA,numWM,numWG))<0) return -1;
				// Sort hits
				if(quickSort(hitsA,0,numHits-1)<0) return -1;
				// Filter hits. Calculte fragments
				if((numFrags=calculateFragments(hitsA,&frags,numHits,atoi(av[3]),atoi(av[4])))<0) return -1;
				free(hitsA); // Free unnecesary space
				// Write frags file
					// Open output stream
					strcat(outputName,&dicMSet[i].name[0]);
					strcat(outputName,"_");
					strcat(outputName,&dicGSet[j].name[0]);
					strcat(outputName,".frags");
					if((fOut = fopen(outputName,"wb"))==NULL){
						fprintf(stderr, "Error opening output fragment file. [%s]\n", outputName);
						return -1;
					}
					strcpy(outputName,av[5]); //Reset name
			
				fwrite(&frags[k],sizeof(frag),numFrags,fOut);
				
				fclose(fOut);

				// Free space
				free(geno);
			}

			free(metag); // Free space
		}

		fclose(dR);
		fclose(dW);
		fclose(dP);
		++i;		
	}// for

	free(dicMSet);
	free(dicGSet);

	return 0;
}