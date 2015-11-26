/* @file frags.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description: this file contains de code to compare and create
 * 		common frags between a genome set and a metagenome using 
 * 		their dictionaries.
 */

#include "frags.h"

int main(int ac, char** av){
	if(ac!=7){
		fprintf(stderr,"USE: %s genomeSetFolder metagenomeFolder similarityThreshold minLength wordLength out\n",av[0]);
		return -1;
	}

	// Variables
	dictionaryG* dicGSet;
	dictionaryM* dicMSet;
	hit *hitsA;
//	FragFile *frags;
	FILE *fOut;
	uint64_t numGenomes,numMetags,numHits,numGHits,numFrags;

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

///////////////////////////////////////////////////////////////////
//fprintf(stdout, "TEST::M=%d::G=%d\n",numMetags,numGenomes);
///////////////////////////////////////////////////////////////////


	// Open output file fragment
	char *outF = av[6];
	strcat(outF,".fr");
	if((fOut = fopen(outF,"wb"))==NULL){
		fprintf(stderr, "Error opening output fragment file. [%s]\n", av[6]);
		return -1;
	}

	// Compare each read with each genome
		//Necessary variables
		uint64_t numWM,numWG=-1;
		int i=0,j;
		FILE *dR, *dW, *dP;
		HE *metag, *geno;

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


		// Special case
		if(numGenomes == 1)
			if((numWG = loadGenome(dicGSet[0],&geno,atoi(av[5])))<0) return -1;


		// Compare each read with each genome
		while(!feof(dR)){
			// Load read
			if((numWM = loadRead(dR,dW,dP,&metag,atoi(av[5])))<0) return -1;
			// Compare with each genome
			for(j=0; j<numGenomes & numWM>0; ++j){
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "WM: %d", numWM);
///////////////////////////////////////////////////////////////////
				// Load genome
				if(numGenomes > 1)
					if((numWG = loadGenome(dicGSet[j],&geno,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tWG: %d", numWG);
///////////////////////////////////////////////////////////////////
				// Calc hits 
				if(numWG > 0){
					// For now only 100% are allowed on hits -> No gaps
					if((numHits = hits(metag,geno,&hitsA,numWM,numWG,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tHits: %d",numHits);
///////////////////////////////////////////////////////////////////
					if(numGenomes > 1)
						// Free space
						free(geno);
					// Sort hits
					if(quickSort(hitsA,0,numHits-1)<0) return -1;
///////////////////////////////////////////////////////////////////
//fwrite(&hitsA[0],sizeof(hit),numHits,fOut);
///////////////////////////////////////////////////////////////////
					// Group hits
					if(numHits>0){
						if((numGHits=groupHits(hitsA,numHits))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tG_Hits: %d",numGHits);
///////////////////////////////////////////////////////////////////
					// Filter hits. Calculte fragments
						if((numFrags=calculateFragments(hitsA,numGHits,atoi(av[3]),atoi(av[4]),fOut))<0) return -1;
						free(hitsA); // Free unnecesary space
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tFrags: %d",numFrags);
///////////////////////////////////////////////////////////////////
					}else
						free(hitsA);
					// Write frags file
	//				fwrite(&frags[0],sizeof(FragFile),numFrags,fOut);
				}
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\n");
///////////////////////////////////////////////////////////////////

			}

			free(metag); // Free space
		}

		if(numGenomes == 1) free(geno);

		fclose(dR);
		fclose(dW);
		fclose(dP);

		++i;		
	}// for
	fclose(fOut);
	free(dicMSet);
	free(dicGSet);

	return 0;
}