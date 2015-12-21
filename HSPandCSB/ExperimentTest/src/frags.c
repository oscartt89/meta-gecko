/* @file frags.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description: this file contains de code to compare and create
 * 		common frags between a genome set and a metagenome using 
 * 		their dictionaries.
 */

#include "frags.h"

///////////////////////////////////////////-------/////////////////
#include <time.h> 
//////////////////////////////////////////-------/////////////////

int main(int ac, char** av){
	if(ac!=7){
		fprintf(stderr,"USE: %s genomeSetFolder metagenomeFolder similarityThreshold minLength wordLength out\n",av[0]);
		return -1;
	}

///////////////////////////////////////////-------/////////////////
time_t rawtime;
struct tm * timeinfo;
time ( &rawtime );
timeinfo = localtime ( &rawtime );
fprintf (stdout,"Init-> %s", asctime(timeinfo));
///////////////////////////////////////////-------/////////////////

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
		fprintf(stderr, "MAIN:: No genome dictionaries detected.\n");
		return -1;
	}

	// Load metagenome set dictionaries
	if((numMetags=readMetagenomeSet(av[2],&dicMSet))<0) return -1;
	else if(numMetags==0){
		fprintf(stderr, "MAIN:: No metagenome dictionaries detected.\n");
		free(dicGSet);
		return -1;
	}

///////////////////////////////////////////////////////////////////
//fprintf(stdout, "TEST::M=%" PRIu64 "::G=%" PRIu64 "\n",numMetags,numGenomes);
///////////////////////////////////////////////////////////////////


	// Open output file fragment
	char *outF = av[6];
	strcat(outF,".fr");
	if((fOut = fopen(outF,"wb"))==NULL){
		fprintf(stderr, "MAIN:: Error opening output fragment file. [%s]\n", av[6]);
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
			fprintf(stderr, "MAIN:: Error opening read dictionary.\n");
			return -1;
		}
		if((dW = fopen(&dicMSet[i].W[0],"rb"))==NULL){
			fprintf(stderr, "MAIN:: Error opening words dictionary.\n");
			return -1;
		}
		if((dP = fopen(&dicMSet[i].P[0],"rb"))==NULL){
			fprintf(stderr, "MAIN:: Error opening locations dictionary.\n");
			return -1;
		}

		if(feof(dR) | feof(dW) | feof(dP)){
			fprintf(stderr, "MAIN:: Any of the metagenome dictionaries are empty.\n");
			return -1;
		}


		// Special case
		if(numGenomes == 1)
			if((numWG = loadGenome(dicGSet[0],&geno,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//int numR = 0;
///////////////////////////////////////////////////////////////////
		// Compare each read with each genome
		while(!feof(dR)){
			// Load read
			if((numWM = loadRead(dR,dW,dP,&metag,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//numR++;
///////////////////////////////////////////////////////////////////
			// Compare with each genome
			for(j=0; j<numGenomes && numWM>0; ++j){
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "R%i", numR);
//fprintf(stdout, "\tWM: %" PRIu64 , numWM);
///////////////////////////////////////////////////////////////////
				// Load genome
				if(numGenomes > 1)
					if((numWG = loadGenome(dicGSet[j],&geno,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tGenome_Loaded");
//fprintf(stdout, "\tWG: %" PRIu64 , numWG);
///////////////////////////////////////////////////////////////////
				// Calc hits 
				if(numWG > 0){
					// For now only 100% are allowed on hits -> No gaps
					if((numHits = hits(metag,geno,&hitsA,numWM,numWG,atoi(av[5])))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tHits");
//fprintf(stdout, "\tH: %" PRIu64 ,numHits);
///////////////////////////////////////////////////////////////////
					if(numGenomes > 1){
						// Free space
						free_HE(geno,numWG);
					}

					if(numHits>0){
						// Sort hits
						quickSort_H(hitsA,0,numHits);
						// Group hits
						if((numGHits=groupHits(hitsA,numHits))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tG_Hits");
//fprintf(stdout, "\tG_H: %" PRIu64 ,numGHits);
///////////////////////////////////////////////////////////////////
					// Filter hits. Calculte fragments
						if((numFrags=calculateFragments(hitsA,numGHits,atoi(av[3]),atoi(av[4]),fOut))<0) return -1;
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\tFrags");
//fprintf(stdout, "\tF: %" PRIu64,numFrags);
///////////////////////////////////////////////////////////////////
					}
						free(hitsA);
				}
///////////////////////////////////////////////////////////////////
//fprintf(stdout, "\n");
///////////////////////////////////////////////////////////////////

			}

			free_HE(metag,numWM); // Free space
		}

		if(numGenomes == 1) free_HE(geno,numWG);

		fclose(dR);
		fclose(dW);
		fclose(dP);

		++i;		
	}// for
	fclose(fOut);
	free(dicMSet);
	free(dicGSet);

///////////////////////////////////////////-------/////////////////
//fprintf (stdout,"Init-> %s", asctime(timeinfo));
time ( &rawtime );
timeinfo = localtime ( &rawtime );
fprintf (stdout,"End-> %s", asctime(timeinfo));
///////////////////////////////////////////-------/////////////////



	return 0;
}