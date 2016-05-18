/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes a process to read and take statistical information
 *    from a metagenome-genome fragments file.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "frags.h"

void readFragment(FragFile*,FILE*);

int main(int ac, char** av){
	// Variables
	FragFile fragment;
	float meanL, meanSi, meanSco, minSi, maxSi;
	uint64_t size, num_frags=0, minL, maxL, minSco, maxSco;
	FILE *fragF, *out;

	// Check arguments
    if(ac!=4){
        fprintf(stderr, "Bad call error.\nUSE: fixeReverseFragments fragFile prefix out\n");
        return -1;
    }

    // Open fragment file
    if((fragF = fopen(av[1],"rb"))==NULL){
		fprintf(stderr, "Error opening frag file\n");
		return -1;
	}

    // Open output file
    if(exists(av[3])){
		if((out = fopen(av[3],"at"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[4]);
			return -1;
		}
	}else{
		if((out = fopen(av[3],"wt"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[4]);
			return -1;
		}
		fprintf(out, "FragFile\tSize\tPrefix\tFragments\tMeanL\tMeanScore\tMeanS\tMaxL\tMinL\tMaxScore\tMinScore\tMaxS\tMinS\n");
	}


	// Read first fragment
    readFragment(&fragment,fragF);
    fprintf(stderr, "T\n");

    if(feof(fragF)) return 0;

    // Take stats
    minL = fragment.length;
    maxL = fragment.length;
    minSco = fragment.score;
    maxSco = fragment.score;
    minSi = fragment.similarity;
    maxSi = fragment.similarity;
    meanL = fragment.length;
    meanSi = fragment.similarity;
    meanSco = fragment.score;
    num_frags++;

    // Read all fragments
    readFragment(&fragment,fragF);
    while(!feof(fragF)){
    	if(minL > fragment.length) minL = fragment.length; 
	    if(maxL < fragment.length) maxL = fragment.length;
	    if(minSco > fragment.score) minSco = fragment.score;
	    if(maxSco < fragment.score) maxSco = fragment.score;
	    if(minSi > fragment.similarity) minSi = fragment.similarity;
	    if(maxSi < fragment.similarity) maxSi = fragment.similarity;
	    meanL = (meanL*num_frags + fragment.length)/(num_frags+1);
	    meanSi = (meanSi*num_frags + fragment.similarity)/(num_frags+1);
	    meanSco = (meanSco*num_frags + fragment.score)/(num_frags+1);
	    num_frags++;

	    readFragment(&fragment,fragF);
    }

    fseek(fragF,0L,SEEK_END);
    size = ftell(fragF);

    fclose(fragF);

    // Write stats
    fprintf(out, "%s\t%"PRIu64"\t%s\t%"PRIu64"\t%f\t%f\t%f\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%f\t%f\n",
            av[1],size,av[2],num_frags,meanL,meanSco,meanSi,maxL,minL,maxSco,minSco,maxSi,minSi);

    fclose(out);
    return 0;
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(FragFile *frag, FILE *f){
    char tmpArray[8];

    if(htons(1)==1){
        //big endian
        if(fread(&frag->diag, sizeof(int64_t), 1, f)!=1){
            if(feof(f))return;
            fprintf(stderr,"Error reading the HSP diagonal");
            return;
        }
        if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xStart");
            return;
        }
        if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yStart");
            return;
        }
        if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xEnd");
            return;
        }
        if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yEnd");
            return;
        }
        if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP length");
            return;
        }
        if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP identities");
            return;
        }
        if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP score");
            return;
        }
        if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP similarity");
            return;
        }
        if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqX");
            return;
        }
        if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqY");
            return;
        }
        if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP block");
            return;
        }
        frag->strand = fgetc(f);
    } else {
        //little endian
        if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
            if(feof(f))return;
            fprintf(stderr,"Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *)(&frag->diag), sizeof(int64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xStart");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->xStart), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yStart");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->yStart), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xEnd");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->xEnd), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yEnd");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->yEnd), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP length");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->length), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP identity");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->ident), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP score");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->score), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(float), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP float");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->similarity), sizeof(float)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqX");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->seqX), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqY");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->seqY), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP block");
            return;
        }
        endianessConversion(tmpArray, (char *)(&frag->block), sizeof(int64_t)); 
        frag->strand = fgetc(f);
    }
}