/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdint.h> // unit64_t ...
#include <inttypes.h> 
#include <string.h> // String functions

// Structures
typedef struct{
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
} FragFile;

// Functions
inline void readFragment(FragFile*,FILE*);
inline void writeFragment(FragFile,FILE*);
int fragmentComparator(FragFile,FragFile);

// MAIN
int main(int ac, char** av){
	// Variables
	FILE *forward,*reverse,*out;
	FragFile forw,reve;

	// Check arguments
	if(ac!=5){
		fprintf(stderr, "Bad call error.\nUSE: combineFragments forwardFrags reverseFrags outputFrags delete\nNote: delete=1 -> TRUE\n");
		return -1;
	}

	// Check names
	if(strcmp(av[1],av[2])==0){
		fprintf(stderr, "Forward and reverse files couldn't be the same.\n");
		return -1;
	}else if(strcmp(av[1],av[3])==0){
		fprintf(stderr, "Forward and final files couldn't be the same.\n");
		return -1;
	}else if(strcmp(av[2],av[3])==0){
		fprintf(stderr, "Reverse and final files couldn't be the same.\n");
		return -1;
	}

	// Open files
	if((forward = fopen(av[1],"rb"))==NULL){
		fprintf(stderr, "Error opening forward fragments file.\n");
		return -1;
	}

	if((reverse = fopen(av[2],"rb"))==NULL){
		fprintf(stderr, "Error opening reverse fragments file.\n");
		return -1;
	}

	if((out = fopen(av[3],"wb"))==NULL){
		fprintf(stderr, "Error opening final fragments file.\n");
		return -1;
	}

	// Load first fragments
	readFragment(&forw,forward);
	readFragment(&reve,reverse);

	// Read,compare and write fragments
	while(!feof(forward) && !feof(reverse)){
		if(fragmentComparator(forw,reve)>=0){ // Write reve
			writeFragment(reve,out);
			readFragment(&reve,reverse);
		}else{ // Write forward
			writeFragment(forw,out);
			readFragment(&forw,forward);
		}
	}

	if(feof(forward)){
		while(!feof(reverse)){
			writeFragment(reve,out);
			readFragment(&reve,reverse);
		}
	}else{ // feof(reverse)
		while(!feof(forward)){
			writeFragment(forw,out);
			readFragment(&forw,forward);
		}
	}

	// Delete if it's necessary
	if(strcmp(av[4],"1")==0){
		remove(av[1]);
		remove(av[2]);
	}

	// Ends OK
	return 0;
}


/*
 */
inline void writeFragment(FragFile frag,FILE *fr){
	fwrite(&frag.diag,sizeof(int64_t),1,fr);
	fwrite(&frag.xStart,sizeof(uint64_t),1,fr);
	fwrite(&frag.yStart,sizeof(uint64_t),1,fr);
	fwrite(&frag.xEnd,sizeof(uint64_t),1,fr);
	fwrite(&frag.yEnd,sizeof(uint64_t),1,fr);
	fwrite(&frag.length,sizeof(uint64_t),1,fr);
	fwrite(&frag.ident,sizeof(uint64_t),1,fr);
	fwrite(&frag.score,sizeof(uint64_t),1,fr);
	fwrite(&frag.similarity,sizeof(float),1,fr);
	fwrite(&frag.seqX,sizeof(uint64_t),1,fr);
	fwrite(&frag.seqY,sizeof(uint64_t),1,fr);
	fwrite(&frag.block,sizeof(int64_t),1,fr);
	fwrite(&frag.strand,sizeof(char),1,fr);
}


/**
 * Function to read a fragment from the specified file
 */
inline void readFragment(FragFile *frag, FILE *f){
	if(fread(&frag->diag, sizeof(int64_t), 1, f)!=1){
		if(feof(f))return;
		fprintf(stderr,"Error reading the HSP diagonal\n");
	}
	if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP xStart\n");
	}
	if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP yStart\n");
	}
	if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP xEnd\n");
	}
	if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP yEnd\n");
	}
	if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP length\n");
	}
	if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP identities\n");
	}
	if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP score\n");
	}
	if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP similarity\n");
	}
	if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP seqX\n");
	}
	if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP seqY\n");
	}
	if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
		fprintf(stderr,"Error reading the HSP block\n");
	}
	frag->strand = fgetc(f);
}


/* This function is used to comapre two fragments. The order of comparation is:
 *     1 - SeqX
 *     2 - SeqY
 *     3 - diag
 *     4 - xStart
 *     5 - length
 *     6 - yStart
 *     7 - block
 *     8 - strand
 *     9 - score
 *  @param f1 fragment to be compared.
 *  @param f2 fragment to be compared.
 *  @return negative number if f1<f2, positive number if f1>f2 or zero if f1=f2.
 */
int fragmentComparator(FragFile f1, FragFile f2){
	if(f1.seqX > f2.seqX) return 1;
	else if(f1.seqX < f2.seqX) return -1;

	if(f1.seqY > f2.seqY) return 1;
	else if(f1.seqY < f2.seqY) return -1;

	if(f1.diag > f2.diag) return 1;
	else if(f1.diag < f2.diag) return -1;

	if(f1.xStart > f2.xStart) return 1;
	else if(f1.xStart < f2.xStart) return -1;

	if(f1.length > f2.length) return 1;
	else if(f1.length < f2.length) return -1;

	if(f1.yStart > f2.yStart) return 1;
	else if(f1.yStart < f2.yStart) return -1;

	if(f1.block > f2.block) return 1;
	else if(f1.block < f2.block) return -1;

	if(f1.strand == 'f' && f2.strand == 'r') return 1;
	else if(f1.strand == 'r' && f2.strand == 'f') return -1;

	if(f1.score > f2.score) return 1;
	else if(f1.score < f2.score) return -1;

	return 0;
}