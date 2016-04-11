/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdint.h> // unit64_t ...
#include <inttypes.h> 
#include <string.h> // String functions
#include <ctype.h>
#include <arpa/inet.h>

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
void readFragment(FragFile*,FILE*);
void writeFragment(FragFile,FILE*);
int fragmentComparator(FragFile,FragFile);
void readSequenceLength(uint64_t*,FILE*);
void writeSequenceLength(uint64_t*,FILE*);
void endianessConversion(char*,char*,int);

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

	uint64_t metagLength,genoLength;
	// Read headers
	readSequenceLength(&metagLength,forward);
	readSequenceLength(&genoLength,forward);

		writeSequenceLength(&metagLength,out);
		writeSequenceLength(&genoLength,out);

	readSequenceLength(&metagLength,reverse);
	readSequenceLength(&genoLength,reverse);


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
void writeFragment(FragFile frag, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//Big endian
		fwrite(&frag.diag, sizeof(int64_t), 1, f);
		fwrite(&frag.xStart, sizeof(uint64_t), 1, f);
		fwrite(&frag.yStart, sizeof(uint64_t), 1, f);
		fwrite(&frag.xEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag.yEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag.length, sizeof(uint64_t), 1, f);
		fwrite(&frag.ident, sizeof(uint64_t), 1, f);
		fwrite(&frag.score, sizeof(uint64_t), 1, f);
		fwrite(&frag.similarity, sizeof(float), 1, f);
		fwrite(&frag.seqX, sizeof(uint64_t), 1, f);
		fwrite(&frag.seqY, sizeof(uint64_t), 1, f);
		fwrite(&frag.block, sizeof(int64_t), 1, f);
		fputc(frag.strand, f);
	} else {
		//Little endian
		endianessConversion((char *)(&frag.diag), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		endianessConversion((char *)(&frag.xStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.yStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.xEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.yEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.length), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.ident), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.score), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.similarity), tmpArray, sizeof(float));
		fwrite(tmpArray, sizeof(float), 1, f);
		endianessConversion((char *)(&frag.seqX), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.seqY), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag.block), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		fputc(frag.strand, f);
	}
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
		}
		if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP xStart");
		}
		if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP yStart");
		}
		if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP xEnd");
		}
		if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP yEnd");
		}
		if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP length");
		}
		if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP identities");
		}
		if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP score");
		}
		if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP similarity");
		}
		if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP seqX");
		}
		if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP seqY");
		}
		if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP block");
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
		}
		endianessConversion(tmpArray, (char *)(&frag->xStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP yStart");
		}
		endianessConversion(tmpArray, (char *)(&frag->yStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP xEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->xEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP yEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->yEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP length");
		}
		endianessConversion(tmpArray, (char *)(&frag->length), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP identity");
		}
		endianessConversion(tmpArray, (char *)(&frag->ident), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP score");
		}
		endianessConversion(tmpArray, (char *)(&frag->score), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(float), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP float");
		}
		endianessConversion(tmpArray, (char *)(&frag->similarity), sizeof(float)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP seqX");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqX), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP seqY");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqY), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
			fprintf(stderr,"Error reading the HSP block");
		}
		endianessConversion(tmpArray, (char *)(&frag->block), sizeof(int64_t)); 
		frag->strand = fgetc(f);
	}
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


/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		if(fread(length, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading sequence length");
		}
	} else {
		//little endian
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			fprintf(stderr,"Error reading sequence length");
		}
		endianessConversion(tmpArray, (char *)length, sizeof(uint64_t));
	}
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		fwrite(length, sizeof(uint64_t), 1, f);
	} else {
		//little endian
		endianessConversion((char *)length, tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
	}
}


/*
 */
void endianessConversion(char *source, char *target, int numberOfBytes){
	int i,j;
	for(i=numberOfBytes-1;i>=0;i--){
		j=numberOfBytes-1-i;
		target[j]=source[i];
	}
}