/*
 *
 * Sintax: ./newFragToBalazsVersion fragsFILE fragsFILE.out 
 * 
 * convert the actual fragment file with the new structure to
 * the version of Balazs program
 *
 *                           Oscar Torreno>oscart@uma.es
 *  -----------------------------------------------Oct 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "comparisonFunctions.h"

#define SCORE   4  // it depends on the score matrix
#define MAXLENGTH  1000

struct Fragmento {
	uint32_t diag; 
	uint32_t xIni; 
	uint32_t yIni; 
	uint32_t xFin; 
	uint32_t yFin; 
	uint32_t length; 
	uint32_t ident; 
	uint32_t score; 
	float similarity; 
	uint32_t seqX; //sequence number in the 'X' file 
	uint32_t seqY; //sequence number in the 'Y' file 
	int32_t block; //synteny block id 
	char strand; //'f' for the forward strain and 'r' for the reverse
};

void convert(struct FragFile in, struct Fragmento *out);

int main(int ac, char** av) {
	FILE* fFrags, *fFragsOld;
	uint64_t n1, n2, nFrags = 0;
	uint32_t n11, n22;
	struct FragFile frag;
	struct Fragmento F;

	if (ac != 3){
		fprintf(stderr,"USE: ./newFragsToBalazsVersion fragsFILE fragsFILE.out\n");
		exit(-1);
	}

	// prepared for multiple files
	if ((fFrags = fopen(av[1], "rb")) == NULL){
		fprintf(stderr,"Opening Frags binary input file\n");
		exit(-1);
	}

	if ((fFragsOld = fopen(av[2], "wb")) == NULL){
		fprintf(stderr,"Opening Frags binary output file\n");
		exit(-1);
	}

	readSequenceLength(&n1, fFrags);
	readSequenceLength(&n2, fFrags);
	n11=(uint32_t)n1;
	n22=(uint32_t)n2;
	fwrite(&n11, sizeof(uint32_t), 1, fFragsOld);
	fwrite(&n22, sizeof(uint32_t), 1, fFragsOld);
	fprintf(stderr, "working with fragsFile=%s SeqX=%" PRIu64 " seqY=%" PRIu64 "\n", av[1], n1,
			n2);

	readFragment(&frag, fFrags);
	while (!feof(fFrags)) {
		convert(frag, &F);
		fwrite(&F, sizeof(struct Fragmento), 1, fFragsOld);
		nFrags++;
		readFragment(&frag, fFrags);
	}
	fprintf(stdout, "nFrags:%" PRIu64 "\n", nFrags);

	fclose(fFrags);
	fclose(fFragsOld);

	return 0;
}

void convert(struct FragFile in, struct Fragmento *out){
	out->diag=(uint32_t)in.diag;
	out->xIni=(uint32_t)in.xStart;
	out->yIni=(uint32_t)in.yStart;
	out->xFin=(uint32_t)in.xEnd;
	out->yFin=(uint32_t)in.yEnd;
	out->length=(uint32_t)in.length;
	out->ident=(uint32_t)in.ident;
	out->score=(uint32_t)in.score;
	out->similarity=(float)in.similarity;
	out->seqX=(uint32_t)in.seqX;
	out->seqY=(uint32_t)in.seqY;
	out->block=(int32_t)in.block;
	out->strand=(char)in.strand;
}
