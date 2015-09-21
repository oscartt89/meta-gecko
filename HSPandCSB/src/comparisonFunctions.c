#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include "structs.h"
#include "commonFunctions.h"

#define MAXBUF 10000000

int readHashEntry(hashentry *h, FILE *f, uint64_t freqThr) {
	do {
		if(fread(h, sizeof(hashentry), 1, f)!=1){
			if(ferror(f))terror("Error reading hash entry from the file");
		}
	} while (!feof(f) && h->num > freqThr);

	if (feof(f))
		return -1;

	return h->num;
}

void loadWordOcurrences(hashentry he, location** pos, FILE** f) {
	// Load word position for he.word
	if (he.num > MAXBUF) {
		*pos = (location*) realloc(pos, he.num * sizeof(location));
		if (pos == NULL)
			terror("memory for position array");
	}
	fseek(*f, he.pos, SEEK_SET);
	if(fread(*pos, sizeof(location), he.num, *f)!=he.num){
		terror("Not possible to read the number of elements described");
	}
}

unsigned long scoreMax(char *seq, char *seq2, uint64_t len, int point) {
	//CHANGE WHEN USING PAM MATRIX
	//Scor=0; for (i=0;ii<len;i++) if (Seq[i]==Seq2[i])scor+=Ptos
	/* TO DELETE WARNINGS*/
	if(seq+1){
	}

	if(seq2+1){
	}
	/**********************/

	return len * point;
}

struct Sequence* LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nStruct,
		int fAst) {
	char c;
	uint64_t lon = 0, k = 0, ns;
	uint64_t lonFinal = 0;
	struct Sequence *sX, *sX2; //sX will be the first elem. sX2 will generate all the structure

	//Initialize
	*n = 0;
	*nStruct = 0;

	//Memory
	ns = 1;
	if ((sX = (struct Sequence*) malloc(sizeof(struct Sequence))) == NULL)
		terror("Memory...");

	if (!fAst)
		while ((c = getc(f)) != '>' && !feof(f))
			; //start seq
	if (feof(f))
		return 0;
	while ((c = getc(f)) == ' ')
		;

	while (k < MAXLID && c != '\n' && c != ' ') {
		if (feof(f))
			return 0;

		sX->ident[k++] = c;
		c = getc(f);
	}

	sX->ident[k] = 0; //end of data.
	while (c != '\n')
		c = getc(f);
	c = getc(f);

	//start list with sX2
	sX2 = sX;
	while (/*c!='*'&&*/!feof(f)) {
		c = toupper(c);
		if (c == '>') {
			fAst = 1;
			sX2->datos[lon++] = '*';
			while (c != '\n') {
				if (feof(f))
					return 0;
				c = getc(f);
			}
			//break;
		}
		if (isupper(c))
			sX2->datos[lon++] = c;
		if (c == '*') {
			sX2->datos[lon++] = c;
		}
		c = getc(f);

		//Check if the length is the end of this struct
		if (lon >= MAXLS) {
			lonFinal += lon;
			lon = 0;
			ns++;
			if ((sX = (struct Sequence*) realloc(sX,
					ns * sizeof(struct Sequence))) == NULL)
				terror("Memory...");
			sX2 = sX + ns - 1;
		}
	}

	if (lon < MAXLS)
		sX2->datos[lon] = 0x00;

	lonFinal += lon;
	*nStruct = ns;
	*n = lonFinal;
	return sX;
}

char getValue(struct Sequence *s, uint64_t pos, int ns) {
	struct Sequence *aux = s;
	int nActual = 1;

	while (pos >= MAXLS) {
		aux++;
		pos -= MAXLS;
		nActual++;
		if (nActual > ns){
			terror("out of sequence.");
		}
	}

	return aux->datos[pos];
}

long getSeqLength(struct Sequence *s, uint64_t start, int ns) {
	int nActual = 1;
	struct Sequence *aux = s;
	while (start >= MAXLS) {
		aux++;
		start -= MAXLS;
		nActual++;
		if (nActual > ns)
			terror("out of sequence.");
	}
	uint64_t s1 = start;
	while (s1 > 0 && aux->datos[s1] != '*') {
		s1--;
	}
	s1++;
	char *tmp = strchr(aux->datos + s1, '*');
	if (tmp == NULL) {
		return strlen(aux->datos) - s1 + 1;
	}
	return tmp - (aux->datos + s1) + 1;
}

void endianessConversion(char *source, char *target, int numberOfBytes){
	int i,j;
	for(i=numberOfBytes-1;i>=0;i--){
		j=numberOfBytes-1-i;
		target[j]=source[i];
	}
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f){
	char tmpArray[8];

	if(htons(1)==1){
		//big endian
		if(fread(&frag->diag, sizeof(int64_t), 1, f)!=1){
			if(feof(f))return;
			terror("Error reading the HSP diagonal");
		}
		if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xStart");
		}
		if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yStart");
		}
		if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xEnd");
		}
		if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yEnd");
		}
		if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP length");
		}
		if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP identities");
		}
		if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP score");
		}
		if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
			terror("Error reading the HSP similarity");
		}
		if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqX");
		}
		if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqY");
		}
		if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
			terror("Error reading the HSP block");
		}
		frag->strand = fgetc(f);
	} else {
		//little endian
		if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
			if(feof(f))return;
			terror("Error reading the HSP diagonal");
		}
		endianessConversion(tmpArray, (char *)(&frag->diag), sizeof(int64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xStart");
		}
		endianessConversion(tmpArray, (char *)(&frag->xStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yStart");
		}
		endianessConversion(tmpArray, (char *)(&frag->yStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->xEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->yEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP length");
		}
		endianessConversion(tmpArray, (char *)(&frag->length), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP identity");
		}
		endianessConversion(tmpArray, (char *)(&frag->ident), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP score");
		}
		endianessConversion(tmpArray, (char *)(&frag->score), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(float), 1, f)!=1){
			terror("Error reading the HSP float");
		}
		endianessConversion(tmpArray, (char *)(&frag->similarity), sizeof(float)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqX");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqX), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqY");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqY), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
			terror("Error reading the HSP block");
		}
		endianessConversion(tmpArray, (char *)(&frag->block), sizeof(int64_t)); 
		frag->strand = fgetc(f);
	}
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//Big endian
		fwrite(&frag->diag, sizeof(int64_t), 1, f);
		fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
		fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
		fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag->length, sizeof(uint64_t), 1, f);
		fwrite(&frag->ident, sizeof(uint64_t), 1, f);
		fwrite(&frag->score, sizeof(uint64_t), 1, f);
		fwrite(&frag->similarity, sizeof(float), 1, f);
		fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
		fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
		fwrite(&frag->block, sizeof(int64_t), 1, f);
		fputc(frag->strand, f);
	} else {
		//Little endian
		endianessConversion((char *)(&frag->diag), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		endianessConversion((char *)(&frag->xStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->yStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->xEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->yEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->length), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->ident), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->score), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->similarity), tmpArray, sizeof(float));
		fwrite(tmpArray, sizeof(float), 1, f);
		endianessConversion((char *)(&frag->seqX), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->seqY), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->block), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		fputc(frag->strand, f);
	}
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		if(fread(length, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
		}
	} else {
		//little endian
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
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

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment(){
	return 9*sizeof(uint64_t)+2*sizeof(int64_t)+1*sizeof(float)+1*sizeof(char);
}

