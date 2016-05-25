/* words.c
 new version of "words" program. Instead of working with the bin file
 this version works over the plain-text sequence (and do not uses the
 seqio library -for managing big data sets-)
 Problems have been detected in the previous version.

 Using this option makes unnecesary to "format" the sequence (moving
 to binary code). Thus, after masking the Low-Complex regions this
 program can be used as next step

 This new uses "./words seq.IN words.OUT
 where seq.IN is a plain-text sequence
 words.OUT is a binary file of "wentry" structures

 NEXT: this program load the full seq into memory. Need to be modified
 to load partial chunks of sequence (not difficult)
 -----------------------------------------------------4.Feb.2012
 ortrelles @ uma.es
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

static int WORD_SIZE=32;
//static int BITS_PER_BASE=2;
static int BYTES_IN_WORD=8;//(int)ceil(WORD_SIZE/8.*BITS_PER_BASE);

void decodify(word*,char*,int);
void shift_seq(char* seq,int wsize);

void shift_word(word * w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}

void main_FILE(char * inFile, char * outFile){
	FILE *f;
	FILE *f2;
	char c;

	if ((f=fopen(inFile,"rt"))==NULL){
	perror("opening sequence file");
	}
	if ((f2=fopen(outFile,"wb"))==NULL) {
		terror("opening OUT sequence Words file");
	}
	
	c=fgetc(f);
	while(c!='\n'){
		c=fgetc(f);
	}

	wentry WE,WErev;
	WE.strand = 'f';
	WErev.strand = 'r';
	WE.seq=0;
	char sequence[4*BYTES_IN_WORD];
	unsigned long index=0;
	unsigned long inEntry=0;
	unsigned long NW=0;
	unsigned long Tot=0;
	unsigned long NoACGT=0;
	unsigned long NoC=0;
	c=fgetc(f);
	while(!feof(f)){
		if (!isupper(toupper(c))){
			if(c=='>'){
				c = fgetc(f);
				while (c != '\n')
					c = fgetc(f);
				WE.seq++;
				inEntry=0;
				index++;
			}
			NoC++;
			c=fgetc(f);
			continue;
		}
		shift_word(&WE.w);
		switch (c) {
			case 'A': inEntry++; break;
			case 'C':
				WE.w.b[BYTES_IN_WORD-1]|=1;
				inEntry++;
				break;
			case 'G':
				WE.w.b[BYTES_IN_WORD-1]|=2;
				inEntry++;
				break;
			case 'T':
				WE.w.b[BYTES_IN_WORD-1]|=3;
				inEntry++;
				break;
			default :
				inEntry=0; NoACGT++; break;
		}
		///////////////////////////// REVERSE /////////////////////////////
		shift_seq(&sequence[0],4*BYTES_IN_WORD);
		sequence[0]=c;
		///////////////////////////// REVERSE /////////////////////////////
		index++;
		Tot++;
		if(inEntry>=(unsigned long)WORD_SIZE){
			WE.pos=index-WORD_SIZE;
			NW++;
			fwrite(&WE,sizeof(wentry),1,f2);
			///////////////////////////// REVERSE /////////////////////////////
			//Copy reverse sequence
			//decodify(&WE.w,&sequence[0],8);
				// Store reverese
				int i;
				for(i=0;i<(4*BYTES_IN_WORD);++i){
					shift_word(&WErev.w);
					switch (sequence[i]) {
						case 'A': break;
						case 'C':
							WErev.w.b[BYTES_IN_WORD-1]|=1;
							break;
						case 'G':
							WErev.w.b[BYTES_IN_WORD-1]|=2;
							break;
						case 'T':
							WErev.w.b[BYTES_IN_WORD-1]|=3;
							break;
						default :
							break; //Maybe an error message here
					}
				}
			//Update values
			WErev.pos=WE.pos;
			WErev.seq=WE.seq;
			NW++;
			fwrite(&WErev,sizeof(wentry),1,f2);
			///////////////////////////// REVERSE /////////////////////////////
		}
		c=fgetc(f);
	}
	//printf("FILE: Create %d Words --(seqLen=%d NoACGT=%d noChar=%d\n",NW,Tot,NoACGT, NoC);
	fclose(f);

}

int main(int ac, char** av){
	if(ac!=3){
		terror("USE: words seqFile.IN words.OUT");
	}
	main_FILE(av[1], av[2]);
	return 0;
}

/* This function is used to decodify binary sequence to char sequence.
 *   @param w is the word instance to be decodified.
 *   @param seq is the char array where the sequence will be stored. Must have 
 *          enough space to store the sequence.
 *   @param wsize is the size of the word.b char array.
 */
void decodify(word* w, char *seq,int wsize){
	char Alf[] = { 'A', 'C', 'G', 'T' };
	int i;
	unsigned char c;
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		seq[4*i] = Alf[(int) c];
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		seq[4*i+1] = Alf[(int) c];
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		seq[4*i+2] = Alf[(int) c];
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		seq[4*i+3] = Alf[(int) c];
	}
}


void shift_seq(char* seq,int wsize){
	int i;
	for(i=wsize-1;i>0;--i)
		seq[i]=seq[i-1];
}