/* leehd read and displays the hash table from disk 
   Syntax: leehd prefixNameOUT
    prefixNameOUT.h2dW  : index of words-Pos-Ocurrences
    prefixNameOUT.h2dP  : positions
    both must be available
    Any char as third argument means "Verbose mode"
    Feb.2012: computes word frequencies	
                              ortrelles@uma.es / Dic.2011
    ---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <inttypes.h>

#include "metag_common.h"


#define PEQ 1001

int quickReadWord(FILE * f, uint32_t * nRep, Word * w, uint64_t * position, uint16_t WORD_LENGTH){

	size_t tot = fread(w->b, 1, WORD_LENGTH/4, f);
	if(tot != WORD_LENGTH/4){
		if(feof(f)){
			return 0;
		}
		terror("Error reading reading word");
	}
	if(fread(position, sizeof(uint64_t), 1, f) != 1) terror("Error reading position");
	if(fread(nRep, sizeof(uint32_t), 1, f) != 1) terror("Error reading repetitions");
	return 0;
}

int main(int ac, char** av){

	char fname[1024], *W;
	W=(char *)malloc(33*sizeof(char));
	FILE *f1, *f2;

	uint16_t WORD_LENGTH;	
	uint32_t nRep;
	uint64_t position;
	Word w;
	LocationEntry lc;

	uint64_t i=0;
	uint64_t totalWords = 1;
	int flagV=0, flagGetchar = 0;

	if(ac<2)terror("USE: leehd  prefixNameOUT [v={1 for quick, 2 for paused}]\n");
	if (ac==3) flagV=1;
	if(atoi(av[2]) == 2) flagGetchar = 1;
	sprintf(fname,"%s.d2hW",av[1]); // Words file (first level of hash table)
	if ((f1 = fopen(fname,"rb"))==NULL) terror("opening prefix.h2dW file");
	sprintf(fname,"%s.d2hP",av[1]); // Positions file
	if ((f2 = fopen(fname,"rb"))==NULL) terror("opening prefix.h2dP file");


	//Read word length
	if(fread(&WORD_LENGTH, sizeof(uint16_t), 1, f1) != 1) terror("Could not read word length");
	fprintf(stdout, "Word length is %"PRIu16"\n", WORD_LENGTH);

        w.b = (unsigned char*) malloc(WORD_LENGTH/4 * sizeof(unsigned char));


	// kick-off
	quickReadWord(f1, &nRep, &w, &position, WORD_LENGTH);
       
        while(!feof(f1)){

             if (flagV) {showWord(&w, W, WORD_LENGTH);fprintf(stdout, "%.32s", W);}
             if (flagV) fprintf(stdout,"  : pos=%-7" PRIu64 " num=%-7" PRIu32 ":",position, nRep);


	     for(i=0;i<nRep;i++){
		if(fread(&lc,sizeof(LocationEntry),1,f2)!=1){
                        terror("Error reading the word occurrences");
		}
               	fprintf(stdout,"(%c, %" PRIu32 ",%" PRIu64 ") ",lc.strand, lc.seq, lc.pos);
                if (i>10) {fprintf(stdout,"...cont"); break;}

	     }
		fprintf(stdout, "\n");
		if(flagGetchar) getchar();
		
	    	if(!feof(f1) && quickReadWord(f1, &nRep, &w, &position, WORD_LENGTH) == 0) totalWords++;
        }
	free(W);

	fclose(f1);
	fclose(f2);
	fprintf(stdout, "Total words read: %"PRIu64"\n", totalWords);
	exit(0);
}

