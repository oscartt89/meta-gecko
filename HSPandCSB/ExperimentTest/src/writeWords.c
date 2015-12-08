#include <stdio.h>
#include <stdint.h> //unit64_t
#include <inttypes.h>


typedef struct {
	//Each letter is stored using 2 bits
	//We have 4 letters per byte and a
	//maximum of 32 in 'b'
    unsigned char b[8];
} word;

typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint16_t num;
} hashentry;

typedef struct {
    // Index of read
    uint32_t readIndex;
    // Position on word dictionary
    uint64_t pos;
    // Number of different kmers
    uint16_t num;
} Read; // "read" is used on dirent.h

void showWord(word* w, char *ws);

int main(int ac, char** av){
	if(ac!=5){
		fprintf(stderr,"USE: %s PD WD RD out \n",av[0]);
		return -1;
	}

	FILE *P,*R,*W, *out;
	Read re;
	hashentry he;
	uint64_t pos;


	// Open file
	if((P = fopen(av[1],"rb"))==NULL){
		fprintf(stderr, "Error opening P file. [%s]\n", av[1]);
		return -1;
	}
	if((W = fopen(av[2],"rb"))==NULL){
		fprintf(stderr, "Error opening W file. [%s]\n", av[2]);
		return -1;
	}
	if((R = fopen(av[3],"rb"))==NULL){
		fprintf(stderr, "Error opening R file. [%s]\n", av[3]);
		return -1;
	}

    if((out = fopen(av[4],"wt"))==NULL){
        fprintf(stderr, "Error opening OUT file. [%s]\n", av[3]);
        return -1;
    }

    fread(&re,sizeof(Read),1,R); // Read first time
    while(!feof(R)){
    	fprintf(stdout, "Read: %" PRIu32 " (%" PRIu16 ")\n", re.readIndex,re.num);
    	fseek(W,re.pos,SEEK_SET);
    	
    	int i;
    	for(i=0;i<re.num;++i){
    		fread(&he,sizeof(hashentry),1,W);
			fseek(P,he.pos,SEEK_SET);

			int j;
			char seq[64];
			for(j=0;j<he.num;++j){
				fread(&pos,sizeof(uint64_t),1,P);
				fprintf(out, "W::\n");
				showWord(&he.w,seq);
				fprintf(out, "\tWord: %s",seq); 
				
				fprintf(out, "\n\tRead:%" PRIu32 "\tPos:%" PRIu64 "\n",re.readIndex, pos);

			}    		
    	}
    	fread(&re,sizeof(Read),1,R);
    }


}



void showWord(word* w, char *ws) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	int i;
	int wsize = 8;
	unsigned char c;
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
	}
}