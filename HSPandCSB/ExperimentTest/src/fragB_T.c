#include <inttypes.h>
#include <stdio.h> //IO stream

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
}FragFile;

int main(int ac, char** av){
	if(ac!=3){
		fprintf(stderr,"USE: %s fragFile out\n",av[0]);
		return -1;
	}

	FILE *out, *in;

	if((in = fopen(av[1],"rb"))==NULL){
			fprintf(stderr, "MAIN:: Error opening frags file. [%s]\n", av[1]);
			return -1;
	}
	if((out = fopen(av[2],"wt"))==NULL){
		fprintf(stderr, "MAIN:: Error opening out file. [%s]\n", av[2]);
		return -1;
	}

	FragFile frag;

	//Read first
	fread(&frag,sizeof(FragFile),1,in);
	while(!feof(in)){
		fprintf(out, "Fragment: X=%" PRIu64 "  Y=%" PRIu64 "\n\tXStart:%" PRIu64 "\tXEnd:%" PRIu64 "\n\tYStart:%" PRIu64 "\tYEnd:%" PRIu64 "\n\tL:%" PRIu64 "\tS:%f\tSc:%" PRIu64 "\n",
			frag.seqX,frag.seqY,frag.xStart,frag.xEnd,frag.yStart,frag.yEnd,frag.length,frag.similarity,frag.score);
		fread(&frag,sizeof(FragFile),1,in);
	}
}