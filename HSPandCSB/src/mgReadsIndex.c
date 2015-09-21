/* mgReadsIndex.c   create an INdex for fasta access to reads contentDraw a significance profile for a 

   Syntax: ./mgReadIdex reads.file

   out: reads.file.IND 
        a tuple {p, Lacum } where p is the position (offset) or ">" 
		symbol (sequenc start position) and Lacum correspond to the 
		accumulated length of previous reads)


 * ------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>


//Struct for reads index tuple
struct rIndex {
    uint64_t pos;      //Start position of read
    uint64_t Lac;      // accumulated length
};


void terror(char *s){
	fprintf(stderr,"\n****ERR ** %s **\n",s);
	exit(-1);
}



int main(int ac, char **av) {

	FILE   *f, *fOut;
	struct rIndex R;
	int  c;
	int rLen, tLen=0, nReads=0;
	char tmp[1000];

	// ----------------------------------

	if (ac!=2) 
       terror("Syntax: ./readsIndex Reads.file");

	
	if ((f=fopen(av[1],"rt"))==NULL)
		terror("opening Reads file");

	sprintf(tmp,"%s.IND",av[1]);
	if ((fOut=fopen(tmp,"wb"))==NULL)
		terror("opening Reads.INDEX file");

	// process reads file
    fprintf(stdout,"In :%s\nOut:%s\n",av[1],tmp);

    R.Lac=0;
    c=fgetc(f);
	while (!feof(f)) {
	    rLen=0;
		
        while(c!='>') c=(char)fgetc(f);
		R.pos=ftell(f) - 1;
        while(c!='\n') c=(char)fgetc(f);
		while(c!='>' && !feof(f)) {
		   c=toupper(c);
		   if (c>64 && c<91) rLen++;
			c=(char)fgetc(f);
		}
	    if (rLen==0) terror("empty read");
		fwrite(&R,sizeof(struct rIndex),1,fOut); 
		R.Lac+=rLen + 1;
		tLen+=rLen;
		nReads++;
		if ((nReads%100000)==0) printf("> %d\n",nReads);
	}

	fclose(f);
	fclose(fOut);
	fprintf(stdout,"Number of Reads = %d\n",nReads);
	fprintf(stdout,"Tot seq.length  = %d\n",tLen);
	return 0;
}


