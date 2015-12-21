#include <stdio.h>
#include <unistd.h> 

int exists(char*);

int main(int ac, char** av){
	if(ac!=6){
		fprintf(stderr,"USE: %s file R L newTime oldTime\n",av[0]);
		return -1;
	}


	FILE* out;
	fprintf(stderr, "%s\n", av[1]);
	if(exists(av[1])){
		if((out = fopen(av[1],"at"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[1]);
			return -1;
		}
	}else if((out = fopen(av[1],"wt"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[1]);
			return -1;
		}


	// NumReads Length TimeNew TimeOld
	fprintf(out, "%s\t%s\t%s\t%s\n", av[2],av[3],av[4],av[5]);
	fclose(out);
	return 0;
}

int exists(char *file){
	if(access(file,F_OK) != (-1)) return 1;
	else return 0;
}