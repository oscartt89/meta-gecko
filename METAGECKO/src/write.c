#include <stdio.h>
#include <unistd.h> 

/*
 * USE: write FileOut arg1 arg2 ...
 */

int exists(char*);

int main(int ac, char** av){
	FILE* out;
	if(exists(av[1])){
		if((out = fopen(av[1],"at"))==NULL){
			fprintf(stderr, "Error opening file[%s]\n", av[1]);
			return -1;
		}
	}else if((out = fopen(av[1],"wt"))==NULL){
		fprintf(stderr, "Error opening file[%s]\n", av[1]);
		return -1;
	}

	int i;
	for(i=2;i<ac;++i){
		fprintf(out, "%s\t", av[i]);
	}
	fprintf(out, "\n");
	fclose(out);
	return 0;
}

int exists(char *file){
	if(access(file,F_OK) != (-1)) return 1;
	else return 0;
}