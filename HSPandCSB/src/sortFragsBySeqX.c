#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksort.h"
#include "comparisonFunctions.h"

int GT(BaseType a1, BaseType a2) {
	if (a1.seqX > a2.seqX)
		return 1;
	else if(a1.seqX < a2.seqX)
		return 0;
	if (a1.diag > a2.diag)
		return 1;

	return 0;
}

int main(int ac, char** av) {
	FILE *fi, *fo;
	char tmp[1024];
	uint64_t nx, ny;
	struct FragFile F;

	if (ac < 4) {
		printf("USE: sortFragsBySeqX <max_size> <num_proc> <input_file> <output_file>\n");
		exit(1);
	}

	sprintf(tmp, "%s.tmp", av[3]);

	if((fi=fopen(av[3],"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}

	readSequenceLength(&nx, fi);
	readSequenceLength(&ny, fi);

	if((fo=fopen(tmp,"wb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}

	readFragment(&F, fi);
	while(!feof(fi)){
		fwrite(&F, sizeof(struct FragFile), 1, fo);
		readFragment(&F, fi);
	}

	fclose(fi);
	fclose(fo);

	mypsort(atoi(av[1]), atoi(av[2]), tmp, av[4]);

	unlink(tmp);

	if((fi=fopen(av[4],"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}

	sprintf(tmp, "%s.tmp", av[4]);
	if((fo=fopen(tmp,"wb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}

	writeSequenceLength(&nx, fo);
	writeSequenceLength(&ny, fo);

	fread(&F, sizeof(struct FragFile), 1, fi);
	while(!feof(fi)){
		writeFragment(&F, fo);
		fread(&F, sizeof(struct FragFile), 1, fi);
	}

	fclose(fi);
	fclose(fo);

	unlink(av[4]);
	rename(tmp, av[4]);

	return 0;
}

