/*

postprocessingFunctions.h


*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

// Extract a piece of seq and stores in folder/seqx.i.fasta
void extractSeq(char* file ,char* name,char* folder, int id,char c,int ini,int fin,int signo);

int validc(char c);
