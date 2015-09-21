/*
Takes fragfiles with absolute coords and it returns a merged file with the frags
with X-coords in relative sorted by read-genome-diagonal.

In seqX we store the read's number and in seqY the genome number.

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "fragmentv2.h"
#include "comparisonFunctions.h"
#include "commonFunctions.h"

#define GIGA 1073741824

struct Files {
   char *fname;
};

int howManyFiles(char *fname);
struct Files* LoadListOfFiles(char*fListName, int *nF);
int validc(char c);
uint64_t getNreads(FILE* fi);
long int* calcOffsetVector(FILE* fi,uint64_t nreads);
void saveFragments(char** av,uint64_t nreads, long int* offset);

struct FragFile* readFragmentsOfRead(char* s,int* nf,uint64_t *xtotal,uint64_t *ytotal, uint64_t nread, long int offset, long int *offsetNext){
        FILE* fe;
        struct FragFile *fs, *more_fs;
	struct FragFile f;
        int n;
	int fragsAllocated=0;

        if((fe=fopen(s,"rb"))==NULL){
                printf("***ERROR Opening input file");
                exit(-1);
        }

        n=0;

	readSequenceLength(xtotal,fe);
	readSequenceLength(ytotal,fe);

        fs=(struct FragFile*)malloc(1000*sizeof(struct FragFile));
        if(fs==NULL){
                printf("****ERROR: Out of memory loading frags\n");
                exit(-1);
        }

	fragsAllocated = 1000;

	//fprintf(stdout, "Going to offset: %ld in file: %s\n", *offset, s);
	fseek(fe, offset, SEEK_SET);

        n=0;
        readFragment(&f, fe);
        while(!feof(fe)&&f.seqX<=nread){
		if(f.seqX==nread){
			memcpy(&fs[n], &f, sizeof(struct FragFile));
			n++;
			if(n>=fragsAllocated){
				fragsAllocated+=1000;
				more_fs = (struct FragFile*)realloc(fs,fragsAllocated*sizeof(struct FragFile));
				if(more_fs != NULL){
					fs=more_fs;
				} else {
					free(fs);
					printf("****ERROR: Out of memory loading frags\n");
					exit(-1);
				}
			}
		}
		*offsetNext=ftell(fe)-sizeofFragment();
        	readFragment(&f, fe);
        }

        *nf=n;
        fclose(fe);
        return fs;
}


int howManyFiles(char *fname){
        FILE *f;
        char line[1024];
        int nL=0;

        if ((f=fopen(fname,"rt"))==NULL){
            fprintf(stdout,"Opening frags file LIST");
	    exit(-1);
	}

        fgets(line,1024,f);
        while(!feof(f)) {
                if (line[0]!='#' && (int)strlen(line)!=0) nL++;
                fgets(line,1024,f);
        }
        fclose(f);
        return nL;
}

// Load to memory a list of files --------------------------------------------------
// datafile format: fileName[tab]nSeq[tab][format][newLINE]

struct Files* LoadListOfFiles(char*fListName, int *nF) {

        FILE *f, *ff;
        struct Files*L=NULL;
        char line[1024];
        int N=0,nFiles;
        uint64_t xnx,xny;

        nFiles = howManyFiles(fListName);
        if ((L=(struct Files*) calloc(nFiles,sizeof(struct Files)))==NULL)
           terror("memory for list of files");

        if ((f=fopen(fListName,"rt"))==NULL)
                terror("opening List of Files");

        fgets(line,1024,f);
        while(!feof(f)) {

                if (line[(int)strlen(line)-1]=='\n') line[(int)strlen(line)-1]=0x00;

                if (line[0]!='#' && (int)strlen(line)>2) {
                        L[N].fname = strdup(line);

                        if ((ff=fopen(L[N].fname,"rt"))==NULL) {
                                //fprintf(stderr,"abriendo -->%s<--\n",L[N].fname);
                                   // rterror("Opening frags file from LIST X X XX");
                                continue;
                        }

                        readSequenceLength(&xnx, ff);
                        readSequenceLength(&xny, ff);

                        N++;
                        fclose(ff);
                }
                fgets(line,1024,f);
        }
        fclose(f);
        (*nF) = nFiles;
        return L;
}


int main(int ac,char** av){
	
	if(ac!=4){
	printf("\n Uso:\n\tmergeMetagenomeFrags fragFilesList.txt seqX.fasta outputPrefixName \n");exit(-1);
	}

	// Open fastas file.
	FILE *fX;
	if((fX=fopen(av[2],"r"))==NULL){
		printf("Error\n");
	}
	
	// Read fasta file and get number of reads
	uint64_t  nreads=0;
	long int* offsetX;
	
	nreads= getNreads(fX);
	offsetX = calcOffsetVector(fX,nreads);
	
	saveFragments(av,nreads,offsetX);
	
return 0;
}


/* programs */

void saveFragments(char** av,uint64_t nreads, long int* offset){

/***********************************************************/
	struct Files *L;
	int nF, outputSuffix=0;
	long int **offsetFrags; 


	uint64_t r, i, j, k, bytesWritten=0;

	char outputFileName[1024];
	FILE *fs;

	// Read fragments
	struct FragFile* f;
	int nf; // Number of fragments
	uint64_t nx,ny,nytotal=0,nFiles;
	nf=0;
	
	L = LoadListOfFiles(av[1], &nF);

	nFiles = nF;

	
	/* Go over fragment*********************************/
	uint64_t tmp=0;

	sprintf(outputFileName, "%s-%03d.frags", av[3], outputSuffix++);
	if((fs=fopen(outputFileName,"wb"))==NULL){
                printf("***ERROR Opening output file %s\n",outputFileName);
                exit(-1);
        }
	
	writeSequenceLength(&nx,fs);
	writeSequenceLength(&nx,fs);

	offsetFrags = (long int **)malloc((nreads+1)*sizeof(long int *));
	if(offsetFrags == NULL){
		printf("***ERROR not enough memory for offset matrix\n");
		exit(-1);
	}
	for(i=0;i<(nreads+1);i++){
		offsetFrags[i] = (long int *)malloc(nFiles*sizeof(long int));
		if(offsetFrags[i] == NULL){
			printf("***ERROR not enough memory for offset matrix\n");
			exit(-1);
		}
	}
	
	for(i=0;i<nreads;i++)
		for(r=0;r<nFiles;r++)
			offsetFrags[i][r]=2*sizeof(uint64_t);

	for(r=0;r<nreads;r++){
		//fprintf(stdout, "Read: %" PRIu64 "\n", r);
		for(i=0;i<nFiles;i++){
			//fprintf(stdout, "Genome index: %" PRIu64 "\tFile name:%s\n", i, L[i].fname);
			f=readFragmentsOfRead(L[i].fname,&nf,&nx,&ny,r,offsetFrags[r][i],&offsetFrags[r+1][i]);
			//fprintf(stdout,"Fragments for read(%" PRIu64 "): %d\n", r, nf);
			if(r==0)nytotal+=ny;
			for(j=0;j<(uint64_t)nf;j++){
				tmp=f[j].seqX; 
				
				//if(((long int)f[j].xStart-offset[tmp])<0){
				//	fprintf(stdout, "ERR** Offset not valid. xStart=%" PRIu64 " offset[%" PRIu64 "]=%ld\n",f[j].xStart,tmp,offset[tmp]);
				//}
				f[j].xStart=(uint64_t)((long int)f[j].xStart-offset[tmp]);
				f[j].xEnd=(uint64_t)((long int)f[j].xEnd-offset[tmp]);
		
				f[j].diag = f[j].xStart-f[j].yStart;

				f[j].seqY = i;
			}
			for(k=0;k<(uint64_t)nf;k++){
				writeFragment(&f[k], fs);
				bytesWritten+=sizeofFragment();
			}
	
			free(f);
		}
		if(bytesWritten>GIGA){
			bytesWritten=0;
			//Update sequence lengths
			rewind(fs);
			writeSequenceLength(&nx,fs);
			writeSequenceLength(&nytotal,fs);
			fclose(fs);
			
			//Create new file
			sprintf(outputFileName, "%s-%03d.frags", av[3], outputSuffix++);
		        if((fs=fopen(outputFileName,"wb"))==NULL){
		              printf("***ERROR Opening output file %s\n",outputFileName);
	        	      exit(-1);
			}
			//Write sequence lengths
			writeSequenceLength(&nx,fs);
			writeSequenceLength(&nytotal,fs);
		}
		
	}
	rewind(fs);
        writeSequenceLength(&nx,fs);
        writeSequenceLength(&nytotal,fs);
        fclose(fs);

}
/********************/
long int* calcOffsetVector(FILE* fi,uint64_t nreads){
	uint64_t i;
	// Create the offset array [nreads] and set to 0
	long int* offset = (long int*)malloc(sizeof(long int)*nreads);
	for(i=0;i<nreads;i++)offset[i]=0;
/*  Now we go over to fasta file calculating the offset value.*/
	int max=10000;
	char linea[max];
	char readName[max];
	char c;
	long int cont=0; // apunta a la posicion absoluta de cada letra
	long int ant=0; // guarda el ultimo valor de 'cont'
	nreads=0;
	fgets(linea,max,fi);
	while(!feof(fi)){
		if(sscanf(linea,">%s\n",readName)==1){
			//printf(">%s\n",readName);
			offset[nreads]=ant;
			c=toupper(fgetc(fi));
			while(!feof(fi)&&c!='>'){
				if(validc(c)){
					//printf("%c",c);
					cont++;
				}
				c=toupper(fgetc(fi));
			}
			//*Debug. 
			//printf("\n%s : %ld\t %ld\n",readName,offset[nreads],cont-offset[nreads]);
			//*/
			fseek(fi, -1, SEEK_CUR);
			nreads++;
		}
		
		fgets(linea,max,fi);
		cont++;// +1 Due to N between two reads
		ant=cont; 
	
		
	}
	fclose(fi);
	return offset;
}
/*********/
uint64_t getNreads(FILE* fi){
	uint64_t nreads=0;
	char linea[10000];
	char readName[1000];
	int max=10000;
	
	
	fgets(linea,max,fi);
	while(!feof(fi)){
		if(sscanf(linea,">%s\n",readName)==1){
			nreads++;
		}
		fgets(linea,max,fi);
	}
	rewind(fi);
	return nreads;
}


int validc(char c){

return ( (c=='A'||c=='C'||c=='G'||c=='T'||c=='N'));
}
