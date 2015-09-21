/*


extractSeq.c


*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

// Extract a piece of seq and stores in folder/seqx.i.fasta
void extractSeq(char* file ,char* newFile,char* folder,int ini,int fin);

int validc(char c);

//av[2]= folder
//av[3]= seqX
int main(int ac,char** av){

	if(ac!=5){
		printf("*****Use: extractSubSeq fasta newfile ini fin \n");
		printf("%d\n",ac);
		exit(-1);
	}
	
	char* file=av[1];
	char* newfile=av[2];
	int ini,fin;
	ini=atoi(av[3]);
	fin=atoi(av[4]);
	
	
	// Open in file
	FILE* seq;
	if((seq=fopen(file,"r"))==NULL){
		printf("***ERROR Opening seq input file %s\n",file);
		exit(-1);
	}
	
	int temp;
	if(ini > fin){
		temp = ini;
		ini= fin;
		fin = temp;
	}
	
	
	FILE* newSeq;
	if((newSeq=fopen(newfile,"wt"))==NULL){
		printf("***ERROR Opening output file %s\n",newfile);
		exit(-1);
	}
	printf("opened: %s\t file new:%s\n",file,newfile);
	printf("%d-%d\n",ini,fin);
	
	// skip first line and copy
	char c;
	char linea[1000];
	fgets(linea,1000,seq);
	fprintf(newSeq,"%s",linea);
	//getchar();
	
	// Go to ini point
	int i=0;
	c=getc(seq);
	while(i<ini && !feof(seq)){
		//printf("primer bucle\n");
		while(!validc(c) && !feof(seq)){
			//printf("not valid\n");
			c=getc(seq);
			
		}
		//printf("%c",c);
		c=getc(seq);
		i++;
	}
	// copy
	int n=0;
	//printf("%d-%d\n",i,fin);
	while(i<=fin && !feof(seq)){
		
		while(!validc(c) && !feof(seq)){
			c=getc(seq);
			//putc(c,newSeq);
		}
		i++;
		//printf("put- %c\n",c);
		if(validc(c)){
			putc(c,newSeq);
			n++;
		}
		
		c=getc(seq);
		
		
	
	}
	
	
	
	
	
	
	fclose(newSeq);
	fclose(seq);
	return 0;
}


/*

validc.c

*/
int validc(char c){

	if( c=='A' || c== 'C' || c== 'G' || c== 'T' || c=='N')return 1;
	return 0;
}
