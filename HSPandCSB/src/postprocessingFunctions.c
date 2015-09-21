/*


Common functions


*/
#include "postprocessingFunctions.h"

//av[2]= folder
//av[3]= seqX
void extractSeq(char* file,char*name,char* folder , int id,char tipo,int ini,int fin,int signo){

	// Open in file
	FILE* seq;
	if((seq=fopen(file,"r"))==NULL){
		printf("***ERROR Opening seq input file %s\n",file);
		exit(-1);
	}
	
	int temp;
	temp=(int)tipo;
	if(ini > fin){
		temp = ini;
		ini= fin;
		fin = temp;
	}
	
	// Open outfile
	char newfile[100];
	sprintf(newfile,"%s/%s-%d.fasta",folder,name,id);
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
	//fprintf(newSeq,"%s",linea);
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
	
	if(n>100){
		if(signo){
			fprintf(newSeq,", 1\n");
		}else{
			fprintf(newSeq,", -1\n");
		}
	}
	
	
	
	
	fclose(newSeq);
	fclose(seq);
}


/*

validc.c

*/
int validc(char c){

	if( c=='A' || c== 'C' || c== 'G' || c== 'T' || c=='N')return 1;
	return 0;
}
