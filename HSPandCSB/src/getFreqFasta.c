/*
Lee una secuencia y extrae la frecuencia relativa de ACGT
Reads a sequence and calcs ACGT frequencies

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <time.h>


int main(int ac,char** av){



	char c;
	
	int i,j,n;

	

	if(ac<2){
	printf("\n Uso:\n\tprobACGT entrada salida \n");exit(-1);
	}

	j=0;
	// Open files.
	FILE *fe,*fs;
	if((fe=fopen(av[1],"r"))==NULL){
		printf("Error\n");
	}
	if((fs=fopen(av[2],"w"))==NULL){
		printf("Error\n");
	}
	
	// Hay que saltar la primera línea e insertarla en el archivo salida
	// first line

	c=getc(fe);
	while(c!='\n') {

		c=getc(fe);
	}



	i=0;
	int A,C,G,T;
	float pA,pC,pG,pT;
	A=C=G=T=0;
	while(!feof(fe)){
		c=getc(fe);

		//printf("%c.%d\n",c,i);

		switch(c){
			case 'A':
        A++; i++;
        break;	
			case 'C':
        C++; i++;
        break;
			case 'G':
        G++; i++;
        break;
			case 'T':
        T++; i++;
        break;
			default:
				break;	
		}
	}
n=i;
pA=(float)A/i;
pC=(float)C/i;
pG=(float)G/i;
pT=(float)T/i;


	
	fclose(fe);

	fprintf(fs,"A\t%f\nC\t%f\nG\t%f\nT\t%f\n",pA,pC,pG,pT);
	fclose(fs);


return 0;


}
