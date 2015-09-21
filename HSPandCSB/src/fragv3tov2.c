/*
 * fragv3tov2.c
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 *	E-mail: jarjonamedina@uma.es
 *
 *	



We change the current fragment structure (v3) into original one (v2)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <dirent.h>
#include <math.h>
#include <inttypes.h>


#include "fragmentv3.h"
#include "fragmentv2.h"
#include "structs.h"
#include "comparisonFunctions.h"




int  main(int ac,char** av){

	if(ac<4){
	printf("\n Uso:\n\tfragv3tov2 entrada.fragv3 salida.fragv2 salida.csv\n");exit(-1);
	}

	
// Read fragments
struct FragFile* f;

//int nf;
//uint64_t xtotal,ytotal;
//int xtotal,ytotal;


		/* Read Fragments */
		Fragmentv3* fv3;
		int nf;
		int xtotal,ytotal;
		
		fv3=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);
	

fv3=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);
printf("nf:%d\n",nf);

f=(struct FragFile*)malloc(sizeof(struct FragFile)*nf);

int i;






for(i=0;i<nf;i++){
	
	f[i].diag=(unsigned long)fv3[i].diag;
    f[i].xStart=(unsigned long)fv3[i].xIni;
    f[i].yStart=(unsigned long)fv3[i].yIni;
    f[i].xEnd=(unsigned long)fv3[i].xFin;
    f[i].yEnd=(unsigned long)fv3[i].yFin;
 
	f[i].length=(unsigned long)fv3[i].length;
    f[i].ident=(unsigned long)fv3[i].ident;
    f[i].score=(unsigned long)fv3[i].score;
    f[i].similarity=fv3[i].similarity;
    f[i].seqX=(unsigned long)fv3[i].seqX; //sequence number in the 'X' file
    f[i].seqY=(unsigned long)fv3[i].seqY; //sequence number in the 'Y' file
    f[i].block=fv3[i].block;          //synteny block id
    f[i].strand=(char)fv3[i].strand;        //'f' for the forward strain and 'r' for the reverse
	
	
	printf("%ld\t%ld\t%ld\t%ld\t%c\t%d\n",fv3[i].xIni,fv3[i].yIni,fv3[i].xFin,fv3[i].yFin,fv3[i].strand,fv3[i].block);

	
}

    
	FILE* fs;
	if((fs=fopen(av[2],"wb"))==NULL){
		printf("***ERROR Opening output file %s\n",av[4]);
		exit(-1);
	}
	FILE* ftxt;
	if((ftxt=fopen(av[3],"w"))==NULL){
		printf("***ERROR Opening output file %s\n",av[5]);
		exit(-1);
	}
	
	
	
	
	uint64_t xtotal64,ytotal64;
	xtotal64 = xtotal;
	ytotal64 = ytotal;
	
	writeSequenceLength(&xtotal64,fs);
	writeSequenceLength(&ytotal64,fs);
	
	for(i=0;i<nf;i++){
		writeFragment(&f[i], fs);
		fprintf(ftxt,"%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%c\t%"PRIu64"\n",f[i].xStart,f[i].yStart,f[i].xEnd,f[i].yEnd,f[i].strand,f[i].block);
		
	}
	fclose(fs);
	fclose(ftxt);

	return 0;

}

