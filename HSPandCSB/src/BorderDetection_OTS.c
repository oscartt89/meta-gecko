/* BorderDetection  

Given tow CSB, it searches end and start areas and look for hits in these areas.
The output are tow histogram of hits in these areas. 

 
 ------------------arjona@uma.es----Oct/2014-------------*/

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "commonFunctions.h"
#include "structs.h"
#include "comparisonFunctions.h"
#include "fragmentv3.h"
#include "JAMfunctions.h"

#define POINT 4
int condicion(hit h,int xIni,int yIni,int xFin,int yFin);

int main(int ac, char **av) {


	if (ac != 7)
		terror("usage: BorderDetectionOTS Fragments.csb.frag Hits-f.File Hits-r.File CBS1 CSB2 output.csv");

	// Get parameters
	int B1,B2;
	B1 = atoi(av[4]);
	B2 = atoi(av[5]);
	
	
	/* Read Fragments */
	Fragmentv3* f;
	int nf;
	int xtotal,ytotal;
		
	f=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);
	/******************/
	int i=0;
	
	/* Debug
	printf("f.xIni\tf.yIni\tf.xFin\tf.yFin\tf.length\tf.strand\tf.score\tf.block\n");
	
	for(i=0;i<nf;i++){
		printf("%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n",(int)f[i].xIni,(int)f[i].yIni,(int)f[i].xFin,(int)f[i].yFin,(int)f[i].length,f[i].strand,(int)f[i].score,(int)f[i].block);
	}
	*/
	


	
	/* Calculating CSB1 and CBS2 as two fragment */
	Fragmentv3 C1,C2;
	resetFragmentv3(&C1);
	resetFragmentv3(&C2);
	// set new fragments.

	C1.block=B1;
	C2.block=B2;
	
	//printf("%d\t%d\t%d\n",B1,B2,nf);
	for(i=0;i<nf;i++){
	
		if((int)f[i].block == B1){
		
			updateFragmentv3(&C1, f[i]);
		}
		
		if((int)f[i].block == B2){
		
			updateFragmentv3(&C2, f[i]);
		}
	}
	
	// Debug
//	printf("%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n",(int)C1.xIni,(int)C1.yIni,(int)C1.xFin,(int)C1.yFin,(int)C1.length,C1.strand,(int)C1.score,(int)C1.block);
//	printf("%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n",(int)C2.xIni,(int)C2.yIni,(int)C2.xFin,(int)C2.yFin,(int)C2.length,C2.strand,(int)C2.score,(int)C2.block);
	/****************************/
	
	/* calculating the 3 areas coordinates*/
	/*  
		Area 1: 
			X: From C1.xIni to C2.xFin  	// A1xIni , A1xFin
			Y: from C1.yIni to C2.yFin		// A1yIni , A1yFin
		Area 2:
			X: From C1.xFin to C2.xIni		// A2xIni , A2xFin
			Y: From C1.yFin to C1.yIni+d	// A2yIni , A2yFin
		Area 3:
			X: From C1.xFin to C2.xIni		// A3xIni , A3xFin
			Y: From C2.yIni-d to C2.yIni	// A3yIni , A3yFin
			
	*/
	int A1xIni,A1xFin,A1yIni,A1yFin,A2xIni,A2xFin,A2yIni,A2yFin,A3xIni,A3xFin,A3yIni,A3yFin;
	int dy,dx;
	
	dx=C2.xIni-C1.xFin;
	dy=C2.yIni-C1.yFin;
	
	A1xIni=C1.xIni;
	A1xFin=C2.xFin; 	
	A1yIni=C1.yIni; 
	A1yFin=C2.yFin;	

	//	
	A2xIni=C1.xFin;
	A2xFin=C2.xIni;		
	A2yIni=C1.yFin;
	A2yFin=C2.yIni;
	
	A3xIni=C1.xFin;
	A3xFin=C2.xIni;		 
	A3yIni=C1.yFin;
	A3yFin=C2.yIni;
	
	if(dy>dx){
		A2yFin = C1.yFin+dx;
		A3yIni = C2.yIni-dx;
	}
	if(dy<=dx){
		A2xFin=C1.xFin+dy;
		A3xIni=C2.xIni-dy;
	}

//	printf("dx: %d\tdy: %d\n",dx,dy);
//	printf("%d\t%d\t%d\t%d\n",A1xIni,A1yIni,A1xFin,A1yFin);
//	printf("%d\t%d\t%d\t%d\n",A2xIni,A2yIni,A2xFin,A2yFin);
//	printf("%d\t%d\t%d\t%d\n",A3xIni,A3yIni,A3xFin,A3yFin);
	
	/**********************************/
			/*********READ HITS **********/
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
	
	FILE *fH;
	hit h;
	int nHits=0;
	
	if ((fH = fopen(av[2], "rb")) == NULL)
		terror("opening HITS file");
	// read Hits
	
	if(fread(&h, sizeof(hit), 1, fH)!=1){
		terror("Empty hits file");
	}
	while (!feof(fH)) {
		nHits++;
		
		fread(&h, sizeof(hit), 1, fH);
		
		if (condicion(h,A1xIni,A1yIni,A1xFin,A1yFin))printf("HIT,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)h.posX,(int)h.posY,(int)h.posX+8,(int)h.posY+8,'f',0,0,0,0,0,0,0,0);

		//printf("HIT,%d,%d,%d,%d,%c\n",(int)h.posX,(int)h.posY,(int)h.posX+8,(int)h.posY+8,'f');
		
	
	}

	fclose(fH);
	/*******************/


	return 0;
}

/****/
int condicion(hit h,int xIni,int yIni,int xFin,int yFin){

	if (((int)h.posX >= xIni) && ((int)h.posX)<=xFin && ((int)h.posY>= yIni) && ((int)h.posY<=yFin))return 1;
	return 0;

}



