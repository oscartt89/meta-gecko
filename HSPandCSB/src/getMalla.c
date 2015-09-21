/*
 * getMalla.c
 *
 *  Created on: 11/11/2014
 *      Author: jarjonamedina
 * 
 
 If three or more fragments are duplicated both X and Y axis, the define a malla.

 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"



  #define MAX_LEVELS 1000000



int overlapFragFile ( Fragmentv3  f, Fragmentv3  g,char c,int sol);
void quickSortX(Fragmentv3 *f, int elements);
void quickSortY(Fragmentv3 *f, int elements,int ytotal);

int MAX(int a, int b);
int MIN(int a, int b);

int main(int ac,char** av){


if(ac<3){
	printf("Use: getMalla fragsv3 fragsFiltr.fil frags.malla \n");
	exit(-1);
}





/* Open file */
FILE* fout;
	if((fout=fopen(av[2],"wb"))==NULL){
		printf("***ERROR Opening input file %s\n",av[2]);
		exit(-1);
	}
/* Open file */
FILE* foutDUP;

	if((foutDUP=fopen(av[3],"wb"))==NULL){
		printf("***ERROR Opening input file %s\n",av[3]);
		exit(-1);
	}



/* Read Fragments */
	Fragmentv3* f;
	int nf,xtotal,ytotal;
	f=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);

int i,j;


// Reset frags
for(i=0;i<nf;i++){
	f[i].gap=i+3;
	f[i].events=0;
	f[i].duplications=0;
}
/*  Empezamos */
for(i=0;i<nf;i++){
	if(f[i].block>=0){
	
			for(j=i+1;j<nf;j++){
				if(f[j].block>=0){
					if( overlapFragFile (f[i], f[j],'x',90)||overlapFragFile (f[i], f[j],'y',90) ){
						
							f[j].gap=MIN(f[i].gap,f[j].gap);
							f[i].events=-1;
							f[j].events=-1;
						
					}
				}
			
			}
		
		
	}
}

/*   
11/11/14
Ahora mismo todos los duplicados estan marcados con n numero negativo.

Proximo paso: Calcular cuantos tipos de bloques duplicados hay, Y numerarlos de nuevo. De esta manera el min(bloque) será el numero de
tipos duplicados o mallas.
Mas tarde en getCSB, los bloques duplicados se podrán unir en CSB duplicados. Tanto en X como en Y tendremos posibilidad de unir boqes 
duplicados. Cuando vayamos a unir bloque duplicados hay q tener en cuenta q podemos coger el que mejor encaje.
*/
/***********************************/

//**********************/

/* Write fragments */ ////score>1
fwrite(&xtotal,sizeof(int),1,fout);
fwrite(&ytotal,sizeof(int),1,fout);

fwrite(&xtotal,sizeof(int),1,foutDUP);
fwrite(&ytotal,sizeof(int),1,foutDUP);

int ndup;
int nfrags;
ndup=nfrags=0;

//

for(i=0;i<nf;i++){

	
	if(f[i].events==0){
		f[i].duplications=0;
	}else{
		f[i].block=-1*f[i].gap;
	}
	
	fwrite(&f[i],sizeof(Fragmentv3),1,fout);
	printf("%d\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%d\t%Le\t%d\t%d\n",f[i].block,f[i].xIni,f[i].yIni,f[i].xFin,f[i].yFin,f[i].strand,f[i].score,f[i].ident,f[i].length,f[i].gap,f[i].pvalue,f[i].duplications,f[i].block);

}


/*********************** END OF FILTER **************/
fclose(fout);
fclose(foutDUP);
return 0;
}

/*******************************************
*****************************************/


int overlapFragFile ( Fragmentv3 f, Fragmentv3 g,char c,int sol){ 

	double d=0.0;
	d=(double)sol;
	d=0.0;
	int plf,plg;
	long int gxS,fxS,gxE,fxE;
	long int gyS,fyS,gyE,fyE;
	
	if(c=='x'){
		gxS=(long int)g.xIni;
		gxE=(long int)g.xFin;
		fxS=(long int)f.xIni;
		fxE=(long int)f.xFin;
		d= MIN(gxE,fxE)-MAX(gxS,fxS);
		if(d<0)d=0;
	}else{
		gyS=MIN((long int)g.yIni,(long int)g.yFin);
		gyE=MAX((long int)g.yFin,(long int)g.yIni);
		fyS=MIN((long int)f.yIni,(long int)f.yFin);
		fyE=MAX((long int)f.yFin,(long int)f.yIni);
		d= MIN(gyE,fyE)-MAX(gyS,fyS);
		if(d<0)d=0;
	}
	
	if(d){
		
		plf = (int)(100*(double)d/(double)f.length);
		plg = (int)(100*(double)d/(double)g.length);
		
		d=MAX((long int)plf,(long int)plg);
//		printf("-- SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");		
		if(d>=sol)return (int)d ;
	}else{
	
//		printf("-- NO SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");
	
		return 0;
	}
	return 0;
}
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7


void quickSortX(Fragmentv3 *f, int elements) {



  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].xIni>=piv.xIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].xIni<=piv.xIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}}



/*****************/
/*****************/

void quickSortY(Fragmentv3 *f, int elements,int ytotal) {

  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
	long ytemp;
	ytemp=ytotal;
// Change Y component

	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			// we have already done this in fragv2tov3.c
			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;

			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	
		}
	}
//************/
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].yIni>=piv.yIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].yIni<=piv.yIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}

// Change Y component

	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	

			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;


		}
	}
//*************/

}
/****************************************/
/*************************************/



int MAX(int a, int b){if (a>=b)return a;if(b>a)return b;return 0;}
int MIN(int a, int b){if (a>=b)return b;if(b>a)return a;return 0;}

