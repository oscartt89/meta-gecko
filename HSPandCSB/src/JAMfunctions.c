/*

JAMfunctions.c


*/

#include "JAMfunctions.h"

// Update fragment
void updateFragmentv3(Fragmentv3* f, Fragmentv3 g){
			f->xIni=MIN(f->xIni,g.xIni);
			f->yIni=MIN(f->yIni,g.yIni);
			
			f->xFin=MAX(f->xFin,g.xFin);
			f->yFin=MAX(f->yFin,g.yFin);
			
			f->strand=g.strand;
			
			f->score+=g.score;
			f->length = f->xFin-f->xIni;
			f->similarity = f->score/(4*f->length);
			f->ident += g.ident;
}
// PrintFragments
void printFragv3(Fragmentv3 f){
	printf("%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n",(int)f.xIni,(int)f.yIni,(int)f.xFin,(int)f.yFin,(int)f.length,f.strand,(int)f.score,(int)f.block);

}
// Reset fragment
void resetFragmentv3(Fragmentv3* f){
	f->xIni=(int)pow(2,8*sizeof(int));
	f->yIni=(int)pow(2,8*sizeof(int));
	f->xFin=0;
	f->yFin=0;
	f->length=f->score=f->similarity=0;

}
// Fragment2 version
// returns the percentage of overlapping
// Devuelve porcentaje de solapamiento
int overlapFragFile ( struct FragFile f, struct FragFile g,char c,int sol){ 

	double d=0.0;
	d=(double)sol;
	d=0.0;
	int plf,plg;
	long int gxS,fxS,gxE,fxE;
	long int gyS,fyS,gyE,fyE;
	
	if(c=='x'){
		gxS=(long int)g.xStart;
		gxE=(long int)g.xEnd;
		fxS=(long int)f.xStart;
		fxE=(long int)f.xEnd;
		d= MIN(gxE,fxE)-MAX(gxS,fxS);
		if(d<0)d=0;
	}else{
		gyS=MIN((long int)g.yStart,(long int)g.yEnd);
		gyE=MAX((long int)g.yEnd,(long int)g.yStart);
		fyS=MIN((long int)f.yStart,(long int)f.yEnd);
		fyE=MAX((long int)f.yEnd,(long int)f.yStart);
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
		return (int)d ;
	}else{
	
//		printf("-- NO SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");
	
		return 0;
	}
}

// Fragment fragmetnv3
// returns the percentage of overlapping
// Devuelve porcentaje de solapamiento
int overlapFragmentv3 ( Fragmentv3 f, Fragmentv3 g,char c,int sol){ 

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
		return (int)d ;
	}else{
	
//		printf("-- NO SOLAPA ---\n");
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)f.xStart,(long int)f.yStart,(long int)f.xEnd,(long int)f.yEnd,(long int)f.length,(long int)f.strand,(long int)f.ident,(long int)f.score,(long int)f.seqX,(long int)f.seqY,(long int)plf);
//		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long int)g.xStart,(long int)g.yStart,(long int)g.xEnd,(long int)g.yEnd,(long int)g.length,(long int)g.strand,(long int)g.ident,(long int)g.score,(long int)f.seqX,(long int)f.seqY,(long int)plg);
//		printf("-------\n");
	
		return 0;
	}
}

/******************************************/
int MIN(long int a, long int b){if (a>=b)return b;else return a;}
int MAX(long int a, long int b){if (a>=b)return a;else return b;}
/******************************************/
