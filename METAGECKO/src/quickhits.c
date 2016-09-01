/*********

File        dictMgMem.c
Author      EPW <estebanpw@uma.es>
Description Computes the in-memory dictionary of kmers for a metagenome, then searches for hits through the target database and
            produces a histogram of sequences hits

USAGE       <metagenome>            The .fasta file of the metagenome
            <genomes_database>      The .fasta database containing the genomes
            <sequence_hits_histogram>   The output binary file containing the amount of hits per sequence
            <ksize[1=8,2=16,3=32]>  The size of the word, e.g. for k=32 use value 3


**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "metag_common.h"


#define STARTING_SEQS 1000
static uint32_t ksizeTables[4] = {0,2,4,8};


int main(int argc, char ** av){
    
    if(argc != 5) terror("USE: quickhits <metagenome> <database> <sequences_hits_histogram> <ksize[1=8,2=16,3=32]>");
    
    
    
    //metagenome to read kmers from, database to find seeds
    FILE * metagenome, * database, * histoSeqs;
    
        
    //Open metagenome 
    metagenome = fopen64(av[1], "rt");
    if(metagenome == NULL) terror("Could not open metagenome");
    
    //Open database
    database = fopen64(av[2], "rt");
    if(database == NULL) terror("Could not open database");

    //Open output file
    histoSeqs = fopen64(av[3], "wb");
    if(histoSeqs == NULL) terror("Could not open output histogram file");

    //Ksize for seeds
    uint32_t ksizeidx = (uint32_t) atoi(av[4]);
    

    //Allocate memory for the sequences histogram
    uint64_t * seqHits = (uint64_t *) calloc(STARTING_SEQS, sizeof(uint64_t));
    if(seqHits == NULL) terror("Could not allocate memory for the sequences histogram");
    
    //Get size of db
    fseeko64(metagenome, 0L, SEEK_END);
    uint64_t totalSize = ftello64(metagenome);
    fseeko64(metagenome, 0L, SEEK_SET);
    
    //Allocate memory, as much as the size of the file for a table that will be sorted
    unsigned char * tableForward = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    unsigned char * tableReverse = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    if(tableForward == NULL) terror("Could not allocate memory for forward table");
    if(tableReverse == NULL) terror("Could not allocate memory for reverse table");
    
    
    //Variables to read kmers
    char c = 'N'; //Char to read characters

    unsigned char b[8], br[8];
    memset(b, 0, ksizeTables[ksizeidx]);
    memset(br, 0, ksizeTables[ksizeidx]);
    
    //Variables to account for positions
    int strandF = 1, strandR = 1;
    uint64_t pos = 0, crrSeqL = 0; //total length, current sequence length
    uint64_t idx = 0; //index for table access
    uint64_t nAlloc = 1; //Times sequence histogram heap was allocated
    
    //Print info
    fprintf(stdout, "[INFO] Computing tree of mers\n");
    
    c = fgetc(metagenome);
    while(!feof(metagenome)){
        // Check if it's a special line
        if (!isupper(toupper(c))) { // Comment, empty or quality (+) line
            if (c == '>') { // Comment line
                c = fgetc(metagenome);
                while (c != '\n') c = fgetc(metagenome); //Avoid comment line

                crrSeqL = 0; // Reset buffered sequence length

                pos++; // Absolute coordinates: add one for the "*"
            }
            c = fgetc(metagenome); // First char of next sequence
            continue;
        }
        
        if (strandF) shift_byte_array_left(b, ksizeTables[ksizeidx]); // Shift bits sequence
        if (strandR) shift_byte_array_right(br, ksizeTables[ksizeidx]); // Shift bits sequence

        // Add new nucleotide
        switch (c) {
            case 'A': // A = 00
                crrSeqL++;
                if (strandR) br[0] |= 192;
                break;
            case 'C': // C = 01
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 1;
                if (strandR) br[0] |= 128;
                crrSeqL++;
                break;
            case 'G': // G = 10
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 2;
                if (strandR) br[0] |= 64;
                crrSeqL++;
                break;
            case 'T': // T = 11
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 3;
                crrSeqL++;
                break;
            default : // Bad formed sequence
                crrSeqL = 0;
                break;
        }
        pos++;
        if (crrSeqL >= (uint64_t) ksizeTables[ksizeidx]*4) { // Full well formed sequence
            if(strandF){
                memcpy(tableForward+idx, b, ksizeTables[ksizeidx]);
            }
            if(strandR){
                memcpy(tableReverse+idx, b, ksizeTables[ksizeidx]);
            }
            idx += ksizeTables[ksizeidx];
            
        }
        c = fgetc(metagenome);

    }
    fprintf(stdout, "[INFO] Sequence of length %"PRIu64" has %"PRIu64" mers of size k=%d\n", pos, pos-ksizeTables[ksizeidx]*4, ksizeTables[ksizeidx]*4);
    fprintf(stdout, "[INFO] Table of seeds used up %"PRIu64" bytes, which are %"PRIu64" MegaBytes\n", idx, idx/(1024*1024));
    
    fprintf(stdout, "[INFO] Sorting forward table\n");
    QuickSortByteArray(tableForward, 0, idx/ksizeTables[ksizeidx] - 1, ksizeTables[ksizeidx]);
    fprintf(stdout, "[INFO] Sorting reverse table\n");
    QuickSortByteArray(tableReverse, 0, idx/ksizeTables[ksizeidx] - 1, ksizeTables[ksizeidx]);


    //printTable(tableForward, idx/ksizeTables[ksizeidx], ksizeTables[ksizeidx], stdout);
    //printTable(tableReverse, idx/ksizeTables[ksizeidx], ksizeTables[ksizeidx], stdout);







    fprintf(stdout, "[INFO] Generating seeds\n");

    memset(b, 0, ksizeTables[ksizeidx]);
    memset(br, 0, ksizeTables[ksizeidx]);

    c = 'N';
    pos = 0, crrSeqL = 0; //total length, current sequence length
    uint64_t seqNum = 0, nMers = idx/ksizeTables[ksizeidx]; // current sequence index, total number of k-mers
    int64_t found; // If the binary search function has found the kmer

    c = fgetc(database);
    while(!feof(database)){
        // Check if it's a special line
        if (!isupper(toupper(c))) { // Comment, empty or quality (+) line
            if (c == '>') { // Comment line
                c = fgetc(database);
                while (c != '\n') c = fgetc(database); //Avoid comment line

                crrSeqL = 0; // Reset buffered sequence length
                seqNum++; // Current sequence index increases
                pos++; // Absolute coordinates: add one for the "*"
            }
            c = fgetc(database); // First char of next sequence
            continue;
        }
        
        if (strandF) shift_byte_array_left(b, ksizeTables[ksizeidx]); // Shift bits sequence
        if (strandR) shift_byte_array_right(br, ksizeTables[ksizeidx]); // Shift bits sequence

        // Add new nucleotide
        switch (c) {
            case 'A': // A = 00
                crrSeqL++;
                if (strandR) br[0] |= 192;
                break;
            case 'C': // C = 01
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 1;
                if (strandR) br[0] |= 128;
                crrSeqL++;
                break;
            case 'G': // G = 10
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 2;
                if (strandR) br[0] |= 64;
                crrSeqL++;
                break;
            case 'T': // T = 11
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 3;
                crrSeqL++;
                break;
            default : // Bad formed sequence
                crrSeqL = 0;
                break;
        }
        pos++;
        if (crrSeqL >= (uint64_t) ksizeTables[ksizeidx]*4) { // Full well formed sequence
            if(strandF){
                //Binary Search
                found = binarySearchByteArray(b, tableForward, ksizeTables[ksizeidx], nMers);
                
                if(found >= 0){
                    if(seqNum > nAlloc*STARTING_SEQS){
                        nAlloc++;
                        seqHits = (uint64_t *) realloc(seqHits, nAlloc*STARTING_SEQS*sizeof(uint64_t));
                        if(seqHits == NULL) terror("Could not reallocate sequence seed counters");
                    }
                    seqHits[seqNum-1]++;
                }
            }
            if(strandR){
                //Binary Search
                found = binarySearchByteArray(br, tableReverse, ksizeTables[ksizeidx], nMers);
                if(found >= 0){
                    if(seqNum > nAlloc*STARTING_SEQS){
                        nAlloc++;
                        seqHits = (uint64_t *) realloc(seqHits, nAlloc*STARTING_SEQS*sizeof(uint64_t));
                        if(seqHits == NULL) terror("Could not reallocate sequence seed counters");
                    }
                    seqHits[seqNum-1]++;
                }
            }
            
        }
        c = fgetc(database);

    }

    fprintf(stdout, "[INFO] Finished generating seeds\n");
    fprintf(stdout, "[INFO] Writing sequences hits histogram\n");

    /*
    //only for display
    uint64_t i;
    for(i=0;i<seqNum;i++){
        fprintf(stdout, "%"PRIu64"\t%"PRIu64"\n", i, seqHits[i]);
    }
    */

    //Store the sequences as binary
    if(seqNum != fwrite(seqHits, sizeof(uint64_t), seqNum, histoSeqs)) terror("Could not write histogram");

    fprintf(stdout, "[INFO] Done. Deallocating memory\n");

    fclose(metagenome);
    fclose(database);
    fclose(histoSeqs);

    free(seqHits);
    free(tableForward);
    free(tableReverse);
    
    
    
    
    return 0;
}