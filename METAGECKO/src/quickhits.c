/*********

File        dictMgMem.c
Author      EPW <estebanpw@uma.es>
Description Computes the in-memory dictionary of kmers for a metagenome, then searches for hits through the target database and
            produces a histogram of sequences hits

USAGE       <metagenome>            The .fasta file of the metagenome
            <genomes_database>      The .fasta database containing the genomes
            <sequence_hits_histogram>   The output binary file containing the amount of hits per sequence
            <ksize[1=8,2=16,3=32]>  The size of the word, e.g. for k=32 use value 3
[OPTIONAL]
            <indexes_number_sorted_table>   [TODO] How many indexes will be stored for reducing the search space in the binary
                                            search function. E.g. if set to 2, it will store indexes for AA, AC, .. CA, CC, up to TG, TT.
                                            The default is 4.


**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "metag_common.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
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
    
    //If the number of indexes was specified
    uint64_t n_sort_idx = 3, n_sort_idx_table_size;
    if(argc == 6) n_sort_idx = asciiToUint64(av[5]);
    n_sort_idx_table_size = 4;
    uint64_t i;
    for(i=0;i<n_sort_idx;i++){
        n_sort_idx_table_size *= 4; //4 to the power of n_sort_idx
    }

    //Allocate memory for the sequences histogram
    uint64_t * seqHits = (uint64_t *) calloc(STARTING_SEQS, sizeof(uint64_t));
    if(seqHits == NULL) terror("Could not allocate memory for the sequences histogram");
    
    //Get size of db
    fseeko64(metagenome, 0L, SEEK_END);
    uint64_t totalSize = ftello64(metagenome);
    fseeko64(metagenome, 0L, SEEK_SET);
    
    //Allocate memory, as much as the size of the file for a table that will be sorted
    unsigned char * tableForward = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    //unsigned char * tableReverse = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    if(tableForward == NULL) terror("Could not allocate memory for forward table");
    //if(tableReverse == NULL) terror("Could not allocate memory for reverse table");
    
    //Indexes table
    uint64_t * idx_sorted_t_from = (uint64_t *) malloc(n_sort_idx_table_size*sizeof(uint64_t));
    if(idx_sorted_t_from == NULL) terror("Could not allocate memory for sorted indexes");

    //Variables to read kmers
    char c = 'N'; //Char to read characters

    unsigned char b[8], br[8];
    memset(b, 0, ksizeTables[ksizeidx]);
    
    //Variables to account for positions
    int strandF = 1, strandR = 1;
    uint64_t pos = 0, crrSeqL = 0; //total length, current sequence length
    uint64_t idx = 0; //index for table access
    uint64_t nAlloc = 1; //Times sequence histogram heap was allocated
    
    //Print info
    fprintf(stdout, "[INFO] Computing table of seeds\n");
    
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
        
        shift_byte_array_left(b, ksizeTables[ksizeidx]); // Shift bits sequence

        // Add new nucleotide
        switch (c) {
            case 'A': // A = 00
                crrSeqL++;
                break;
            case 'C': // C = 01
                b[ksizeTables[ksizeidx] - 1] |= 1;
                crrSeqL++;
                break;
            case 'G': // G = 10
                b[ksizeTables[ksizeidx] - 1] |= 2;
                crrSeqL++;
                break;
            case 'T': // T = 11
                b[ksizeTables[ksizeidx] - 1] |= 3;
                crrSeqL++;
                break;
            default : // Bad formed sequence
                crrSeqL = 0;
                break;
        }
        pos++;
        if (crrSeqL >= (uint64_t) ksizeTables[ksizeidx]*4) { // Full well formed sequence
            memcpy(tableForward+idx, b, ksizeTables[ksizeidx]);
            /*
            if(strandF){
                memcpy(tableForward+idx, b, ksizeTables[ksizeidx]);
            }
            if(strandR){
                memcpy(tableReverse+idx, b, ksizeTables[ksizeidx]);
            }
            */
            idx += ksizeTables[ksizeidx];
            
        }
        c = fgetc(metagenome);

    }
    fprintf(stdout, "[INFO] Sequence of length %"PRIu64" has %"PRIu64" mers of size k=%d\n", pos, pos-ksizeTables[ksizeidx]*4, ksizeTables[ksizeidx]*4);
    fprintf(stdout, "[INFO] Table of seeds used up %"PRIu64" bytes, which are %"PRIu64" MegaBytes\n", idx, idx/(1024*1024));
    
    fprintf(stdout, "[INFO] Sorting forward table\n");
    QuickSortByteArray(tableForward, 0, idx/ksizeTables[ksizeidx] - 1, ksizeTables[ksizeidx]);

    fprintf(stdout, "[INFO] Computing sorted table index\n");

    //b_comp will hild the current kmer we are looking for to index
    unsigned char b_comp = 0;
    uint64_t k = 0; //k is the index to go through the table
    i = 1;
    idx_sorted_t_from[0] = 0; //First kmer starts on zero, independtly of which one it is


    while(k < idx && i < n_sort_idx_table_size - 1){ //Go through the forward table

        if(b_comp < (unsigned char)tableForward[k]){
            idx_sorted_t_from[i] = (k)/ksizeTables[ksizeidx];
            b_comp += 1;
            i++;
        }else if(b_comp == (unsigned char)tableForward[k]){
            idx_sorted_t_from[i] = (k)/ksizeTables[ksizeidx];
            k+= ksizeTables[ksizeidx];
        }else{
            k+= ksizeTables[ksizeidx];
        }
        
    }
    //If some are missing at the end, just fill with the max
    i++;
    while(i < n_sort_idx_table_size){
        idx_sorted_t_from[i] = (idx/ksizeTables[ksizeidx]) - 1;
        i++;
    }

    /*
    for(i=0;i<25;i++){
        printf("Printing kmer %"PRIu64" -> %"PRIu64", %"PRIu64"\n", i, idx_sorted_t_from[i], idx_sorted_t_from[i+1]);    
    }
    */

    //fprintf(stdout, "[INFO] Sorting reverse table\n");
    //QuickSortByteArray(tableReverse, 0, idx/ksizeTables[ksizeidx] - 1, ksizeTables[ksizeidx]);

    
    printTable(tableForward, idx/ksizeTables[ksizeidx], ksizeTables[ksizeidx], stdout);
    /*
    for(i=0;i<n_sort_idx_table_size;i++){
        fprintf(stdout, "%"PRIu64"\n", idx_sorted_t_from[i]);
    }
    */
    //printTable(tableReverse, idx/ksizeTables[ksizeidx], ksizeTables[ksizeidx], stdout);







    fprintf(stdout, "[INFO] Generating seeds\n");

    memset(b, 0, ksizeTables[ksizeidx]);
    memset(br, 0, ksizeTables[ksizeidx]);

    c = 'N';
    pos = 0, crrSeqL = 0; //total length, current sequence length
    uint64_t seqNum = 0; // current sequence index
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
                found = binarySearchByteArray(b, tableForward, ksizeTables[ksizeidx], idx_sorted_t_from[b[0]], idx_sorted_t_from[b[0]+1]);
                //found = binarySearchByteArray(b, tableForward, ksizeTables[ksizeidx], 0, idx/ksizeTables[ksizeidx]);
                
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
                found = binarySearchByteArray(br, tableForward, ksizeTables[ksizeidx], idx_sorted_t_from[br[0]], idx_sorted_t_from[br[0]+1]);
                //found = binarySearchByteArray(br, tableForward, ksizeTables[ksizeidx], 0, idx/ksizeTables[ksizeidx]);
                
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
    free(idx_sorted_t_from);
    //free(tableReverse);
    
    
    
    
    return 0;
}