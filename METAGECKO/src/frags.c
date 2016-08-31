/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes the workflow of GECKO for create
 *    metagenome dictionaries.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "frags.h"
#include <omp.h>

/* This main contains the workflow to find hits, filter and extend it to generate
 * fragments over a length and similarity wanted. There are two ways to invoke the
 * the program:
 *   @use frag metagDict metagFile genoDic genoFile out minS minL f/r prefix
 *   @use frag metagDict metagFile genoDic genoFile out minS minL f/r prefix startIndx
 * Where the parameters used are:
 *   @param metagDict the string with the basename of the metagenome dictionary to be
 *          used. Must contain the relative/absolute path to the file if it's not in
 *          the invocation folder. Note: the base name is <basename>.(d2hP/W)
 *   @param metagFile is the FASTA (or any otther extension) metagenome file used to
 *          generate the metagenome dictionary given.
 *   @param genoDict the string with the basename of the genome dictionary to be
 *          used. Must contain the relative/absolute path to the file if it's not in
 *          the invocation folder. Note: the base name is <basename>.(d2hP/W)
 *   @param genoFile is the FASTA (or any otther extension) genome file used to
 *          generate the genome dictionary given.
 *   @param out is the basename of the file where fragmetns will be stored.
 *   @param minS is the minimum similarity that must have a fragment to be stored.
 *   @param minL is the minimum length that must have a fragment to be stored.
 *   @param f/r a char that indicates the direction of the genome dictionary forward (f)
 *          or reverse (r).
 *   @param prefix is the subsequence length that will be used of the dictionaries words.
 *          Note: length = prefix * 4. If length is bigger than dictionary words length an
 *          error will be launched.
 * The result of the program will be the following:
 *   @file <out>.frags file with all fragments generated that satisfies the similarity and
 *         length thresholds.
 * The program will write some messages in default output stream that shows the action that
 * is being done by the program each moment.
 * Warning: if any error is launched, the program automatically stops, if it happens is
 * normal that the following files appear in your output folder:
 *   @file <out>.hindx is an auxiliary file with information about some buffer used in the 
 *         program process.
 *   @file <out>.hts is an auxilary file with information about the seed generated during 
 *         the comparisson process.
 */
int main(int ac, char **av) {
    // Check arguments
    if (ac != 10 && ac != 9) {
        fprintf(stderr,
                "Bad call error.\nUSE: frags metagDic metagFile genoDic genoFile out minS minL prefix\nUSE: frags metagDic metagFile genoDic genoFile out minS minL prefix startIndx\n");
        return -1;
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, "\tFrags: Starting fragments program.\n");
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Necessary variables
    char *fname; // File names handler
    int prefixSize; // Word size


    // Memory for file names handler
    if ((fname = (char *) malloc(sizeof(char) * MAX_FILE_LENGTH)) == NULL) {
        fprintf(stderr, "Error allocating memory for file names handler.\n");
        return -1;
    }


    // Check arguments
    strcpy(fname, av[1]);
    if (!exists(strcat(fname, ".d2hP"))) { // Check metagenome dict
        fprintf(stderr, "Error:: couldn't find metagenome dictionary.\n");
        return -1;
    }
    if (!exists(av[2])) { // Check metagenome file
        fprintf(stderr, "Error:: Metagenome file specified doesn't exists\n");
        return -1;
    }
    strcpy(fname, av[3]);
    if (!exists(strcat(fname, ".d2hP"))) { // Check genome dict
        fprintf(stderr, "Error:: couldn't find genome dictionary.\n");
        return -1;
    }
    if (!exists(av[4])) { // Check metagenome dict
        fprintf(stderr, "Error:: Genome file specified doesn't exists\n");
        return -1;
    }
    if (is_float(av[6])) { // Check similarity threshold
        fprintf(stderr, "Error:: Similarity threshold specified isn't a float number\n");
        return -1;
    } else if (atof(av[6]) < 0 || atof(av[6]) > 100) {
        fprintf(stderr, "Error:: Similarity threshold specified isn't contained in range [0,100]\n");
        return -1;
    }
    if (!is_int(av[7])) { // Check length threshold
        fprintf(stderr, "Error:: Length threshold specified isn't a number.\n");
        return -1;
    } else if (atoi(av[7]) < 0) {
        fprintf(stderr, "Error:: Similarity threshold must be positive.\n");
        return -1;
    }
    if (!is_int(av[8])) { // Check prefix
        fprintf(stderr, "Error:: Prefix specified isn't a number.\n");
        return -1;
    } else if (atoi(av[8]) < 1) {
        fprintf(stderr, "Error:: Prefix must be >1.\n");
        return -1;
    }


    // Variables
    FILE *mW, *mP, *gW, *gP; // Dictionaries
    FILE *hIndx, *hts; // Intermediate files
    Hit *buffer;
    uint64_t hitsInBuffer = 0, genomeLength, nStructs, metagenomeLength;
    uint16_t mWL, gWL;
    uint16_t BytesGenoWord = 8, BytesMetagWord;
    uint64_t buffersWritten = 0; // Init global variable (frags.h)1
    
    prefixSize = atoi(av[8]); // Prefix array length
    bool removeIntermediataFiles = false; // Internal variable to delete intermediate files

    // Allocate necessary memory
    // Memory for buffer
    if ((buffer = (Hit *) malloc(sizeof(Hit) * MAX_BUFF)) == NULL) {
        fprintf(stderr, "Error allocating memory for hits buffer.\n");
        return -1;
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, "\tFrags: Opening/creating necessary files.");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Open current necessary files
    // Open metagenome positions file
    strcpy(fname, av[1]);
    if ((mP = fopen(strcat(fname, ".d2hP"), "rb")) == NULL) {
        fprintf(stderr, "Error opening metagenome positions dictionaries.\n");
        return -1;
    }

    // Open metagenome words file
    strcpy(fname, av[1]);
    if ((mW = fopen(strcat(fname, ".d2hW"), "rb")) == NULL) {
        fprintf(stderr, "Error opening metagenome words dictionaries.\n");
        return -1;
    }
    // Read words header = WordLength
    // Check
    if (fread(&mWL, sizeof(uint16_t), 1, mW) != 1) {
        fprintf(stderr, "Error, couldn't find word length on metagenome dictionary.\n");
        return -1;
    } else if (mWL % 4 != 0) {
        fprintf(stderr, "Error, word length of metagenome dictionary isn't a 4 multiple.\n");
        return -1;
    } else
        BytesMetagWord = mWL / 4;

    // Open genome postions file
    strcpy(fname, av[3]);
    if ((gP = fopen(strcat(fname, ".d2hP"), "rb")) == NULL) {
        fprintf(stderr, "Error opening genome positions dictionaries.\n");
        return -1;
    }

    // Open genome words file
    strcpy(fname, av[3]);
    if ((gW = fopen(strcat(fname, ".d2hW"), "rb")) == NULL) {
        fprintf(stderr, "Error opening genome words dictionaries.\n");
        return -1;
    }

    // Read words header = WordLength
    // Check
    if (fread(&gWL, sizeof(uint16_t), 1, gW) != 1) {
        fprintf(stderr, "Error, couldn't find word length on genome dictionary.\n");
        return -1;
    } else if (gWL % 4 != 0) {
        fprintf(stderr, "Error, word length of genome dictionary isn't a 4 multiple.\n");
        return -1;
    } else
        BytesGenoWord = gWL / 4;

    if (BytesMetagWord * 4 < prefixSize || BytesGenoWord * 4 < prefixSize) {
        fprintf(stderr, "Error: prefix is too long.\n");
        return -1;
    }

    //Open intermediate files
    strcpy(fname, av[5]); // Copy outDic name
    if ((hIndx = fopen(strcat(fname, ".hindx"), "wb")) == NULL) {
        fprintf(stderr, "Error opening buffer index file.\n");
        return -1;
    }

    // Open hits repo
    strcpy(fname, av[5]);
    if ((hts = fopen(strcat(fname, ".hts"), "wb")) == NULL) {
        fprintf(stderr, "Error opening hits repository.\n");
        return -1;
    }

    

#pragma omp parallel num_threads(2)
    {
#pragma omp sections
        {
#pragma omp section
            {
                /////////////////////////// CHECKPOINT ///////////////////////////
                fprintf(    stdout, "\n\tFrags: Loading sequences of genome.\n");
                fflush(stdout);
                /////////////////////////// CHECKPOINT ///////////////////////////

                // Load sequences
                getSeqDBLength(av[4], &genomeLength, &nStructs);
                fprintf(stdout, "{REPORT} Genome DB has length %"PRIu64"\n", genomeLength);

                /////////////////////////// CHECKPOINT ///////////////////////////
                fprintf(stdout, "\t(Loaded) sequences of genome\n");
                /////////////////////////// CHECKPOINT ///////////////////////////
            }
#pragma omp section
            {
                /////////////////////////// CHECKPOINT ///////////////////////////
                fprintf(stdout, "\tFrags: Loading sequences of metagenome.\n");
                fflush(stdout);
                /////////////////////////// CHECKPOINT ///////////////////////////

                LoadMetagenome(av[2], &metagenomeLength);
                fprintf(stdout, "{REPORT} Metagenome has length %"PRIu64"\n", metagenomeLength);

                /////////////////////////// CHECKPOINT ///////////////////////////
                fprintf(stdout, "\t(Loaded) sequences of metagenome\n");
                /////////////////////////// CHECKPOINT ///////////////////////////
            }
        }
    }
    
    // Search hits
    // Prepare necessary variables
    HashEntry we[2]; // [0]-> Metagenome [1]-> Genome
    uint64_t lastFirstHit = (uint64_t) ftello(gW);
    bool firstmatch = true;
    int cmp;

    // Take memory
    if ((we[0].seq = (unsigned char *) malloc(sizeof(unsigned char) * BytesMetagWord)) == NULL) {
        fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
        return -1;
    } else we[0].WB = BytesMetagWord;

    if ((we[1].seq = (unsigned char *) malloc(sizeof(unsigned char) * BytesGenoWord)) == NULL) {
        fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
        return -1;
    } else we[1].WB = BytesGenoWord;

    // Read first entrances
    if (readHashEntrance(&we[0], mW, BytesMetagWord) < 0) return -1;

    if (readHashEntrance(&we[1], gW, BytesGenoWord) < 0) return -1;

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, "\t(Done)\n");
    fprintf(stdout, "\tFrags: Generating seeds.");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////


    // Search
    if (prefixSize == BytesMetagWord * 4 && prefixSize == BytesGenoWord * 4) {

        while (!feof(mW) && !feof(gW)) {
            if ((cmp = wordcmp(we[0].seq, we[1].seq, prefixSize)) == 0) { // Hit
                if (generateHits(buffer, we[0], we[1], mP, gP, hIndx, hts, &hitsInBuffer, prefixSize,
                                 &buffersWritten, genomeLength, metagenomeLength) < 0)
                    return -1;
                fflush(stdout);
            }

            //getchar();

            // Load next word
            if (cmp <= 0) { // New metagenome word is necessary
                if (!feof(mW)) if (readHashEntrance(&we[0], mW, BytesMetagWord) < 0) return -1;
            }
            if (cmp >= 0) { // New genome word is necessary
                if (!feof(gW)) if (readHashEntrance(&we[1], gW, BytesGenoWord) < 0) return -1;
            }
        }
    } else {
        while (!feof(mW)) {
            // Check hit
            if ((cmp = wordcmp(we[0].seq, we[1].seq, prefixSize)) == 0) { // Hit
                if (generateHits(buffer, we[0], we[1], mP, gP, hIndx, hts, &hitsInBuffer, prefixSize,
                                 &buffersWritten, genomeLength, metagenomeLength) < 0)
                    return -1;
                if (firstmatch) {
                    lastFirstHit = ((uint64_t) ftello(gW) - size_of_HashEntry(BytesGenoWord));
                    firstmatch = false;
                }
            }
            // Check if could be more
            if (cmp >= 0) { // Could be more
                // Load next genome word
                readHashEntrance(&we[1], gW, BytesGenoWord);
                if (feof(gW)) { // End of genome file
                    // Load next metagenome word
                    if (readHashEntrance(&we[0], mW, BytesMetagWord) < 0) return -1;
                    // Reset values and come back at dict
                    firstmatch = true;
                    fseeko(gW, lastFirstHit, SEEK_SET); // Reset geno dict
                    readHashEntrance(&we[1], gW, BytesGenoWord);
                }
            } else if (cmp < 0) { // No more matches, take next metag word
                // Load next metagenome word
                if (readHashEntrance(&we[0], mW, BytesMetagWord) < 0) return -1;
                // Reset values and come back at dict
                firstmatch = true;
                fseeko(gW, lastFirstHit, SEEK_SET); // Reset geno dict
                readHashEntrance(&we[1], gW, BytesGenoWord);
            }
        }
    }

    //If there are NO hits at all
    if (hitsInBuffer == 0 && buffersWritten == 0) {
        /////////////////////////// CHECKPOINT ///////////////////////////
        fprintf(stdout, "\tFrags: No match found.\n");
        fprintf(stdout, "\tFrags: Closing the program.\n");
        /////////////////////////// CHECKPOINT ///////////////////////////

        // Free auxiliar buffers
        free(we[0].seq);
        free(we[1].seq);
        free(buffer);

        // Close files
        fclose(mW);
        fclose(gW);
        fclose(mP);
        fclose(gP);
        fclose(hIndx);
        fclose(hts);

        // Remove intermediate files
        if (removeIntermediataFiles) {
            strcpy(fname, av[5]);
            remove(strcat(fname, ".hts"));
            strcpy(fname, av[5]);
            remove(strcat(fname, ".hindx"));
        }

        free(fname);

        // End program
        return 0;
    }

    //If more than one buffer has been written, write the last one
    if (hitsInBuffer > 0 && buffersWritten > 0) {
        writeHitsBuff(buffer, hIndx, hts, hitsInBuffer, prefixSize, &buffersWritten);
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, " (Generated)\n");
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Write hits if there is only one buffer

    if (hitsInBuffer > 0 && buffersWritten <= 0) { // Only one buffer
        // Sort buffer
        //quicksort_H(buffer, 0, hitsInBuffer - 1);
        writeHitsBuff(buffer, hIndx, hts, hitsInBuffer, prefixSize, &buffersWritten);
        
        // Close unnecesary files
        fclose(mW);
        fclose(gW);
        fclose(mP);
        fclose(gP);
        fclose(hIndx);
        fclose(hts);

        // Free unnecesary variables
        free(we[0].seq);
        free(we[1].seq);
    }

    // Free malloc blocks
    free(fname);


    // Everything finished OK
    return 0;
}
