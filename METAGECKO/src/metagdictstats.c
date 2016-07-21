/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes a process to read and take statistical information
 *    from a metagenome dictionary file.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "metag_common.h"

int main(int ac, char **av) {
    // Variables
    uint16_t wl;
    unsigned int BYTES_IN_WORD;
    long read;
    uint32_t aux32, minReps, maxReps;
    uint64_t numWords = 0, aux64;
    float meanReps = 0;
    uint64_t PSize, WSize;
    FILE *dict, *stats_file;
    char *fname;

    // Check
    if (ac != 3) {
        fprintf(stderr, "Bad call error.\nUSE: metagdictstats dictName outputFile\n");
        return -1;
    }

    // Memory for file names handler
    if ((fname = (char *) malloc(sizeof(char) * 1024)) == NULL) {
        fprintf(stderr, "Error allocating memory for file names handler.\n");
        return -1;
    }

    // Open metagenome file
    strcpy(fname, av[1]); // Copy outDic name
    if ((dict = fopen(strcat(fname, ".d2hW"), "rb")) == NULL) {
        fprintf(stderr, "Error opening metagenome dictionary file.\n");
        return -1;
    }

    //Reading Word length
    if (fread(&wl, sizeof(uint16_t), 1, dict) != 1) {
        fprintf(stderr, "Error reading word length from header\n");
        return -1;
    }

    BYTES_IN_WORD = wl / 4;

    // Open output stats file
    if (exists(av[2])) {
        if ((stats_file = fopen(av[2], "at")) == NULL) {
            fprintf(stderr, "Error opening file[%s]\n", av[4]);
            return -1;
        }
    } else {
        if ((stats_file = fopen(av[2], "wt")) == NULL) {
            fprintf(stderr, "Error opening file[%s]\n", av[4]);
            return -1;
        }
        fprintf(stats_file, "Metagenome\tWSize\tPSize\tWordLength\tNumWords\tMaxRep\tMinRep\tMeanRep\n");
    }

    fprintf(stdout, "Working with kmer length: %"
    PRIu16
    "\n", wl);
    unsigned char kmer[BYTES_IN_WORD];

    // Start to read
    // Read first kmer
    if (fread(kmer, sizeof(unsigned char), BYTES_IN_WORD, dict) != BYTES_IN_WORD) {
        fprintf(stderr, "Error reading first kmer sequence.\n");
        return -1;
    }
    if (fread(&aux64, sizeof(uint64_t), 1, dict) != 1) {
        fprintf(stderr, "Error reading first kmer position on PFile.\n");
        return -1;
    }
    if (fread(&aux32, sizeof(uint32_t), 1, dict) != 1) {
        fprintf(stderr, "Error reading first kmer number of repetitions\n");
        return -1;
    }

    // Init stats
    minReps = aux32;
    maxReps = aux32;

    char W[33];
    word ww;
    memcpy(ww.b,kmer,BYTES_IN_WORD);
    Word w;
    w.b = (unsigned char *) malloc(BYTES_IN_WORD*sizeof(unsigned char));
    memcpy(w.b, ww.b, BYTES_IN_WORD);
    showWord(&w, W, BYTES_IN_WORD*4);
    fprintf(stdout, "%s ", W);

    fprintf(stdout, "OFFSET: %" PRIu64 " REPETITIONS: %" PRIu32 "\n", aux64, aux32);

    while (!feof(dict)) {
        meanReps = (meanReps * numWords + aux32) / (numWords + 1);
        numWords++;

        if (aux32 < minReps) minReps = aux32;
        if (aux32 > maxReps) maxReps = aux32;

        if ((read = fread(kmer, sizeof(unsigned char), BYTES_IN_WORD, dict)) != BYTES_IN_WORD) {
            return -1;
        }

        char W[33];
        word ww;
        memcpy(ww.b,&kmer,BYTES_IN_WORD);

	w.b = (unsigned char *) malloc(BYTES_IN_WORD*sizeof(unsigned char));
	memcpy(w.b, ww.b, BYTES_IN_WORD);

        showWord(&w, W, BYTES_IN_WORD*4);
        fprintf(stdout, "%s ", W);

        if (fread(&aux64, sizeof(uint64_t), 1, dict) != 1) {
            fprintf(stderr, "Error reading kmer position on PFile.\n");
            return -1;
        }
        if (fread(&aux32, sizeof(uint32_t), 1, dict) != 1) {
            fprintf(stderr, "Error reading kmer number of repetitions\n");
            return -1;
        }

        fprintf(stdout, "OFFSET: %" PRIu64 " REPETITIONS: %" PRIu32 "\n", aux64, aux32);
    }

    // Take W size
    fseek(dict, 0L, SEEK_END); // go to end
    WSize = ftell(dict) / 1000000; // MB

    // Close unnecessary varaibles
    fclose(dict);

    // Open PDict file
    strcpy(fname, av[1]); // Copy outDic name
    if ((dict = fopen(strcat(fname, ".d2hP"), "rb")) == NULL) {
        fprintf(stderr, "Error opening P metagenome dictionary file.\n");
        return -1;
    }

    // Take P size
    fseek(dict, 0L, SEEK_END); // go to end
    PSize = ftell(dict) / 1000000; // MB

    // Close unnecessary varaibles
    fclose(dict);
    free(fname);

    // Print info
    fprintf(stats_file, "%s\t%"
    PRIu64
    "\t%"
    PRIu64
    "\t%"
    PRIu16
    "\t%"
    PRIu64
    "\t%"
    PRIu32
    "\t%"
    PRIu32
    "\t%f\n",
            av[1], WSize, PSize, wl, numWords, maxReps, minReps, meanReps);

    fclose(stats_file);

    // End
    return 0;
}
