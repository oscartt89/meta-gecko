#include <stdio.h>
#include <stdlib.h>
#include "metag_common.h"

/* Buffered reader to read char by char
 *  @param s: Error mesage to print to stderr before exiting
 */


void terror(char *s) {
	fprintf(stderr, "ERR**** %s ****\n", s);
	exit(-1);
}


/* Buffered reader to read char by char
 *  @param buffer: An already allocated buffer of chars
 *  @param pos: Current position of last char read in buffer (should be initialized to READBUF+1
 *  @param read: Number of chars read at iteration
 *	@param f: File descriptor from which to read
 *	@return: The read char
 */

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

/* Translates an unsigned char word of ACTG letters compressed in 2 bits to a readable char
 *  @param w: Word containing the unsigned char bits word
 *  @param ws: Char to receive readable translation
 *  @param WORD_LENGTH: Length in bits of the word to translate
 */

void showWord(Word *w, char *ws, uint16_t WORD_LENGTH) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	int i;
	int wsize = WORD_LENGTH/4;
	unsigned char c;
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
	}
	ws[32]='\0';
}


/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
 
int wordcmp(unsigned char *w1, unsigned char *w2, int n) {
    int i = 0, limit;

    if (n % 4 != 0) {
        w1[n / 4] = w1[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w1[n / 4] = w1[n / 4] << (2 * (3 - ((n - 1) % 4)));
        w2[n / 4] = w2[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w2[n / 4] = w2[n / 4] << (2 * (3 - ((n - 1) % 4)));
        limit = (n / 4) + 1;
    } else {
        limit = n / 4;
    }

    for (i = 0; i < limit; i++) {
        if (w1[i] < w2[i]) return -1;
        if (w1[i] > w2[i]) return +1;
    }
    return 0;
}


/* This function is used to check if a file exists or not.
 *  @param file is a string with the absolute/relative path to the file.
 *  @return a positive number if the file exists and the program have access
 *          or zero in other cases.
 */
int exists(char *file) {
    if (access(file, F_OK) != (-1)) return 1;
    else return 0;
}

/* This function is used to check if a string given is an integer.
 *  @param str is the string to be checked.
 *  @return a positive number if it's an integer or zero in other cases. 
 */
int is_int(char const *str) {
    int integer = atoi(str); // Return the first integer found on the string
    char str2[1024];
    sprintf((char *) &str2, "%d", integer); // int -> string
    int isInteger = strcmp(str2, str) == 0; // Check if there are equals => String==Integer
    return isInteger;
}


/* This function is used to check if a string given is a float.
 *  @param str is the string to be checked.
 *  @return a positive number if it's an float or zero in other cases. 
 */
int is_float(char const *str) {
    float number = atof(str); // Return the first float found on the string
    char str2[1024];
    sprintf((char *) &str2, "%f", number); // float -> string
    int isFloat = strcmp(str2, str) == 0; // Check if there are equals => String==Float
    return isFloat;
}
