
/* 	This function compares two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param bytes_to_check: length of BOTH arrays in bytes
 *  @retun a positive number if w1>w2, a negative number if w1<w2 and zero if they are equal.
 */
 
int wordcmpbytearray(const unsigned char *w1, const unsigned char *w2, int bytes_to_check) {
	int i;
	for (i=0;i<bytes_to_check;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}

/* 	This function is used to shift left bits in a unsigned char array
 *	@b: char array representing the word using compressed 2-bit characters
 */
 
void shift_byte_array_left(unsigned char * b, unsigned int bytes_to_check) {
    unsigned int i;
    for (i = 0; i < bytes_to_check - 1; i++) {
        b[i] <<= 2;
        b[i] |= (b[i + 1] >> 6);
    }
    b[bytes_to_check - 1] <<= 2;
}

/* 	This function is used to shift right bits in a unsigned char array
 *	@b: char array representing the word using compressed 2-bit characters
 */
void shift_byte_array_right(unsigned char * b, unsigned int bytes_to_check) {
    unsigned int i;
    for (i = bytes_to_check - 1; i > 0; i--) {
        b[i] >>= 2;
        b[i] |= (b[i - 1] << 6);
    }
    b[i] >>= 2;
}