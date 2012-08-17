#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define PACK_T uint64_t
#define PACK_SIZE (sizeof(PACK_T) * 8)

unsigned long long cuberoot(unsigned long long n)
{
    unsigned long long ret = pow(n + 0.5, 1.0/3);
    if (n < 100000000000001ULL)
        return ret;
    if (n >= 18446724184312856125ULL)
        return 2642245ULL;
    if (ret * ret * ret > n) {
        ret--;
        while (ret * ret * ret > n)
            ret--;
        return ret;
    }
    while ((ret + 1) * (ret + 1) * (ret + 1) <= n)
        ret++;
    return ret;
}

/* The inverse function of the pyramidal numbers is well
bounded below by the cube root function.*/
unsigned long long invpyr(unsigned long long n)
{
    return cuberoot(n*6)+1;
}

/* Generate the p'th pyramidal number.*/
unsigned long long pyramid(unsigned long long p)
{
    return (p*(p*p - 1))/6;
}

/* Extended shift by pos:
    This function extends bitwise left shift to a type of size 2*sizeof(PACK_T).
    This function returns two PACK_T integers where retval[1] is ch shifted
    left by pos and retval[0] is the result of the overflow off the left of
    ch << pos.*/
PACK_T* ext_shift(PACK_T ch, uint8_t pos)
{
    static PACK_T retval[2];
    retval[0] = ch << pos;
    if (pos==0)
        retval[1] = 0;
    else 
        retval[1] = ch >> (PACK_SIZE-pos);
    return retval;
}

/* Mark a particular bit in the array set. */
void mark(PACK_T *set, unsigned long long pos)
{
    uint8_t setbit = pos % PACK_SIZE;
    unsigned long long int offset = pos / PACK_SIZE;
    set[offset] |= 1ULL<<setbit;
    return;
}

/* Test a particulary bit in the array set. */
PACK_T test(const PACK_T *set, unsigned long long pos)
{
    uint8_t testbit = pos % PACK_SIZE;
    unsigned long int offset = pos / PACK_SIZE;
    return set[offset] & 1ULL<<testbit;
}

/* Shifts the entire array source by pos positions and set the bits of target
    which correspond to the shifted 1 bits of set. set and target are offset by
    segoffset*PACK_SIZE bits and bits which are 1 in target are preserved.
    Example: PACT_T = unit8_t
             length = 2;
             set[0] = 0b00001010;
             set[1] = 0b11111111;
             target[0] = 0b11111111;
             target[1] = 0b00000000;
             shift = 12;
             segoffset = 2;
             Resultant target[0] = 0b11111111;
             target[1] = 0b11110000;
             */
void shift_and_mark(PACK_T *source, PACK_T *target, unsigned long long length, unsigned long long shift, unsigned long long segoffset)
{
    uint8_t bitoff = shift % PACK_SIZE;
    unsigned long long offset = shift / PACK_SIZE;
    PACK_T *shifted;
    unsigned long long i;
    
    for (i = 0; i < length; i++) {
        if(i + offset >= length + segoffset)
            break;
        shifted = ext_shift(source[i], bitoff);
        if(i + offset >= segoffset)
            target[(i + offset) - segoffset] |= shifted[0];
        if(i + offset + 1 >= length + segoffset)
            break;
        if(i + offset + 1 >= segoffset)
            target[(i + offset + 1) - segoffset] |= shifted[1];
    }
    return;
}

/* A fast method for doing population counts for bits. Equivalent to gcc's
   __builtin_popcountll instruction. Needs to be changed if PACK_T is not 
   uint64_t. */
uint64_t pop4(PACK_T *elements)
{
    uint64_t x, y, u, v;
    x = elements[0];
    y = elements[1];
    u = elements[2];
    v = elements[3];
    enum { m1 = 0x5555555555555555, 
         m2 = 0x3333333333333333, 
         m3 = 0x0F0F0F0F0F0F0F0F, 
         m4 = 0x000000FF000000FF };
    
    x = x - ((x >> 1) & m1);
    y = y - ((y >> 1) & m1);
    u = u - ((u >> 1) & m1);
    v = v - ((v >> 1) & m1);
    x = (x & m2) + ((x >> 2) & m2);
    y = (y & m2) + ((y >> 2) & m2);
    u = (u & m2) + ((u >> 2) & m2);
    v = (v & m2) + ((v >> 2) & m2);
    x = x + y; 
    u = u + v; 
    x = (x & m3) + ((x >> 4) & m3);
    u = (u & m3) + ((u >> 4) & m3);
    x = x + u; 
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x & m4; 
    x = x + (x >> 32);
    return x & 0x00000000000001FF;
}

/* Count the number of set bits in set. */
unsigned long long count_set(PACK_T *set, unsigned long long length)
{
    unsigned long long int i,j;
    j = 0;
    for(i = 0; i < length; i += 4)
        j += pop4(set+i);
    return j;
}

/* Save the results of set in a file.*/
void save(char *fname, PACK_T *set, unsigned long long length)
{
    FILE *ptr_file;
    unsigned long long i;
    ptr_file = fopen(fname, "wb");
    if (!ptr_file)
        printf("Unable to open file!");
    else
        fwrite(set, sizeof(PACK_T), length, ptr_file);
    fclose(ptr_file);
    return;
}

void genn1(PACK_T *set, unsigned long long length, unsigned long long offset)
{
    unsigned long long i;
    unsigned long long numelements = length * PACK_SIZE;
    unsigned long long maxp = invpyr(offset+numelements) + 1;
    for(i=2;i<=maxp;i++) {
        if (pyramid(i) < offset + 1)
            continue;
        if (pyramid(i) > offset + numelements)
            break;
        mark(set, pyramid(i) - offset - 1);
    }
    return;
}

void genn2(PACK_T *set, unsigned long long length, unsigned long long offset)
{
    unsigned long long i,j,tmp;
    unsigned long long numelements = length * PACK_SIZE;
    unsigned long long maxp = invpyr(offset+numelements) + 1;
    for(i=2;i<=maxp;i++) {
        if (pyramid(i) > offset + numelements)
            break;
        for(j=2; j <= maxp; j++) {
            tmp = pyramid(i) + pyramid(j);
            if (tmp < offset + 1)
                continue;
            if (tmp > offset + numelements) 
                break;
            mark(set, tmp - offset - 1);
        }
    }
    return;
}

void gennext(PACK_T *srcset, PACK_T *dstset, unsigned long long length, unsigned long long offset)
{
    unsigned long long i,j,tmp;
    unsigned long long numelements = length * PACK_SIZE;
    unsigned long long maxp = invpyr(offset+numelements) + 1;
    for(i=2; i < maxp; i++) {
        if (pyramid(i) + numelements + 1 < offset)
            continue;
        if (pyramid(i) > numelements + offset)
            break;
        shift_and_mark(srcset, dstset, length, pyramid(i), offset/PACK_SIZE);
    }
}

/*
    This program aims to be an aid for studying the Pollock Conjecture.
    Integers are represented as positions in a bit array consisting of elements
    of type PACK_T. In practice PACK_T should be the largest representable 
    unsinged integer type on any given architecture. 
    Pyramidal (more correctly tetrahedral) numbers are a count of the number of
    balls in a 3-D triangular pyramid .
    Overall strategy is to first mark the integers which are pyramidal
    then the integers which are the sum of two pyramidals. For sums
    of pyramidals above 2 the summations will be performed on batches
    using bit shifts in order to realize a speedup roughly proportional
    to the size of the largest machine integer. For example with a machine
    integer of size 64-bits the speedup should be roughly 64 times that of
    using integer multiplication.
*/
int main(int argc, char *argv[])
{
    unsigned long long numints,length,numsegs, segelements, seglength, maxp, offset;
    segelements = 1000000000ULL;
    seglength = segelements / PACK_SIZE;
    /* By default assume we want 1,000,000,000 digits checked*/
    if (argc == 1) {
        /*bit 0 in our packing represents 1, bit 1 represents 2 etc*/ 
        numsegs = 1;
        maxp = invpyr(segelements) + 1;
    }
    else if ((argc == 2) && (atoll(argv[1]) > 0)) {
        numsegs = atoll(argv[1]);
        maxp = invpyr(numsegs*segelements) + 1;
    }
    else if ((argc == 3) && (atoll(argv[1]) > 0) && (atoll(argv[2]) % 256 == 0)) {
        numsegs = atoll(argv[1]);
        segelements = atoll(argv[2]);
        seglength = segelements / PACK_SIZE;
        maxp = invpyr(numsegs*segelements) + 1;
    }
    else {
        printf("Usage is %s blocknum: Where blocknum is the identifier of a billion\
        integer block to check Pollock's conjecture on.\n\n\
        \nExample: %s 1\n\nComputes the relevant statistics for the first billion integers.\n\
        Note that the statistics are for a given block of integers only.\n\
        Alternate usage is %s blocknum blocklen: Where blocklen is the length of each block\n\
        and blocklen must be divisible by 256.\n", argv[0], argv[0], argv[0]);
        return 0;
    }
    offset = (numsegs - 1) * segelements;
    
    /* numberset is a large bit field representing the set of numbers representable
    as the sum of fewer than 5 pyramidal numbers
    */
    PACK_T *numberset = calloc(seglength, sizeof(PACK_T));
    PACK_T *tmpset = calloc(seglength, sizeof(PACK_T));
    if ((numberset == NULL) | (tmpset == NULL)) {
        printf("Failed to allocate memory. \n");
        return -1;
    }
    
    unsigned long long int i, j, tmp, n1, n2, n3, n4, n5, curroffset, prev, prev2;

    genn1(tmpset, seglength, offset);
    n1 = count_set(tmpset, seglength);
    printf("N1 number in segment %llu: %llu\n", numsegs, n1);
    fflush(stdout);
    
    /* 2 */
    genn2(tmpset, seglength, offset);
    n2 = count_set(tmpset, seglength) - n1;
    printf("N2 number in segment %llu: %llu\n", numsegs, n2);
    fflush(stdout);
    
    /* 3 */
    for(i=0; i < numsegs*2/3 + 1; i++) {
        curroffset = i * segelements;
        memset(numberset, 0, seglength * sizeof(PACK_T));
        genn1(numberset, seglength, curroffset);
        genn2(numberset, seglength, curroffset);
        gennext(numberset, tmpset, seglength, offset - curroffset);
    }
    n3 = count_set(tmpset, seglength) - n2 - n1;
    printf("N3 number in segment %llu: %llu\n", numsegs, n3);
    fflush(stdout);
    
    /* 4 */
    memset(numberset, 0, seglength * sizeof(PACK_T));
    genn1(numberset, seglength, offset+segelements);
    genn2(numberset, seglength, offset+segelements);
    gennext(tmpset, numberset, seglength, segelements);
    n4 = count_set(numberset, seglength);
    printf("N4 coverage of the next segment using fast algorithm %llu\n", n4);
    fflush(stdout);
    
    free(numberset);
    free(tmpset);
    return 0;
}
