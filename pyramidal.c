#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define PACK_T uint64_t
#define PACK_SIZE (sizeof(PACK_T) * 8)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
void ext_shift(PACK_T ch, uint8_t pos, PACK_T retval[2])
{
    retval[0] = ch << pos;
    if (pos==0)
        retval[1] = 0;
    else 
        retval[1] = ch >> (PACK_SIZE-pos);
    return;
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

uint32_t popcnt(const uint64_t* buf, int len) {
  uint64_t cnt[4];
  for (int i = 0; i < 4; ++i) {
    cnt[i] = 0;
  }

  for (int i = 0; i < len; i+=4) {
    __asm__(
        "popcnt %4, %4  \n\t"
        "add %4, %0     \n\t"
        "popcnt %5, %5  \n\t"
        "add %5, %1     \n\t"
        "popcnt %6, %6  \n\t"
        "add %6, %2     \n\t"
        "popcnt %7, %7  \n\t"
        "add %7, %3     \n\t" // +r means input/output, r means intput
        : "+r" (cnt[0]), "+r" (cnt[1]), "+r" (cnt[2]), "+r" (cnt[3])
        : "r"  (buf[i]), "r"  (buf[i+1]), "r"  (buf[i+2]), "r"  (buf[i+3]));
  }
  return cnt[0] + cnt[1] + cnt[2] + cnt[3];
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
    PACK_T shifted[2];
    unsigned long long i;
    long long off_the_end;
    off_the_end = length + segoffset - offset - 1;
    if (off_the_end < 0) printf("Failure!!!!!\n");
    long long stop = MIN(length-1, off_the_end);
    const unsigned long long access_offset = offset - segoffset;
    const long long start_offset = MAX(0, (long long)segoffset - (long long)offset);
    
    if (start_offset > 0) {
        i = start_offset - 1;
        ext_shift(source[i], bitoff, shifted);
        target[i + 1 + access_offset] |= shifted[1];
    }

    for (i = start_offset; i < stop; i++) {
        ext_shift(source[i], bitoff, shifted);
        target[i + access_offset] |= shifted[0];
        target[i + 1 + access_offset] |= shifted[1];
    }
    if (i + access_offset > 0) {
        i = stop;
        ext_shift(source[i], bitoff, shifted);
        target[i + access_offset] |= shifted[0];
    }
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

unsigned long long pyr(unsigned long long p)
{
    return p*(p+1)*(p+2)/6;
}

void genn3(PACK_T *set, unsigned long long length, unsigned long long offset)
{
    unsigned long long i, j, k, pyri, pyrj, pyrk;
    unsigned long long numelements = length * PACK_SIZE;
    unsigned long long maxi = offset + numelements;
    unsigned long long maxp = invpyr(maxi);
    for (i = 1; i <= invpyr(maxi/2); ++i) {
        pyri = pyr(i);
        for (j = i; j <= invpyr(maxi/2); ++j) {
            pyrj = pyr(j);
            unsigned long long mink = j, pij = pyri + pyrj;
            if (pij < offset) mink = invpyr(offset - pij);
            if (mink>j) mink--;
            for (k = mink; k < maxp; ++k) {
                pyrk = pyr(k);
                unsigned long long sum = pij + pyrk;
                if (sum > maxi) break;
                if (sum < offset + 1) continue;
                mark(set, sum - offset - 1);
            }
        }
    }
}

void gennext(PACK_T *srcset, PACK_T *dstset, unsigned long long length, unsigned long long offset)
{
    unsigned long long i, pyr;
    unsigned long long numelements = length * PACK_SIZE;
    unsigned long long maxp = invpyr(offset+numelements) + 1;
    unsigned long long deltap = invpyr(offset);
    for(i=2; i < maxp; i++) {
	pyr = pyramid(i);
        if (pyr + numelements + 1 < offset)
            continue;
        if (pyr > numelements + offset)
            break;
        shift_and_mark(srcset, dstset, length, pyr, offset/PACK_SIZE);
        if (i>=deltap) {
            if ((i-deltap) % 50 == 0) {
                unsigned long long cnt = popcnt(dstset, length);
                if (cnt == numelements) {
                    break;
                }
            }
        }
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
    unsigned long long numsegs, segelements, seglength, offset;
    segelements = 1000000000ULL;
    seglength = segelements / PACK_SIZE;
    /* By default assume we want 1,000,000,000 digits checked*/
    if (argc == 1) {
        /*bit 0 in our packing represents 1, bit 1 represents 2 etc*/ 
        numsegs = 1;
    }
    else if ((argc == 2) && (atoll(argv[1]) > 0)) {
        numsegs = atoll(argv[1]);
    }
    else if ((argc == 3) && (atoll(argv[1]) > 0) && (atoll(argv[2]) % 256 == 0)) {
        numsegs = atoll(argv[1]);
        segelements = atoll(argv[2]);
        seglength = segelements / PACK_SIZE;
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
    
    unsigned long long int n1, n2, n3, n4;

    genn1(tmpset, seglength, offset);
    n1 = popcnt(tmpset, seglength);
    printf("N1 number in segment %llu: %llu\n", numsegs, n1);
    fflush(stdout);
    
    /* 2 */
    genn2(tmpset, seglength, offset);
    n2 = popcnt(tmpset, seglength) - n1;
    printf("N2 number in segment %llu: %llu\n", numsegs, n2);
    fflush(stdout);
    
    /* 3 */
    memcpy(numberset, tmpset, seglength * sizeof(PACK_T));
    genn3(numberset, seglength, offset);
    n3 = popcnt(numberset, seglength) - n2 - n1;
    printf("N3 number in segment %llu: %llu\n", numsegs, n3);
    fflush(stdout);
    
    /* 4 */
    memset(tmpset, 0, seglength * sizeof(PACK_T));
    genn1(tmpset, seglength, offset+segelements);
    genn2(tmpset, seglength, offset+segelements);
    gennext(numberset, tmpset, seglength, segelements);
    n4 = popcnt(tmpset, seglength);
    printf("N4 coverage of the next segment using fast algorithm %llu\n", n4);
    fflush(stdout);
    
    free(numberset);
    free(tmpset);
    return 0;
}
