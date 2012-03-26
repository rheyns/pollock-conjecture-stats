#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define PACK_T uint64_t
#define PACK_SIZE (sizeof(PACK_T) * 8)

unsigned long long pyramid(unsigned long long p)
{
    return (p*(p*p - 1))/6;
}

PACK_T* ext_shift(PACK_T ch, uint8_t pos)
{
    static PACK_T retval[2];
    retval[0] = ch << pos;
    if (pos==0) retval[1] = 0;
    else retval[1] = ch >> (PACK_SIZE-pos);
    return retval;
}

void mark(PACK_T *set, unsigned long long pos)
{
    uint8_t bitsize = PACK_SIZE;
    uint8_t setbit = pos % bitsize;
    unsigned long long int offset = pos / bitsize;
    set[offset] |= 1ULL<<setbit;
    return;
}

PACK_T test(const PACK_T *set, unsigned long long pos)
{
    uint8_t bitsize = PACK_SIZE;
    uint8_t testbit = pos % bitsize;
    unsigned long int offset = pos / bitsize;
//    printf("testbit %d, offset %ld, s 0x%llX, b 0x%llX\n", testbit, offset, set[offset], 1LL<<testbit);
    return set[offset] & 1ULL<<testbit;
}

void shift_and_mark(PACK_T *source, PACK_T *target, unsigned long long length, unsigned long long pos)
{
    uint8_t bitsize = PACK_SIZE;
    uint8_t bitoff = pos % bitsize;
    unsigned long long offset = pos / bitsize;
    PACK_T *shifted;
    unsigned long long i;
    
    for (i = 0; i < length; i++)
    {
        if(i + offset >= length) break;
        shifted = ext_shift(source[i], bitoff);
        target[i + offset] |= shifted[0];
        if(i + offset + 1 >= length) break;
        target[i + offset + 1] |= shifted[1];
    }
    return;
}

void pr(const PACK_T *set, unsigned long long length)
{
    unsigned long long i;
    for(i=0;i<length*PACK_SIZE;i++)
    {
        if(test(set, i)) printf("X, %lld\n", i);
       // else printf("O, %ld\n", i);
    }
}

void save(char *fname, PACK_T *set, unsigned long long length)
{/*
    FILE *ptr_file;
    unsigned long long i;
    ptr_file = fopen(fname, "wb");
    if (!ptr_file)
    {
        printf("Unable to open file!");
    }
    else
    {
        fwrite(set, PACK_SIZE, length, ptr_file);
    }
    fclose(ptr_file);*/
    return;
}

/*
    Overall strategy is to first mark the integers which are pyramidal
    then the integers which are the sum of two pyramidals. For sums
    of pyramidals above 3 the summations will be performed on batches
    using bit shifts in order to realize a speedup roughly proportional
    to the size of the largest machine integer. For example with a machine
    integer of size 64-bits the speedup should be roughly 64 times that of
    using integer multiplication.
*/
int main(int argc, char *argv[])
{
    unsigned long long numints,length, maxp;
    if (argc != 2) /* By default assume we want 1,000,000 digits checked*/
    {
        /*We add one since indexing is based on the cardinal value of numbers, in other
        words bit 0 in our packing always remains empty since the first pyramidal number is one.
        Therefor the cardinality for the set of 1000000 numbers is 1000001.*/ 
        numints = 1000000+1;
        length = numints / PACK_SIZE;
        if(numints % PACK_SIZE != 0) length++;
        maxp = (long long) floor(6.0*numints)+1ULL;
    }
    else if ((argc == 2) && (atoll(argv[1]) != 0))
    {
        /*We add one since indexing is based on the cardinal value of numbers, in other
        words bit 0 in our packing always remains empty since the first pyramidal number is one.
        Therefor the cardinality for argv[1] numbers is argv[1]+1*/ 
        numints = atoll(argv[1]) + 1;
        length = numints / PACK_SIZE;
        if (numints % PACK_SIZE != 0) length++;
        maxp = (long long) floor(6.0*numints)+1ULL;  
    }
    else
    {
        printf("Usage is %s length: Where length is the number of integers to check.\n", argv[0]);
    }
    /* numberset is a large bit field representing the set of numbers representable
    as the sum of fewer than 5 pyramidal numbers
    */

    printf("length %lld\n", length);
    PACK_T *numberset = calloc(length, sizeof(PACK_T));
    PACK_T *tmpset = calloc(length, sizeof(PACK_T));
    if (numberset == NULL | tmpset == NULL)
    {
        printf("Failed to allocate memory. \n");
        return -1;
    }
    
    unsigned long long int i, j, tmp, n1, n2, n3, n4, n5;
    for(i=2;i<=maxp;i++)
    {
        if (pyramid(i)>=length*PACK_SIZE) break;
        mark(numberset, pyramid(i));
    }
    printf("First stage done.\n");
    
    j = 0;
    for(i = 0; i < numints; i++)
    {
        if(test(numberset,i)) j++;
        if (i%1000000000ULL == 0) printf("N1 %lld : %lld\n",i,j);
    }
    n1 = j;
    printf("Number in N1: %lld\n", n1);
    
    save("N1.bin", numberset, length);
    
    /* 2 */
    for(i=2;i<=maxp;i++)
    {
        if (pyramid(i) >= numints) break;
        for(j=2; j <= maxp; j++)
        {
            tmp = pyramid(i) + pyramid(j);
            if (tmp >= numints) break;
            mark(numberset, tmp);
        }
    }
    
    printf("Second stage done.\n");
    
    j = 0;
    for(i=0;i < numints; i++)
    {
        if(test(numberset,i)) j++;
        if (i%1000000000ULL == 0) printf("N2 %lld : %lld\n",i,j);
    }
    n2 = j - n1;
    printf("Number in N2: %lld\n", n2);
    
    save("N2N1.bin", numberset, length);
    
    /* 3 */
    for(i=2; i < maxp; i++)
    {
        if (pyramid(i) >= numints) break;
        shift_and_mark(numberset, tmpset, length, pyramid(i));
    }
    

    for(i=0; i<length; i++)
    {
        numberset[i] |= tmpset[i];
    }
    printf("Third stage done.\n");
    
    j = 0;
    for(i=0; i < numints; i++)
    {
        if(test(numberset,i)) j++;
        if (i%1000000000ULL == 0) printf("N3 %lld : %lld\n",i,j);
    }
    n3 = j - n2 - n1;
    printf("Number in N3: %lld\n", n3);
    
    save("N3N2N1.bin", numberset, length);
    
    /* 4 */
    for(i=2; i < maxp; i++)
    {
        if (pyramid(i) >= numints) break;
        shift_and_mark(numberset, tmpset, length, pyramid(i));
    }
    
    for(i=0; i<length; i++)
    {
        numberset[i] |= tmpset[i];
    }
    
    printf("Fourth stage done.\n");
    
    j = 0;
    for(i=0; i< numints; i++)
    {
        if(test(numberset,i)) j++;
        if (i%1000000000ULL == 0) printf("N4 %lld : %lld\n",i,j);
    }
    n4 = j - n3 - n2 - n1;
    printf("Number in N4: %lld\n", n4);
    
    save("N4N3N2N1.bin", numberset, length);
    
    /* Count the unmarked numbers, should be 241*/
    j = 0;
    for(i=1; i < numints; i++)
    {
        if(!test(numberset, i)) j++;
    }
    n5 = j;
    
    printf("Number of integers requiring 5 pyramidal integers, %lld\n", n5);
    
    free(numberset);
    free(tmpset);
    return 0;
}
