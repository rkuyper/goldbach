#include "include/mcbsp.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Perform the sieve of Eratosthenes (sequentially).
//Input:  integer n, integer array numlist of length n+1 initialized to
//        all ones.
//Output: numlist is modified such that numlist[i]==1 if i is prime
//        and numlist[i]==0 if i is not prime.
void seqsieve (int n, int *numlist) {
    int sqrt_n,i,j;
    //Perform the sieve.
    sqrt_n = (int) sqrt(n);
    numlist[0] = 0;
    numlist[1] = 0;
    for (i=2;i<=sqrt_n;i++)
        if (numlist[i])
            for (j=i*i;j <= n; j += i)
                numlist[j] = 0;
}

int P;
//Block distribution: determine the size of the block on processor s.
//for an array of length n.
//Input:  integer b for the block size, integer s for the processor number,
//        integer n for the array length.
//Output: the size of the block on processor s.
int block_size (int b, int s, int n) {
    if ((s+1)*b <= n)
        return b;
    else if (n > s*b)
        return n - s*b;
    else
        return 0;
}

//Block distribution: determine the global index of the local index 0.
//Input:  integer b for the block size, integer s for the processor number,
//        integer n for the array length.
//Output: the global index of the local index 0; we return -1 if the block has
//        length 0.
int block_init (int b, int s, int n) {
    if (s*b < n)
        return s*b;
    else
        return -1;
}

//Perform the sieve of Eratosthenes in parallel, using a block distribution.
//At the end of the function, the array numlist is distributed over the
//processors in a block distribution. We have, for global indices i, that at
//the end numlist[i]==1 if i is prime and numlist[i]==0 if i is not prime.
//We first run a sequential sieve on processor 0 up to sqrt(n), broadcast the
//results and then perform the crossouts in parallel.
//Input: none
//Output: none
void bspsieve()
{
    int *numlist,n,i,j,p,s,b,n_loc,n_init,offset,sqrt_n;
    int *prime_list,prime_count,cur_prime_sq;
    bsp_begin(P);
    p = bsp_nprocs();
    s = bsp_pid();
    if (s == 0) {
        printf("What is n?\n");
        scanf("%d", &n);
        if (n < 1)
            bsp_abort("n needs to be a positive integer.\n");
        if (n < p*p)
            bsp_abort("n should be at least the square  of p.\n");
    }
    //Register and broadcast n
    bsp_push_reg(&n,sizeof(int));
    bsp_sync();
    bsp_get(0,&n,0,&n,sizeof(int));
    bsp_sync();
    //Determine the block size b
    if ((n+1) % p == 0)
        b = (n+1) / p;
    else
        b = (n+1) / p + 1;
    //Determine the local block size and the global index of the local index 0
    n_loc = block_size(b,s,n+1);
    n_init = block_init(b,s,n+1);
    //Prepare the list that will contain the prime numbers <= sqrt(n), i.e. the
    //result of the first sequential sieve. Observe that there are less than
    //sqrt(n) prime numbers <= sqrt(n), so we will use this as length of
    //the array.
    sqrt_n = (int) sqrt(n);
    prime_list = malloc(sqrt_n*sizeof(int));
    bsp_push_reg(&prime_count,sizeof(int));
    bsp_push_reg(prime_list,sqrt_n*sizeof(int));
    bsp_sync();
    //Initialize numlist.
    numlist = malloc(n_loc*sizeof(int));
    for (i=0;i<n_loc;i++)
        numlist[i]=1;
    //Perform the first sequential sieve up to sqrt(n) and broadcast the
    //results.
    if (s == 0) {
        seqsieve(sqrt_n,numlist);
        //Put the prime numbers up to sqrt(n) in a list
        prime_count = 0;
        for (i=0;i <= sqrt_n; i++)
            if (numlist[i]) {
                prime_list[prime_count] = i;
                prime_count++;
            }
        //Put the results in all processors.
        for (i=0; i < p; i++) {
            bsp_put(i,&prime_count,&prime_count,0,sizeof(int));
            bsp_put(i,prime_list,prime_list,0,prime_count*sizeof(int));
        }
    }
    bsp_sync();
    //Perform the parallel sieve using the result from the first sequential
    //sieve.
    for (i=0; i < prime_count; i++) {
        //We start crossing out at the square of the prime number.
        cur_prime_sq = prime_list[i]*prime_list[i];
        //The square of the current prime is before our initial index.
        if (cur_prime_sq < n_init) {
            //Calculate the local index of the first multiple of the current
            //prime in our block
            offset = n_init % prime_list[i];
            if (offset == 0)
                j = 0;
            else
                j = prime_list[i] - offset;
        }
        //The square of the current prime is in or beyond our block.
        else
            j = cur_prime_sq - n_init;
        //Perform the cross-outs
        for (;j < n_loc; j += prime_list[i])
            numlist[j] = 0;
    }
    bsp_end();
}

int main (int argc, char * argv[]) {
    bsp_init(bspsieve, argc, argv);

    // Sequential part
    printf("How many processors?\n");
    scanf("%d",&P);
    if (P > bsp_nprocs()) {
        printf("Sorry, not enough available.\n");
        exit(1);
    }

    // Parallel part
    bspsieve();

    // Sequential part
    exit(0);
}