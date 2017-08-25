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

//First runs the sieve of Eratosthenes in parallel, and then checks in parallel
//whether every even number 2 < k <= n is the sum of two primes. See documentation for
//exact algorithm.
//Input: none
//Output: none
void bspgoldbach()
{
    int *numlist,n,i,j,p,s,b,n_loc,n_init,offset,sqrt_n;
    int *prime_list,prime_count,cur_prime_sq;
    int *full_numlist,sqrt_p,*row_primes,*column_primes,prow,pcolumn;
    int b2,sum,counterex,endpoint;
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

    //Prepare the list that will contain the prime numbers on our own processor.
    //Observe that there are less than b local prime numbers, so we will use
    //this as length of the array.
    sqrt_n = (int) sqrt(n);
    prime_list = malloc(b*sizeof(int));
    bsp_push_reg(&prime_count,sizeof(int));
    bsp_push_reg(prime_list,b*sizeof(int));

    //Prepare the arrays that will contain the prime numbers in our row and
    //column to check the Goldbach conjecture.
    sqrt_p = (int) sqrt(p);
    row_primes = malloc(sqrt_p*b*sizeof(int));
    column_primes = malloc(sqrt_p*b*sizeof(int));
    for (i=0;i<sqrt_p*b;i++) {
        row_primes[i]=0;
        column_primes[i]=0;
    }
    bsp_push_reg(row_primes,sqrt_p*b*sizeof(int));
    bsp_push_reg(column_primes,sqrt_p*b*sizeof(int));

    //Initialize numlist.
    full_numlist = malloc((n+1)*sizeof(int));
    numlist = &full_numlist[n_init];
    bsp_push_reg(numlist,n_loc*sizeof(int));
    for (i=0;i<n_loc;i++)
        numlist[i]=1;
    bsp_sync();

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

    //We only need the prime number 2 to write 4 as the sum of two primes;
    //every other sum of 2 plus a different prime is odd. However, 2 plus a
    //different prime could again be a prime (which cannot happen for odd
    //primes, since their sum is even). To avoid improper behaviour we declare 2
    //to be non-prime for now.
    if (s == 0)
        numlist[2]=0;

    //Determine our processor row and column
    prow = s/sqrt_p;
    pcolumn = s%sqrt_p;

    //Determine the prime numbers on our processor
    prime_count=0;
    for (i=0;i<n_loc;i++)
    if (numlist[i]) {
        prime_list[prime_count]=n_init + i;
        prime_count++;
    }

    //Broadcast the prime numbers to the correct processors.
    //Broadcast to the correct row
    for (i=prow*sqrt_p;i < (prow+1)*sqrt_p; i++)
        bsp_put(i,prime_list,row_primes,pcolumn*b*sizeof(int),
            prime_count*sizeof(int));

    //Broadcast to the correct column
    //It might seem off that we start at s/sqrt_p; this is because we combine
    //sqrt_p parts of the original block distribution into one new part.
    for (i=s/sqrt_p;i < p;i += sqrt_p)
        bsp_put(i,prime_list,column_primes,pcolumn*b*sizeof(int),
            prime_count*sizeof(int));
    bsp_sync();

    //Add our primes
    //Since summation is commutative we only use half of one of the lists if we
    //are not on the diagonal (i.e. prow != pcolumn), and on the diagonal we
    //make sure not to add primes twice.
    b2 = b*sqrt_p;
    if (prow < pcolumn)
        endpoint = b2/2;
    else
        endpoint = b2;
    for (i=0;i<b2;i++)
        if (row_primes[i] != 0) {
            if (pcolumn < prow)
                j = b2/2;
            else if (pcolumn == prow)
                j = i;
            else j=0;
            for (j=0;j<b2;j++)
                if (column_primes[j] != 0) {
                    sum = row_primes[i]+column_primes[j];
                    if (sum <= n)
                        full_numlist[sum]=2;
                }
        }

    //Put the result back in the original block distribution used for the sieve.
    for (i=0;i<=n;i++)
        if (full_numlist[i] == 2)
            bsp_put(i/b,&full_numlist[i],numlist,(i%b)*sizeof(int),sizeof(int));
    bsp_push_reg(&counterex,sizeof(int));
    bsp_sync();

    //Determine the index to start checking from. Recall that we marked 2 as
    //non-prime, therefore 4=2+2 is an exception and we start at 6.
    if (n_init < 6)
        i = 6;
    else if (n_init % 2 == 0)
        i = 0;
    else
        i = 1;

    //Check if every even number >=6 was marked
    for (;i<n_loc && numlist[i]==2;i += 2) {}
    counterex = 0;
    if (i < n_loc) {
        counterex = n_init + i;
        bsp_put(0,&counterex,&counterex,0,sizeof(int));
    }
    bsp_sync();
    if (s == 0) {
        if (counterex)
            printf("Goldbach Conjecture is false, counterexample: %d\n",
                counterex);
        else
            printf("Goldbach Conjecture verified up to %d\n",n);
    }
bsp_end();
}

int main (int argc, char * argv[]) {
    bsp_init(bspgoldbach, argc, argv);
    int sqrt_P;

    // Sequential part
    printf("How many processors?\n");
    scanf("%d",&P);
    if (P > bsp_nprocs()) {
        printf("Sorry, not enough available.\n");
        exit(1);
    }
    sqrt_P = (int) sqrt(P);
    if (sqrt_P*sqrt_P != P) {
        printf("The number of processors should be a square.\n");
        exit(1);
    }
    // Parallel part
    bspgoldbach();

    // Sequential part
    exit(0);
}