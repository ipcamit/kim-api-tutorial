/* main.c : behaves like LAMMPS*/
#include "driver.h"
#include <stdio.h>

int main(void)
{
    void* ctx = NULL;
    compute_fn fn = make_driver(5, &ctx);   /* factory gives fn ptr + ctx */
    int y = fn(ctx, 3);                     /* call through adapter */
    printf("%d\n", y);                      /* --> 8 */
    return 0;
}

// g++ -std=c++17 -c driver.cpp 
// gcc  -std=c11   -c call.c    
// g++ driver.o call.o -o demo  
// ./demo                       
