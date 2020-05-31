//
//  BaransGreensFn.h
//  
//
//  Created by Baran Serajelahi  on 2019-06-21.
//
// August 1, 2019. Next part of the code.

#define DEBG 0
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;} // for gausj
#define TINY 1.0e-20 // for ludcmp()

//Function declerations included in Secomb et al main.cpp
/*
 void input(void);
 void analyzenet(void);
 void picturenetwork(float *nodvar, float *segvar, const char fname[]);
 void greens(void);
 void contour(const char fname[]);
 void histogram(const char fname[]);
 void setuparrays0();
 void setuparrays1(int nseg, int nnod);
 void setuparrays2(int nnv, int nnt);
 void cmgui(float *segvar);
 void postgreens(void);
 */

/* The following is a list of unique headers from all of the .cpp files supplied by Secomb et al.
 
 These are the standard header includes in Baran codes using the green function
 given by Secomb.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <calloc.h>
// #include <stddef.h>
#include "nrutil.h"
#include <string>// I was able to add this in when I changed the file names from .c to .cpp
// I needed to change file names to .cpp from .c so that I could avoid the following error: Compiler error: “initializer element is not a compile-time constant”
// The headers below give fatal errors

/*
#include <string>
#include <Windows.h>
#include <cstdio>
#include <malloc.h>
*/



