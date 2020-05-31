/*
 Date: 24 June 2019.
 This project is to get the Timothy Secomb codes working on Goldman group
 Mac computers, general use Middlesex linux computers.
 This project has following authors:
 Baran Serajelahi.
 Sanjay R. Kharche.
 Daniel Goldman.

 The original codes were obtained from following link:
 https://physiology.arizona.edu/sites/default/files/version4.0.zip_.txt
*/

// this is the header file where I declare and define all my user functions.
#include "GreensFn.h"
#include "myVariables.c"
// Function declarations (all of these functions appear in sequence at the end of this file
void WriteExnodeHeader(FILE *exnode);
void WriteExelemHeader(FILE *exelem);
float rtflsp(float(*func)(float), float x1, float x2, float xacc);
float rtbis(float(*func)(float), float x1, float x2, float xacc);
float func(float x);
float bloodconc(float p, float h);
float bloodconcp(float p, float h);
void blood(float c, float hem, float *p, float *pp);
void convect(int isp);
void gaussj(double **a, int n, double **b, int m);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
double bicgstab(double **a, double *b, double *x, int n, double eps, int itmax);
float *eval(int slsegdiv, float req, float *x);
void tissrate(int nsp, float *p, float *mtiss, float *mptiss);
void postgreens(void);
void picturenetwork(float *nodvar, float *segvar, const char []);
void cmgui(float *segvar);
void contr_lines(FILE *ofp, int m, int n, float scalefac, int nl,
                 float xmin, float xmax, float ymin, float ymax, float *cl, float **zv);
void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,
                 float xmin, float xmax, float ymin, float ymax, float *cl, float **zv,
                 int showscale, int lowcolor, int hatch, int plotcontour);
void contour(const char fname[]);
void histogram(const char fname[]);
void initgreens();
void putrank(void);

//-----------------------------------------------------------------------------/
int main( int argc, char *argv [] )
{

// these are the main function variables.

    //-----------------------------------------------------------------------------/
    // This is where I begin working on input.cpp
    /*
     The original input file was called SoluteParams.dat.
     This file is for read only at this stage.
     
     We made changes to the file format as follows:
     1. All text was removed;
     2. No numerical values were added or removed;
     3. The sequence of all numerical values in SoluteParams.dat was strictly preserved; and
     4. The new version of our file is called simpleSoluteParams.
     */
    ifp = fopen("simpleSoluteParams", "r");
    if(ifp!=NULL) printf("simpleSoluteParams opened successfully.\n");
    else { printf("could not open params file, see L85 in main. exiting.\n"); exit(-1); }
    
    
    //alternative to the above 3 lines
    //if(ifp = fopen("simpleSoluteParams", "r")) printf("simpleSoluteParams opened successfully.\n");
    //else { printf("could not open params file, see L85 in main. exiting.\n"); exit(-1); }
    
/*
Run greens is an integer that has values 1, 2, or 3. // rungreens is not used in the code anywhere - Baran
g0method is an integer used to choose the method for cal    lating g0.
linmethod is also an integer, it determines the technique used to solved the linear system of equations. linmethod = 3 for bi-conjugate stab method.
*/
    fscanf(ifp, "%d %d %d", &rungreens, &g0method, &linmethod);
//-----------------------------------------------------------------------------/
    if (g0method != 1 && g0method != 2){
        printf("Error: simpleSoluteParams, g0method must be 1 or 2. Exiting. \n"); exit(-54);
}
    if (linmethod < 1 || linmethod > 3){
        printf("*** simpleSoluteParams, linmethod must be 1, 2 or 3. Exiting. \n"); exit(-54);
}
    fscanf(ifp, "%d %d %d", &nmaxvessel, &nmaxtissue, &nmax);
    fscanf(ifp, "%f", &errfac);
    fscanf(ifp, "%f", &lowflowcrit);
    fscanf(ifp, "%f", &p50);
    fscanf(ifp, "%f", &fn);
    fscanf(ifp, "%f", &cs);
    fscanf(ifp, "%f", &alphab);
    fscanf(ifp, "%f", &q0fac);
    fscanf(ifp, "%d %d", &nsp, &ntissparams);
    

    if(DEBG==1){
        printf("%d %d %d\n",rungreens, g0method, linmethod);
        printf("%d %d %d\n",nmaxvessel, nmaxtissue, nmax);
        printf("%f\n", errfac);
        printf("%f\n", lowflowcrit);
        printf("%f\n", p50);
        printf("%f\n", fn);
        printf("%f\n", cs);
        printf("%f\n", alphab);
        printf("%f\n", q0fac);
        printf("%d %d\n", nsp, ntissparams);
    }

//-----------------------------------------------------------------------------/
    // declare memory for permsolute as a 1D array of size nsp
    // the code below replaces: permsolute = ivector(1, nsp);
    permsolute = (int *) calloc(nsp+1,sizeof(int ));
    // I will use pm3 policy to explicitly initialize all my variables and arrays. - BS
    for(i=1;i<=nsp;i++) permsolute[i] = -34;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++) printf("%d %d %d\n", nsp, i, permsolute[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for diffsolute as a 1D array of size nsp.
    // the code below replaces: diffsolute = ivector(1, nsp);
    diffsolute = (int *) calloc(nsp+1,sizeof(int ));
    // I will use pm3 policy to explicitly initialize all my variables and arrays. - BS
    for(i=1;i<=nsp;i++) diffsolute[i] = -33;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++) printf("%d %d %d\n", nsp, i, diffsolute[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for oxygen as a 1D array of size nsp.
    // the code below replaces: oxygen = ivector(1, nsp);
    oxygen = (int *) calloc(nsp+1,sizeof(int ));
    // I will use pm3 policy to explicitly initialize all my variables and arrays. - BS
    for(i=1;i<=nsp;i++) oxygen[i] = -35;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++) printf("%d %d %d\n", nsp, i, oxygen[i]);
    }
    
    //-----------------------------------------------------------------------------/
    // This "vector" is not done using nrutils, there was something that did not work - BS
    // declare memory for pref as a 1D array of size nsp
    // the code below replaces pref = vector(1,nsp);
    pref = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pref[i] = -234.04;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, pref[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for diff as a 1D array of size nsp
    diff = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) diff[i] = 1223.01;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, diff[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for g0 as a 1D array of size nsp
    g0 = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) g0[i] = 14.01;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, g0[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for g0fac as a 1D array of size nsp
    g0fac = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) g0fac[i] = 45.0123;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, g0fac[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for tissparam as a 2D array of size ntissparamsxnsp
    // the following code replaces:tissparam = matrix(1, ntissparams, 1, nsp);
    // nsp is the number of reacting species (it is 7 for the Tumor1998 data)
    // ntissparams in the number of tissue parameters for each reacting species (it is 3 for the Tumor1998 data) - BS
    // In my code ntissparams will be the first iterator (rows), and nsp will be the second iterator (columns).
    tissparam = (float **) calloc(ntissparams+1, sizeof(float *));
    for(i=1;i<=ntissparams;i++)
    {
        tissparam[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=ntissparams;i++)
        for(j=1;j<=nsp;j++)
            tissparam[i][j] = -520.909;
    
    if(DEBG==1) {
        for(i=1;i<=ntissparams;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf \n", ntissparams, i, nsp, j, tissparam[i][j]);
    }
    //-----------------------------------------------------------------------------/
    // read data from simpleSoluteParams into arrays
    for (isp = 1; isp <= nsp; isp++) {
        fscanf(ifp, "%d %d %d", &permsolute[isp], &diffsolute[isp], &oxygen[isp]);
        if (diffsolute[isp] != 0 && diffsolute[isp] != 1) printf("*** Error: simpleSoluteParams, diffsolute[isp] must be 0 or 1. \n");// 1 if the solute is diffuses
        if (oxygen[isp] != 0 && oxygen[isp] != 1) printf("*** Error: simpleSoluteParams, oxygen[isp] must be 0 or 1. \n");// the solute is oxygen or is isn't
        if (oxygen[isp] == 1) permsolute[isp] = 1; //oxygen is permeable
        if (permsolute[isp] == 1) diffsolute[isp] = 1;  //permeable solutes must be diffusible
        fscanf(ifp, "%f", &pref[isp]); //typical maximum value
        fscanf(ifp, "%f", &diff[isp]); //D*alpha (for oxygen) otherwise, this is tissue diffusivity
        diff[isp] = diff[isp] * 1.e8; //convert diffusivity from cm^2/s to microns^2/s
        for (i = 1; i <= ntissparams; i++) {
            fscanf(ifp, "%f", &tissparam[i][isp]);
        }
        fscanf(ifp, "%f", &g0[isp]);// initial estimate used for tissue consumption
        fscanf(ifp, "%f", &g0fac[isp]);
    }
    
    if(DEBG==1) {
        for (isp = 1; isp <= nsp; isp++) {
            printf("%d %d %d %d %d \n", nsp, isp, permsolute[isp], diffsolute[isp], oxygen[isp]);
        }
        for (isp = 1; isp <= nsp; isp++) {
            printf("%d %d %f \n", nsp, isp, pref[isp]);
        }
        for (isp = 1; isp <= nsp; isp++) {
            printf("%d %d %f \n", nsp, isp, diff[isp]);
        }
        for (isp = 1; isp <= nsp; isp++) {
            printf("%d %d %f \n", nsp, isp, g0[isp]);
        }
        for (isp = 1; isp <= nsp; isp++) {
            printf("%d %d %f \n", nsp, isp, g0fac[isp]);
        }
    }
    
    if(DEBG==1) {
        for (isp = 1; isp <= nsp; isp++)
            for(i = 1; i<= ntissparams;i++)
        printf("%d %d %f \n", isp, i, tissparam[i][isp]);
    }
    fclose(ifp);
    
    //parameters for blood.  TWS January 2012
    plow = 0.1*p50;
    phigh = 5.*p50;
    clowfac = cs * (1.0 - 1.0 / (1.0 + pow((plow / p50), fn)));// Hill equation
    chighfac = cs * (1.0 - 1.0 / (1.0 + pow((phigh / p50), fn)));// Hill equation
    pphighfac = cs * fn / p50 * pow(phigh / p50, (fn - 1)) / SQR(1. + pow(phigh / p50, fn));
    if(DEBG==1){
        printf("plow = %5.40f, phigh = %5.40f, clowfac = %5.40f, chighfac = %5.40f, pphighfac = %5.40f\n", plow, phigh, clowfac, chighfac, pphighfac);
    }
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    //The section below is used to do multiple runs with variations in some parameters.
    nvaryparams = 0;    //default if no VaryParams.dat is present
    nruns = 1;
    //-----------------------------------------------------------------------------/
    // declare nsp slots in memory for solutefac.
    solutefac = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) solutefac[i] = 12.01;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, solutefac[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare nsp slots in memory for intravascfac.
    intravascfac = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) intravascfac[i] = 13.02;
    
    if(DEBG==1){
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf \n", nsp, i, intravascfac[i]);
    }
    
    for (isp = 1; isp <= nsp; isp++) {
        solutefac[isp] = 1.;
        intravascfac[isp] = 1.; // This is a factor that scales the resistance of the vessel walls wrt solute #isp, here it is initialized to 1. If there is no myVarayParams file, it will remain as 1.
    }
    
    if(DEBG==1){
        for (isp = 1; isp<= nsp; isp++){
            printf("%d %d %lf %lf \n", nsp, isp, solutefac[isp], intravascfac[isp]);
           
        }
    }

    if (ifp = fopen("myVaryParams", "r")){
        fscanf(ifp, "%i", &nvaryparams);
        if (nvaryparams) {
            //ivaryparams = imatrix(1, nvaryparams, 1, 3);
            ivaryparams = (int **) calloc(nvaryparams+1, sizeof(int *));
            for(i=1;i<=nvaryparams;i++) {
                ivaryparams[i] = (int *) calloc(3+1, sizeof(int ));
            }
            for(i=1;i<=nvaryparams;i++)
                for(j=1;j<=3;j++)
                    ivaryparams[i][j] = 45;
            
            if(DEBG==1) {
                for(i=1;i<=nvaryparams;i++)
                    for(j=1;j<=3;j++)
                        printf("%d %d %d %d %i\n", nvaryparams, i, 3, j, ivaryparams[i][j]);
            }
            //Up to 3 parameters can be varied. Allowable parameters to vary:
            //1: q0fac, 2: solutefac[isp], 3: diff[isp], 4: intravascfac[isp], 5: tissparam[i][isp]
            //values stored in ivaryparams[i][1], with i and or isp stored in ivaryparams[i][2 or 3]
            //-----------------------------------------------------------------------------/
            // Set up paramtyp(1, nvaryparams) as 1D array of integers
            paramtyp = (int *) calloc(nvaryparams+1, sizeof(int ));
            for(i=1;i<=nvaryparams;i++) paramtyp[i] = 12;
            
            if(DEBG==1) {
                for(i=1;i<=nvaryparams;i++)
                    printf("%d %d %i\n", nvaryparams, i, paramtyp[i]);
            }
            //-----------------------------------------------------------------------------/
            // Set up paramind(1,nvaryparams,1,3) as a 2D array of integers
            paramind = (int **) calloc(nvaryparams+1, sizeof(int *));
            for(i=1;i<=nvaryparams;i++) {
                paramind[i] = (int *) calloc(3+1, sizeof(int ));
            }
            for(i=1;i<=nvaryparams;i++)
                for(j=1;j<=3;j++)
                    paramind[i][j] = 13;
            
            if(DEBG==1) {
                for(i=1;i<=nvaryparams;i++)
                    for(j=1;j<=3;j++)
                        printf("%d %d %d %d %i\n", nvaryparams, i , 3, j, paramind[i][j]);
            }
            //-----------------------------------------------------------------------------/
            // This section was changed to avoid string matching - Baran
            for (i = 1; i <= nvaryparams; i++) {
                for (j = 1; j <= 3; j++) ivaryparams[i][j] = 0;
                fscanf(ifp, "%i %i %i", &paramtyp[i], &paramind[i][1], &paramind[i][2]);
                if (paramtyp[i] == 1) {
                    ivaryparams[i][1] = 1;
                    printf("Variable parameter %i: q0fac\n", i);
                }
                else if (paramtyp[i] == 2) {
                    ivaryparams[i][1] = 2;
                    ivaryparams[i][2] = paramind[i][1];
                    printf("Variable parameter %i: solutefac[%i]\n", i, ivaryparams[i][2]);
                }
                else if (paramtyp[i] == 3) {
                    ivaryparams[i][1] = 3;
                    ivaryparams[i][2] = paramind[i][1];
                    printf("Variable parameter %i: diff[%i]\n", i, ivaryparams[i][2]);
                }
                else if (paramtyp[i] == 4) {
                    ivaryparams[i][1] = 4;
                    ivaryparams[i][2] = paramind[i][1];
                    printf("Variable parameter %i: intravascfac[%i]\n", i, ivaryparams[i][2]);
                }
                else if (paramtyp[i] == 5) {
                    ivaryparams[i][1] = 5;
                    ivaryparams[i][2] = paramind[i][1];
                    ivaryparams[i][3] = paramind[i][2];
                    printf("Variable parameter %i: tissparams[%i][%i]\n", i, ivaryparams[i][2], ivaryparams[i][3]);
                }
            }
            
            if(DEBG==1) {
                for(i=1;i<=nvaryparams;i++)
                    printf("%d %d %d %d\n", nvaryparams, i, paramind[i][1], paramind[i][2]);
            }
            
            if(DEBG==1) {
                for(i=1;i<=nvaryparams;i++)
                    printf("%d %d %i %i %i\n", nvaryparams, i ,paramtyp[i], ivaryparams[i][2], ivaryparams[i][3]);
                
            }
            
            fscanf(ifp, "%i", &nruns);
            
            if(DEBG==1) {
                printf("%i\n",nruns);
            }
            //-----------------------------------------------------------------------------/
            //paramvalue = matrix(1, nruns, 1, nvaryparams);
            paramvalue = (float **) calloc(nruns+1, sizeof(float *));
            for(i=1;i<=nruns;i++) {
                paramvalue[i] = (float *) calloc(nvaryparams+1, sizeof(float ));
            }
            for(i=1;i<=nruns;i++)
                for(j=1;j<=nvaryparams;j++)
                    paramvalue[i][j] = 09.09;
            
            if(DEBG==1) {
                for(i=1;i<=nruns;i++)
                    for(j=1;j<=nvaryparams;j++)
                        printf("%d %d %d %d %f\n", nruns, i , nvaryparams, j, paramvalue[i][j]);
            }
            
            for (i = 1; i <= nruns; i++) {
                fscanf(ifp, "%i", &idum);
                for (j = 1; j <= nvaryparams; j++) fscanf(ifp, "%lf", &paramvalue[i][j]);
                fscanf(ifp, "%*[^\n]");
            }
        }
        fclose(ifp);
           
    }
    
    if(DEBG==1){
        printf("%d \n", nvaryparams);

    }
    
    //-----------------------------------------------------------------------------/
    // intravascular oxygen resistance data.  Assume zero unless specified in data file.
    // the code below replaces: resisdiam = matrix(1, 20, 1, nsp);
    // declare memory for resisdiam as a 2D array of size 20xnsp
    resisdiam = (float **) calloc(20+1, sizeof(float *));
    for(i=1;i<=20;i++)
    {
        resisdiam[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    
    for(i=1;i<=20;i++)
        for(j=1;j<=nsp;j++)
            resisdiam[i][j] = 10.09;
    
        if(DEBG==1){
                for(i=1;i<=20;i++)
                    for(j=1;j<=nsp;j++)
                        printf("%d %d %d %d %lf \n", 20, i, nsp, j, resisdiam[i][j]);
            }
    // the code below replaces: resis = matrix(1, 20, 1, nsp);
    // declare memory for resis as a 2D array of size 20xnsp
    resis = (float **) calloc(20+1, sizeof(float *));
    for(i=1;i<=20;i++)
    {
        resis[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    
    for(i=1;i<=20;i++)
        for(j=1;j<=nsp;j++)
            resis[i][j] = -1.07;
    
    if(DEBG==1){
        for(i=1;i<=20;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf \n", 20, i, nsp, j, resis[i][j]);
    }
    
    // declare memory for nresis as a 1D array of size nsp
    // the code below replaces: nresis = ivector(1, nsp);
    nresis = (int *) calloc(nsp+1,sizeof(int ));
    for(isp=1;isp<=nsp;isp++) nresis[i] = -4;
        
    if(DEBG==1){
            for(isp=1;isp<=nsp;isp++)
                printf("%d %d %d \n", nsp, i, nresis[i]);
        }
    
    ifp = fopen("simpleIntravascRes", "r");// This file was simplified from the Secomb group file following the same rules as were followed for simpleSoluteParams. - Baran
    if(ifp!=NULL) printf ("simpleIntravascRes opened successfully. \n");
    else { printf("could not open simpleIntravascRes file. exiting. \n"); exit(-1); }
    
    for (isp = 1; isp <= nsp; isp++) {
        fscanf(ifp, "%i", &nresis[isp]);
        if (nresis[isp] > 20) printf("Error: too many points in IntravascRes, nresis = %i > 20 \n", nresis[isp]);
        if (nresis[isp] > 0) {
            for (i = 1; i <= nresis[isp]; i++) fscanf(ifp, "%f %f", &resisdiam[i][isp], &resis[i][isp]);
        }
    }
    fclose(ifp);
    
    if(DEBG==1) {
        for (isp=1;isp<=nsp;isp++) printf("%d %d %d \n", nsp, isp, nresis[isp]);
    }
    
    if(DEBG==1) {
        for (isp=1;isp<=nsp;isp++)
            for(i=1;i<=nresis[isp];i++)
                printf("%d %d %d %d %f %f \n", nsp, isp, nresis[isp], i, resisdiam[i][isp], resis[i][isp]);
            }
    //-----------------------------------------------------------------------------/
    //network data file
    ifp = fopen("simpleNetwork", "r");// This file was simplified from the Secomb group file following the same rules as were followed for SoluteParams. See the comment beging at line 26. - Baran
    if (ifp!=NULL) printf("simpleNetwork opened successfully. \n");
    else { printf("could not open simpleNetwork file. exiting. \n"); exit(-1); }
    //dimensions of box in microns; vertex must be at origin
    fscanf(ifp, "%f %f %f", &alx, &aly, &alz);
    fscanf(ifp, "%i %i %i", &mxx, &myy, &mzz);
    fscanf(ifp, "%f", &lb);//Every tissue point is within lb of some vessel in the network.
    fscanf(ifp, "%f", &maxl);
    fscanf(ifp, "%i", &nodsegm);
    //number of segments in vessel network
    fscanf(ifp, "%i", &nseg);
    
    if(DEBG==1) {
        printf("%f %f %f \n", alx, aly, alz);
        printf("%d %d %d \n", mxx, myy, mzz);
        printf("%f \n", lb);
        printf("%f \n", maxl);
        printf("%d \n", nodsegm);
        printf("%d\n", nseg);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for segname as a 1D array of size nseg. For each segment it holds the number that represents that segment. - Baran
    // the code below replaces: segname = ivector(1, nseg);
    segname = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) segname[i] = 39;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %d\n", nseg, i, segname[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for segtyp as a 1D array of size nseg. For each segment if holds the type of that segment. The type affects how things will be calculated for that segment down the line. As and example, it could be used to have different permeabiliy to solutes in different segments types. - Baran
    // the code below replaces: segtyp = ivector(1, nseg);
    segtyp = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) segtyp[i] = 38;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %d\n", nseg, i, segtyp[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for segnodname as a 2D array of size 2xnseg. For each segment it holds the numbers that indicate the pair of nodes that is spanned by the segment. - BS
    // the code below replaces: segnodname = imatrix(1, 2, 1, nseg);
    // I had to change the declaration of segnodname from int to double in "myVariables.c" - BS
    segnodname = (float **) calloc(2+1, sizeof(float *));
    for(i=1;i<=2;i++)
    {
        segnodname[i] = (float *) calloc(nseg+1, sizeof(float ));
    }
    for(i=1;i<=2;i++)
        for(j=1;j<=nseg;j++)
            segnodname[i][j] = 6;
    
    if(DEBG==1) {
        for(i=1;i<=2;i++)
            for(j=1;j<=nseg;j++)
                printf("%d %d %d %d %lf \n", 2, i, nseg, j, segnodname[i][j]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for diam as a 1D array of size nseg
    // the code below replaces: diam = vector(1, nseg);
    diam = (float *) calloc(nseg+1, sizeof(float ));// The +1 is needed so that there will be nseg+1 memory locations allocated; I need nseg+1 because the first one is left unused.
    for (i=1;i<=nseg;i++) diam[i] = 13.033;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
                printf("%d %d %lf \n", nseg, i, diam[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for qdata as a 1D array of size nseg
    // the code below replaces: qdata = vector(1, nseg);
    qdata = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) qdata[i] = 12.02;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf \n", nseg,i, qdata[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for hd as a 1D array of size nseg
    // the code below replaces: hd = vector(1, nseg);
    hd = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) hd[i] = 24.8;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf \n", nseg, i, hd[i]);
    }
    //-----------------------------------------------------------------------------/
    // read in segment connectivity data etc...
    // qdata is relative flow??
    for (iseg = 1; iseg <= nseg; iseg++)
        fscanf(ifp, "%i %i %f %f %f %f %f",&segname[iseg], &segtyp[iseg], &segnodname[1][iseg], &segnodname[2][iseg], &diam[iseg], &qdata[iseg], &hd[iseg]);
    
    if(DEBG==1) {
        for(iseg=1;iseg<=nseg;iseg++)
            printf("%d %d %i %i %f %f %f %f %f \n", nseg, iseg, segname[iseg], segtyp[iseg], segnodname[1][iseg], segnodname[2][iseg], diam[iseg], qdata[iseg], hd[iseg]);
    }
    //-----------------------------------------------------------------------------/
    // read in the total number of nodes
    fscanf(ifp, "%i", &nnod);
    
    if(DEBG==1) {
        printf("%i \n", nnod);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for nodname as a 1D array of integers of size nnod
    //  the code below replaces: nodname = ivector(1, nnod);
    nodname = (int *) calloc(nnod+1, sizeof(int ));
    for(i=1;i<=nnod;i++) nodname[i] = 12;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %d\n", nnod, i, nodname[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for cnode and a 2D array of size 3xnnod
    // the code below replaces: cnode = matrix(1, 3, 1, nnod);
    cnode  = (float **) calloc(3+1, sizeof(float *));
    for(i=1;i<=3;i++)
    {
        cnode[i] = (float *) calloc(nnod+1, sizeof(float ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nnod;j++)
            cnode[i][j] = -13.0900;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nnod;j++)
                printf("%d %d %d %d %lf \n", 3, i, nnod, j, cnode[i][j]);
    }
    //-----------------------------------------------------------------------------/
    // read in node location data
    for (i = 1; i <= nnod; i++) fscanf(ifp, "%i %f %f %f%*[^\n]", &nodname[i], &cnode[1][i], &cnode[2][i], &cnode[3][i]);
    
    if(DEBG==1
       ) {
        for (i = 1; i <= nnod; i++) printf("%d %d %i %f %f %f \n", nnod, i, nodname[i], cnode[1][i], cnode[2][i], cnode[3][i]);
    }
    //-----------------------------------------------------------------------------/
    // number of boundary nodes
    fscanf(ifp, "%i ", &nnodbc);
    
    if(DEBG==1) {
        printf("%i \n", nnodbc);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for bcnodname as a 1D array of integers of size nnodbc
    //  the code below replaces: bcnodname = ivector(1, nnodbc);
    bcnodname = (int *) calloc(nnodbc+1, sizeof(int ));
    for(i=1;i<=nnodbc;i++) bcnodname[i] = 19;
    
    if(DEBG==1) {
        for(i=1;i<=nnodbc;i++)
            printf("%d %d %d\n", nnodbc, i, bcnodname[i]);
    }
    
    //-----------------------------------------------------------------------------/
    // declare memory for bcnod as a 1D array of integers of size nnodbc
    //  the code below replaces: bcnod = ivector(1, nnodbc);
    bcnod = (int *) calloc(nnodbc+1, sizeof(int ));
    for(i=1;i<=nnodbc;i++) bcnod[i] = 11;
    
    if(DEBG==1) {
        for(i=1;i<=nnodbc;i++)
            printf("%d %d %d\n", nnodbc, i, bcnod[i]);
    }
    // this one looks like it is unused in the Secomb group code
    //-----------------------------------------------------------------------------/
    // declare memory for bctyp as a 1D array of integers of size nnodbc
    //  the code below replaces: bctyp = ivector(1, nnodbc);
    bctyp = (int *) calloc(nnodbc+1, sizeof(int ));
    for(i=1;i<=nnodbc;i++) bctyp[i] = 18;
    
    if(DEBG==1) {
        for(i=1;i<=nnodbc;i++)
            printf("%d %d %d\n", nnodbc, i, bctyp[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for bcprfl as a 1D array of size nnodbc
    // the code below replaces: bcprfl = vector(1, nnodbc);
    bcprfl = (float *) calloc(nnodbc+1, sizeof(float ));
    for (i=1;i<=nnodbc;i++) bcprfl[i] = -12.008;
    
    if(DEBG==1) {
        for (i=1;i<=nnodbc;i++) printf("%d %d %lf \n", nnodbc,i, bcprfl[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for bchd as a 1D array of size nnodbc
    // the code below replaces: bchd = vector(1, nnodbc);
    bchd = (float *) calloc(nnodbc+1, sizeof(float ));
    for (i=1;i<=nnodbc;i++) bchd[i] = -14.008;
    
    if(DEBG==1) {
        for (i=1;i<=nnodbc;i++) printf("%d %d %lf \n", nnodbc, i, bchd[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for bcp as a 2D array of size nnodbcxnsp
    // the code below replaces: bcp = matrix(1, nnodbc, 1, nsp);
    bcp = (float **) calloc(nnodbc+1, sizeof(float *));
    for(i=1;i<=nnodbc;i++) {
        bcp[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnodbc;i++)
        for(j=1;j<=nsp;j++)
            bcp[i][j] = 90.81;
    
    if(DEBG==1) {
        for(i=1;i<=nnodbc;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf \n",nnodbc, nsp, i, j, bcp[i][j]);
    }
    //-----------------------------------------------------------------------------/
    // read in boundry node data
    for (i = 1; i <= nnodbc; i++) {
        fscanf(ifp, "%d %d %f %f", &bcnodname[i], &bctyp[i], &bcprfl[i], &bchd[i]);
        for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) fscanf(ifp, "%f", &bcp[i][isp]);
        fscanf(ifp, "%*[^\n]");    //ignore any 'extra' solutes in data file
    }
    
    if(DEBG==1) {
        for (i = 1; i <= nnodbc; i++) {
            printf("%d %d %d %d %lf %lf \n", nnodbc, i, bcnodname[i], bctyp[i], bcprfl[i], bchd[i]);
            for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) printf("%d %d %d %d %lf \n", nnodbc, i, nsp, isp, bcp[i][isp]);
        }
    }
    fclose(ifp);
    //-----------------------------------------------------------------------------/
    //v = total box volume, vol = volume represented by each tissue point; req = radius of equivalent sphere
    v = alx*aly*alz;
    vol = v / (mxx*myy*mzz);// total volume divided by number of tissue points,
    printf("vol = %5.5f\n", vol);
    if (mzz == 1) req = pow(vol*1. / alz / pi1, 0.5);//2d version (if the number of tissue points in the z direction is 1; in the planar case) do this
    else req = pow(vol*0.75 / pi1, 0.333333);
    
    if(DEBG==1) {
        printf("%f\n", req);
    }
    //-----------------------------------------------------------------------------/
    //Read parameters for slice on which P is computed for contour plot
    // declare memory for nl as a 1D array of integers of size nsp
    // the code below replaces: nl = ivector(1, nsp);
    nl = (int *) calloc(nsp+1, sizeof(int ));
    for(i=1;i<=nsp;i++) nl[i] = 15;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %d\n", nsp, i, nl[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for xsl0 as a 1D array of size 3
    // the code below replaces: xsl0 = vector(1, 3);
    xsl0 = (float *) calloc(3+1, sizeof(float ));
    for (i=1;i<=3;i++) xsl0[i] = 90.07;
    
    if(DEBG==1) {
        for (i=1;i<=3;i++) printf("%d %d %lf \n", 3, i, xsl0[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for xsl1 as a 1D array of size 3
    // the code below replaces: xsl1 = vector(1, 3);
    xsl1 = (float *) calloc(3+1, sizeof(float ));
    for (i=1;i<=3;i++) xsl1[i] = 809.1202;
    
    if(DEBG==1) {
        for (i=1;i<=3;i++) printf("%d %d %lf \n", 3, i, xsl1[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for xsl2 as a 1D array of size 3
    // the code below replaces: xsl2 = vector(1, 3);
    xsl2 = (float *) calloc(3+1, sizeof(float ));
    for (i=1;i<=3;i++) xsl2[i] = 7404.023;
    
    if(DEBG==1) {
        for (i=1;i<=3;i++) printf("%d %d %lf \n", 3, i, xsl2[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for clmin as a 1D array of size ns
    // the code below replaces: clmin = vector(1, nsp);
    clmin = (float *) calloc(nsp+1, sizeof(float ));
    for (i=1;i<=nsp;i++) clmin[i] = -607.01;
    
    if(DEBG==1) {
        for (i=1;i<=nsp;i++) printf("%d %d %lf \n", nsp, i, clmin[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for clint as a 1D array of size nsp
    // the code below replaces: clint = vector(1, nsp);
    clint = (float *) calloc(nsp+1, sizeof(float ));
    for (i=1;i<=nsp;i++) clint[i] = 90.01;
    
    if(DEBG==1) {
        for (i=1;i<=nsp;i++) printf("%d %d %lf \n", nsp, i, clint[i]);
    }
    //-----------------------------------------------------------------------------/
    // read in contour parameters
    ifp = fopen("simpleContourParams", "r");/// This file was simplified from the Secomb group file following the same rules as were followed for simpleSoluteParams. - Baran
    if (ifp!=NULL) printf("simpleContourParams opened successfully. \n");
    else { printf("could not open simpleContourParams file. exiting. \n"); exit(-1); }
    
    fscanf(ifp, "%f %f %f %i", &xsl0[1], &xsl0[2], &xsl0[3], &slsegdiv);
    fscanf(ifp, "%f %f %f %i", &xsl1[1], &xsl1[2], &xsl1[3], &nsl1);
    fscanf(ifp, "%f %f %f %i", &xsl2[1], &xsl2[2], &xsl2[3], &nsl2);
    nlmax = 1;
    
    if(DEBG==1) {
        printf("%f %f %f %i \n", xsl0[1], xsl0[2], xsl0[3], slsegdiv);
        printf("%f %f %f %i \n", xsl1[1], xsl1[2], xsl1[3], nsl1);
        printf("%f %f %f %i \n", xsl2[1], xsl2[2], xsl2[3], nsl2);
    }
    
    for (isp = 1; isp <= nsp; isp++) {
        fscanf(ifp, "%f %f %d", &clmin[isp], &clint[isp], &nl[isp]);
        if (nl[isp] > nlmax) nlmax = nl[isp];
    }
    
    if(DEBG==1) {
        for(isp=1;isp<=nsp;isp++)
        printf("%d %d %f %f %d \n", nsp, isp, clmin[isp], clint[isp], nl[isp]);
    }
    fclose(ifp);
    //-----------------------------------------------------------------------------/
    xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));// farthest point from the origin of the planar slice in the x direction
    ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));// farthest point from the origin of the planar slice in the y direction
    scalefac = FMIN(500. / xmax, 700. / ymax);//updated April 2010
    //-----------------------------------------------------------------------------/
    // declare memory for cl as an 1D array of size nlmax
    // the code below replaces cl = vector(1, nlmax);
    cl = (float *) calloc(nlmax+1, sizeof(float ));
    for (i=1;i<=nlmax;i++) cl[i] = -12.09;
    
    if(DEBG==1) {
        for (i=1;i<=nlmax;i++) printf("%d %d %lf \n", nlmax, i, cl[i]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for zv as a 2D array of size nsl1xnsl2
    // the code below replaces: zv = matrix(1, nsl1, 1, nsl2);
    zv = (float **) calloc(nsl1+1, sizeof(float *));
    for (i=1;i<=nsl1;i++) {
        zv[i] = (float *) calloc(nsl2+1, sizeof(float ));
        }
    for(i=1;i<=nsl1;i++)
        for(j=1;j<=nsl2;j++)
            zv[i][j] = 19.08;
    
    if(DEBG==1) {
        for(i=1;i<=nsl1;i++)
            for(j=1;j<=nsl2;j++)
                printf("%d %d %d %d %lf \n", nsl1, i, nsl2, j, zv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    // declare memory for psl as a 3D array of size nsl1xnsl2xnsp
    // the code below replaces: psl = f3tensor(1, nsl1, 1, nsl2, 1, nsp);
    psl = (float ***) calloc(nsl1+1, sizeof(float **));
    for (i=1;i<=nsl1;i++) {
        psl[i] = (float **) calloc(nsl2+1, sizeof(float *));
    }
    for(i=1;i<=nsl1;i++) {
        for(j=1;j<=nsl2;j++) {
            psl[i][j] = (float  *) calloc(nsp+1, sizeof(float ));
         }
    }
    for(i=1;i<=nsl1;i++)
        for(j=1;j<=nsl2;j++)
            for(k=1;k<=nsp;k++)
                psl[i][j][k] = -1040.24020;
    
    if(DEBG==1) {
        for(i=1;i<=nsl1;i++)
            for(j=1;j<=nsl2;j++)
                for(k=1;k<=nsp;k++)
                    printf("%d %d %d %d %d %d %lf \n", nsl1, i, nsl2, j, nsp, k, psl[i][j][k]);
    }
    //-----------------------------------------------------------------------------/
    //Read parameters to run postgreens
    npostgreensparams = 0;
    if (ifp = fopen("myPostGreensParams", "r")) {
        fgets(bb, max, ifp);
        fscanf(ifp, "%i", &npostgreensout);
        fscanf(ifp, "%i", &npostgreensparams);
        if (npostgreensout)
        //-----------------------------------------------------------------------------/
        //postgreensout = vector(1, npostgreensout);
        postgreensout = (float *) calloc(npostgreensout+1, sizeof(float ));
        for(i=1;i<=npostgreensout;i++) postgreensout[i] = 12.09;
            
        if(DEBG==1) {
            for(i=1;i<=npostgreensout;i++) printf("%d %d %lf\n", npostgreensout, i, postgreensout[i]);
            }
        //-----------------------------------------------------------------------------/
        if (npostgreensparams) {
        //-----------------------------------------------------------------------------/
        //postgreensparams = vector(1, npostgreensparams);
        postgreensparams = (float *) calloc(npostgreensparams+1, sizeof(float ));
        for(i=1;i<=npostgreensparams;i++) postgreensparams[i] = 12.099;
            
        if(DEBG==1) {
            for(i=1;i<=npostgreensparams;i++) printf("%d %d %lf\n", npostgreensparams, i, postgreensparams[i]);
            }
        //-----------------------------------------------------------------------------/
            for (i = 1; i <= npostgreensparams; i++) fscanf(ifp, "%f%*[^\n]", &postgreensparams[i]);
        }
        if(DEBG==1){
            printf("Printed\n");
        }
    }
    //-----------------------------------------------------------------------------/
    // these two lines appear in the Secomb group code after input() and before setuparrays0()
    is2d = 0; //set to 1 for 2d version, 0 otherwise
    if (mzz == 1) is2d = 1; //assumes 2d version if all tissue points lie in one z-plane
    
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on setuparrays0()
    //errvesselcount = ivector(1, nsp);
    errvesselcount = (int *) calloc(nsp+1, sizeof(int ));
    for(i=1;i<=nsp;i++) errvesselcount[i] = 12;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %i\n", nsp, i, errvesselcount[i]);
    }
    //-----------------------------------------------------------------------------/
    //errtissuecount = ivector(1, nsp);
    errtissuecount = (int *) calloc(nsp+1, sizeof(int ));
    for(i=1;i<=nsp;i++) errtissuecount[i] = 13;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %i\n", nsp, i, errtissuecount[i]);
    }
    //-----------------------------------------------------------------------------/
    //imaxerrvessel = ivector(1, nsp);
    imaxerrvessel = (int *) calloc(nsp+1, sizeof(int ));
    for(i=1;i<=nsp;i++) imaxerrvessel[i] = 14;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %i\n", nsp, i, imaxerrvessel[i]);
    }
    //-----------------------------------------------------------------------------/
    //imaxerrtissue = ivector(1, nsp);
    imaxerrtissue = (int *) calloc(nsp+1, sizeof(int ));
    for(i=1;i<=nsp;i++) imaxerrtissue[i] = 15;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %i\n", nsp, i, imaxerrtissue[i]);
    }
    //-----------------------------------------------------------------------------/
    //nbou = i3tensor(1, mxx, 1, myy, 1, mzz);
    nbou = (int ***) calloc(mxx+1, sizeof(int **));
    for(i=1;i<=mxx;i++) {
        nbou[i] = (int **) calloc(myy+1, sizeof(int *));
    }
    for(i=1;i<=mxx;i++) {
        for(j=1;j<=myy;j++) {
            nbou[i][j] = (int *) calloc(mzz+1, sizeof(int ));
        }
    }
    for(i=1;i<=mxx;i++)
        for(j=1;j<=myy;j++)
            for(k=1;k<=mzz;k++)
                nbou[i][j][k] = 16;
    
    if(DEBG==1) {
        for(i=1;i<=mxx;i++)
            for(j=1;j<=myy;j++)
                for(k=1;k<=mzz;k++)
                    printf("%d %d %d %d %d %d %i\n", mxx, i, myy, j, mzz, k, nbou[i][j][k]);
    }
    //-----------------------------------------------------------------------------/
    pmin = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmin[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %lf\n", nsp, i, pmin[i]);
    }
    //-----------------------------------------------------------------------------/
    pmax = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.02;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    pmean = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmean[i] = -18.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmean[i]);
    }
    //-----------------------------------------------------------------------------/
    mtiss = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    mptiss = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) mptiss[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, mptiss[i]);
    }
    //-----------------------------------------------------------------------------/
    g0old = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) g0old[i] = -21.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, g0old[i]);
    }
    //-----------------------------------------------------------------------------/
    ptt = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) ptt[i] = -22.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, ptt[i]);
    }
    //-----------------------------------------------------------------------------/
    ptpt = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) ptpt[i] = -23.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, ptpt[i]);
    }
    //-----------------------------------------------------------------------------/
    qtsum = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) qtsum[i] = -24.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, qtsum[i]);
    }
    //-----------------------------------------------------------------------------/
    qvsum = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) qvsum[i] = -25.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, qvsum[i]);
    }
    //-----------------------------------------------------------------------------/
    errvessel = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) errvessel[i] = -26.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, errvessel[i]);
    }
    //-----------------------------------------------------------------------------/
    errtissue = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) errtissue[i] = -27.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, errtissue[i]);
    }
    //-----------------------------------------------------------------------------/
    dqvsumdg0 = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    dqtsumdg0 = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.02;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    g0facnew = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    pinit = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.03;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    p = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.401;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    epsvessel = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    epstissue = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.11;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    eps = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -17.21;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    qvfac = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -15.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    mptissref = (float *) calloc(nsp+1, sizeof(float ));
    for(i=1;i<=nsp;i++) pmax[i] = -12.01;
    
    if(DEBG==1) {
        for(i=1;i<=nsp;i++)
            printf("%d %d %f\n", nsp, i, pmax[i]);
    }
    
    //-----------------------------------------------------------------------------/
    x = (float *) calloc(3+1, sizeof(float ));
    for(i=1;i<=3;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            printf("%d %d %f\n", 3, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    y = (float *) calloc(3+1, sizeof(float ));
    for(i=1;i<=3;i++) pmax[i] = -17.02;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            printf("%d %d %f\n", 3, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    ss = (float *) calloc(3+1, sizeof(float ));
    for(i=1;i<=3;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            printf("%d %d %f\n", 3, i, pmax[i]);
    }
    
    //-----------------------------------------------------------------------------/
    axt = (float *) calloc(mxx+1, sizeof(float ));
    for(i=1;i<=mxx;i++) pmax[i] = -37.01;
    
    if(DEBG==1) {
        for(i=1;i<=mxx;i++)
            printf("%d %d %f\n", mxx, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    ayt = (float *) calloc(myy+1, sizeof(float ));
    for(i=1;i<=myy;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=myy;i++)
            printf("%d %d %f\n", myy, i, pmax[i]);
    }
    //-----------------------------------------------------------------------------/
    azt = (float *) calloc(mzz+1, sizeof(float ));
    for(i=1;i<=mzz;i++) pmax[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=mzz;i++)
            printf("%d %d %f\n", mzz, i, pmax[i]);
    }
    
    //-----------------------------------------------------------------------------/
    //dtt = f3tensor(1, mxx, 1, myy, 1, mzz);
    dtt = (float ***) calloc(mxx+1, sizeof(float **));
    for (i=1;i<=mxx;i++) {
        dtt[i] = (float **) calloc(myy+1, sizeof(float *));
    }
    for(i=1;i<=mxx;i++) {
        for(j=1;j<=myy;j++) {
            dtt[i][j] = (float *) calloc(mzz+1, sizeof(float ));
        }
    }
    for(i=1;i<=mxx;i++)
        for(j=1;j<=myy;j++)
            for(k=1;k<=mzz;k++)
                dtt[i][j][k] = -3304405.09;
    
    if(DEBG==1) {
        for(i=1;i<=mxx;i++)
            for(j=1;j<=myy;j++)
                for(k=1;k<=mzz;k++)
                    printf("%d %d %d %d %d %d %lf\n", mxx, i, myy, j, mzz, k, dtt[i][j][k]);
    }
    
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on setuparrays1(nseg, nnod)
    //nspoint = ivector(1, nseg);
    nspoint = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) nspoint[i] = 1;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %i\n", nseg, i, nspoint[i]);
    }
    //-----------------------------------------------------------------------------/
    //istart = ivector(1, nseg);
    istart = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) istart[i] = 2;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %i\n", nseg, i, istart[i]);
    }
    //-----------------------------------------------------------------------------/
    //lowflow = ivector(1, nseg);
    lowflow = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) lowflow[i] = 3;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %i\n", nseg, i, lowflow[i]);
    }
    //-----------------------------------------------------------------------------/
    //nk = ivector(1, nnod); //not nseg - error fixed 20 April 2010
    nk = (int *) calloc(nnod+1, sizeof(int ));
    for(i=1;i<=nnod;i++) nk[i] = 4;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %i\n", nnod, i, nk[i]);
    }
    //-----------------------------------------------------------------------------/
    //ista = ivector(1, nseg);
    ista = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) ista[i] = 5;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %i\n", nseg, i, ista[i]);
    }
    //-----------------------------------------------------------------------------/
    //iend = ivector(1, nseg);
    iend = (int *) calloc(nseg+1, sizeof(int ));
    for(i=1;i<=nseg;i++) iend[i] = 6;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %i\n", nseg, i, iend[i]);
    }
    //-----------------------------------------------------------------------------/
    //nodrank = ivector(1, nnod);
    nodrank = (int *) calloc(nnod+1, sizeof(int ));
    for(i=1;i<=nnod;i++) nodrank[i] = 0;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %i\n", nnod, i, nodrank[i]);
    }
    //-----------------------------------------------------------------------------/
    //nodout = ivector(1, nnod);
    nodout = (int *) calloc(nnod+1, sizeof(int ));
    for(i=1;i<=nnod;i++) nodout[i] = 8;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %i\n", nnod, i, nodout[i]);
    }
    //-----------------------------------------------------------------------------/
    //nodtyp = ivector(1, nnod);
    nodtyp = (int *) calloc(nnod+1, sizeof(int ));
    for(i=1;i<=nnod;i++) nodtyp[i] = 9;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %i\n", nnod, i, nodtyp[i]);
    }
    //-----------------------------------------------------------------------------/
    //nodnod = imatrix(1, nodsegm, 1, nnod);
    nodnod = (int **) calloc(nodsegm+1, sizeof(int *));
    for(i=1;i<=nodsegm;i++) {
        nodnod[i] = (int *) calloc(nnod+1, sizeof(int ));
    }
    for(i=1;i<=nodsegm;i++)
        for(j=1;j<=nnod;j++)
            nodnod[i][j] = 10;
    
    if(DEBG==1) {
        for(i=1;i<=nodsegm;i++)
            for(j=1;j<=nnod;j++)
                printf("%d %d %d %d %i\n", nodsegm, i, nnod, j, nodnod[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //nodseg = imatrix(1, nodsegm, 1, nnod);
    nodseg = (int **) calloc(nodsegm+1, sizeof(int *));
    for(i=1;i<=nodsegm;i++) {
        nodseg[i] = (int *) calloc(nnod+1, sizeof(int ));
    }
    for(i=1;i<=nodsegm;i++)
        for(j=1;j<=nnod;j++)
            nodseg[i][j] = 11;
    
    if(DEBG==1) {
        for(i=1;i<=nodsegm;i++)
            for(j=1;j<=nnod;j++)
                printf("%d %d %d %d %i\n", nodsegm, i, nnod, j, nodseg[i][j]);
    }
    
    //-----------------------------------------------------------------------------/
    //segvar = vector(1, nseg);
    segvar = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) segvar[i] = 901.000;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, segvar[i]);
    }
    //-----------------------------------------------------------------------------/
    //q = vector(1, nseg);//added Novamber 2016
    q = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) q[i] = 901.000;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, q[i]);
    }
    //-----------------------------------------------------------------------------/
    //qq = vector(1, nseg);
    qq = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) qq[i] = 901.000;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, qq[i]);
    }
    //-----------------------------------------------------------------------------/
    //oxflux = vector(1, nodsegm);
    oxflux = (float *) calloc(nodsegm+1, sizeof(float ));
    for(i=1;i<=nodsegm;i++) oxflux[i] = 901.000;
    
    if(DEBG==1) {
        for(i=1;i<=nodsegm;i++)
            printf("%d %d %lf\n", nodsegm, i, oxflux[i]);
    }
    //-----------------------------------------------------------------------------/
    //rseg = vector(1, nseg);
    rseg = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) rseg[i] = 901.900;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, rseg[i]);
    }
    //-----------------------------------------------------------------------------/
    //cbar = vector(1, nseg);
    cbar = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) cbar[i] = 901.900;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, cbar[i]);
    }
    //-----------------------------------------------------------------------------/
    //lseg = vector(1, nseg);
    lseg = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) lseg[i] = 901.900;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, lseg[i]);
    }
    //-----------------------------------------------------------------------------/
    //ds = vector(1, nseg);
    ds = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) ds[i] = 901.9000;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, ds[i]);
    }
    //-----------------------------------------------------------------------------/
    //nodvar = vector(1, nnod);
    nodvar = (float *) calloc(nnod+1, sizeof(float ));
    for(i=1;i<=nnod;i++) nodvar[i] = 901.90000;
    
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            printf("%d %d %lf\n", nnod, i, nodvar[i]);
    }
    //-----------------------------------------------------------------------------/
    //segc = vector(1, nseg);
    segc = (float *) calloc(nseg+1, sizeof(float ));
    for(i=1;i<=nseg;i++) segc[i] = 901.90000;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            printf("%d %d %lf\n", nseg, i, segc[i]);
    }
    //-----------------------------------------------------------------------------/
    //gamma1 = matrix(1, nseg, 1, nsp);
    gamma1 = (float **) calloc(nseg+1, sizeof(float *));
    for(i=1;i<=nseg;i++) {
        gamma1[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nseg;i++)
        for(j=1;j<=nsp;j++)
            gamma1[i][j] = -201.01;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nseg, i, nsp, j, gamma1[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //qvseg = matrix(1, nseg, 1, nsp);
    qvseg = (float **) calloc(nseg+1, sizeof(float *));
    for(i=1;i<=nseg;i++) {
        qvseg[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nseg;i++)
        for(j=1;j<=nsp;j++)
            qvseg[i][j] = -201.02;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nseg, i, nsp, j, qvseg[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pvseg = matrix(1, nseg, 1, nsp);
    pvseg = (float **) calloc(nseg+1, sizeof(float *));
    for(i=1;i<=nseg;i++) {
        pvseg[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nseg;i++)
        for(j=1;j<=nsp;j++)
            pvseg[i][j] = -201.03;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nseg, i, nsp, j, pvseg[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pevseg = matrix(1, nseg, 1, nsp);
    pevseg = (float **) calloc(nseg+1, sizeof(float *));
    for(i=1;i<=nseg;i++) {
        pevseg[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nseg;i++)
        for(j=1;j<=nsp;j++)
            pevseg[i][j] = -201.04;
    
    if(DEBG==1) {
        for(i=1;i<=nseg;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nseg, i, nsp, j, pevseg[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //start = matrix(1, 3, 1, nseg);
    start = (float **) calloc(3+1, sizeof(float *));
    for(i=1;i<=3;i++) {
        start[i] = (float *) calloc(nseg+1, sizeof(float ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nseg;j++)
            start[i][j] = -201.05;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nseg;j++)
                printf("%d %d %d %d %lf\n", 3, i, nseg, j, start[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //scos = matrix(1, 3, 1, nseg);
    scos = (float **) calloc(3+1, sizeof(float *));
    for(i=1;i<=3;i++) {
        scos[i] = (float *) calloc(nseg+1, sizeof(float ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nseg;j++)
            scos[i][j] = -201.06;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nseg;j++)
                printf("%d %d %d %d %lf\n", 3, i, nseg, j, scos[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //end = matrix(1, 3, 1, nseg);
    end = (float **) calloc(3+1, sizeof(float *));
    for(i=1;i<=3;i++) {
        end[i] = (float *) calloc(nseg+1, sizeof(float ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nseg;j++)
            end[i][j] = -201.07;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nseg;j++)
                printf("%d %d %d %d %lf\n", 3, i, nseg, j, end[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //rsta = f3tensor(1, 3, 1, 16, 1, nseg);
    rsta = (float ***) calloc(3+1, sizeof(float **));
    for(i=1;i<=3;i++) {
        rsta[i] = (float **) calloc(16+1, sizeof (float *));
    }
    for(i=1;i<=3;i++) {
        for(j=1;j<=16;j++) {
            rsta[i][j] = (float *) calloc(nseg+1, sizeof(float ));
        }
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=16;j++)
            for(k=1;k<=nseg;k++)
                rsta[i][j][k] = 90.033339;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=16;j++)
                for(k=1;k<=nseg;k++)
                    printf("%d %d %d %d %d %d %lf\n", 3, i, 16 , j , nseg, k, rsta[i][j][k]);
    }
    //-----------------------------------------------------------------------------/
    //rend = f3tensor(1, 3, 1, 16, 1, nseg);
    rend = (float ***) calloc(3+1, sizeof(float **));
    for(i=1;i<=3;i++) {
        rend[i] = (float **) calloc(16+1, sizeof (float *));
    }
    for(i=1;i<=3;i++) {
        for(j=1;j<=16;j++) {
            rend[i][j] = (float *) calloc(nseg+1, sizeof(float ));
        }
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=16;j++)
            for(k=1;k<=nseg;k++)
                rend[i][j][k] = 90.033340;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=16;j++)
                for(k=1;k<=nseg;k++)
                    printf("%d %d %d %d %d %d %lf\n", 3, i, 16 , j , nseg, k, rend[i][j][k]);
    }
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on analyzenet();
    // Find node numbers corresponding to segment nodes - for all segments
    for(inod=1; inod<=nnod; inod++) nodtyp[inod] = 0;
    for(iseg=1; iseg<=nseg; iseg++)    {
        for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[1][iseg]){
            ista[iseg] = inod;
            goto foundit1;
        }
        printf("*** Error: No matching node found for nodname %d\n", segnodname[1][iseg]);
    foundit1:;
        for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[2][iseg]){
            iend[iseg] = inod;
            goto foundit2;
        }
        printf("*** Error: No matching node found for nodname %d\n", segnodname[2][iseg]);
    foundit2:;
    }
    
    if(DEBG==1) {
        for(iseg=1;iseg<=nseg;iseg++)
            printf("%d %d %d %d %lf %lf\n", nseg, iseg , segnodname[1][iseg], segnodname[2][iseg], ista[iseg], iend[iseg]);
    }
    
    if(DEBG==1) {
        for(iseg=1;iseg<=nseg;iseg++) if(segnodname[1][iseg] == ista[iseg] && segnodname[2][iseg])
            printf("node assignments for segment %i are OK!\n", iseg );
        else
            printf("*****Error");
    }
    
    //Setup nodtyp, nodseg and nodnod
    //goes segment by segment, looking at parent and child node of the segment
    //each time a node appears its node type (nodtyp) is increased by one
    //at the end, nodtyp is the number of segments that meet at that node
    //if the nodtyp for any node exceeds a user set maximum (nodsegm), error messages are printed
    //finally the arrays nodseg and nodnod are setup
    for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
        inod1 = ista[iseg];
        inod2 = iend[iseg];
        nodtyp[inod1] = nodtyp[inod1] + 1;
        nodtyp[inod2] = nodtyp[inod2] + 1;
        if(nodtyp[inod1] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod1);
        if(nodtyp[inod2] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod2);
        nodseg[nodtyp[inod1]][inod1] = iseg;
        nodseg[nodtyp[inod2]][inod2] = iseg;
        nodnod[nodtyp[inod1]][inod1] = inod2;
        nodnod[nodtyp[inod2]][inod2] = inod1;
    }
    
    if(DEBG==1) {
        for(iseg=1;iseg<=nseg;iseg++)
            printf("%d %d %d %d\n", nseg, iseg, nodtyp[ista[iseg]], nodtyp[iend[iseg]]);
    }
    if(DEBG==1) {
        for(i=1;i<=nnod;i++)
            for(j=1;j<=nodtyp[i];j++)
                printf("%d %d %d %d %d %d\n", nnod, i, nodtyp[i], j, nodseg[j][i], nodnod[j][i]);
    }
    
    //Find node numbers corresponding to boundary nodes
    for (inodbc=1; inodbc<=nnodbc; inodbc++){
        for(inod=1; inod<=nnod; inod++) if(nodname[inod] == bcnodname[inodbc]){
            bcnod[inodbc] = inod;
            if(nodtyp[inod] != 1) printf("*** Error: Boundary node %i is not a 1-segment node\n", inod);
            goto foundit;
        }
        printf("*** Error: No matching node found for nodname %i\n", bcnodname[inodbc]);
    foundit:;
    }
    
    if(DEBG==1) {
        for(inodbc=1;inodbc<=nnodbc;inodbc++) if(bcnod[inodbc] == bcnodname[inodbc])
             printf("node name of boundry node %i copied sucessfully!\n", inodbc);
    }
    //Calculate total inflow to network. TWS July 2018
    totalq = 0.;
    for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
        inod = bcnod[inodbc];    //boundary node
        iseg = nodseg[1][inod];
        // the following two cases are inflow, the complementary cases (q[iseg] > 0. && inod == iend[iseg]) and (q[iseg] < 0. && inod == ista[iseg]) are outflow
        if (qdata[iseg] > 0. && inod == ista[iseg]) totalq += qdata[iseg]*q0fac; //changed from q[iseg] - Baran
        if (qdata[iseg] < 0. && inod == iend[iseg]) totalq -= qdata[iseg]*q0fac; //changed from q[iseg] - Baran
       
        if(DEBG==1) {
            printf("%d %d %f %f %f\n", nnodbc, inodbc, totalq, q[iseg], qdata[iseg]);
        }
    }
    // calculate length of each segment
    // start[k][iseg] = coordinates of starting point of segment i
    for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
        rseg[iseg] = diam[iseg]/2.0;
        lseg[iseg] = 0.;
        for(k=1; k<=3; k++){
            start[k][iseg] = cnode[k][ista[iseg]];
            end[k][iseg] = cnode[k][iend[iseg]];
            ss[k] = end[k][iseg] - start[k][iseg];
            lseg[iseg] += SQR(ss[k]);//lseg[iseg] = lseg[iseg] + SQR(ss[k])
        }
        lseg[iseg] = sqrt(lseg[iseg]);
        if(lseg[iseg] == 0.) printf("*** Error: segment %i has zero length\n",segname[iseg]);//added May 2010, modified July 2017
        for(j=1; j<=3; j++)    scos[j][iseg] = ss[j]/lseg[iseg];//cos = adj/hyp
        //find points on vessel cylinders (used for tissue region bounds) ??
        for(k=0; k<=15; k++){
            t = k*2*pi1/16.0;
            sintheta = sqrt(1.0 - SQR(scos[3][iseg]));// sin^2 + cos^2 = 1
            if(sintheta > 0.0001){
                cosfi = -scos[2][iseg]/sintheta;
                sinfi = scos[1][iseg]/sintheta;
            }
            else {
                cosfi = 1.;
                sinfi = 0.;
            }
            // Start and end discretization of lumens circles being assigned here - Baran(g)
            rsta[1][k+1][iseg] = rseg[iseg]*cosfi*cos(t)-rseg[iseg]*scos[3][iseg]*sinfi*sin(t) + start[1][iseg];
            rsta[2][k+1][iseg] = rseg[iseg]*sinfi*cos(t)+rseg[iseg]*scos[3][iseg]*cosfi*sin(t) + start[2][iseg];
            rsta[3][k+1][iseg] = rseg[iseg]*sintheta*sin(t) + start[3][iseg];
            rend[1][k+1][iseg] = rseg[iseg]*cosfi*cos(t)-rseg[iseg]*scos[3][iseg]*sinfi*sin(t) + end[1][iseg];
            rend[2][k+1][iseg] = rseg[iseg]*sinfi*cos(t)+rseg[iseg]*scos[3][iseg]*cosfi*sin(t) + end[2][iseg];
            rend[3][k+1][iseg] = rseg[iseg]*sintheta*sin(t) + end[3][iseg];
        }
    }
    
    if(DEBG==1){
        for(iseg=1;iseg<=nseg;iseg++) printf("%d %d %d %f\n", nseg, iseg,3, lseg[iseg]);
    }
    if(DEBG==1){
        for(iseg=1;iseg<=nseg;iseg++) for(j=1;j<=3;j++) printf("%d %d %d %d %f\n", nseg, iseg,3,j, scos[j][iseg]);
    }
    if(DEBG==1){
        for(iseg=1;iseg<=nseg;iseg++) for(j=1;j<=3;j++) printf("%d %d %d %d %f\n", nseg, iseg,3,j, sqrt(1.0 - SQR(scos[j][iseg])));
    }
    //subdivide segments into small elements as needed
    //compute coordinates of source points
    //nnv = total number of elements
    nnv = 0;
    for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
        m = lseg[iseg]/maxl + 1; //1seg[iseg]/maxl gets coverted to the nearest integer since m is declared as an integer. The +1 is so that m will be 1 rather than 0 when lseg is small
        nspoint[iseg] = m; //number of source points per segment is equal to the number of subsegments
        istart[iseg] = nnv + 1;// the name of the parent node of each segment when numbered taking into account the new nodes that are included after dividing into subsegments
        ds[iseg] = lseg[iseg]/m;// length of subsegments
        nnv += m;
    }
    printf("Total vessel points = %i\n", nnv);
    //compute coordinates of tissue points
    delx = alx/mxx;
    for(i=1; i<=mxx; i++) axt[i] = (i-0.5f)*delx;
    dely = aly/myy;
    for(i=1; i<=myy; i++) ayt[i] = (i-0.5f)*dely;
    delz = alz/mzz;
    for(i=1; i<=mzz; i++) azt[i] = (i-0.5f)*delz;
    //create array of tissue points inside tissue domain boundaries using method 1 or 2
    // method can be 1, 2, 3 or 5
    method = 1; // - Baran
    //nnt = outboun(1);
    //-----------------------------------------------------------------------------/
    // outboun.cpp
    // ----------------------------------------------------
    // method = 1: finds the smallest convex region inside the cuboid
    // which excludes tissue node points that have a distance to the nearest
    // vessel greater than a value specified by the user (lb).  Any point
    // outside a region between two planes containing the required points is
    // excluded.  This is repeated for multiple plane orientations, determined by am.
    // ----------------------------------------------------
    // method = 2: finds all tissue points within a distance lb of the vessels, but
    // does not make a convex region.  Fills in 'holes' in interior, whether 2D or 3D.
    // ----------------------------------------------------
    // method = 3: as in method 1, but remove a strip at top and bottom of rectangular domain
    // ----------------------------------------------------
    // method = 5: retinal angiogenesis model. As in method 2, but exclude regions r < IR and r > OR
    // ----------------------------------------------------
    // Output:    nnt, total tissue node points inside the region.
    // nbou > 1 if inside region, value gives tissue point index
    // TWS 2010
    // Version 3.0, May 17, 2011.
    //-----------------------------------------------------------------/
    if (method == 1 || method == 3) {
        for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) nbou[i][j][k] = 1;
        am = 6;
        for (a = 0; a <= am; a++) for (b = -am; b <= am; b++) for (c = -am; c <= am; c++) {
            if (a != 0 || b != 0 || c != 0) {
                aa = 0.;
                bbp = 0.;
                cc = 0.;
                if (a != 0) aa = 1.0 / a;
                if (b != 0) bbp = 1.0 / b;
                if (c != 0) cc = 1.0 / c;
                abc = sqrt(SQR(aa) + SQR(bbp) + SQR(cc));
                ddmin = 1.e8;
                ddmax = -1.e8;
                for (iseg = 1; iseg <= nseg; iseg++)    if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
                    dd = (aa*start[1][iseg] + bbp * start[2][iseg] + cc * start[3][iseg]) / abc;
                    ddmin = FMIN(dd - rseg[iseg] - lb, ddmin);
                    ddmax = FMAX(dd + rseg[iseg] + lb, ddmax);
                    dd = (aa*end[1][iseg] + bbp * end[2][iseg] + cc * end[3][iseg]) / abc;
                    ddmin = FMIN(dd - rseg[iseg] - lb, ddmin);
                    ddmax = FMAX(dd + rseg[iseg] + lb, ddmax);
                }
                for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
                    t = (aa*axt[i] + bbp * ayt[j] + cc * azt[k]) / abc;
                    if (t > ddmax || t < ddmin) nbou[i][j][k] = 0;
                }
            }
        }
    }
    if (method == 2 || method == 5) {
        for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
            nbou[i][j][k] = 0;
            for (jseg = 1; jseg <= nseg; jseg++) if (segtyp[jseg] == 4 || segtyp[jseg] == 5) {
                x11 = (axt[i] - start[1][jseg])*scos[2][jseg] - (ayt[j] - start[2][jseg])*scos[1][jseg];
                x22 = (ayt[j] - start[2][jseg])*scos[3][jseg] - (azt[k] - start[3][jseg])*scos[2][jseg];
                x33 = (azt[k] - start[3][jseg])*scos[1][jseg] - (axt[i] - start[1][jseg])*scos[3][jseg];
                disp2 = SQR(x11) + SQR(x22) + SQR(x33);
                ds2 = SQR(axt[i] - start[1][jseg]) + SQR(ayt[j] - start[2][jseg]) + SQR(azt[k] - start[3][jseg]);
                de2 = SQR(axt[i] - end[1][jseg]) + SQR(ayt[j] - end[2][jseg]) + SQR(azt[k] - end[3][jseg]);
                if (FMAX(ds2, de2) - disp2 > SQR(lseg[jseg])) d = sqrt(FMIN(ds2, de2)) - rseg[jseg];
                else d = sqrt(disp2) - rseg[jseg];
                if (d < lb) nbou[i][j][k] = 1;
            }
        }
        for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) if (nbou[i][j][k] == 0) {
            imin = mxx;
            imax = 1;
            for (ii = 1; ii <= mxx; ii++) if (nbou[ii][j][k] == 1) {
                imin = IMIN(ii, imin);
                imax = IMAX(ii, imax);
            }
            jmin = myy;
            jmax = 1;
            for (jj = 1; jj <= myy; jj++) if (nbou[i][jj][k] == 1) {
                jmin = IMIN(jj, jmin);
                jmax = IMAX(jj, jmax);
            }
            kmin = mzz;
            kmax = 1;
            for (kk = 1; kk <= mzz; kk++) if (nbou[i][j][kk] == 1) {
                kmin = IMIN(kk, kmin);
                kmax = IMAX(kk, kmax);
            }
            if (i >= imin && i <= imax) if (j >= jmin && j <= jmax) if (k >= kmin && k <= kmax) nbou[i][j][k] = 1;
        }
        if (method == 5) {
            for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
                d = SQR(axt[i] - alx / 2.) + SQR(ayt[j] - aly / 2.);
                if (d < SQR(IR)) nbou[i][j][k] = 0;          //exclude region near origin
                if (d > SQR(OR)) nbou[i][j][k] = 0;          //exclude region at outer edge
            }
        }
        
    }
    //For 2D AV parallel vessels, set up exclusion of tisspoints 100 microns from the open "ends" as a boundary to stabiize the networks
    if (method == 3) {
        for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) if (ayt[j] >= (aly - 250.0)) nbou[i][j][k] = 0;
        for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) if (ayt[j] <= 250.0) nbou[i][j][k] = 0;
    }
    nnt = 0;
    for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) if (nbou[i][j][k] == 1) {
        nnt++;
        nbou[i][j][k] = nnt;
    }
    printf("Tissue node points used in simulation  = %i\n", nnt);
    printf("V of simulation region/V of cubic region = %f\n", nnt / (1.0*mxx*myy*mzz));
    
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on setuparrays2(nnv, nnt);
    //tissfix = imatrix(1, nnt, 1, nsp);//added September 2010
    tissfix = (int **) calloc(nnt+1, sizeof(int *));
    for(i=1;i<=nnt;i++) {
        tissfix[i] = (int *) calloc(nsp+1, sizeof(int ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            tissfix[i][j] = 34;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %d\n", nnt, i, nsp, j, tissfix[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //tisserr = matrix(1, nnt, 1, nsp);
    tisserr = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        tisserr[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            tisserr[i][j] = 35.06;
            
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, tisserr[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //dmtissdp = matrix(1, nnt, 1, nsp);
    dmtissdp = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        dmtissdp[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            dmtissdp[i][j] = -35.07;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, dmtissdp[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //mainseg = ivector(1, nnv);
    mainseg = (int *) calloc(nnv+1, sizeof(int ));
    for(i=1;i<=nnv;i++)
        mainseg[i] = -35;
        
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            printf("%d %d %d\n", nnv, i, mainseg[i]);
        }
    //-----------------------------------------------------------------------------/
    //indx = ivector(1, nnv + 1);        //added March 2010
    indx = (int *) calloc(nnv+1+1, sizeof(int ));
    for(i=1;i<=nnv;i++)
        indx[i] = -36;// this probably is an iterator index somewhere - SK
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            printf("%d %d %d\n", nnv, i, indx[i]);
    }
    //-----------------------------------------------------------------------------/
    //tisspoints = imatrix(1, 3, 1, nnt);
    tisspoints = (int **) calloc(3+1, sizeof(int *));
    for(i=1;i<=3;i++) {
        tisspoints[i] = (int *) calloc(nnt+1, sizeof(int ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nnt;j++)
            tisspoints[i][j] = -37;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nnt;j++)
                printf("%d %d %d %d %d\n", 3, i ,nnt, j, tisspoints[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //dtmin = vector(1, nnt);//added July 2011
    dtmin = (float *) calloc(nnt+1, sizeof(float ));
    for(i=1;i<=nnt;i++)
        dtmin[i] = -38.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            printf("%d %d %lf\n", nnt, i, dtmin[i]);
    }
    //-----------------------------------------------------------------------------/
    //rhstiss = dvector(1, nnt);    //added April 2016
    rhstiss = (float *) calloc(nnt+1, sizeof(float ));
    for(i=1;i<=nnt;i++)
        rhstiss[i] = -42.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            printf("%d %d %lf\n", nnt, i, rhstiss[i]);
    }
    //-----------------------------------------------------------------------------/
    //matxtiss = dvector(1, nnt);
    matxtiss = (float *) calloc(nnt+1, sizeof(float ));
    for(i=1;i<=nnt;i++)
        matxtiss[i] = -43.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            printf("%d %d %lf\n", nnt, i, matxtiss[i]);
    }
    
    //-----------------------------------------------------------------------------/
    //rhs = vector(1, nnv);            //added March 2010
    rhs = (float *) calloc(nnv+1, sizeof(float ));
    for(i=1;i<=nnv;i++)
        rhs[i] = -39.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            printf("%d %d %lf\n", nnv, i, rhs[i]);
    }
    
    //-----------------------------------------------------------------------------/
    //qvtemp = vector(1, nnv);
    qvtemp = (float *) calloc(nnv+1, sizeof(float ));
    for(i=1;i<=nnv;i++)
        qvtemp[i] = -40.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            printf("%d %d %lf\n", nnv, i, qvtemp[i]);
    }
    //-----------------------------------------------------------------------------/
    //sumal = vector(1, nnv);            //added August 2010
    sumal = (float *) calloc(nnv+1, sizeof(float ));
    for(i=1;i<=nnv;i++)
        sumal[i] = -41.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            printf("%d %d %lf\n", nnv, i, sumal[i]);
    }
    
    //-----------------------------------------------------------------------------/
    //qv = matrix(1, nnv, 1, nsp);
    qv = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        qv[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            qv[i][j] = -40004.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, qv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pv = matrix(1, nnv, 1, nsp);
    pv = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        pv[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            pv[i][j] = -40005.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, pv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pev = matrix(1, nnv, 1, nsp);
    pev = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        pev[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            pev[i][j] = -40006.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, pev[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pvt = matrix(1, nnv, 1, nsp);
    pvt = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        pvt[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            pvt[i][j] = -40007.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, pvt[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pvprev = matrix(1, nnv, 1, nsp);
    pvprev = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        pvprev[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            pvprev[i][j] = -40008.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, pvprev[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //qvprev = matrix(1, nnv, 1, nsp);
    qvprev = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        qvprev[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            qvprev[i][j] = -40009.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, qvprev[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //cv = matrix(1, nnv, 1, nsp);
    cv = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        cv[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            cv[i][j] = -40001.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, cv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //dcdp = matrix(1, nnv, 1, nsp);
    dcdp = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        dcdp[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            dcdp[i][j] = -40002.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, dcdp[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //cv0 = matrix(1, nnv, 1, nsp);
    cv0 = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        cv0[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            cv0[i][j] = -40003.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, cv0[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //conv0 = matrix(1, nnv, 1, nsp);
    conv0 = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        conv0[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nsp;j++)
            conv0[i][j] = -40004.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nsp, j, conv0[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //gvv = matrix(1, nnv, 1, nnv);
    gvv = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        gvv[i] = (float *) calloc(nnv+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nnv;j++)
            gvv[i][j] = -40005.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nnv;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nnv, j, gvv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //qt = matrix(1, nnt, 1, nsp);
    qt = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        qt[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            qt[i][j] = -40006.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, qt[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //pt = matrix(1, nnt, 1, nsp);
    pt = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        pt[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            pt[i][j] = -40007.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, pt[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //ptprev = matrix(1, nnt, 1, nsp);
    ptprev = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        ptprev[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            ptprev[i][j] = -40007.08;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, ptprev[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //ptv = matrix(1, nnt, 1, nsp);
    ptv = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        ptv[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            ptv[i][j] = -40008.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, ptv[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //qcoeff1 = matrix(1, nnt, 1, nsp);
    qcoeff1 = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        qcoeff1[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            qcoeff1[i][j] = -40009.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, qcoeff1[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //qtp = matrix(1, nnt, 1, nsp);    //added April 2016
    qtp = (float **) calloc(nnt+1, sizeof(float *));
    for(i=1;i<=nnt;i++) {
        qtp[i] = (float *) calloc(nsp+1, sizeof(float ));
    }
    for(i=1;i<=nnt;i++)
        for(j=1;j<=nsp;j++)
            qtp[i][j] = -40088.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnt;i++)
            for(j=1;j<=nsp;j++)
                printf("%d %d %d %d %lf\n", nnt, i ,nsp, j, qtp[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //ax = matrix(1, 3, 1, nnv);
    ax = (float **) calloc(3+1, sizeof(float *));
    for(i=1;i<=3;i++) {
        ax[i] = (float *) calloc(nnv+1, sizeof(float ));
    }
    for(i=1;i<=3;i++)
        for(j=1;j<=nnv;j++)
            ax[i][j] = -40098.09;
    
    if(DEBG==1) {
        for(i=1;i<=3;i++)
            for(j=1;j<=nnv;j++)
                printf("%d %d %d %d %lf\n", 3, i ,nnv, j, ax[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //al = matrix(1, nnv, 1, nnv);        //August 2010
    al = (float **) calloc(nnv+1, sizeof(float *));
    for(i=1;i<=nnv;i++) {
        al[i] = (float *) calloc(nnv+1, sizeof(float ));
    }
    for(i=1;i<=nnv;i++)
        for(j=1;j<=nnv;j++)
            al[i][j] = -40118.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv;i++)
            for(j=1;j<=nnv;j++)
                printf("%d %d %d %d %lf\n", nnv, i ,nnv, j, al[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //mat = dmatrix(1, nnv + 1, 1, nnv + 1);    //March 2010
    mat = (double **) calloc(nnv+1+1, sizeof(double *));
    for(i=1;i<=nnv+1;i++) {
        mat[i] = (double *) calloc(nnv+1+1, sizeof(double ));
    }
    for(i=1;i<=nnv+1;i++)
        for(j=1;j<=nnv+1;j++)
            mat[i][j] = -40118.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv+1;i++)
            for(j=1;j<=nnv+1;j++)
                printf("%d %d %d %d %lf\n", nnv+1, i ,nnv+1, j, mat[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //rhsg = dmatrix(1, nnv + 1, 1, 2);    //March 2010
    rhsg = (double **) calloc(nnv+1+1, sizeof(double *));
    for(i=1;i<=nnv+1;i++) {
        rhsg[i] = (double *) calloc(2+1, sizeof(double ));
    }
    for(i=1;i<=nnv+1;i++)
        for(j=1;j<=2+1;j++)
            rhsg[i][j] = -40118.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv+1;i++)
            for(j=1;j<=2+1;j++)
                printf("%d %d %d %d %lf\n", nnv+1, i ,2+1, j, rhsg[i][j]);
    }
    //-----------------------------------------------------------------------------/
    //rhsl = dvector(1, nnv + 1);        //March 2010
    rhsl = (double *) calloc(nnv+1+1, sizeof(double ));
    for(i=1;i<=nnv+1;i++)
        rhsl[i] = -340.09;
    
    if(DEBG==1) {
        for(i=1;i<=nnv+1;i++)
            printf("%d %d %lf\n", nnv+1, i, rhsl[i]);
    }
    //-----------------------------------------------------------------------------/
    //matx = dvector(1, nnv + 1);        //March 2010
    matx = (double *) calloc(nnv+1+1, sizeof(double ));
    for(i=1;i<=nnv+1;i++)
        matx[i] = -340.10;
    
    if(DEBG==1) {
        for(i=1;i<=nnv+1;i++)
            printf("%d %d %lf\n", nnv+1, i, matx[i]);
    }
    // the following two lines appear right after setuparrays2 in the Secomb code
    for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
    for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];
    
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on picturenetwork(nodvar, segvar, "Current\\NetNodesSegs.ps");

    //fname[80] = "netNodesSegs.ps"; I removed this line and just put the string into fopen directly because there were problems opening the file.
    xmin = 0.;
    xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));//repeated
    ymin = 0.;
    ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));//repeated
    
    cosp = matrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {    //set up matrix of direction cosines
        cosp[1][i] = (xsl1[i] - xsl0[i]) / xmax;
        cosp[2][i] = (xsl2[i] - xsl0[i]) / ymax;
    }
    cosp[3][1] = cosp[1][2] * cosp[2][3] - cosp[1][3] * cosp[2][2];
    cosp[3][2] = cosp[1][3] * cosp[2][1] - cosp[1][1] * cosp[2][3];
    cosp[3][3] = cosp[1][1] * cosp[2][2] - cosp[1][2] * cosp[2][1];
    //Determine range of z values
    zmin = 1.e6;
    zmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        zcoord = 0.;
        for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cosp[3][i];
        zmin = FMIN(zmin, zcoord - 1.);
        zmax = FMAX(zmax, zcoord + 1.);
    }

    picfac = FMIN(500. / xmax, 700. / ymax);//updated April 2010
    ofp = fopen("netNodesSegs", "w");
    if(ofp!=NULL) printf("netNodesSegs opened succesfully.\n");
    else printf("netNodesSegs could not be opened\n");

    fprintf(ofp, "%%!PS-Adobe-2.0\n");
    fprintf(ofp, "%%%%Pages: 1\n");
    fprintf(ofp, "%%%%EndComments\n");
    fprintf(ofp, "%%%%Page: 1 1\n");
    fprintf(ofp, "/mx {%g mul 50 add} def\n", picfac);
    fprintf(ofp, "/my {%g mul 50 add} def\n", picfac);//updated April 2010
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/sl {setlinewidth} def\n");
    fprintf(ofp, "/sc {setrgbcolor} def\n");
    fprintf(ofp, "/s {stroke} def\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "8 scalefont\n");
    fprintf(ofp, "setfont\n");
    
    fprintf(ofp, "newpath\n");
    fprintf(ofp, "%g mx %g my m\n", xmin, ymin);
    fprintf(ofp, "%g mx %g my l\n", xmax, ymin);
    fprintf(ofp, "%g mx %g my l\n", xmax, ymax);
    fprintf(ofp, "%g mx %g my l\n", xmin, ymax);
    fprintf(ofp, "closepath\n");
    fprintf(ofp, "stroke\n");
    //show tissue points
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) {
        showpoint = 0;
        for (k = 1; k <= mzz; k++) if (nbou[i][j][k] > 0) showpoint = 1;
        if (showpoint == 1) fprintf(ofp, "%g mx %g my m (.) show\n", axt[i], ayt[j]);
        for (k = 1; k <= mzz; k++) if (nbou[i][j][k] == 1022) fprintf(ofp, "%g mx %g my m (X) show\n", axt[i], ayt[j]);
    }
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "4 scalefont\n");
    fprintf(ofp, "setfont\n");
    //plot vessels according to segvar in order from bottom to top according to z-coordinate
    xzmin = 1.e6;
    xzmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        xzmin = FMIN(xzmin, segvar[iseg]);
        xzmax = FMAX(xzmax, segvar[iseg]);
    }
    for (ilevel = 1; ilevel <= nlevel; ilevel++) {
        zbottom = zmin + (ilevel - 1)*(zmax - zmin) / nlevel;
        ztop = zmin + ilevel * (zmax - zmin) / nlevel;
        for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
            zcoord = 0.;
            for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cosp[3][i];
            if (zcoord >= zbottom && zcoord < ztop) {
                if (xzmin != xzmax) xz = (segvar[iseg] - xzmin) / (xzmax - xzmin);
                else xz = 0.75;
                blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
                green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
                red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
                fprintf(ofp, "%f %f %f sc\n", red, green, blue);
                fprintf(ofp, "%g sl\n", picfac*diam[iseg] * diamfac);//line widths scaled up by diamfac
                xs = 0.;
                ys = 0.;
                for (i = 1; i <= 3; i++) {
                    xs += (cnode[i][ista[iseg]] - xsl0[i])*cosp[1][i];
                    ys += (cnode[i][ista[iseg]] - xsl0[i])*cosp[2][i];
                }
                fprintf(ofp, "%g mx %g my m ", xs, ys);
                xs = 0.;
                ys = 0.;
                for (i = 1; i <= 3; i++) {
                    xs += (cnode[i][iend[iseg]] - xsl0[i])*cosp[1][i];
                    ys += (cnode[i][iend[iseg]] - xsl0[i])*cosp[2][i];
                }
                fprintf(ofp, "%g mx %g my l s \n", xs, ys);
            }
        }
    }

    //label nodes in black
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for (inod = 1; inod <= nnod; inod++) {
        xs = 0.;
        ys = 0.;
        for (i = 1; i <= 3; i++) {
            xs += (cnode[i][inod] - xsl0[i])*cosp[1][i];
            ys += (cnode[i][inod] - xsl0[i])*cosp[2][i];
        }
        //comment out next two lines to remove node numbers
        fprintf(ofp, "%g mx %g my m ", xs + 0.5 / picfac, ys);
        fprintf(ofp, "(%g) show\n", nodvar[inod]);
    }
    //label segments in blue
    fprintf(ofp, "0 0 1 setrgbcolor\n");//blue
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        xs = 0.;
        ys = 0.;
        for (i = 1; i <= 3; i++) {
            xs += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*cosp[1][i];
            ys += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*cosp[2][i];
        }
        //comment out next two lines to remove segment numbers
        fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac, ys);
        fprintf(ofp, "(%g) show\n", segvar[iseg]);
    }
    //create a color bar
    float cp;
    float cbbox = 15.; //size of boxes
    float cbx = 560; //origin of color bar
    float cby = 100;//origin of color bar
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "8 scalefont\n");
    fprintf(ofp, "setfont\n");
    for (k = 0; k <= 10; k++) {
        xz = k * 0.1;
        cp = xzmin + (xzmax - xzmin)*xz;
        blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
        green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
        red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
        fprintf(ofp, "%f %f %f setrgbcolor\n", red, green, blue);
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx, cby + k * cbbox, cbx + cbbox, cby + k * cbbox, cbx + cbbox, cby + (k + 1)*cbbox, cbx, cby + (k + 1)*cbbox);
        if (k > 0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n", cbx + cbbox * 1.1, cby + cbbox * (k - 0.1), cp);
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx, cby, cbx + cbbox, cby, cbx + cbbox, cby + cbbox * 11, cbx, cby + cbbox * 11);
    fprintf(ofp, "showpage\n");
    fclose(ofp);
    free_matrix(cosp, 1, 3, 1, 3);
    
    for (iseg = 1; iseg <= nseg; iseg++)
        segvar[iseg] = log(fabs(qdata[iseg]));
    
    cmgui(segvar);  //<------------------------------------
    
    ofp = fopen("CurrentSummary.out", "w");
    //print headings for summary output file
    fprintf(ofp, "imain kmain ");
    for (j = 1; j <= nvaryparams; j++) {
        switch (ivaryparams[j][1]) {
            case 1:
            {
                fprintf(ofp, "   q0fac    ");
                break;
            }
            case 2:
            {
                fprintf(ofp, " solutefac[%i]", ivaryparams[j][2]);
                break;
            }
            case 3:
            {
                fprintf(ofp, " diff[%i]     ", ivaryparams[j][2]);
                break;
            }
            case 4:
            {
                fprintf(ofp, " intravascfac[%i]", ivaryparams[j][2]);
                break;
            }
            case 5:
            {
                fprintf(ofp, " tissparam[%i][%i]", ivaryparams[j][2], ivaryparams[j][3]);
                break;
            }
        }
    }
    for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "  pmean[%i]  ", isp);
    for (j = 1; j <= npostgreensout; j++) fprintf(ofp, " postgreens[%i]", j);
    fprintf(ofp, "\n");
    //printf("nruns==%f\n", nruns);
    
    //The following loop allows running a series of cases with varying parameters
    for (imain = 1; imain <= nruns; imain++) {
        sprintf(numstr, "%03i", imain);    //need 3-digit frame number for file name. November 2016
        for (j = 1; j <= nvaryparams; j++) {
            switch (ivaryparams[j][1])
            {
                case 1:
                {
                    q0fac = paramvalue[imain][j];
                    break;
                }
                case 2:
                {
                    isp = ivaryparams[j][2];    //updated November 2016
                    if (isp <= nsp) solutefac[isp] = paramvalue[imain][j];
                    break;
                }
                case 3:
                {
                    isp = ivaryparams[j][2];
                    if (isp <= nsp) diff[isp] = paramvalue[imain][j];
                    break;
                }
                case 4:
                {
                    isp = ivaryparams[j][2];
                    if (isp <= nsp) intravascfac[isp] = paramvalue[imain][j];
                    break;
                }
                case 5:
                {
                    isp = ivaryparams[j][3];
                    if (isp <= nsp) tissparam[ivaryparams[j][2]][isp] = paramvalue[imain][j];
                    break;
                }
            }
        }
    }
    
    //-----------------------------------------------------------------------------/
    //-----------------------------------------------------------------------------/
    // This is where I begin working on greens();
    //-----------------------------------------------------------------------------/
    // Green's function approach for multiple reacting species.
    // T.W. Secomb July 2007 - based on Greens.f by R. Hsu.
    // See http://www.physiology.arizona.edu/people/secomb/greens.html
     
    // Variable array sizes and bounds, using Numerical Recipes utilities.
    // Tissue-vessel and tissue-tissue matrices computed on the fly to save memory.
    // No non-dimensionalization.  Lengths, diameters in microns, times in s.
    // Flows in nanoliters/min
    // Oxygen concentrations in cm^3/cm^3
    // Consumption rates in cm^3/cm^3/s - changed in V2
    // Partial pressures in mmHg
    // Diffusivities in cm2/s, converted to micron^2/s
     
    // Special parameters for oxygen
    // p50, fn --- parameters in the Hill equation
    // cs --- red blood cell oxygen binding capacity in cm^3 O2/cm^3
    // alphab --- average solubility in blood in cm^3 O2/cm^3/mmHg
    // gamma1 --- intravascular resistance, varies with vessel diameter, in mmHg.cm.s/cm^3 O2
     
    // Main variables:
    // gvv --- Green's function matrix for vessels
    // mat --- matrix for vessel strengths
    // al --- matrix giving dependence of vessel convective fluxes on source strengths
    // lseg --- segment length
    // ds --- subsegment length
    // qv --- oxygen efflux from subsegment
    // pv --- PO2 in the subsegment
    // cv --- oxygen concentration in the subsegment
    // qt --- tissue source strength
    // pvt --- PO2 on vessels due to source terms in tissue
    // ptv --- PO2 in tissue due to source terms on vessels
    // q --- flow rate, qq = abs(q)
     
    // Version 2.0, May 2010.
    // With 9/08 updates.  New vessel-vesel interaction coefficient. January 2009
    // With alternate terms for 2D version.  May 2009
    // With g0 computed as part of linear system, for permeable solutes.  March 2010
    // g0method = 1:  include g0 in linear system to be solved - fastest method *****
    // g0method = 2:  theoretical estimates of dqsum/dg0 - was used in Version 1
    // For impermeable solutes, method 2 is always used
    // With choice of Gauss-Jordan, LU or biconjugate gradient linear solvers.  March 2010
    // linmethod = 1:  Gaussian elimination - was used in Version 1
    // linmethod = 2:  LU decomposition
    // linmethod = 3:  biconjugate gradient (iterative) - fastest method *****
    // Does not require that species 1 is oxygen. Allows for non-diffusing solutes.  April 2010
    // Creates log file. April 2010.
    // During tissue loop, scales vessel sources so that qvsum = qtsum, for faster convergence.  April 2010.
    // Modified for compatibility with Mac XCode compiler.  April 2010.
    // Includes intravascular resistance for all solutes.  May 2010.
    // Includes non-diffusible solutes.  May 2010.
     
    // Version 3.0, May 17, 2011.
    // Uses convect instead of genalpha and genalphahd.
    // This gives improved convergence if hematocrit is non-uniform in network
     
    // Version 4.0, March 1, 2018.
    // Includes cmgui to generate files needed for visualization using CMGUI
    // In soluteparams.dat, the parameter for total inflow is replaced by a flow factor
    // Includes capability to run multiple cases, controlled by VaryParams.dat
    // Includes capability to compute dependent variables (e.g. cell survival) after
    // greens has run, as specified by PostGreensParams.dat
    // Problem-specific c code (used as #include files) is needed as follows:
    // tissrate.cpp.txt specifies dependence of reaction rates on concentrations
    // postgreens.cpp.dat specifies computation of dependent variables
    // Produces a summary output file from multiple runs (summary.out)
    // Output files (TissueLevels.dat, TissueSources.dat, VesselLevels.dat, VesselSources.dat)
    // are formatted differently, with one line per tissue point or vessel point
    // Previous combined output file (GreensRes.out) is no longer generated
    // Capability to restart a run has been removed (was not functional)
    // All output files and a copy of input files are placed in a folder called "Current"
    // Algorithm for non-diffusible solutes is improved by using reaction rates based on
    // updated levels of other solutes at each iteration of Newton method
    //-----------------------------------------------------------------------------/
    // Some parameters for the solution methods
    w2d = alz;  //needed for 2d
    r2d = sqrt(SQR(alx) + SQR(aly) + SQR(alz));// calculates the radius of a sphere centered at the origin and having the point (alx,aly,alz) on its surface
    if(DEBG==1) {
        printf("%f %f %f\n", alx, SQR(alx), r2d);
    }
    lam3d = 0.2;//    lam3d = 0.2;    //Underrelax iteration of tissue levels
    lam2d = 0.1;    //Use this one for 2D permeable solutes only.  May 2010
    if (is2d) lam = lam2d;
    else lam = lam3d;
    int bicgstabit = 2000; //parameter for biconjugate gradient method.  March 2010
    double bicgstaberr = 0.0001; //parameter for biconjugate gradient method.  March 2010
    
    clock_t tstart, tfinish, tstart1, tfinish1; // Used to time the speed of the calculation
   
    //setup mainseg (must be done after setuparrays2)
    //Identifies the segment that each vessel point belongs to
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5)
        for (i = 0; i < nspoint[iseg]; i++) mainseg[istart[iseg] + i] = iseg;
    //identify vessel points
    for (i = 1; i <= nnv; i++) {
        iseg = mainseg[i];
        isn = i - istart[iseg]; // each vessel point is assigned a number based on it's position along the segment that it is a part of. The numbering starts from zero (0)
        for (j = 1; j <= 3; j++) ax[j][i] = start[j][iseg] + scos[j][iseg] * ds[iseg] * (isn + 0.5);// the x,y,z coordinates are determined for each vessel point
    }
    //index tissue points
    for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
        nt = nbou[i][j][k];
        if (nt > 0) {
            tisspoints[1][nt] = i;
            tisspoints[2][nt] = j;
            tisspoints[3][nt] = k;
        }
    }
    //calculate the distance of tissue points to the nearest vessel.
    dtave = 0.;
    for (itp = 1; itp <= nnt; itp++) {
        i = tisspoints[1][itp];
        j = tisspoints[2][itp];
        k = tisspoints[3][itp];
        dtmin[itp] = 1.e6;
        for (jseg = 1; jseg <= nseg; jseg++) if (segtyp[jseg] == 4 || segtyp[jseg] == 5) {
            x11 = (axt[i] - start[1][jseg])*scos[2][jseg] - (ayt[j] - start[2][jseg])*scos[1][jseg];
            x22 = (ayt[j] - start[2][jseg])*scos[3][jseg] - (azt[k] - start[3][jseg])*scos[2][jseg];
            x33 = (azt[k] - start[3][jseg])*scos[1][jseg] - (axt[i] - start[1][jseg])*scos[3][jseg];
            disp2 = SQR(x11) + SQR(x22) + SQR(x33);
            //axt[i] is the x-coordinate of itp, since i is the x iterator in tisspoints.
            //start[1][jseg] is the x-coordinate of the starting node (one with the lower number) of  jseg.
            ds2 = SQR(axt[i] - start[1][jseg]) + SQR(ayt[j] - start[2][jseg]) + SQR(azt[k] - start[3][jseg]);
            de2 = SQR(axt[i] - end[1][jseg]) + SQR(ayt[j] - end[2][jseg]) + SQR(azt[k] - end[3][jseg]);
            if (FMAX(ds2, de2) - disp2 > SQR(lseg[jseg])) d = FMAX(sqrt(FMIN(ds2, de2)) - rseg[jseg], 0.);
            else d = FMAX(sqrt(disp2) - rseg[jseg], 0.);
            if (d < dtmin[itp]) dtmin[itp] = d;
        }
        dtave += dtmin[itp];
    }
    dtave = dtave / nnt; // the distance between tissue points and the nearest vessel, averaged over all of the tissue points
    vdom = nnt * vol; // The volume of the domain, calculated by multiplying the volume of each tissue mesh region by the total number of tissue points
    tlength = 0.; // holds the total length of all vessels in the network
    tlengthq = 0.;//added 8/09 // The sum of the absolute values of the flow over all segments (the flow per length is weighted by the length in the calculation to get it)
    tlengthqhd = 0.;//added 8/10 // Total volume of RBCs flowing through the network
    for (iseg = 1; iseg <= nseg; iseg++) {
        
        q[iseg] = qdata[iseg] * q0fac;//scaling by q0fac. November 2016
        qq[iseg] = fabs(q[iseg]);
        tlength += lseg[iseg];//moved to greens November 2016
        tlengthq += lseg[iseg] * qq[iseg];//added 8/09
        tlengthqhd += lseg[iseg] * qq[iseg] * hd[iseg];//added 8/10
    }
    if(DEBG==1) for (iseg = 1; iseg <= nseg; iseg++) {
        printf("%f %f %f %f %f\n", q[iseg], qq[iseg], tlength, tlengthq, tlengthqhd);
    }
    den = sqrt(vdom / tlength); //reported, but not used in the later code
    if (greensverbose) {
        printf("The volume of the tissue domain is %7.5f \n", vdom);
        printf("Average distance from tissue node to the nearest vessel = %f\n", dtave);
        printf("Sqrt(Tissue Volume/vessel length) = %f\n", den);
        printf("Capillary density = %8.1f /mm2\n", tlength / vdom * 1.e6); // the factor 1.e6 converts 1/microns^2 to 1/mm^2
        printf("q0fac = %f\n", q0fac);
        printf("Perfusion = %f cm3/cm3/min (uncorrected for path length effect)\n", totalq*q0fac / vdom * 1.e6);    //added q0fac November 2016. The 1.e6 factor is 1.e-6/1.e-12, to convert the units of flow from nL/min to cm^3/min and convert the units of volume from micrometers^3 to cm^3.
        //printf("The total blood flow to the region is %f\n", totalq);
        pathlength = 0.;
        if (q0fac > 0.) for (iseg = 1; iseg <= nseg; iseg++) pathlength += fabs(q[iseg])*lseg[iseg] / totalq / q0fac;    //added q0fac November 2016
        printf("Flow-weighted path length = %f micron\n", pathlength);
    }
    
    //Calculate intravascular or wall transport resistance. Zero unless specified in simpleIntravascRes.
    //If not oxygen, assume value from data is 1/(wall permeability in um/s)
    //nresis[isp] specifies the number of different vessel diameters for which there is a wall resistance specified in simpleIntravascRes - Baran
    // If nresis[isp] = 0, that means that no resistance is specified because solute #isp does not permeate the vessel wall - Baran
    // If nresis[isp] > 1 and if the diameter of a segment falls between two of the values for which resistances are specified, the resistance is caluculated by a linear interpolation - Baran
    // The resistances for the current network, that result for each segment and each solute, are stored in gamma1[iseg][isp]
    for (isp = 1; isp <= nsp; isp++) {
        for (iseg = 1; iseg <= nseg; iseg++) {
            gamma1[iseg][isp] = 0.;
            if (nresis[isp] != 0) {
                gamma1[iseg][isp] = resis[1][isp];// assign the first value that appears in simpleIntravascRes as a possible value for resistance (regardless of the segment diameter), this here to take care of the case when there is only one possibility. The next bit of code assigns the resistance based on the segment diameter for solutes for which there is a list of diameters and their resistances provided. - Baran
                for (j = 2; j <= nresis[isp]; j++) if (diam[iseg] <= resisdiam[j][isp] && diam[iseg] > resisdiam[j - 1][isp])
                    gamma1[iseg][isp] = resis[j - 1][isp] + (resis[j][isp] - resis[j - 1][isp])
                    *(diam[iseg] - resisdiam[j - 1][isp]) / (resisdiam[j][isp] - resisdiam[j - 1][isp]);
                if (diam[iseg] > resisdiam[nresis[isp]][isp]) gamma1[iseg][isp] = resis[nresis[isp]][isp];
                if (oxygen[isp] != 1) gamma1[iseg][isp] = gamma1[iseg][isp] / pi1 / diam[iseg];    //was /2./pi1/diam[isp]; fixed July 2013
                gamma1[iseg][isp] = gamma1[iseg][isp] * intravascfac[isp];        //scaling from VaryParams.dat
            }
        }
    }
    if(DEBG==1) {
        for(isp=1;isp<=nsp;isp++) printf("intravasfac[%i] = %f\n", isp, intravascfac[isp]);
        }
    
    //vessel ~ vessel matrix elements gvv
    //Uses empirical fit to results from elliptical integral form for diagonal elements, updated 2009
    //if center of one segment lies within the other segment, calculate gvv as for self-interaction term
    //based on larger radius and larger length. (Jan. 08)
    for (i = 1; i <= nnv; i++) {
        iseg = mainseg[i];
        for (j = 1; j <= nnv; j++) {
            jseg = mainseg[j];
            //this section modified to give better behavior for 'hairpin' structures
            dist = sqrt(SQR(ax[1][j] - ax[1][i]) + SQR(ax[2][j] - ax[2][i]) + SQR(ax[3][j] - ax[3][i])); // the distance between the vessel points
            if (dist < FMAX(sqrt(ds[iseg] * rseg[iseg]), sqrt(ds[jseg] * rseg[jseg]))) {
                dsmax = FMAX(ds[iseg], ds[jseg]);
                rsegmax = FMAX(rseg[iseg], rseg[jseg]);
                //Version 2.0 of 3-D interaction coefficients for close or coincident segments.  See Sep. 2009 notes.
                gvarfac = 0.6*exp(-0.45*dsmax / rsegmax);
                //for distinct vessels close together, make distance rsegmax in following calculation, to improve convergence. TWS2011
                if (iseg != jseg) dist = rsegmax;
                gvv[i][j] = (1.298 / (1. + 0.297*powf(dsmax / rsegmax, 0.838)) - gvarfac * SQR(dist / rsegmax))*fac / rsegmax; //  fac = 1. / 4. / pi1
                //for 2D version, additional terms give effect of boundaries (reflection method)
                if (is2d) {
                    gvv[i][j] -= fac * 2. / w2d * 0.926*SQR(1. - 1. / (1. + 0.36*dsmax / w2d))*powf(1. + dsmax / w2d, 0.27);
                    gvv[i][j] += fac * 2. / w2d * (log(r2d / w2d + 0.5 + 0.27 / r2d * w2d) - 0.117);
                }
            }
            else {
                if (is2d) gvv[i][j] = fac * 2. / w2d * log(r2d / dist);
                else gvv[i][j] = fac / dist; // gvv[i][j] = 1/4pi|(x-x*)|,  notice that there are no G_2 or G_3 terms here that would be required if the boundary condition were net no flux (G_2 needed) or pointwise no flux (G_2 and G_3 needed)
            }
        }
    }
    // tissue ~ vessel, vessel ~ tissue: compute matrix elements gvt, gtv on fly as needed
    // tissue ~ tissue: construct matrix of distances from a corner node
    //diagonal elements of tissue ~ tissue matrix gtt
    if (is2d) gtt = fac / w2d * (2.*log(r2d / req) + 0.5);
    else gtt = 1.2*fac / req;
    for (jx = 1; jx <= mxx; jx++) for (jy = 1; jy <= myy; jy++) for (jz = 1; jz <= mzz; jz++) {
        dist = sqrt(SQR(axt[1] - axt[jx]) + SQR(ayt[1] - ayt[jy]) + SQR(azt[1] - azt[jz]));
        if (jx*jy*jz != 1) {
            if (is2d) dtt[jx][jy][jz] = fac * 2. / w2d * log(r2d / dist);
            else dtt[jx][jy][jz] = fac / dist;
        }
        else dtt[jx][jy][jz] = gtt;
    }
    //detect and label vessels with very low q*hd - updated April 2010 - test if oxygen is one of the solutes
    qhdcrit = 0.;
    for (isp = 1; isp <= nsp; isp++) if (oxygen[isp] == 1) qhdcrit = lowflowcrit * tissparam[1][isp];
    for (iseg = 1; iseg <= nseg; iseg++) {
        lowflow[iseg] = 0;
        if (qq[iseg] * (hd[iseg] + 0.01) < qhdcrit) lowflow[iseg] = 1;//Added 0.01 September 2011 to allow for high flow, zero hematocrit channels
    }
    /****************************************************************************
     initgreens - initial tissue source strengths, given uniform solute field, initial g0
     initial vessel source strengths based on uniform efflux rate from all vessels
     use this version only if values are not available from previous call to greens
     TWS Jan 08
     Version 2.0, May 1, 2010.
     Version 3.0, May 17, 2011.
     "Initial estimates of tissue oxygen levels and vessel source
     strengths are defined, and the tissue source strengths φj are
     computed using Eq. (18)."
    *******************************************************************************/
    initgreens();
    /************************************************************************
     putrank - generate list of nodes in order of flow direction
     nodrank --- if nodrank[i] < nodrank[j], node j is not upstream of node i
     So node j is upstream of node i (node j feeds solute to node i) only if nodrank[i] >= nodrank[j] - Baran
     The node ranks are probably used to implement equation 20 from the paper - Baran
     Version 2.0, May 1, 2010.
     Version 3.0, May 17, 2011.
     "For a series of
     segments forming an unbranched vessel, αi j takes the value
     1 if segment j is upstream of segment i , 12
     if i = j , and
     0 otherwise. For segments forming a branched network,
     the values of αij depend on the partition of oxygen fluxes
     at diverging bifurcations and the summation of fluxes at
     converging bifurcations."
     *************************************************************************/
    putrank();
   
    //-----------------------------------------------------------------------------/
    if(DEBG==1){
        for(i=1;i<=nnod;i++) printf("%d %d %d\n", nnod, i, nodrank[i]);
    }
    //greens() continues here
    //    for(isp=1; isp<=nsp; isp++){//for testing purposes only
    //        convect(isp);
    //        testconvect(isp);
    //    }

    //create log file
    ofp1 = fopen("Current\\GreensLog.txt", "w");
    fprintf(ofp1, "GreensLog.txt\n");
    fclose(ofp1);
    tstart = clock();
    //-------------------- start of main loop -----------------------
    for (kmain = 1; kmain <= nmax; kmain++) {
        tstart1 = clock();
        if (greensverbose) printf("\n----- kmain = %i -----\n", kmain);
        else printf(" %i", kmain);
        for (isp = 1; isp <= nsp; isp++) {
            for (itp = 1; itp <= nnt; itp++)    ptprev[itp][isp] = pt[itp][isp];
            if (permsolute[isp] == 1) for (i = 1; i <= nnv; i++) pvprev[i][isp] = pv[i][isp];
            g0old[isp] = g0[isp];
        }
        //-------------------- start of vessel loop -----------------------
        //compute contribution pvt from tissue source strengths qt
        for (i = 1; i <= nnv; i++) {
            for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
                if (g0method == 1 && permsolute[isp] == 1) pvt[i][isp] = 0.;
                else pvt[i][isp] = g0[isp];
            }
            for (itp = 1; itp <= nnt; itp++) {
                dist = sqrt(SQR(ax[1][i] - axt[tisspoints[1][itp]])
                            + SQR(ax[2][i] - ayt[tisspoints[2][itp]])
                            + SQR(ax[3][i] - azt[tisspoints[3][itp]]));
                if (dist <= req) {
                    if (is2d) gvt = fac / w2d * (2.*log(r2d / req) + 1. - SQR(dist / req));
                    else gvt = fac * (1.5 - 0.5*SQR(dist / req)) / req;
                }
                else {
                    if (is2d) gvt = fac * 2. / w2d * log(r2d / dist);
                    else gvt = fac / dist;
                }
                for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) pvt[i][isp] += gvt / diff[isp] * qt[itp][isp];//May 2010 permsolute
            }
        }
        //compute blood solute levels and PO2
        for (kvessel = 1; kvessel <= nmaxvessel; kvessel++) {
            convflagv = 1;
            for (isp = 1; isp <= nsp; isp++) {// for each solute (isp)
                qvsum[isp] = 0.;
                dqvsumdg0[isp] = 0.;
                if (permsolute[isp] == 1) {// only for the permeating solutes
                    ineg = 0;
                    ihigh = 0;
                    convect(isp);// Recursively compute theive flux of solute isp throughout the network and use it to compute the al[i][j] (see equation 20 of the paper) as well as the concentrations of solute in each segment.  - Baran
                    for (i = 1; i <= nnv; i++) {// for each vessel segment
                        iseg = mainseg[i];// get the associating segment index
                        qvprev[i][isp] = qv[i][isp];// save the previous source strengths
                        if (oxygen[isp] == 1) {// only for oxygen
                            if (lowflow[iseg] != 1) {//only do this if not a lowflow segment.
                                if (cv[i][isp] < 0.) {// if the oxygen concentration is negative
                                    ineg++;
                                    if (ineg == 1 && greensverbose) printf("*** Warning: cblood is negative -%i", segname[iseg]);
                                    if (ineg > 1 && greensverbose) printf("-%i", segname[iseg]);
                                }
                                if (cv[i][isp] > bloodconc(150., hd[iseg])) {// if the oxygen concentration is high
                                    ihigh++;
                                    if (ihigh == 1 && greensverbose) printf("*** Warning: cblood is high +%i", segname[iseg]);
                                    if (ihigh > 1 && greensverbose) printf("+%i", segname[iseg]);
                                }
                                blood(cv[i][isp], hd[iseg], &pv[i][isp], &dcdp[i][isp]);
                                if(DEBG==0){
                                    printf("nnv = %d, i = %d, isp = %d\n", nnv, i, isp);
                                    printf("The concentration of O2 is = %5.40f\n", cv[i][isp]);
                                    printf("The partial pressure of O2 in the blood is = %5.40f\n", pv[i][isp]);
                                    printf("The partial pressure of O2 in the blood is = %5.40f\n", p);
                                }
                            }
                        }
                        else {
                            pv[i][isp] = cv[i][isp];// for solutes that are not oxygen
                            dcdp[i][isp] = 1.; // for solutes that are not oxygen
                        }
                    }
                    if (ineg > 0 || ihigh > 0) if (greensverbose) printf("\n");
                    //generate linear system to be solved
                    for (i = 1; i <= nnv; i++) {// for each vessel
                        iseg = mainseg[i]; // get the associated segment index
                        rhs[i] = pv[i][isp] - pvt[i][isp];
                        for (j = 1; j <= nnv; j++) {
                            jseg = mainseg[j];
                            mat[i][j] = gvv[i][j] / diff[isp] + al[i][j] / dcdp[i][isp] / qq[iseg] / flowfac;
                            if (i == j) mat[i][j] += gamma1[iseg][isp] / ds[iseg];
                            rhs[i] += al[i][j] * qv[j][isp] / dcdp[i][isp] / qq[iseg] / flowfac;
                            if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1) {  //low q*hd
                                if (i == j) mat[i][j] = 1.;
                                else mat[i][j] = 0.;
                            }
                        }
                        if (oxygen[isp] == 1 && lowflow[iseg] == 1) rhs[i] = qvprev[i][isp];  //low q*hd
                    }
                    
                    //solve system of linear algebraic equations: Sum mat[i][j]*qv[j]= rhs[i]
                    if (g0method == 1) {
                        for (i = 1; i <= nnv; i++) {
                            if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1) mat[i][nnv + 1] = 0.;  //low q*hd
                            else mat[i][nnv + 1] = 1.;
                            mat[nnv + 1][i] = 1.;
                            mat[nnv + 1][nnv + 1] = 0.;
                        }
                        if (linmethod == 1) {
                            for (i = 1; i <= nnv; i++) rhsg[i][1] = rhs[i];
                            rhsg[nnv + 1][1] = -qtsum[isp];
                            gaussj(mat, nnv + 1, rhsg, 1);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = rhsg[i][1];
                                qvsum[isp] += qv[i][isp];
                            }
                            g0[isp] = rhsg[nnv + 1][1];
                        }
                        if (linmethod == 2) {
                            ludcmp(mat, nnv + 1, indx, &dd);
                            for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
                            rhsl[nnv + 1] = -qtsum[isp];
                            lubksb(mat, nnv + 1, indx, rhsl);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = rhsl[i];
                                qvsum[isp] += qv[i][isp];
                            }
                            g0[isp] = rhsl[nnv + 1];
                        }
                        if (linmethod == 3) {
                            for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
                            rhsl[nnv + 1] = -qtsum[isp];
                            for (i = 1; i <= nnv; i++) matx[i] = qv[i][isp];
                            matx[nnv + 1] = g0[isp];
                            bicgstab(mat, rhsl, matx, nnv + 1, bicgstaberr, bicgstabit);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = matx[i];
                                qvsum[isp] += qv[i][isp];
                            }
                            g0[isp] = matx[nnv + 1];
                        }
                    }
                    if (g0method == 2) {
                        if (linmethod == 1) {
                            for (i = 1; i <= nnv; i++) {
                                rhsg[i][1] = rhs[i];
                                rhsg[i][2] = -1.;
                            }
                            gaussj(mat, nnv, rhsg, 2);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = rhsg[i][1];
                                qvsum[isp] += qv[i][isp];
                                if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsg[i][2];
                            }
                        }
                        if (linmethod == 2) {
                            ludcmp(mat, nnv, indx, &dd);
                            for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
                            lubksb(mat, nnv, indx, rhsl);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = rhsl[i];
                                qvsum[isp] += qv[i][isp];
                            }
                            for (i = 1; i <= nnv; i++) rhsl[i] = -1.;
                            lubksb(mat, nnv, indx, rhsl);
                            for (i = 1; i <= nnv; i++)
                                if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsl[i];
                        }
                        if (linmethod == 3) {
                            for (i = 1; i <= nnv; i++) {
                                rhsl[i] = rhs[i];
                                matx[i] = qv[i][isp];
                            }
                            bicgstab(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
                            for (i = 1; i <= nnv; i++) {
                                qv[i][isp] = matx[i];
                                qvsum[isp] += qv[i][isp];
                            }
                            for (i = 1; i <= nnv; i++) {
                                rhsl[i] = -1.;
                                matx[i] = 0.;
                            }
                            bicgstab(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
                            for (i = 1; i <= nnv; i++)
                                if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += matx[i];
                        }
                    }
                    //for low q*hd segments, calculate efflux based on change in extravascular oxygen level
                    //save values in qvtemp to avoid influence on eval, update qv, underrelax - July 2008
                    for (i = 1; i <= nnv; i++) {
                        iseg = mainseg[i];
                        if (oxygen[isp] == 1 && lowflow[iseg] == 1) {
                            for (j = 1; j <= 3; j++) x[j] = ax[j][i] - 0.5*scos[j][iseg] * ds[iseg];
                            p = eval(1, req, x);
                            p[isp] = FMAX(p[isp], 0.);
                            pv[i][isp] = p[isp] / 2.;   //added August 2009
                            qvtemp[i] = q[iseg] * flowfac*bloodconc(p[isp], hd[iseg]); //q here (not qq) April 2008
                            for (j = 1; j <= 3; j++) x[j] = ax[j][i] + 0.5*scos[j][iseg] * ds[iseg];
                            p = eval(1, req, x);
                            p[isp] = FMAX(p[isp], 0.);
                            pv[i][isp] += p[isp] / 2.;  //added August 2009
                            qvtemp[i] -= q[iseg] * flowfac*bloodconc(p[isp], hd[iseg]);
                            cv[i][isp] = bloodconc(pv[i][isp], hd[iseg]);    //added October 2012
                        }
                    }
                    for (i = 1; i <= nnv; i++) if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1)
                        qv[i][isp] = 0.5*qvtemp[i] + 0.5*qvprev[i][isp];    //underrelax
                    errvessel[isp] = 0.;
                    imaxerr = 0;
                    errvesselcount[isp] = 0;
                    for (i = 1; i <= nnv; i++) {
                        dif = qv[i][isp] - qvprev[i][isp];
                        //If qv is large, use relative rather than absolute error
                        if (qv[i][isp] != 0.) dif = dif * FMIN(1., epsvessel[isp] / errfac / fabs(qv[i][isp]));
                        if (fabs(dif) >= errvessel[isp]) {    //>= added Dec. 2011
                            imaxerrvessel[isp] = mainseg[i];
                            errvessel[isp] = fabs(dif);
                        }
                        if (fabs(dif) > epsvessel[isp]) errvesselcount[isp]++;
                    }
                    if (greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
                    if (greensverbose) printf("Solute %i: kvessel = %i, errvessel_q = %f, imaxerr = %i, g0 = %f\n",
                                              isp, kvessel, errvessel[isp], imaxerrvessel[isp], g0[isp]);
                    if (errvesselcount[isp] > 0) convflagv = 0;
                }
            }
            if (convflagv) goto vesselconv;
        }
        for (isp = 1; isp <= nsp; isp++) if (errvesselcount[isp] > 0)
            if (!greensverbose) printf("*** Warning: solute %i, %i vessel source strengths not converged\n",
                                      isp, errvesselcount[isp]);
    vesselconv:;
        //--------------------  end of vessel loop --------------------
        //--------------------  start of tissue loop --------------------
        //Compute tissue source strengths iteratively by successive relaxation: updated qt values are immediately used.
        //Continually scales up qv values so that their sum equals updated sum of qt values. (required for consistency, see equation (19) in the paper - Baran)
        //contribution ptv from vessel source strengths qv
        for (itp = 1; itp <= nnt; itp++) {
            for (isp = 1; isp <= nsp; isp++) ptv[itp][isp] = 0.;
            for (i = 1; i <= nnv; i++) {
                dist = sqrt(SQR(ax[1][i] - axt[tisspoints[1][itp]])
                            + SQR(ax[2][i] - ayt[tisspoints[2][itp]]) + SQR(ax[3][i] - azt[tisspoints[3][itp]]));
                if (dist <= req) {
                    if (is2d) gtv = fac / w2d * (2.*log(r2d / req) + 1. - SQR(dist / req));
                    else gtv = fac * (1.5 - 0.5*SQR(dist / req)) / req;
                }
                else {
                    if (is2d) gtv = fac * 2. / w2d * log(r2d / dist);
                    else gtv = fac / dist;
                }
                for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) ptv[itp][isp] += gtv / diff[isp] * qv[i][isp];
            }
        }
        for (isp = 1; isp <= nsp; isp++) qvfac[isp] = 1.;
        for (ktissue = 1; ktissue <= nmaxtissue; ktissue++) {
            //Scale all qv, qvsum and ptv values so that qvsum = qtsum.
            
            for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
                qvfac[isp] = -qtsum[isp] / qvsum[isp];
                if (fabs(qvfac[isp]) > 2.) qvfac[isp] = 1.;  //avoid extreme values
                if (fabs(qvfac[isp]) < 0.5) qvfac[isp] = 1.;  //avoid extreme values
            }
            
            convflagt = 1;
            for (isp = 1; isp <= nsp; isp++) {
                qtsum[isp] = 0;
                errtissue[isp] = 0.;
                dqtsumdg0[isp] = 0.;
                errtissuecount[isp] = 0; //added June 2009
            }
            //----------------------------------------------------/
            
            for (itp = 1; itp <= nnt; itp++) {    //contribution ptt from tissue source strengths qt
                ix = tisspoints[1][itp];
                iy = tisspoints[2][itp];
                iz = tisspoints[3][itp];
                for (isp = 1; isp <= nsp; isp++)    ptt[isp] = 0.;//all solutes
                for (jtp = 1; jtp <= nnt; jtp++) {
                    jx = tisspoints[1][jtp];
                    jy = tisspoints[2][jtp];
                    jz = tisspoints[3][jtp];
                    ixdiff = abs(ix - jx) + 1;
                    iydiff = abs(iy - jy) + 1;
                    izdiff = abs(iz - jz) + 1;
                    for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) ptt[isp] += dtt[ixdiff][iydiff][izdiff] * qt[jtp][isp];
                }
                for (isp = 1; isp <= nsp; isp++) ptpt[isp] = pt[itp][isp];
                for (isp = 1; isp <= nsp; isp++) {
                    if (diffsolute[isp]) {        //diffusible solutes
                        pt[itp][isp] = (1. - lam)*pt[itp][isp]
                        + lam * (ptv[itp][isp] * qvfac[isp] + g0[isp] + ptt[isp] / diff[isp]);//underrelaxation
                    }
                    else {    //non-diffusible - use Newton method to solve for pt.
                        tissrate(nsp, ptpt, mtiss, mptiss);    //update tissrates - April 2015
                        if (mptiss[isp] == 0.) printf("*** Error: mptiss[%i] = 0 at tissue point %i\n", isp, itp);
                        else pt[itp][isp] -= mtiss[isp] / mptiss[isp];
                    }
                    ptpt[isp] = pt[itp][isp];
                }
                tissrate(nsp, ptpt, mtiss, mptiss);    //update tissrates
                for (isp = 1; isp <= nsp; isp++) {  //replace qt with value based on updated pt - all solutes
                    dif = mtiss[isp] * vol - qt[itp][isp];
                    qt[itp][isp] = mtiss[isp] * vol;
                    qtsum[isp] += qt[itp][isp];
                    if (diffsolute[isp] == 1) dqtsumdg0[isp] += mptiss[isp] * vol;
                    if (fabs(dif) > errtissue[isp]) {
                        errtissue[isp] = fabs(dif);
                        imaxerrtissue[isp] = itp;
                    }
                    if (fabs(dif) > epstissue[isp]) errtissuecount[isp]++;
                }
            }
            
            for (isp = 1; isp <= nsp; isp++) {
                if (greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp] * qvfac[isp]);//May 2010
                if (greensverbose) printf("Solute %i: ktissue = %i, errtissue_q = %f, imaxerr = %i, g0 = %f\n",
                                          isp, ktissue, errtissue[isp], imaxerrtissue[isp], g0[isp]);
                if (errtissuecount[isp] > 0) convflagt = 0;
            }
            if (kmain > 1 && convflagt) goto tissueconv;  //force full number of iterations when kmain = 1.
        }
        for (isp = 1; isp <= nsp; isp++) if (errtissuecount[isp] > 0)
            if (!greensverbose) printf("*** Warning: solute %i, %i tissue source strengths not converged\n", isp, errtissuecount[isp]);
    tissueconv:;
        //Print log file.  April 2010
        ofp1 = fopen("Current\\GreensLog.txt", "a");
        kvessel = IMIN(kvessel, nmaxvessel);
        ktissue = IMIN(ktissue, nmaxtissue);
        fprintf(ofp1, "\n----- kmain = %i, kvessel = %i, ktissue = %i -----\n", kmain, kvessel, ktissue);
        for (isp = 1; isp <= nsp; isp++) {
            if (diffsolute[isp] == 1) fprintf(ofp1, "Solute %i: qtsum = %f, qvsum = %f, g0 = %f\n",
                                              isp, qtsum[isp], qvsum[isp] * qvfac[isp], g0[isp]);
            if (permsolute[isp] == 1) fprintf(ofp1, "Solute %i: errvessel_q = %f, imaxerr = %i\n",
                                              isp, errvessel[isp], segname[imaxerrvessel[isp]]);
            if (diffsolute[isp] == 1) fprintf(ofp1, "Solute %i: errtissue_q = %f, imaxerr = %i\n",
                                              isp, errtissue[isp], imaxerrtissue[isp]);
        }
        fclose(ofp1);
        //--------------------  end of tissue loop --------------------
        //Update g0.  If permsolute[isp] != 1, always use method 2.
        //Method 2 is based on derivative wrt g0 - new version September 2009 - automatic estimation of g0fac
        for (isp = 1; isp <= nsp; isp++) g0facnew[isp] = 0.;
        for (itp = 1; itp <= nnt; itp++) {
            for (isp = 1; isp <= nsp; isp++) ptpt[isp] = pt[itp][isp];
            tissrate(nsp, ptpt, mtiss, mptiss);
            ix = tisspoints[1][itp];
            iy = tisspoints[2][itp];
            iz = tisspoints[3][itp];
            for (jtp = 1; jtp <= nnt; jtp++) {
                jx = tisspoints[1][jtp];
                jy = tisspoints[2][jtp];
                jz = tisspoints[3][jtp];
                ixdiff = abs(ix - jx) + 1;
                iydiff = abs(iy - jy) + 1;
                izdiff = abs(iz - jz) + 1;
                for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) g0facnew[isp] += dtt[ixdiff][iydiff][izdiff] / diff[isp] * mptiss[isp] * vol;
            }
        }
        for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1 && (g0method == 2 || permsolute[isp] == 0)) {
            g0facnew[isp] = 1. / (1. - g0facnew[isp] / nnt);
            dqsumdg0 = FMIN(dqvsumdg0[isp], 0.) + FMIN(dqtsumdg0[isp], 0.)*g0facnew[isp];
            if (fabs(dqsumdg0) > 1.e-6) {
                dif = (qvsum[isp] + qtsum[isp]) / dqsumdg0 * g0fac[isp];//This g0fac should normally be 1.0.
                g0[isp] -= dif;
            }
        }
        //Convergence based on changes in pv, pt and g0.  Express relative to eps[isp].
        convflag = 1;
        if (greensverbose) printf("\n");
        for (isp = 1; isp <= nsp; isp++) {
            err = 0.;
            imaxerr = 0;
            if (permsolute[isp] == 1) for (i = 1; i <= nnv; i++) {
                dif = fabs(pv[i][isp] - pvprev[i][isp]) / eps[isp];
                if (dif > err) {
                    imaxerr = mainseg[i];
                    err = dif;
                }
            }
            errvessel[isp] = err;
            imaxerrvessel[isp] = imaxerr;
            err = 0.;
            imaxerr = 0;
            for (itp = 1; itp <= nnt; itp++) {
                dif = fabs(pt[itp][isp] - ptprev[itp][isp]) / eps[isp];
                if (dif > err) {
                    imaxerr = itp;
                    err = dif;
                }
            }
            errtissue[isp] = err;
            imaxerrtissue[isp] = imaxerr;
            if (errvessel[isp] > err) {
                imaxerr = imaxerrvessel[isp];
                err = errvessel[isp];
            }
            else imaxerr = -imaxerr;
            dif = fabs(g0[isp] - g0old[isp]) / eps[isp];
            if (dif > err) {
                imaxerr = 0;
                err = dif;
            }
            if (greensverbose) printf("Solute %i: err = %f, imaxerr = %i (- for tissue point)\n", isp, err, imaxerr);
            if (greensverbose && imaxerr > 0) if (lowflow[imaxerr]) printf("Solute %i: max error is at a low-flow segment\n", isp);
            if (err > 1.0) convflag = 0; //changed from 1.0 - Baran
        }
        //Print log file - April 2010
        ofp1 = fopen("Current\\GreensLog.txt", "a");
        for (isp = 1; isp <= nsp; isp++) {
            if (permsolute[isp] == 1) fprintf(ofp1, "Solute %i: errvessel_p = %f, imaxerr = %i\n",
                                              isp, errvessel[isp], segname[imaxerrvessel[isp]]);
            fprintf(ofp1, "Solute %i: errtissue_p = %f, imaxerr = %i\n",
                    isp, errtissue[isp], imaxerrtissue[isp]);
        }
        fclose(ofp1);
        tfinish1 = clock();
        duration = (float)(tfinish1 - tstart1) / CLOCKS_PER_SEC;
        if (greensverbose) printf("\nkmain = %i, %2.3f seconds for step\n", kmain, duration);
        if (convflag && convflagv && convflagt) goto mainconv;
    }
    printf("\n convflag = %i, convflagv = %i, convflagt = %i", convflag, convflagv, convflagt);
    printf("\n*** Warning: tissue or vessel solute levels not converged");

    mainconv:;
    //--------------------  end of main loop --------------------
    tfinish = clock();
    
    duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
    printf("\n%i iterations, %2.1f seconds for main loop\n", kmain, duration);
    ofp1 = fopen("Current\\GreensLog.txt", "a");
    fprintf(ofp1, "\n%i iterations, %2.1f seconds for main loop\n", kmain, duration);
    fclose(ofp1);
    //Scale all qv values so that qvsum = qtsum.  April 2010.
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
        qvsum[isp] *= qvfac[isp];
        for (i = 1; i <= nnv; i++) qv[i][isp] *= qvfac[isp];
    }
    //general output file
    strcpy(fname, "Current\\GreensRes");
    strcat(fname, numstr);
    strcat(fname, ".out");
    ofp = fopen(fname, "w");
    fprintf(ofp, "%i %i %i %i %i %i\n", nnv, nseg, mxx, myy, mzz, nnt);
    fprintf(ofp, "Scaling factor for flows q0fac = %f\n", q0fac);
    for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "g0[%i] = %f\n", isp, g0[isp]);
    fprintf(ofp, "\n");
    //extravascular solute levels
    for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
        fprintf(ofp, "\nSolute %i\n", isp);
        fprintf(ofp, "Segment");
        if (permsolute[isp] == 1) fprintf(ofp, "Efflux Pvessel Ptissue Cvessel");
        fprintf(ofp, "\n");
        for (i = 1; i <= nnv; i++) {
            if(DEBG==1) {
                printf("isp = %d, i = %d, pv = %5.40f\n", isp, i, pv[i][isp]);
            }
            pev[i][isp] = pvt[i][isp];
            if (g0method == 1 && permsolute[isp] == 1) pev[i][isp] += g0[isp];
            if (permsolute[isp] == 1) for (j = 1; j <= nnv; j++) pev[i][isp] += gvv[i][j] * qv[j][isp] / diff[isp];
        }
        for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
            qvseg[iseg][isp] = 0.;
            pevseg[iseg][isp] = 0.;
            pvseg[iseg][isp] = 0.;
        }
        for (i = 1; i <= nnv; i++) {
            iseg = mainseg[i];
            if (permsolute[isp] == 1) fprintf(ofp, "%4i %4i %10.10f %10.40f %10.10f %10.10f\n",
                                              i, iseg, qv[i][isp], pv[i][isp], pev[i][isp], cv[i][isp]);
            qvseg[iseg][isp] += qv[i][isp];
            pevseg[iseg][isp] += pev[i][isp] / nspoint[iseg];
            pvseg[iseg][isp] += pv[i][isp] / nspoint[iseg];
        }
        fprintf(ofp, "Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
    }
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) {
        fprintf(ofp, "Solute %i: segment length pvseg pevseg qvseg gamma\n", isp);
        for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5)
            fprintf(ofp, "%4i %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                    segname[iseg], lseg[iseg], pvseg[iseg][isp], pevseg[iseg][isp], qvseg[iseg][isp], gamma1[iseg][isp]);
    }
    fclose(ofp);
    
    //Vessel levels for all vessel points
    strcpy(fname, "Current\\VesselLevels");
    strcat(fname, numstr);
    strcat(fname, ".out");
    ofp = fopen(fname, "w");
    fprintf(ofp, "Vessel levels\n");
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) {
        pmax[isp] = -1.e8;
        pmean[isp] = 0.;
        pmin[isp] = 1.e8;
        for (i = 1; i <= nnv; i++) {
            pmean[isp] += pv[i][isp];
            pmax[isp] = FMAX(pv[i][isp], pmax[isp]);
            pmin[isp] = FMIN(pv[i][isp], pmin[isp]);
        }
        pmean[isp] = pmean[isp] / nnv;
        fprintf(ofp, "   Solute %i  ", isp);
    }
    fprintf(ofp, "\npmean\n");
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp])    fprintf(ofp, "%12f ", pmean[isp]);
    fprintf(ofp, "\npmin\n");
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pmin[isp]);
    fprintf(ofp, "\npmax\n");
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pmax[isp]);
    fprintf(ofp, "\nvalues\n");
    for (i = 1; i <= nnv; i++) {
        for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pv[i][isp]);
        fprintf(ofp, "\n");
    }
    fclose(ofp);
    //Vessel source strengths for all vessel points
    strcpy(fname, "Current\\VesselSources");
    strcat(fname, numstr);
    strcat(fname, ".out");
    ofp = fopen(fname, "w");
    fprintf(ofp, "Vessel sources\n");
    for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "   Solute %i  ", isp);
    fprintf(ofp, "\n");
    for (i = 1; i <= nnv; i++) {
        for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", qv[i][isp]);
        fprintf(ofp, "\n");
    }
    fclose(ofp);
    
    //Tissue levels for all tissue points
    strcpy(fname, "Current\\TissueLevels");
    strcat(fname, numstr);
    strcat(fname, ".txt");
    ofp = fopen(fname, "w");
    fprintf(ofp, "Tissue levels\n");
    for (isp = 1; isp <= nsp; isp++) {
        pmax[isp] = -1.e8;
        pmean[isp] = 0.;
        pmin[isp] = 1.e8;
        for (itp = 1; itp <= nnt; itp++) {
            pmean[isp] += pt[itp][isp];
            pmax[isp] = FMAX(pt[itp][isp], pmax[isp]);
            pmin[isp] = FMIN(pt[itp][isp], pmin[isp]);
        }
        pmean[isp] = pmean[isp] / nnt;
        fprintf(ofp, "   Solute %i  ", isp);
    }
    fprintf(ofp, "\npmean\n");
    for (isp = 1; isp <= nsp; isp++)    fprintf(ofp, "%12f ", pmean[isp]);
    fprintf(ofp, "\npmin\n");
    for (isp = 1; isp <= nsp; isp++)    fprintf(ofp, "%12f ", pmin[isp]);
    fprintf(ofp, "\npmax\n");
    for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pmax[isp]);
    fprintf(ofp, "\nvalues\n");
    for (itp = 1; itp <= nnt; itp++) {
        for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pt[itp][isp]);
        fprintf(ofp, "\n");
    }
    fclose(ofp);
    
    //Tissue levels for all tissue points 
    strcpy(fname, "TissueLevels");
    strcat(fname, numstr);
    strcat(fname, ".txt");
    ofp = fopen(fname, "w");
    for (isp = 1; isp <= nsp; isp++) {
        pmax[isp] = -1.e8;
        pmean[isp] = 0.;
        pmin[isp] = 1.e8;
        for (itp = 1; itp <= nnt; itp++) {
            pmean[isp] += pt[itp][isp];
            pmax[isp] = FMAX(pt[itp][isp], pmax[isp]);
            pmin[isp] = FMIN(pt[itp][isp], pmin[isp]);
        }
        pmean[isp] = pmean[isp] / nnt;
    }
    fprintf(ofp, "%12f\n ", pmean[1]);
    fprintf(ofp, "%12f\n ", pmin[1]);
    fprintf(ofp, "%12f\n ", pmax[1]);
    for (itp = 1; itp <= nnt; itp++) {
        fprintf(ofp, "%12f\n ", pt[itp][1]);

    }
    fclose(ofp);
    
    //Tissue source strengths for all tissue points
    strcpy(fname, "Current\\TissueSources");
    strcat(fname, numstr);
    strcat(fname, ".out");
    ofp = fopen(fname, "w");
    fprintf(ofp, "Tissue sources\n");
    for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "   Solute %i  ", isp);
    fprintf(ofp, "\n");
    for (itp = 1; itp <= nnt; itp++) {
        for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", qt[itp][isp]);
        fprintf(ofp, "\n");
    }
    fclose(ofp);
    
    fprintf(ofp, "%4i  %4i  ", imain, kmain);
    for (j = 1; j <= nvaryparams; j++) fprintf(ofp, "%12f ", paramvalue[imain][j]);
    for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pmean[isp]);
    
    if (npostgreensparams) postgreens();
    
    if (npostgreensout) for (j = 1; j <= npostgreensout; j++) fprintf(ofp, "%12f ", postgreensout[j]);
    fprintf(ofp, "\n");
    
    for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = pvseg[iseg][1];
    for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];
    
    strcpy(fname, "Current\\NetNodesOxygen");
    strcat(fname, numstr);
    strcat(fname, ".ps");
    picturenetwork(nodvar, segvar, fname);
    
    cmgui(segvar);
    
    strcpy(fname, "Current\\Contour");
    strcat(fname, numstr);
    strcat(fname, ".ps");
    contour(fname);
    
    strcpy(fname, "Current\\Histogram");
    strcat(fname, numstr);
    strcat(fname, ".out");
    histogram(fname);
    
    fclose(ofp);
    
    printf("Success! \n");
    return 0;
}

    


// the two functions WriteExnodeHeader and WriteExelemHeader are needed for creating a postscript file inside main();
void WriteExnodeHeader(FILE *exnode) // Write initial section of .exnode file
{
    //    fprintf(exnode, "Region: /vessels\n");
    fprintf(exnode, "Group name: vessels\n");
    fprintf(exnode, " #Fields=3\n");
    fprintf(exnode, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  x.  Value index=1, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  y.  Value index=2, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  z.  Value index=3, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exnode, "  1.  Value index=4, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  1.  Value index=5, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  2.  Value index=6, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  3.  Value index=7, #Derivatives=0, #Versions=1\n");
}

void WriteExelemHeader(FILE *exelem)  // Write initial section of .exelem file
{
    //   fprintf(exelem, "Region: /vessels\n");
    fprintf(exelem, "Group name: vessels\n");
    fprintf(exelem, " Shape.  Dimension=1\n");
    fprintf(exelem, " #Scale factor sets= 1\n");
    fprintf(exelem, "  l.Lagrange, #Scale factors= 2\n");
    fprintf(exelem, " #Nodes= 2\n #Fields=3\n");
    fprintf(exelem, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");
    fprintf(exelem, "   2.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");
    fprintf(exelem, "   3.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");
}
void convect(int isp)
{
    extern int nnodbc, nseg, nnodfl, nodsegm, nsp, nnv, nsegfl;
    extern int *bcnod, *nodtyp, *nodrank, *nodout, *segtyp, *permsolute, *oxygen, *nspoint, *istart;
    extern int *segname, *nodname, **nodseg, *mainseg;
    extern float *bifpar, *hd, *qq, *q, *bchd, *diam, **pv, **bcp, *ds, *segc, **cv, **qv, flowfac;
    extern float **al;
    extern float *solutefac;    //April 2015
    
    int i, j, k, ii, jj, inod, iseg, jseg, in, isegk, nodt, nin, nout, ineg, ihigh;
    float fluxsumin, pb, pp;
    float sumin, sumout, hdsumin, hdsumout;
    
    int *isegkk; //added June 2013 to check for errors in segment sequence
    isegkk = ivector(1, nseg);
    for (iseg = 1; iseg <= nseg; iseg++) isegkk[iseg] = 0;
    
    isegk = 0;    //number of segments processed
    for (i = 1; i <= nnv; i++) for (j = 1; j <= nnv; j++) al[i][j] = 0.;
    for (iseg = 1; iseg <= nseg; iseg++) {segc[iseg] = 0.;    //added June 2013
        if(DEBG==1) printf("%f %f\n", segc[iseg], qq[iseg]);
    }

    //set convective fluxes in segments connected to inflow boundary nodes using boundry values specified in simpleSoluteParams.
    //segc is the convective flux of solute
    for (j = 1; j <= nnodbc; j++) {//for each boundry node
        inod = bcnod[j];
        if(DEBG==1) printf("%d\n", nodout[inod]);
        if (nodout[inod] == 1) {//only for inflow nodes
            iseg = nodseg[1][inod];
            //values modified according to VaryParams.dat, April 2015
            if(DEBG==1) printf("%f %f %f %f\n", bcp[j][isp], solutefac[isp], hd[iseg], qq[iseg]);
            if (oxygen[isp] == 1) segc[iseg] = bloodconc(bcp[j][isp] * solutefac[isp], hd[iseg])*qq[iseg] * flowfac;
            else segc[iseg] = bcp[j][isp] * solutefac[isp] * qq[iseg] * flowfac;
            isegkk[iseg] = 1;
        }
         if(DEBG==1) printf("%f %f %f %f\n", segc[iseg], qq[iseg], bcp[j][isp], bloodconc(bcp[j][isp] * solutefac[isp], hd[iseg]));
    }
    ineg = 0;
    ihigh = 0;
    for (in = 1; in <= nnodfl; in++) {    //scan nodes in downstream order
        inod = nodrank[in];
        nodt = nodtyp[inod];
        nout = nodout[inod];
        nin = nodt - nout;
        if (nodt > 1) {    //don't do for network boundary nodes
            sumin = 0.;
            hdsumin = 0.;
            fluxsumin = 0.;
            for (ii = nout + 1; ii <= nodt; ii++) { //inflows
                iseg = nodseg[ii][inod];
                if (isegkk[iseg] == 0) printf("*** Error: wrong segment sequence in convect, segment %i ***\n", iseg);
                sumin += qq[iseg] * flowfac;//Nano litres per minute of inflow (flowfac converts the units to ?)
                hdsumin += qq[iseg] * flowfac*hd[iseg];//Nano litres per minute inflow of RBCs
                fluxsumin += segc[iseg];//Convective flux of solute in
            }
            //calculate solute level going into node
            if(DEBG==1){
                printf("%f %f \n", pb, pp);
            }
            if (oxygen[isp] == 1) blood(fluxsumin / sumin, hdsumin / sumin, &pb, &pp);//gives the partial pressure of oxygen in the blood given the concentration (fluxsumin/sumin), hematocrit (hdsumin/sumin)
            else pb = fluxsumin / sumin;
            //assign solute levels going out of node
            sumout = 0.;
            hdsumout = 0.;
            for (ii = 1; ii <= nout; ii++) {        //outflows
                iseg = nodseg[ii][inod];
                isegkk[iseg] = 1;
                sumout += qq[iseg] * flowfac;    //To check conservation of flow
                hdsumout += qq[iseg] * flowfac*hd[iseg];//To check conservation of hematocrit
                if (oxygen[isp] == 1) segc[iseg] = bloodconc(pb, hd[iseg])*qq[iseg] * flowfac;
                else segc[iseg] = pb * qq[iseg] * flowfac;
                if (q[iseg] >= 0.) i = istart[iseg];
                else i = istart[iseg] + nspoint[iseg] - 1;
                for (jj = nout + 1; jj <= nodt; jj++) {    //inflows
                    jseg = nodseg[jj][inod];
                    if (q[jseg] >= 0.) j = istart[jseg] + nspoint[jseg] - 1;
                    else j = istart[jseg];
                    al[i][j] = qq[iseg] * flowfac / sumin;
                    // al[i][j] for consecutive segments i and j is the proportion of the flow from segment j, that is directed to segment i.
                    //calculate alpha values across node by taking the ratio of the flow down each output segment i with the total flow into the nod under study
                    if (oxygen[isp] == 1 && nout > 1)    //if nout=1, hd[iseg] = hdsumin/sumin
                        al[i][j] *= bloodconcp(pb, hd[iseg]) / bloodconcp(pb, hdsumin / sumin);
                    for (k = 1; k <= nnv; k++) al[i][k] += al[i][j] * al[j][k];//calculate other alpha values
                }
            }
            if (sumin + sumout != 0.) if (fabs(sumin - sumout) / (sumin + sumout) > 0.01)
                printf("*** Error: Flow conservation violation at node %i\n", nodname[inod]);
            if (hdsumin + hdsumout != 0.) if (fabs(hdsumin - hdsumout) / (hdsumin + hdsumout) > 0.01)
                printf("*** Error: Hematocrit conservation violation at node %i\n", nodname[inod]);
        }
        //subsegments of outflow segments - convective fluxes and alpha matrix, including network boundary nodes
        for (ii = 1; ii <= nout; ii++) {        //outflows
            iseg = nodseg[ii][inod];
            for (jj = 0; jj < nspoint[iseg]; jj++) {    //convective fluxes
                if (q[iseg] >= 0.) i = istart[iseg] + jj;
                else i = istart[iseg] + nspoint[iseg] - jj - 1;
                segc[iseg] -= qv[i][isp] / 2.;
                cv[i][isp] = segc[iseg] / qq[iseg] / flowfac;
                if(DEBG==1) {
                    printf("nsp = %d, isp = %d, qv = %10.20f, segc = %10.20f, DCv = %10.30f, Blood Oxygen Conc = %30.30f cm^3 O2/cm^3 blood, Cb = %f cm^3 O2/cm^3 blood\n" , nsp, isp, qv[i][isp] , segc[iseg],qv[i][isp]/2/qq[iseg]/flowfac, cv[i][isp], cs);
                }
                segc[iseg] -= qv[i][isp] / 2.;
                if(DEBG==1) printf("%f %f\n", segc[iseg], qq[iseg]);
               
            }
            for (jj = 1; jj < nspoint[iseg]; jj++) {    //alpha matrix
                if (q[iseg] >= 0.) {
                    j = istart[iseg] + jj - 1;
                    i = istart[iseg] + jj;
                }
                else {
                    j = istart[iseg] + nspoint[iseg] - jj;
                    i = istart[iseg] + nspoint[iseg] - jj - 1;
                }
                al[i][j] = 1.;
                for (k = 1; k <= nnv; k++) al[i][k] += al[i][j] * al[j][k];
            }
        }
        isegk += nout;
    }
    if (isegk != nsegfl) printf("*** Error in convect, %i of %i segments processed\n", isegk, nseg);
    for (i = 1; i <= nnv; i++) al[i][i] = 0.5;
    free_ivector(isegkk, 1, nseg);
}
void blood(float c, float h, float *p, float *pp)
{
    
    extern float fn, alphab, p50, cs, cext, hext;
    extern float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
    float pl, ph;
    float clow, chigh, pphigh;
    
    if (h < 1.e-6) {    //changed from 1.e-3, January 2012
        *p = c / alphab;
        *pp = alphab;
        printf("h < 1.e-6\n");
        return;
    }
    //changed for better behavior in severe hypoxia. 2 Oct 08
    if (c < 0.) {//Only the term that discribes how O2 is dissolved in the blood is included in this case.
        *p = c / alphab;
        *pp = alphab;
        printf("c < 0. \n");
        return;
    }
    clow = clowfac * h + alphab * plow;
    if (c < clow) {
        *p = c * plow / clow;
        *pp = clow / plow; //A linear approximation of the PO2 vs C_O2 relationship, clow/plow is rise over run.
        printf("c < clow\n");
        return;
    }
    chigh = chighfac * h + alphab * phigh;
    if (c < chigh) {// It looks like alphab*p (the amount of oxygen dissolved in the blood) is ignored when the concentration is higher than clow. -Baran
        if (c / h / cs < 1.) {
            ph = pow((c / h / cs) / (1.0 - c / h / cs), 1. / fn)*p50;
            pl = 0.;    //0. here to be sure to bracket the root.  June 2009
            printf("c/h/cs < 1.\n");
        }
        else {
            ph = phigh;
            pl = plow;//changed August 2010
            printf("c/h/cs >= 1.\n");
        }
        hext = h;
        cext = c;
        if(DEBG==1) { printf(" ph = %5.10f, pl = %5.10f\n", ph, pl);}
            
        *p = rtflsp(func, pl, ph, 0.001f);//lower tolerance, August 2010
        //if false position algorithm fails, use bisection!  June 2009.
        if (*p < 0.) {*p = rtbis(func, pl, ph, 0.001f);
            printf("rtbis is used to find roots because rtflsp failed");
        }
        *pp = cs * h*fn / p50 * pow(*p / p50, (fn - 1)) / SQR(1. + pow(*p / p50, fn)) + alphab;
        printf("c < chigh\n");
        return;
    }
    pphigh = pphighfac * h + alphab;
    *p = phigh + (c - chigh) / pphigh;
    *pp = pphigh;
    printf("c > chigh\n");
    return;
}
//Using the false position method, find the root of a iftion func known to lie between x1 and x2.
//The root, returned as rtflsp, is refined until its accuracy is +- xacc.

float rtflsp(float(*func)(float), float x1, float x2, float xacc)
{
    int j;
    int maxit = 30;
    float fl, fh, xl, xh, swap, dx, del, f, rtf;
    fl = (*func)(x1);
    fh = (*func)(x2);
    if (fl*fh > 0.0) {
        printf("\n*** Warning: Root must be bracketed in rtflsp %f %f\n", x1, x2);
        return -1.0;
    }
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xl = x2;
        xh = x1;
        swap = fl;
        fl = fh;
        fh = swap;
    }
    dx = xh - xl;
    for (j = 1; j <= maxit; j++) {
        rtf = xl + dx * fl / (fl - fh);
        f = (*func)(rtf);
        if (f < 0) {
            del = xl - rtf;
            xl = rtf;
            fl = f;
        }
        else {
            del = xh - rtf;
            xh = rtf;
            fh = f;
        }
        dx = xh - xl;
        if (fabs(del) < xacc || f == 0.0) return rtf;
    }
    //    printf("\n*** Warning: Maximum number of iterations exceeded in rtflsp %f %f\n", x1,x2);
    return -1.0;
}

float rtbis(float(*func)(float), float x1, float x2, float xacc)
{
    int j;
    int maxit = 30;
    float dx, f, fmid, xmid, rtb;
    f = (*func)(x1);
    fmid = (*func)(x2);
    if (f*fmid >= 0.0) {
        printf("*** Error: Root must be bracketed in rtbis %f %f\n", x1, x2);
        return 0.0;
    }
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    for (j = 1; j <= maxit; j++) {
        fmid = (*func)(xmid = rtb + (dx *= 0.5));
        if (fmid <= 0.) rtb = xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    printf("*** Error: Maximum number of iterations exceeded in rtbis %f %f\n", x1, x2);
    return 0.0;
}

float func(float p)
{
    extern float cext, hext;
    float y;
    y = bloodconc(p, hext) - cext;
    return y;
}

float bloodconc(float p, float h)
{
    extern float fn, p50, alphab, cs;
    extern float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
    float cv;
    if (p < 0.) cv = alphab * p;  //This gives more consistent behavior when p<0
    else if (p < plow) cv = clowfac * h*p / plow + alphab * p;//The Hill equation is approximated by a linear equation, where clowfac*h/plow aproximates the derivative on the Hill relationship at low oxygen pressures (0.1*p50)
    else if (p < phigh) cv = cs * h*(1.0 - 1.0 / (1.0 + pow((p / p50), fn))) + alphab * p;//Hill equation + dissolved fraction.
    else cv = (chighfac + (p - phigh)*pphighfac)*h + alphab * p;//The Hill equation is approximated by it's linear part when oxygen pressures are high (5*p50)
    return cv;
}
float bloodconcp(float p, float h)
{
    extern float fn, p50, alphab, cs;
    extern float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
    float cvp;
    if (p < 0.) cvp = alphab;  //This gives more consistent behavior when p<0
    else if (p < plow) cvp = clowfac * h / plow + alphab;
    else if (p < phigh) cvp = cs * h*fn / p50 * pow(p / p50, (fn - 1)) / SQR(1. + pow(p / p50, fn)) + alphab;
    else cvp = pphighfac * h + alphab;
    return cvp;
}
void gaussj(double **a, int n, double **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int i, icol, irow, j, k, l, ll;
    double big, dum, pivinv;
    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for (j = 1; j <= n; j++) ipiv[j] = 0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++) if (ipiv[j] != 1) for (k = 1; k <= n; k++) if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
                big = fabs(a[j][k]);
                irow = j;
                icol = k;
            }
        }
        else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
                for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
                    }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
        for (ll = 1; ll <= n; ll++) if (ll != icol) {
            dum = a[ll][icol];
            a[ll][icol] = 0.0;
            for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
            for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
        }
    }
    for (l = n; l >= 1; l--) if (indxr[l] != indxc[l]) for (k = 1; k <= n; k++) SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
}
void ludcmp(double **a, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;
    vv = dvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++) if ((temp = fabs(a[i][j])) > big) big = temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n) {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv, 1, n);
}

void lubksb(double **a, int n, int *indx, double *b)
{
    int i, ii = 0, ip, j;
    double sum;
    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii) for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}
double bicgstab(double **a, double *b, double *x, int n, double eps, int itmax)
{
    double lu, lunew, beta, delta, er, gamma1, t1, t2, err;
    double *r, *rs, *v, *s, *t, *p;
    int i, j, kk;
    r = dvector(1, n);
    rs = dvector(1, n);
    v = dvector(1, n);
    s = dvector(1, n);
    t = dvector(1, n);
    p = dvector(1, n);
    lu = 0.;
    for (i = 1; i <= n; i++) {
        r[i] = 0.;
        for (j = 1; j <= n; j++)    r[i] += a[i][j] * x[j];
        r[i] -= b[i];
        p[i] = r[i];
        rs[i] = 1.;
        lu += r[i] * rs[i];
    }
    kk = 1;
    do
    {
        t1 = 0.;
        for (i = 1; i <= n; i++) {
            v[i] = 0.;
            for (j = 1; j <= n; j++) v[i] += a[i][j] * p[j];
            t1 += v[i] * rs[i];
        }
        delta = -lu / t1;
        for (i = 1; i <= n; i++) s[i] = r[i] + delta * v[i];
        for (i = 1; i <= n; i++) {
            t[i] = 0.;
            for (j = 1; j <= n; j++) t[i] += a[i][j] * s[j];
        }
        t1 = 0.;
        t2 = 0.;
        for (i = 1; i <= n; i++) {
            t1 += s[i] * t[i];
            t2 += t[i] * t[i];
        }
        gamma1 = -t1 / t2;
        err = 0.;
        lunew = 0.;
        for (i = 1; i <= n; i++) {
            r[i] = s[i] + gamma1 * t[i];
            er = delta * p[i] + gamma1 * s[i];
            x[i] += er;
            if (fabs(er) > err) err = fabs(er);
            lunew += r[i] * rs[i];
        }
        beta = lunew * delta / (lu*gamma1);
        lu = lunew;
        for (i = 1; i <= n; i++) p[i] = r[i] + beta * (p[i] + gamma1 * v[i]);
        kk += 1;
    }
    while (kk < itmax && err > eps);
    free_dvector(r, 1, n);
    free_dvector(rs, 1, n);
    free_dvector(v, 1, n);
    free_dvector(s, 1, n);
    free_dvector(t, 1, n);
    free_dvector(p, 1, n);
    if (err > eps) printf("*** Warning: linear solution using BICGSTB not converged, err = %e\n", err);
    return err;
}

float *eval(int slsegdiv, float req, float *x)
{
    extern int mxx, myy, mzz, nnt, nnv, nseg, nnod, nsp;
    extern int *mainseg, **tisspoints, *permsolute, *diffsolute, ***nbou;
    extern int is2d; //needed for 2d version
    extern float fac;
    extern float *axt, *ayt, *azt, *ds, *diff, *g0, *y, *p;
    extern float **qt, **start, **scos, **qv, **ax, **pt;
    extern float w2d, r2d; //needed for 2d version
    
    float dist2, gtt, gtv, lamx, lamy, lamz, r2d2 = SQR(r2d), req2 = SQR(req);
    int i, j, k, ii, jj, kk, iseg, itp, isp;
    
    //initialize to g0
    for (isp = 1; isp <= nsp; isp++) p[isp] = g0[isp];
    //add contributions from tissue sources
    for (itp = 1; itp <= nnt; itp++) {
        dist2 = SQR(x[1] - axt[tisspoints[1][itp]])
        + SQR(x[2] - ayt[tisspoints[2][itp]])
        + SQR(x[3] - azt[tisspoints[3][itp]]);
        if (dist2 <= req2) {
            if (is2d) gtt = fac / w2d * (log(r2d2 / req2) + 1. - dist2 / req2);
            else gtt = fac * (1.5 - 0.5*dist2 / req2) / req;
        }
        else {
            if (is2d) gtt = fac / w2d * log(r2d2 / dist2);
            else gtt = fac / sqrt(dist2);
        }
        for (isp = 1; isp <= nsp; isp++)    if (diffsolute[isp]) p[isp] += gtt / diff[isp] * qt[itp][isp];
    }
    //add contributions from vessel sources.  Subdivide subsegments.
    //Note that vessel point is at midpoint of subsegment
    for (i = 1; i <= nnv; i++) {
        iseg = mainseg[i];
        for (k = 1; k <= slsegdiv; k++) {
            for (j = 1; j <= 3; j++)    y[j] = ax[j][i] + scos[j][iseg] * ds[iseg] * (-0.5 + (k - 0.5) / slsegdiv);
            dist2 = SQR(x[1] - y[1]) + SQR(x[2] - y[2]) + SQR(x[3] - y[3]);
            if (dist2 <= req2) {
                if (is2d) gtv = fac / w2d * (log(r2d2 / req2) + 1. - dist2 / req2);
                else gtv = fac * (1.5 - 0.5*dist2 / req2) / req;
            }
            else {
                if (is2d) gtv = fac / w2d * log(r2d2 / dist2);
                else gtv = fac / sqrt(dist2);
            }
            for (isp = 1; isp <= nsp; isp++)    if (permsolute[isp]) p[isp] += gtv / diff[isp] * qv[i][isp] / slsegdiv;
        }
    }
    //for non-diffusible solute, calculate by interpolating values at tissue points.  May 2010.
    for (isp = 1; isp <= nsp; isp++)    if (diffsolute[isp] == 0) {
        i = 0;
        while (x[1] > axt[i + 1] && i < mxx) i++;
        if (i == 0) lamx = 1.;
        else if (i == mxx) lamx = 0.;
        else lamx = (x[1] - axt[i]) / (axt[i + 1] - axt[i]);
        j = 0;
        while (x[2] > ayt[j + 1] && j < myy) j++;
        if (j == 0) lamy = 1.;
        else if (j == myy) lamy = 0.;
        else lamy = (x[2] - ayt[j]) / (ayt[j + 1] - ayt[j]);
        k = 0;
        while (x[3] > azt[k + 1] && k < mzz) k++;
        if (k == 0) lamz = 1.;
        else if (k == mzz) lamz = 0.;
        else lamz = (x[3] - azt[k]) / (azt[k + 1] - azt[k]);
        p[isp] = 0.;
        for (ii = 0; ii <= 1; ii++) for (jj = 0; jj <= 1; jj++) for (kk = 0; kk <= 1; kk++)
            if (i + ii >= 1 && i + ii <= mxx && j + jj >= 1 && j + jj <= myy && k + kk >= 1 && k + kk <= mzz) {
                itp = nbou[i + ii][j + jj][k + kk];
                if (itp != 0) p[isp] += (1. - lamx + ii * (2.*lamx - 1))*(1. - lamy + jj * (2.*lamy - 1))
                    *(1. - lamz + kk * (2.*lamz - 1))*pt[itp][isp];
            }
    }
    return p;
}
void tissrate(int nsp, float *p, float *mtiss, float *mptiss)
{
    extern int *oxygen;
    extern float **tissparam;
    
    int isp;
    
#include "tissrate.cpp.txt"
}
//-----------------------------------------------------------------/
// postgreens.cpp - analyzes results from greens
// Uses parameters from PostGreensParams.dat
// Includes problem-specific code from postgreens.cpp.dat
// Example of usage: compute survival fraction of cells from drug concentration
// TWS, May 2015
//-----------------------------------------------------------------/
void postgreens(void)
{
    extern int max, nsp, nnt, npostgreensparams, npostgreensout;
    extern float **pt, *dtmin, *postgreensparams, *postgreensout;
    extern char numstr[6];
    char fname[80];
    FILE *ofp;
    
    int    i, isp, itp;
    
    strcpy(fname, "Current\\PostGreens");
    strcat(fname, numstr);
    strcat(fname, ".out");
    ofp = fopen(fname, "w");
    //-----------------------------------------------------------------/
#include "postgreens.cpp.dat"
    //-----------------------------------------------------------------/
    fclose(ofp);
}
//-----------------------------------------------------------------/
// picturenetwork.cpp - project network on z = 0 plane.  TWS Dec. 07
// Uses parameters from CountourParams.dat
// Labels nodes with nodevar and segments with segvar (must be float).
// Generates a postscript file.
// Version 2.0, May 1, 2010.  Added abbrevations of setlinewidth and stroke, May 09.
// Added visualization of tissue points, April 2010.
// Version 3.0, May 17, 2011.  Plots in z-order, for better results with 3D networks.
 //-----------------------------------------------------------------/
#define _CRT_SECURE_NO_DEPRECATE
void picturenetwork(float *nodvar, float *segvar, const char [])
{
    extern int max, mxx, myy, mzz, nseg, nnod;
    extern int *segtyp, *ista, *iend;
    extern int ***nbou;//added April 2010
    extern float *axt, *ayt;//added April 2010
    extern float *diam, **cnode, *xsl0, *xsl1, *xsl2;;
    int i, j, k, iseg, inod, showpoint, ilevel, nlevel = 100;
    float xmin, xmax, ymin, ymax, xs, ys, picfac, red, green, blue, xz, xzmin, xzmax;
    float diamfac = 1., zcoord, zbottom, ztop, zmin, zmax;
    float **cos;
    
    FILE *ofp;
    
    xmin = 0.;
    xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));
    ymin = 0.;
    ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));
    
    cos = matrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {    //set up matrix of direction cosines
        cos[1][i] = (xsl1[i] - xsl0[i]) / xmax;
        cos[2][i] = (xsl2[i] - xsl0[i]) / ymax;
    }
    cos[3][1] = cos[1][2] * cos[2][3] - cos[1][3] * cos[2][2];
    cos[3][2] = cos[1][3] * cos[2][1] - cos[1][1] * cos[2][3];
    cos[3][3] = cos[1][1] * cos[2][2] - cos[1][2] * cos[2][1];
    //Determine range of z values
    zmin = 1.e6;
    zmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        zcoord = 0.;
        for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
        zmin = FMIN(zmin, zcoord - 1.);
        zmax = FMAX(zmax, zcoord + 1.);
    }
    
    picfac = FMIN(500. / xmax, 700. / ymax);//updated April 2010
    ofp = fopen(fname, "w");
    fprintf(ofp, "%%!PS-Adobe-2.0\n");
    fprintf(ofp, "%%%%Pages: 1\n");
    fprintf(ofp, "%%%%EndComments\n");
    fprintf(ofp, "%%%%Page: 1 1\n");
    fprintf(ofp, "/mx {%g mul 50 add} def\n", picfac);
    fprintf(ofp, "/my {%g mul 50 add} def\n", picfac);//updated April 2010
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/sl {setlinewidth} def\n");
    fprintf(ofp, "/sc {setrgbcolor} def\n");
    fprintf(ofp, "/s {stroke} def\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "8 scalefont\n");
    fprintf(ofp, "setfont\n");
    
    fprintf(ofp, "newpath\n");
    fprintf(ofp, "%g mx %g my m\n", xmin, ymin);
    fprintf(ofp, "%g mx %g my l\n", xmax, ymin);
    fprintf(ofp, "%g mx %g my l\n", xmax, ymax);
    fprintf(ofp, "%g mx %g my l\n", xmin, ymax);
    fprintf(ofp, "closepath\n");
    fprintf(ofp, "stroke\n");
    //show tissue points
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) {
        showpoint = 0;
        for (k = 1; k <= mzz; k++) if (nbou[i][j][k] > 0) showpoint = 1;
        if (showpoint == 1) fprintf(ofp, "%g mx %g my m (.) show\n", axt[i], ayt[j]);
        for (k = 1; k <= mzz; k++) if (nbou[i][j][k] == 1022) fprintf(ofp, "%g mx %g my m (X) show\n", axt[i], ayt[j]);
    }
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "4 scalefont\n");
    fprintf(ofp, "setfont\n");
    //plot vessels according to segvar in order from bottom to top according to z-coordinate
    xzmin = 1.e6;
    xzmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        xzmin = FMIN(xzmin, segvar[iseg]);
        xzmax = FMAX(xzmax, segvar[iseg]);
    }
    for (ilevel = 1; ilevel <= nlevel; ilevel++) {
        zbottom = zmin + (ilevel - 1)*(zmax - zmin) / nlevel;
        ztop = zmin + ilevel * (zmax - zmin) / nlevel;
        for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
            zcoord = 0.;
            for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
            if (zcoord >= zbottom && zcoord < ztop) {
                if (xzmin != xzmax) xz = (segvar[iseg] - xzmin) / (xzmax - xzmin);
                else xz = 0.75;
                blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
                green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
                red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
                fprintf(ofp, "%f %f %f sc\n", red, green, blue);
                fprintf(ofp, "%g sl\n", picfac*diam[iseg] * diamfac);//line widths scaled up by diamfac
                xs = 0.;
                ys = 0.;
                for (i = 1; i <= 3; i++) {
                    xs += (cnode[i][ista[iseg]] - xsl0[i])*cos[1][i];
                    ys += (cnode[i][ista[iseg]] - xsl0[i])*cos[2][i];
                }
                fprintf(ofp, "%g mx %g my m ", xs, ys);
                xs = 0.;
                ys = 0.;
                for (i = 1; i <= 3; i++) {
                    xs += (cnode[i][iend[iseg]] - xsl0[i])*cos[1][i];
                    ys += (cnode[i][iend[iseg]] - xsl0[i])*cos[2][i];
                }
                fprintf(ofp, "%g mx %g my l s \n", xs, ys);
            }
        }
    }
    
    //label nodes in black
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for (inod = 1; inod <= nnod; inod++) {
        xs = 0.;
        ys = 0.;
        for (i = 1; i <= 3; i++) {
            xs += (cnode[i][inod] - xsl0[i])*cos[1][i];
            ys += (cnode[i][inod] - xsl0[i])*cos[2][i];
        }
        //comment out next two lines to remove node numbers
        fprintf(ofp, "%g mx %g my m ", xs + 0.5 / picfac, ys);
        fprintf(ofp, "(%g) show\n", nodvar[inod]);
    }
    //label segments in blue
    fprintf(ofp, "0 0 1 setrgbcolor\n");//blue
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        xs = 0.;
        ys = 0.;
        for (i = 1; i <= 3; i++) {
            xs += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*cos[1][i];
            ys += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*cos[2][i];
        }
        //comment out next two lines to remove segment numbers
        fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac, ys);
        fprintf(ofp, "(%g) show\n", segvar[iseg]);
    }
    //create a color bar
    float c;
    float cbbox = 15.; //size of boxes
    float cbx = 560; //origin of color bar
    float cby = 100;//origin of color bar
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "8 scalefont\n");
    fprintf(ofp, "setfont\n");
    for (k = 0; k <= 10; k++) {
        xz = k * 0.1;
        c = xzmin + (xzmax - xzmin)*xz;
        blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
        green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
        red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
        fprintf(ofp, "%f %f %f setrgbcolor\n", red, green, blue);
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx, cby + k * cbbox, cbx + cbbox, cby + k * cbbox, cbx + cbbox, cby + (k + 1)*cbbox, cbx, cby + (k + 1)*cbbox);
        if (k > 0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n", cbx + cbbox * 1.1, cby + cbbox * (k - 0.1), c);
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx, cby, cbx + cbbox, cby, cbx + cbbox, cby + cbbox * 11, cbx, cby + cbbox * 11);
    fprintf(ofp, "showpage\n");
    fclose(ofp);
    free_matrix(cos, 1, 3, 1, 3);
}
//-----------------------------------------------------------------/
// cmgui - for Greens.  TWS December 2011
// Based on code provided by Gib Bogle
// Produces greens.exelem and greens.exnode
// To visualize, use CMGUI (e.g. cmgui_install.msi) at:
// http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/
// Start CMGUI, File > Open > Com File
// Navigate to the directory containing the data files and select greens.com.txt (as listed below).
// At the bottom of the window that pops up, click the All button.  Change the view with the Layout,
// Up and Front selections on the left hand side of the Graphics window.  Rotate by holding the left
// button down and moving the mouse, pan with the middle button and zoom with the right button.
// Use Graphics > Scene editor to change the display.  Note also the Spectrum editor.
// ++++++++++++++++++greens.com.txt+++++++++++++++++++++++++++++++++
// # Create a material in addition to the default.
// gfx cre mat gold ambient 1 0.7 0 diffuse 1 0.7 0 specular 0.5 0.5 0.5 shininess 0.8
 
// gfx create spectrum jet
// gfx modify spectrum jet clear overwrite_colour
// gfx modify spectrum jet linear range 0 1 red   colour_range 0 1 ambient diffuse component 1
// gfx modify spectrum jet linear range 0 1 green colour_range 0 1 ambient diffuse component 2
// gfx modify spectrum jet linear range 0 1 blue  colour_range 0 1 ambient diffuse component 3
 
// # Read in the reticular mesh (group vessels) and hide the axes.
// gfx read nodes greens001.exnode
// gfx read elements greens001.exelem
 
// # The radius of the vessel is stored in component 1 of field
// # 'vessel_radius', defined over the elements in the vessels group.
 
// # Destroy the default lines.
// gfx modify g_element vessels lines delete
 
// gfx destroy node all
// gfx modify g_element vessels general clear;
// gfx modify g_element vessels cylinders coordinate coordinates tessellation default local circle_discretization 12 radius_scalar vessel_radius scale_factor 1 native_discretization NONE data node_colour spectrum jet
// gfx modify g_element vessels node_points coordinate coordinates local glyph sphere general size "0*0*0" centre 0,0,0 font default orientation vessel_radius scale_factors "2*2*2" data node_colour spectrum jet
 
// # Open the graphics window and turn on perspective (if desired).
// gfx cre win 1
// gfx mod win 1 view perspective
//-----------------------------------------------------------------/
void cmgui(float *segvar)
{
    extern int max, mxx, myy, mzz, nseg, nnod;
    extern int *segtyp, *ista, *iend;
    extern float *diam, **cnode;
    extern char numstr[6];
    
    int iseg, is, in;
    float red, green, blue, xz, xzmin, xzmax;
    char fname[80];
    FILE *exelem, *exnode;
    
    strcpy(fname, "Current\\greens");
    strcat(fname, numstr);
    strcat(fname, ".exelem");
    exelem = fopen(fname, "w");
    strcpy(fname, "Current\\greens");
    strcat(fname, numstr);
    strcat(fname, ".exnode");
    exnode = fopen(fname, "w");
    WriteExelemHeader(exelem);
    WriteExnodeHeader(exnode);
    
    xzmin = 1.e6;
    xzmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        xzmin = FMIN(xzmin, segvar[iseg]);
        xzmax = FMAX(xzmax, segvar[iseg]);
    }
    
    is = 0;
    in = 0;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        if (xzmin != xzmax) xz = (segvar[iseg] - xzmin) / (xzmax - xzmin);
        else xz = 0.75;
        blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
        green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
        red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
        is++;
        in++;
        //  write to elements file
        fprintf(exelem, "Element: %d 0 0\n", is);
        fprintf(exelem, "  Nodes: %d %d\n", in, in + 1);
        fprintf(exelem, "  Scale factors: 1 1\n");
        //  write to nodes file
        fprintf(exnode, "Node: %d\n", in);
        fprintf(exnode, "%6.1f %6.1f %6.1f\n", cnode[1][ista[iseg]], cnode[2][ista[iseg]], cnode[3][ista[iseg]]);
        fprintf(exnode, "%6.2f\n", diam[iseg] / 2.);
        fprintf(exnode, "%6.2f %6.2f %6.2f\n", red, green, blue);
        in++;
        fprintf(exnode, "Node: %d\n", in);
        fprintf(exnode, "%6.1f %6.1f %6.1f\n", cnode[1][iend[iseg]], cnode[2][iend[iseg]], cnode[3][iend[iseg]]);
        fprintf(exnode, "%6.2f\n", diam[iseg] / 2.);
        fprintf(exnode, "%6.2f %6.2f %6.2f\n", red, green, blue);
    }
    fclose(exelem);
    fclose(exnode);
}
//-----------------------------------------------------------------/
// contr_lines - generates contour lines and labels with contour heights
// Variables -
// m n nl:  dimensions of array, no. of contour levels
// scalefac: determines size of plot
// xmin xmax ymin ymax:  boundaries of box
// cl(nl):  array of contour levels
// zv(m,n):  array of heights
// Output to a postscript file.
// TWS, November 1989. Converted to C September 2007.  Revised February 2009.
//-----------------------------------------------------------------/
void contr_lines(FILE *ofp, int m, int n, double scalefac, int nl,
                 float xmin, float xmax, float ymin, float ymax, float *cl, float **zv)
{
    int i, j, k, iwsp, in, in1, nsq, nsq1, isq, isq1;
    int **nwsp;
    float c, d, di1, dj1, cx, cy;
    float *xv, *yv;
    float **wsp;
    const int ncmax = 10000;
    xv = vector(1, m);
    yv = vector(1, n);
    wsp = matrix(1, ncmax, 1, 2);
    nwsp = imatrix(1, ncmax, 1, 2);
    cx = 50; //origin of contour plot
    cy = 50; //origin of contour plot
    for (i = 1; i <= m; i++) xv[i] = xmin + (i - 1)*(xmax - xmin) / (m - 1);
    for (j = 1; j <= n; j++) yv[j] = ymin + (j - 1)*(ymax - ymin) / (n - 1);
    fprintf(ofp, "/mx {%g sub %g mul %g add} def\n", xmin, scalefac, cx);
    fprintf(ofp, "/my {%g sub %g mul %g add} def\n", ymin, scalefac, cy);
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "newpath\n");
    fprintf(ofp, "%g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l\n",
            xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax);
    fprintf(ofp, "closepath stroke\n");
    fprintf(ofp, "/Times-Roman findfont 8 scalefont setfont\n");
    //Begin contour generation
    for (k = 1; k <= nl; k++) {
        iwsp = 0;
        c = cl[k];
        for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) {
            //Look at vertical segments of mesh to detect crossing points.  Label with index numbers of adjacent two squares
            d = zv[i][j] - c;
            if (j != n) {
                dj1 = zv[i][j + 1] - c;
                if ((d >= 0 && dj1 < 0) || (d < 0 && dj1 >= 0)) {
                    iwsp += 1;
                    if (iwsp >= ncmax) printf("*** Error: ncmax too small in contr\n");
                    wsp[iwsp][1] = xv[i];
                    wsp[iwsp][2] = (d*yv[j + 1] - dj1 * yv[j]) / (d - dj1);
                    if (i == 1) nwsp[iwsp][1] = 0; else nwsp[iwsp][1] = i - 1 + (j - 1)*m;
                    if (i == m) nwsp[iwsp][2] = 0; else nwsp[iwsp][2] = i + (j - 1)*m;
                }
            }
            //Look at horizontal segments
            if (i != m) {
                di1 = zv[i + 1][j] - c;
                if ((d >= 0 && di1 < 0) || (d < 0 && di1 >= 0)) {
                    iwsp += 1;
                    wsp[iwsp][1] = (d*xv[i + 1] - di1 * xv[i]) / (d - di1);
                    wsp[iwsp][2] = yv[j];
                    if (j == 1) nwsp[iwsp][1] = 0; else nwsp[iwsp][1] = i + (j - 2)*m;
                    if (j == n) nwsp[iwsp][2] = 0; else nwsp[iwsp][2] = i + (j - 1)*m;
                }
            }
        }
        //Work through crossing points in sequence.  If nwsp(1 or 2) is not zero, find another point with matching value.
        //Join them and set nwsp values to zero at both.  Continue around loops.
        for (in = 1; in <= iwsp; in++) for (isq = 1; isq <= 2; isq++) {
            nsq = nwsp[in][isq];
            if (nsq != 0) {
                fprintf(ofp, "newpath\n");
                fprintf(ofp, "%g mx %g my m\n", wsp[in][1], wsp[in][2]);
                fprintf(ofp, "(%g) show\n", c);  //label contour
                fprintf(ofp, "%g mx %g my m\n", wsp[in][1], wsp[in][2]); //return to starting point
                nwsp[in][isq] = 0;
                do {
                    for (in1 = 1; in1 <= iwsp; in1++) for (isq1 = 1; isq1 <= 2; isq1++) {
                        nsq1 = nwsp[in1][isq1];
                        if (nsq1 == nsq && nsq != 0) {
                            fprintf(ofp, "%g mx %g my l\n", wsp[in1][1], wsp[in1][2]);
                            nwsp[in1][isq1] = 0;
                            nsq = nwsp[in1][3 - isq1];
                            nwsp[in1][3 - isq1] = 0;
                        }
                    }
                } while (nsq != 0);
                fprintf(ofp, "stroke\n");
            }
        }
    }
}
//-----------------------------------------------------------------/
// contr_shade - generates contour lines with colored shading between contours
// Variables -
// m n nl:  dimensions of array, no. of contour levels
// scalefac: determines size of plot
// xmin xmax ymin ymax:  boundaries of box
// cl(nl):  array of contour levels
// zv(m,n):  array of heights
// Output to a postscript file.
// Version including color shading of region, TWS February 2009.  Updated June 2009.
// Rewritten November 2010 to combine polygons with the same color before writing to
// the PostScript file.  This yields much smaller files.  They key is that each region
// is always enlcosed within a single closed curve, even if it contains 'holes', so
// that it is colored correctly.  When polygons are combined, only consecutive sides
// merged, so that the curve is not split.  Annular regions contain a cut.
// s-polygon: polygon to be shaded within a given rectangle ('square')
// r-polygon: polygon to be shaded consisting of multiple s-polygons
// flagr - region found
// flags - square found
// flaga - adjacent square found
// Outline of method:
// ____________________________________________
// label entire region with lowest color.
// subdivide region into squares
// for contour k
// label all squares type 0
// set number of r-polygons nrp = 0
// do
// set flagr = 0
// search for a square of type 0 with >= 1 corner above contour level
// if found
// increment nrp
// do
// set flags = 0
// search for a square of type -nrp
// if found
// do
// set flaga = 0
// search current then neighboring squares for a square of type = 0
// with >= 1 corner in common with current square above contour level
// if found
// set flagr = 1
// set flags = 1
// set flaga = 1
// define s-polygon for found square
// combine s-polygon with existing r-polygon
// label found square type -nrp
// else label current square type nrp
// while flaga = 1
// endif
// while flags = 1
// endif
// while flagr = 1
// end do
//-----------------------------------------------------------------/
void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl, float xmin,
                 float xmax, float ymin, float ymax, float *cl, float **zv,
                 int showscale, int lowcolor, int hatch, int plotcontour)
{
    fprintf(ofp, "%%!PS\n");
    int i, j, k, iwsp, flagr, flags, flaga, nrp, nspoly, nrpoly, irpoly, ispoly;
    int iin, iint, jjn, jjnt, iseed, jseed, match1, match2, dnrpoly, ispoly1, irpoly1, iipoly;
    int irpolystart, irpolyend, ispolystart, ispolyend;
    int in, in1, in2, in3, in4, inh;
    int *ii, *jj, **corners, **sqtype;
    const int ncmax = 100000, nrpolym = 10000;
    double xz, xv12, yv12, xv23, yv23, xv34, yv34, xv41, yv41;
    double cx, cy, cbx, cby, cbbox, polytolx, polytoly;
    double *xv, *yv, *red, *blue, *green, *dd, *xvv, *yvv, *xspoly, *yspoly, *xrpoly, *yrpoly;
    double **wsp;
    
    cx = 50; //origin of contour plot
    cy = 50;//origin of contour plot
    cbbox = 15.; //size of boxes
    cbx = 560; //origin of color bar
    cby = 100;//origin of color bar
    polytolx = 1.e-6*(xmax - xmin) / (m - 1);
    polytoly = 1.e-6*(ymax - ymin) / (n - 1);
    
    xv = dvector(1, m);
    yv = dvector(1, n);
    xvv = dvector(1, 4);
    yvv = dvector(1, 4);
    red = dvector(0, nl);
    green = dvector(0, nl);
    blue = dvector(0, nl);
    wsp = dmatrix(1, ncmax, 1, 4);
    corners = imatrix(1, 4, 1, 2);
    sqtype = imatrix(1, m - 1, 1, n - 1);
    xspoly = dvector(1, 6);
    yspoly = dvector(1, 6);
    xrpoly = dvector(1, nrpolym);
    yrpoly = dvector(1, nrpolym);
    dd = dvector(1, 4);
    ii = ivector(1, 4);
    jj = ivector(1, 4);
    corners[1][1] = 0;
    corners[1][2] = 0;
    corners[2][1] = 1;
    corners[2][2] = 0;
    corners[3][1] = 1;
    corners[3][2] = 1;
    corners[4][1] = 0;
    corners[4][2] = 1;
    
    for (i = 1; i <= m; i++) xv[i] = xmin + (i - 1)*(xmax - xmin) / (m - 1);
    for (j = 1; j <= n; j++) yv[j] = ymin + (j - 1)*(ymax - ymin) / (n - 1);
    fprintf(ofp, "/mx {%g sub %g mul %g add} def\n", xmin, scalefac, cx);
    fprintf(ofp, "/my {%g sub %g mul %g add} def\n", ymin, scalefac, cy);
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/s {stroke} def\n");
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "8 scalefont\n");
    fprintf(ofp, "setfont\n");
    if (hatch > 0) {
        fprintf(ofp, "/diagonals1\n");
        fprintf(ofp, "{ newpath\n");
        fprintf(ofp, "-700 10 550\n");
        fprintf(ofp, "{ 0 moveto\n");
        fprintf(ofp, "700 700 rlineto } for\n");
        fprintf(ofp, "stroke } def\n");
        
        fprintf(ofp, "/diagonals2\n");
        fprintf(ofp, "{ newpath\n");
        fprintf(ofp, "-700 10 550\n");
        fprintf(ofp, "{ 700 moveto\n");
        fprintf(ofp, "700 -700 rlineto } for\n");
        fprintf(ofp, "stroke } def\n");
        
        fprintf(ofp, "/verticals\n");
        fprintf(ofp, "{ newpath\n");
        fprintf(ofp, "50 10 550\n");
        fprintf(ofp, "{ 0 moveto\n");
        fprintf(ofp, "0 700 rlineto } for\n");
        fprintf(ofp, "stroke } def\n");
        
        fprintf(ofp, "/horizontals\n");
        fprintf(ofp, "{ newpath\n");
        fprintf(ofp, "100 10 800\n");
        fprintf(ofp, "{ 0 exch moveto\n");
        fprintf(ofp, "500 0 rlineto } for\n");
        fprintf(ofp, "stroke } def\n");
    }
    //Set up colors using Matlab 'jet' scheme
    for (k = 0; k <= nl; k++) {
        xz = float(k) / float(nl);
        blue[k] = FMIN(FMAX(1.5 - 4 * fabs(xz - 0.25), 0.), 1.);
        green[k] = FMIN(FMAX(1.5 - 4 * fabs(xz - 0.5), 0.), 1.);
        red[k] = FMIN(FMAX(1.5 - 4 * fabs(xz - 0.75), 0.), 1.);
    }
    //Color whole region with lowest color
    if (lowcolor > 0) {
        fprintf(ofp, "%f %f %f setrgbcolor\n", red[0], green[0], blue[0]);
        fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cf\n",
                xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax);
    }
    //Analyze each rectangle separately. Overwrite lower colors
    iwsp = 0;
    for (k = 1; k <= nl; k++) {
        fprintf(ofp, "%f %f %f setrgbcolor\n", red[k], green[k], blue[k]);
        for (i = 1; i < m; i++) for (j = 1; j < n; j++) sqtype[i][j] = 0;
        nrp = 0;
        do {
            //search for a square of type 0 with >= 1 corner above contour level
            flagr = 1;
            for (i = 1; i < m; i++) for (j = 1; j < n; j++) if (sqtype[i][j] == 0) for (in = 1; in <= 4; in++)
                if (zv[i + corners[in][1]][j + corners[in][2]] > cl[k]) goto foundit;
            flagr = 0;
        foundit:;
            if (flagr == 1) {
                nrp++;
                iseed = i;
                jseed = j;        //'seed' region nrp
                nrpoly = 0;        //initialize r-polygon
                do {
                    //search for a square of type -nrp
                    flags = 1;
                    if (i == iseed && j == jseed && sqtype[i][j] == 0) goto foundit1;
                    for (i = 1; i < m; i++) for (j = 1; j < n; j++) if (sqtype[i][j] == -nrp) goto foundit1;
                    flags = 0;
                foundit1:;
                    if (flags == 1) {
                        do {
                            //search current then neighboring squares for square type 0, with a corner in common with current square above contour level
                            flaga = 1;
                            for (inh = 0; inh <= 4; inh++) {
                                iin = i;
                                jjn = j;
                                if (inh == 1) iin++;
                                if (inh == 2) jjn++;
                                if (inh == 3) iin--;
                                if (inh == 4) jjn--;
                                if (iin > 0 && iin < m && jjn > 0 && jjn < n) if (sqtype[iin][jjn] == 0) {
                                    for (in = 1; in <= 4; in++) {
                                        iint = iin + corners[in][1];
                                        jjnt = jjn + corners[in][2];
                                        if ((iint == i || iint == i + 1) && (jjnt == j || jjnt == j + 1)
                                            && (zv[iint][jjnt] > cl[k])) goto foundit2;
                                    }
                                }
                            }
                            flaga = 0;
                        foundit2:;
                            if (flaga == 1) {
                                sqtype[iin][jjn] = -nrp;
                                //define s-polygon for found square (iin,jjn) - may have up to 6 sides
                                in1 = in;                //this is always in the region
                                in2 = in1 % 4 + 1;
                                in3 = in2 % 4 + 1;
                                in4 = in3 % 4 + 1;
                                for (in = 1; in <= 4; in++) {
                                    ii[in] = iin + corners[in][1];
                                    jj[in] = jjn + corners[in][2];
                                    dd[in] = zv[ii[in]][jj[in]] - cl[k];
                                    xvv[in] = xv[ii[in]];
                                    yvv[in] = yv[jj[in]];
                                }
                                if (dd[in1] != dd[in2]) {
                                    xv12 = (dd[in1] * xv[ii[in2]] - dd[in2] * xv[ii[in1]]) / (dd[in1] - dd[in2]);
                                    yv12 = (dd[in1] * yv[jj[in2]] - dd[in2] * yv[jj[in1]]) / (dd[in1] - dd[in2]);
                                }
                                if (dd[in2] != dd[in3]) {
                                    xv23 = (dd[in2] * xv[ii[in3]] - dd[in3] * xv[ii[in2]]) / (dd[in2] - dd[in3]);
                                    yv23 = (dd[in2] * yv[jj[in3]] - dd[in3] * yv[jj[in2]]) / (dd[in2] - dd[in3]);
                                }
                                if (dd[in3] != dd[in4]) {
                                    xv34 = (dd[in3] * xv[ii[in4]] - dd[in4] * xv[ii[in3]]) / (dd[in3] - dd[in4]);
                                    yv34 = (dd[in3] * yv[jj[in4]] - dd[in4] * yv[jj[in3]]) / (dd[in3] - dd[in4]);
                                }
                                if (dd[in4] != dd[in1]) {
                                    xv41 = (dd[in4] * xv[ii[in1]] - dd[in1] * xv[ii[in4]]) / (dd[in4] - dd[in1]);
                                    yv41 = (dd[in4] * yv[jj[in1]] - dd[in1] * yv[jj[in4]]) / (dd[in4] - dd[in1]);
                                }
                                xspoly[1] = xvv[in1];
                                yspoly[1] = yvv[in1];
                                if (dd[in2] > 0) {                        //corners 1,2 are this color
                                    xspoly[2] = xvv[in2];
                                    yspoly[2] = yvv[in2];
                                    if (dd[in3] > 0) {                    //corners 1,2,3 are this color
                                        xspoly[3] = xvv[in3];
                                        yspoly[3] = yvv[in3];
                                        if (dd[in4] > 0) {                //corners 1,2,3,4 are this color
                                            xspoly[4] = xvv[in4];
                                            yspoly[4] = yvv[in4];
                                            nspoly = 4;
                                        }
                                        else {                            //corners 1,2,3,not 4 are this color
                                            xspoly[4] = xv34;
                                            yspoly[4] = yv34;
                                            xspoly[5] = xv41;
                                            yspoly[5] = yv41;
                                            nspoly = 5;
                                            iwsp++;
                                            wsp[iwsp][1] = xv34;
                                            wsp[iwsp][2] = yv34;
                                            wsp[iwsp][3] = xv41;
                                            wsp[iwsp][4] = yv41;
                                        }
                                    }
                                    else {                                //corners 1,2,not 3 are this color
                                        xspoly[3] = xv23;
                                        yspoly[3] = yv23;
                                        iwsp++;
                                        wsp[iwsp][1] = xv23;
                                        wsp[iwsp][2] = yv23;
                                        if (dd[in4] > 0) {                //corners 1,2,not 3,4 are this color
                                            xspoly[4] = xv34;
                                            yspoly[4] = yv34;
                                            xspoly[5] = xvv[in4];
                                            yspoly[5] = yvv[in4];
                                            nspoly = 5;
                                            wsp[iwsp][3] = xv34;
                                            wsp[iwsp][4] = yv34;
                                        }
                                        else {                            //corners 1,2,not 3,not 4 are this color
                                            xspoly[4] = xv41;
                                            yspoly[4] = yv41;
                                            nspoly = 4;
                                            wsp[iwsp][3] = xv41;
                                            wsp[iwsp][4] = yv41;
                                        }
                                    }
                                }
                                else {                                    //corners 1,not 2 are this color
                                    xspoly[2] = xv12;
                                    yspoly[2] = yv12;
                                    iwsp++;
                                    wsp[iwsp][1] = xv12;
                                    wsp[iwsp][2] = yv12;
                                    if (dd[in3] > 0) {                    //corners 1,not 2,3 are this color
                                        xspoly[3] = xv23;
                                        yspoly[3] = yv23;
                                        xspoly[4] = xvv[in3];
                                        yspoly[4] = yvv[in3];
                                        wsp[iwsp][3] = xv23;
                                        wsp[iwsp][4] = yv23;
                                        if (dd[in4] > 0) {                //corners 1,not 2,3,4 are this color
                                            xspoly[5] = xvv[in4];
                                            yspoly[5] = yvv[in4];
                                            nspoly = 5;
                                        }
                                        else {                            //corners 1,not 2,3,not 4 are this color
                                            xspoly[5] = xv34;
                                            yspoly[5] = yv34;
                                            xspoly[6] = xv41;
                                            yspoly[6] = yv41;
                                            nspoly = 6;
                                            iwsp++;
                                            wsp[iwsp][1] = xv34;
                                            wsp[iwsp][2] = yv34;
                                            wsp[iwsp][3] = xv41;
                                            wsp[iwsp][4] = yv41;
                                        }
                                    }
                                    else {                                //corners 1,not 2,not 3 are this color
                                        if (dd[in4] > 0) {                //corners 1,not 2,not 3,4 are this color
                                            xspoly[3] = xv34;
                                            yspoly[3] = yv34;
                                            xspoly[4] = xvv[in4];
                                            yspoly[4] = yvv[in4];
                                            nspoly = 4;
                                            wsp[iwsp][3] = xv34;
                                            wsp[iwsp][4] = yv34;
                                        }
                                        else {                            //corners 1,not 2,not 3,not 4 are this color
                                            xspoly[3] = xv41;
                                            yspoly[3] = yv41;
                                            nspoly = 3;
                                            wsp[iwsp][3] = xv41;
                                            wsp[iwsp][4] = yv41;
                                        }
                                    }
                                }
                                if (iwsp > ncmax) printf("*** Error: too many contour points.  Increase ncmax ***\n");
                                //combine s-polygon with existing r-polygon, eliminating redundant segments
                                if (nrpoly == 0) {                        //initiate r-polygon
                                    for (ispoly = 1; ispoly <= nspoly; ispoly++) {
                                        xrpoly[ispoly] = xspoly[ispoly];
                                        yrpoly[ispoly] = yspoly[ispoly];
                                    }
                                    nrpoly = nspoly;
                                }
                                else {                    //search r-polygon and s-polygon for one side that matches
                                    for (irpoly = nrpoly; irpoly >= 1; irpoly--) for (ispoly = 1; ispoly <= nspoly; ispoly++) {
                                        ispoly1 = ispoly % nspoly + 1;
                                        irpoly1 = irpoly % nrpoly + 1;
                                        if ((fabs(xrpoly[irpoly] - xspoly[ispoly1]) < polytolx)
                                            && (fabs(yrpoly[irpoly] - yspoly[ispoly1]) < polytoly)
                                            && (fabs(xrpoly[irpoly1] - xspoly[ispoly]) < polytolx)
                                            && (fabs(yrpoly[irpoly1] - yspoly[ispoly]) < polytoly)) goto foundit3;
                                    }
                                    printf("*** Error: matching segment not found ***\n");
                                foundit3:;
                                    match1 = 0;
                                    for (iipoly = 2; iipoly < nspoly; iipoly++) {    //search for further matches nearby
                                        irpoly1 = (irpoly - iipoly + nrpoly) % nrpoly + 1;
                                        ispoly1 = (ispoly + iipoly - 1) % nspoly + 1;
                                        if ((fabs(xrpoly[irpoly1] - xspoly[ispoly1]) < polytolx)
                                            && (fabs(yrpoly[irpoly1] - yspoly[ispoly1]) < polytoly)) match1++;
                                        else goto nomatch1;
                                    }
                                nomatch1:;
                                    match2 = 0;
                                    for (iipoly = 2; iipoly < nspoly; iipoly++) {    //search other way for further matches nearby
                                        irpoly1 = (irpoly + iipoly - 1) % nrpoly + 1;
                                        ispoly1 = (ispoly - iipoly + nspoly) % nspoly + 1;
                                        if ((fabs(xrpoly[irpoly1] - xspoly[ispoly1]) < polytolx)
                                            && (fabs(yrpoly[irpoly1] - yspoly[ispoly1]) < polytoly)) match2++;
                                        else goto nomatch2;
                                    }
                                nomatch2:;
                                    ispolystart = (ispoly + match1) % nspoly + 1;                //first node of s-polygon to include
                                    ispolyend = (ispoly - match2 - 1 + nspoly) % nspoly + 1;    //last node of s-polygon to include
                                    irpolystart = (irpoly + match2) % nrpoly + 1;                //first node of s-polygon to include
                                    irpolyend = (irpoly - match1 - 1 + nrpoly) % nrpoly + 1;    //last node of s-polygon to include
                                    if ((fabs(xrpoly[irpolystart] - xspoly[ispolyend]) > polytolx)
                                        || (fabs(yrpoly[irpolystart] - yspoly[ispolyend]) > polytoly))
                                        printf("*** Error: r and s nodes do not match ****\n");
                                    if ((fabs(xrpoly[irpolyend] - xspoly[ispolystart]) > polytolx)
                                        || (fabs(yrpoly[irpolyend] - yspoly[ispolystart]) > polytoly))
                                        printf("*** Error: r and s nodes do not match ****\n");
                                    dnrpoly = nspoly - 2 - 2 * match1 - 2 * match2;
                                    if (irpolystart > irpolyend) {
                                        if (dnrpoly > 0) {        //expand the arrays xrpoly and yrpoly
                                            for (iipoly = nrpoly; iipoly >= irpolystart; iipoly--) {
                                                xrpoly[iipoly + dnrpoly] = xrpoly[iipoly];
                                                yrpoly[iipoly + dnrpoly] = yrpoly[iipoly];
                                            }
                                        }
                                        if (dnrpoly < 0) {        //contract the arrays xrpoly and yrpoly
                                            for (iipoly = irpolystart; iipoly <= nrpoly; iipoly++) {
                                                xrpoly[iipoly + dnrpoly] = xrpoly[iipoly];
                                                yrpoly[iipoly + dnrpoly] = yrpoly[iipoly];
                                            }
                                        }
                                    }
                                    nrpoly += dnrpoly;
                                    if (nrpoly > nrpolym) printf("*** Error: too many polygon points.  Increase nrpolym ***\n");
                                    if (nrpoly < irpolyend) for (iipoly = irpolyend; iipoly > nrpoly; iipoly--) {    //otherwise these values get lost!
                                        xrpoly[iipoly - nrpoly] = xrpoly[iipoly];
                                        yrpoly[iipoly - nrpoly] = yrpoly[iipoly];
                                    }
                                    for (iipoly = 1; iipoly <= nspoly - 2 - match1 - match2; iipoly++) {
                                        irpoly1 = (irpolyend + iipoly - 1) % nrpoly + 1;
                                        ispoly1 = (ispolystart + iipoly - 1) % nspoly + 1;
                                        xrpoly[irpoly1] = xspoly[ispoly1];
                                        yrpoly[irpoly1] = yspoly[ispoly1];
                                    }
                                }
                            }
                            else sqtype[i][j] = nrp;
                        }
                        while (flaga == 1);
                    }
                }
                while (flags == 1);
                //display polygons after combining
                if (hatch > 0) fprintf(ofp, "/clippingpath {\n");
                fprintf(ofp, "n %g mx %g my m ", xrpoly[1], yrpoly[1]);
                for (irpoly = 2; irpoly <= nrpoly; irpoly++) {
                    fprintf(ofp, "%g mx %g my l ", xrpoly[irpoly], yrpoly[irpoly]);
                    if (irpoly % 5 == 1) fprintf(ofp, "\n");
                }
                if (hatch > 0) {
                    fprintf(ofp, "} def\n");
                    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
                    fprintf(ofp, "gsave clippingpath clip\n");
                    if (hatch == 1) fprintf(ofp, "diagonals1 grestore\n");
                    if (hatch == 2) fprintf(ofp, "diagonals2 grestore\n");
                }
                else fprintf(ofp, "cf\n");
            }
        }
        while (flagr == 1);
    }
    //outline contours
    if (plotcontour > 0) {
        fprintf(ofp, "0 0 0 setrgbcolor\n");//black
        for (in = 1; in <= iwsp; in++) fprintf(ofp, "n %g mx %g my m %g mx %g my l s\n",
                                               wsp[in][1], wsp[in][2], wsp[in][3], wsp[in][4]);
    }
    //draw a box
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cs\n",
            xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax);
    //create a color bar
    if (showscale > 0) {
        for (k = 0; k <= nl; k++) {
            fprintf(ofp, "%f %f %f setrgbcolor\n", red[k], green[k], blue[k]);
            fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                    cbx, cby + k * cbbox, cbx + cbbox, cby + k * cbbox, cbx + cbbox, cby + (k + 1)*cbbox, cbx, cby + (k + 1)*cbbox);
            if (k > 0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n", cbx + cbbox * 1.1, cby + cbbox * (k - 0.1), cl[k]);
        }
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
                cbx, cby, cbx + cbbox, cby, cbx + cbbox, cby + cbbox * (nl + 1), cbx, cby + cbbox * (nl + 1));
    }
    free_dvector(xv, 1, m);
    free_dvector(yv, 1, n);
    free_dvector(xvv, 1, 4);
    free_dvector(yvv, 1, 4);
    free_dvector(red, 0, nl);
    free_dvector(green, 0, nl);
    free_dvector(blue, 0, nl);
    free_dmatrix(wsp, 1, ncmax, 1, 4);
    free_imatrix(corners, 1, 4, 1, 2);
    free_imatrix(sqtype, 1, m - 1, 1, n - 1);
    free_dvector(xspoly, 1, 6);
    free_dvector(yspoly, 1, 6);
    free_dvector(xrpoly, 1, nrpolym);
    free_dvector(yrpoly, 1, nrpolym);
    free_dvector(dd, 1, 4);
    free_ivector(ii, 1, 4);
    free_ivector(jj, 1, 4);
}

/**********************************************************
 contour.cpp - generate data for contour plot.  TWS Dec. 07
 Version 3.0, May 17, 2011.
 Produces a single postscript file with a page for each solute
 ***********************************************************/
void contour(const char fname[])
{
    extern int max, nsp, nseg, *segtyp, *nl, *ista, *iend;
    extern int slsegdiv, nsl1, nsl2;
    extern float pi1, req, scalefac;
    extern float *x, *p, *diam, **cnode, **pvseg;
    extern float *xsl0, *xsl1, *xsl2, *clmin, *clint, *cl, **zv, ***psl;
    
    int i, iseg, isp, isl1, isl2, ilevel, nlevel = 100;
    float xmin, ymin, xmax, ymax, xs1, ys1, xs2, ys2, **cos;
    float red, green, blue, xz, xzmin, xzmax;
    float diamfac = 1., zcoord, zbottom, ztop, zmin, zmax;
    
    
    FILE *ofp;
    
    printf("Generating data for contour plots...");
    
    xmin = 0.;
    xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));
    ymin = 0.;
    ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));
    cos = matrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {    //set up matrix of direction cosines
        cos[1][i] = (xsl1[i] - xsl0[i]) / xmax;
        cos[2][i] = (xsl2[i] - xsl0[i]) / ymax;
    }
    cos[3][1] = cos[1][2] * cos[2][3] - cos[1][3] * cos[2][2];
    cos[3][2] = cos[1][3] * cos[2][1] - cos[1][1] * cos[2][3];
    cos[3][3] = cos[1][1] * cos[2][2] - cos[1][2] * cos[2][1];
    
    //Determine range of z values
    zmin = 1.e6;
    zmax = -1.e6;
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        zcoord = 0.;
        for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
        zmin = FMIN(zmin, zcoord - 1.);
        zmax = FMAX(zmax, zcoord + 1.);
    }
    
    //Calculate P on a planar slice through the region, specified by three corners and number of points along each edge
    //Subdivide vessel segments by slsegdiv for smoother plots
    for (isl1 = 1; isl1 <= nsl1; isl1++)    for (isl2 = 1; isl2 <= nsl2; isl2++) {
        for (i = 1; i <= 3; i++) x[i] = xsl0[i] + (isl1 - 1)*(xsl1[i] - xsl0[i]) / (nsl1 - 1) + (isl2 - 1)*(xsl2[i] - xsl0[i]) / (nsl2 - 1);
        p = eval(slsegdiv, req, x);
        for (isp = 1; isp <= nsp; isp++) psl[isl1][isl2][isp] = p[isp];
    }
    xmin = 0.;
    ymin = 0.;
    ofp = fopen(fname, "w");
    fprintf(ofp, "%%!PS-Adobe-2.0\n");
    fprintf(ofp, "%%%%Pages: %i\n", nsp);
    fprintf(ofp, "%%%%EndComments\n");
    for (isp = 1; isp <= nsp; isp++) {
        fprintf(ofp, "%%%%Page: %i %i\n", isp, isp);
        for (isl1 = 1; isl1 <= nsl1; isl1++)    for (isl2 = 1; isl2 <= nsl2; isl2++) zv[isl1][isl2] = psl[isl1][isl2][isp];
        for (i = 1; i <= nl[isp]; i++) cl[i] = clmin[isp] + (i - 1)*clint[isp];
        
        if (isp == 1) {
            cl[1] = 1.;        //override contour levels to give contours at p = 1, 2 and 5
            cl[2] = 2.;
            cl[3] = 5.;
        }
        
        
        contr_shade(ofp, nsl1, nsl2, scalefac, nl[isp], xmin, xmax, ymin, ymax, cl, zv, 1, 1, 0, 0);
        //        contr_lines(ofp,nsl1,nsl2,scalefac,nl[isp],xmin,xmax,ymin,ymax,cl,zv);
        
        fprintf(ofp, "/sl {setlinewidth} def\n");
        fprintf(ofp, "/sc {setrgbcolor} def\n");
        fprintf(ofp, "/s {stroke} def\n");
        fprintf(ofp, "1 setlinecap\n");
        
        //Plot projection of network in contour plane
        //plot vessels according to pvseg in order from bottom to top according to z-coordinate
        xzmin = clmin[isp];
        xzmax = clmin[isp] + (nl[isp] - 1)*clint[isp];
        for (ilevel = 1; ilevel <= nlevel; ilevel++) {
            zbottom = zmin + (ilevel - 1)*(zmax - zmin) / nlevel;
            ztop = zmin + ilevel * (zmax - zmin) / nlevel;
            for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
                zcoord = 0.;
                for (i = 1; i <= 3; i++)    zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
                if (zcoord >= zbottom && zcoord < ztop) {
                    if (xzmin != xzmax) xz = (pvseg[iseg][isp] - xzmin) / (xzmax - xzmin);
                    else xz = 0.75;
                    blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
                    green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
                    red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
                    
                    xs1 = 0.;
                    ys1 = 0.;
                    xs2 = 0.;
                    ys2 = 0.;
                    for (i = 1; i <= 3; i++) {
                        xs1 += (cnode[i][ista[iseg]] - xsl0[i])*cos[1][i];
                        ys1 += (cnode[i][ista[iseg]] - xsl0[i])*cos[2][i];
                        xs2 += (cnode[i][iend[iseg]] - xsl0[i])*cos[1][i];
                        ys2 += (cnode[i][iend[iseg]] - xsl0[i])*cos[2][i];
                    }
                    fprintf(ofp, "0 0 0 sc\n");    //Plot vessels slightly larger in black to outline
                    fprintf(ofp, "%g sl\n", scalefac*diam[iseg] * diamfac + 2.);
                    fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1, ys1, xs2, ys2);
                    fprintf(ofp, "%f %f %f sc\n", red, green, blue);
                    fprintf(ofp, "%g sl\n", scalefac*diam[iseg] * diamfac);//line widths scaled up by diamfac
                    fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1, ys1, xs2, ys2);
                }
            }
        }
        fprintf(ofp, "0 0 0 setrgbcolor\n");        //black
        fprintf(ofp, "/Times-Roman findfont 12 scalefont setfont\n");
        fprintf(ofp, "50 30 moveto\n");
        fprintf(ofp, "(Solute %i) show\n", isp);  //show solute number
        //create a scale bar
        float barlength = 100;
        if (xmax > 500.) barlength = 200.;
        if (xmax > 1500.) barlength = 500.;
        fprintf(ofp, "%g sl\n", 2.);
        fprintf(ofp, "%g %g m %g %g l stroke\n", 120., 25., 120. + scalefac * barlength, 25.);
        fprintf(ofp, "%g %g m (%g mm) show\n", 120., 30., barlength / 1000.);
        fprintf(ofp, "showpage\n");
    }
    fclose(ofp);
    printf("done\n");
}
//-----------------------------------------------------------------/
// Evaluate histograms of solute levels.  TWS December 07.
// Version 2.0, May 1, 2010.
// Version 3.0, May 17,2011.
//-----------------------------------------------------------------/
void histogram(const char fname[])
{
    extern int mxx, myy, mzz, nnt, nnv, nseg, nnod, nsp;
    extern float *pmin, *pmax, *pmean, *pref;
    extern float **pt;
    int i, j, itp, isp, nctop, ncbot, nstat, nstatm = 101;
    float step, dev;
    float *po2samp, *stat, *cumu, *mstat;
    FILE *ofp;
    po2samp = vector(1, nstatm);
    stat = vector(1, nstatm);
    cumu = vector(1, nstatm);
    mstat = vector(1, nstatm);
    
    ofp = fopen(fname, "w");
    for (isp = 1; isp <= nsp; isp++) {
        step = pref[isp] / 100;
        nctop = pmax[isp] / step + 1.;
        if (nctop > 100) nctop = 100;
        ncbot = pmin[isp] / step;
        if (ncbot <= 2) ncbot = 0;
        nstat = nctop - ncbot + 1;
        dev = 0.;
        if (nstat > nstatm) printf("*** nstatm too small in histogram\n");
        for (i = 1; i <= nstat; i++) {
            po2samp[i] = step * (i - 1 + ncbot);
            mstat[i] = 0;
        }
        for (itp = 1; itp <= nnt; itp++) {
            dev = dev + SQR(pmean[isp] - pt[itp][isp]);
            for (j = 1; j <= nstat; j++)    if (pt[itp][isp] <= po2samp[j]) {
                mstat[j] = mstat[j] + 1;
                goto binned;
            }
        binned:;
        }
        dev = sqrt(dev / nnt);
        for (i = 1; i <= nstat; i++) stat[i] = mstat[i] * 100. / nnt;
        cumu[1] = stat[1];
        for (i = 2; i <= nstat; i++) cumu[i] = cumu[i - 1] + stat[i];
        
        fprintf(ofp, "Histogram data for solute %i\n", isp);
        fprintf(ofp, "value  %% cumul. %%\n");
        for (i = 1; i <= nstat; i++) fprintf(ofp, "%g %7.2f %7.2f\n", po2samp[i], stat[i], cumu[i]);
        fprintf(ofp, "Mean = %f deviation = %f min = %f max  = %f\n", pmean[isp], dev, pmin[isp], pmax[isp]);
    }
    fclose(ofp);
    free_vector(po2samp, 1, nstatm);
    free_vector(stat, 1, nstatm);
    free_vector(cumu, 1, nstatm);
    free_vector(mstat, 1, nstatm);
}
/****************************************************************************
 initgreens - initial tissue source strengths, given uniform solute field, initial g0
 initial vessel source strengths based on uniform efflux rate from all vessels
 use this version only if values are not available from previous call to greens
 TWS Jan 08
 Version 2.0, May 1, 2010.
 Version 3.0, May 17, 2011.
 *******************************************************************************/
void initgreens()
{
    extern int nnt, nnv, nsp;
    extern int *lowflow, *mainseg, *permsolute;
    extern int *oxygen, *diffsolute; //added April 2010
    
    extern float vol, errfac, tlength;
    extern float tlengthq, tlengthqhd;//added 8/09
    
    extern float *mtiss, *mptiss, *epsvessel, *epstissue, *eps, *errvessel, *errtissue, *pinit, *p;
    extern float *g0, *qtsum, *pref, *ds, *qq, *hd;
    extern float **qt, **qv, **pt, **pv, **tissparam, *mptissref;
    
    int isp, i, itp;
    
    for (isp = 1; isp <= nsp; isp++)    pinit[isp] = g0[isp];    //initial values based on g0
    if(DEBG==1){
        for(i=1;i<=nsp;i++) printf("%d %d %f %f %f\n", nsp,i, pinit[i], mtiss[i], mptiss[i]);
    }
    tissrate(nsp, pinit, mtiss, mptiss);
    if(DEBG==1){
        for(i=1;i<=nsp;i++) printf("%d %d %f %f %f\n", nsp,i, pinit[i], mtiss[i], mptiss[i]);
    }
    for (isp = 1; isp <= nsp; isp++) {
        mptissref[isp] = mptiss[isp];
        qtsum[isp] = 0.;
        for (itp = 1; itp <= nnt; itp++) {
            qt[itp][isp] = mtiss[isp] * vol;// solute field assumed uniform
            pt[itp][isp] = pinit[isp];
            qtsum[isp] += qt[itp][isp];
        }
        for (i = 1; i <= nnv; i++) {
            qv[i][isp] = 0.;
            pv[i][isp] = 0.;
            if (permsolute[isp] == 1) {
                if (oxygen[isp] == 1) {
                    qv[i][isp] = -qtsum[isp] * ds[mainseg[i]] * qq[mainseg[i]] * hd[mainseg[i]] / tlengthqhd;//modified 8/10
                    if (lowflow[mainseg[i]] == 1) qv[i][isp] = 0.;  //low q*hd
                }
                else qv[i][isp] = -qtsum[isp] * ds[mainseg[i]] * qq[mainseg[i]] / tlengthq;//modified 8/09
                pv[i][isp] = pinit[isp];
            }
        }
        if(DEBG==1) for (isp = 1; isp <= nsp; isp++) for (i = 1; i <= nnv; i++)
        {
            printf("%d %d %d %f %f %f %f %f %f\n", nnv, i, mainseg[i], tlengthq, hd[mainseg[i]], qq[mainseg[i]], ds[mainseg[i]], qtsum[isp], qv[i][isp]);
        }
    }
    //set error bounds, proportional to errfac, based on pref separately for each solute, updated April 2015
    for (isp = 1; isp <= nsp; isp++) pinit[isp] = 0.;
    for (isp = 1; isp <= nsp; isp++) {
        pinit[isp] = pref[isp];
        if(DEBG==1){
            for(i=1;i<=nsp;i++) printf("%d %d %f %f %f\n", nsp,i, pinit[i], mtiss[i], mptiss[i]);
        }
        tissrate(nsp, pinit, mtiss, mptiss);
        epstissue[isp] = FMAX(fabs(mtiss[isp])*vol*errfac, 0.001);    //error bound for tissue source strengths
        epsvessel[isp] = nnt * epstissue[isp] / nnv;    //error bound for vessel source strengths
        eps[isp] = pref[isp] * errfac;    //error bound for concentrations or PO2
        pinit[isp] = 0.;
    }
}
/************************************************************************
 putrank - generate list of nodes in order of flow direction
 nodrank --- if nodrank[i] < nodrank[j], node j is not upstream of node i
 So node j is upstream of node i (node j feeds solute to node i) only if nodrank[i] >= nodrank[j] - Baran
 The node ranks are probably used to implement equation 20 from the paper - Baran
 Version 2.0, May 1, 2010.
 Version 3.0, May 17, 2011.
 *************************************************************************/

void putrank(void)
{
    extern int nseg, nnod, nnodfl, nnodbc, nsegfl;
    extern int *nodrank, *nodtyp, *nodout, *segtyp, *nk, *ista, *iend;
    extern int **nodseg, **nodnod;
    extern float *q;
    
    int inod, j, jseg, nod1, nod2, iseg, flag;
    
    //construct node table; count outputs from node; output nodes precede input nodes
    for (inod = 1; inod <= nnod; inod++) {
        nodtyp[inod] = 0;
        nodout[inod] = 0;
    }
    nsegfl = 0;    //added TWS 2010
    if(DEBG==1) {
        for(i=1;i<=nseg;i++) printf("%d %d %f\n", nseg, i, q[i]);
    }

    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        if (q[iseg] >= 0) {
            nod1 = ista[iseg];
            nod2 = iend[iseg];
        }
        else {
            nod1 = iend[iseg];
            nod2 = ista[iseg];
        }
        nodtyp[nod1]++;// nodtype is the number of inflows + outflows from the node
        nodseg[nodtyp[nod1]][nod1] = iseg;
        nodnod[nodtyp[nod1]][nod1] = nod2;
        nodout[nod1]++;//nodout is the number of outflows from the node
        nsegfl++;
    }
    for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
        if (q[iseg] >= 0) {
            nod1 = ista[iseg];
            nod2 = iend[iseg];
        }
        else {
            nod1 = iend[iseg];
            nod2 = ista[iseg];
        }
        nodtyp[nod2]++;
        nodseg[nodtyp[nod2]][nod2] = iseg;
        nodnod[nodtyp[nod2]][nod2] = nod1;
    }
    for (inod = 1; inod <= nnod; inod++)    if (nodtyp[inod] == 0) printf("***Warning: Node %i is not related to any segment\n", inod);
    //assign low ranks to inflow nodes
    nnodfl = 0;
    for (inod = 1; inod <= nnod; inod++) {
        nk[inod] = 0;
        if (nodtyp[inod] == 1 && nodout[inod] == 1) {
            nnodfl++;
            nk[inod] = 1;
            nodrank[nnodfl] = inod;
        }
    }
    if(DEBG==1) for (inod=1;inod<=nnod;inod++) {
        printf("%d %d %d %d %d\n", nnod, inod, nk[inod],nodrank[1],nodrank[2]);
    }
    //assign increasing ranks to downstream connected nodes
    //j runs over the number of segments that feed each node - Baran
    flag = 1;
    while (flag == 1) {
        flag = 0;
        for (inod = 1; inod <= nnod; inod++)    if (nk[inod] == 0 && nodtyp[inod] > 0) {//downstream from the boundary (since nk[inod]==0 and connected because nodtyp[inod]!=0)
            for (j = nodout[inod] + 1; j <= nodtyp[inod]; j++) {//Inflows
                jseg = nodseg[j][inod];
                if (inod == iend[jseg] && (nk[ista[jseg]] == 0 || q[jseg] <= 0.)) goto skipnode;
                if (inod == ista[jseg] && (nk[iend[jseg]] == 0 || q[jseg] >= 0.)) goto skipnode;
            }
            nnodfl++;
            nk[inod] = 1;
            nodrank[nnodfl] = inod;
            flag = 1;
        skipnode:;
        }
    }
    if(DEBG==1) for (inod=1;inod<=nnod;inod++) {
        printf("%d %d %d %d %d\n", nnod, inod, nk[inod],nodrank[1], nodrank[2]);
    }
}


#undef SWAP
