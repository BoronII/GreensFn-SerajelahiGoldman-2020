// Variable declarations. All variables must be local, and passed as argument. No variables should be global as done in the Secomb code.

    
    int max = 100;
	int nmaxvessel, nmaxtissue, nmax, rungreens, g0method, linmethod;
    int mxx, myy, mzz, nnt, nseg, nnod, nnodfl, nnv, nsp, nnodbc, nodsegm, nsegfl, kmain;
    int slsegdiv, nsl1, nsl2;
    int is2d; //needed for 2d version
    int nvaryparams, nruns, ntissparams, npostgreensparams, npostgreensout;    //needed for varying parameters, postgreens
    
    float  *dtmin;//added July 2011
    int *mainseg, *permsolute, *nodrank, *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp, *lowflow;
    int *nodname, *segname, *segtyp, *nspoint, *istart, *nl, *nk, *indx, *ista, *iend;
    int *errvesselcount, *errtissuecount;
    int *imaxerrvessel, *imaxerrtissue, *nresis;  //added April 2010
    int *oxygen, *diffsolute; //added April 2010
    int **nodseg, **tisspoints, **nodnod;
    float **segnodname; //change from int to double Aug 2019 - Baran
    int ***nbou;
    int **tissfix;    //added September 2010
    float  **tisserr, **dmtissdp, *mptissref;//September 2010;
    int **ivaryparams;    //added April 2015
    
    //double  gtt;    //added September 2010
    float  fn, c, alphab, p50, cs, cext, hext, req, q0fac, totalq, flowfac = 1.e6 / 60.;
    float  plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
    float  pi1 = atan(1.)*4., fac = 1. / 4. / pi1;
    float  lb, maxl, v, vol, vdom, errfac, tlength, alx, aly, alz, lowflowcrit;
    float  tlengthq, tlengthqhd;//added 8/09
    float  xmax, ymax, scalefac;
    float  w2d, r2d; //needed for 2d version
    
    float  *axt, *ayt, *azt, *ds, *diff, *pmin, *pmax, *pmean, *pref, *g0, *g0fac, *g0facnew, *sumal;
    float  *rseg, *q, *qdata, *qq, *hd, *oxflux, *segc, *bcprfl, *bchd, *nodvar, *segvar, *qvtemp, *qvfac;//added qdata November 2016
    float   **start, **scos, **ax, **resisdiam, **resis, **bcp; //added April 2010
    float   **qv, **qt, **pv, **pev, **pt;
    float  **qvseg, **pvseg, **pevseg;
    float  **paramvalue, *solutefac, *intravascfac, *postgreensparams, *postgreensout;    //added April 2015

    float  *x, *y, *lseg, *ss, *cbar, *mtiss, *mptiss, *dqvsumdg0, *dqtsumdg0;
    float  *epsvessel, *epstissue, *eps, *errvessel, *errtissue, *pinit, *p;
    float  *rhs, *rhstest, *g0old, *ptt, *ptpt, *qtsum, *qvsum;
    float  **pvt, **pvprev, **qvprev, **cv, **dcdp, **tissparam;
    float  **ptprev, **ptv, **gamma1, **qcoeff1, **cv0, **conv0;
    float  **gvv, **end, **al;
    float  ***rsta, ***rend, ***dtt;
    float  *xsl0, *xsl1, *xsl2, *clmin, *clint, *cl, **zv, ***psl;
    float  **qtp;
    float *rhstiss, *matxtiss;
    double **mat, **rhsg, *rhsl, *matx;
    
    char numstr[6];

    int i, j, k, iseg, isp, nlmax, idum;
    FILE *ifp;
    char bb[100];// used to hold the first line of SoluteParam.dat in the Secomb groups code, unused so far in this code
    int method;
    //int outboun(int method);
    int inod, inod1, inod2, inodbc, m;
    float sintheta, sinfi, cosfi, t, delx, dely, delz;

    int jseg, a, b, am;
    int ii, jj, kk, ll, mm, imin, imax, jmin, jmax, kmin, kmax;
    double aa, bbp, cc, abc, ddmin, ddmax, dd;// renamed bb here to bb1 since bb is declared a char earlier, I will have to change all instances of bb in the outboun() section of the code to bb1
    float ds2, de2, disp2, d, x11, x22, x33;
    float OR = 1000., IR = 200.;

    char fname[80];
    int showpoint, ilevel, nlevel = 100;
    float xmin, ymin, xs, ys, picfac, red, green, blue, xz, xzmin, xzmax;
    float diamfac = 1., zcoord, zbottom, ztop, zmin, zmax;
    float **cosp;

    FILE *ofp, *ofp1;



    float *diam, **cnode;

    int is, in;
    FILE *exelem, *exnode;

    int *paramtyp, **paramind, imain;

    int ix, iy, iz, jx, jy, jz, nt, ineg, ihigh, imaxerr;
    int ixdiff, iydiff, izdiff, isn, ktissue, kvessel, itp, jtp, convflag, convflagt, convflagv;
    int greensverbose = 1;

    float duration, rsegmax, dsmax, gvarfac;
    float gtt, gtv, gvt, dist, dtave, den, dqsumdg0;
    float dif, err, qhdcrit, pathlength;
    float lam, lam3d, lam2d;

    float pcr,m0,phi,kinstabA,kfA,krA,kAH,kinstabH,kfH,krH,kHP,kinstabM,kfM,krM,kMP,po2;

    int nod1, nod2,flag;

    float ***tissueLevels;

    float PO2max; // added Jan 27th, used to determine distance after which the strength of a source can be ignore.
    float cutoff_dist; // added Jan 27th, used to determine distance after which the strength of a source can be ignore.

    int **vt, **tv;
    int *size_vt, *size_tv;


