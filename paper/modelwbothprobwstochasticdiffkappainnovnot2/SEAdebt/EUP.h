//// USE OF OPTIONS ////
#define OMP 1  // (0 == no omp, 1 == omp on VFI, 2 == omp on params)
#define nbthread 2
#define CRS 0 // (1 = sobol, 0 =  no sobol)
// Choose which type of Borrowing constraint to use (Buera type (1) or CD type (2))
#define BRC 1
#define indexPM 0
#define indexPMwrite 0
#define debtadj 1
#define SEA 1
double taxdebtfixed = 0.011306767508575;
////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <limits>
//#include <iostream>
//#include <iomanip>
#include <assert.h>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <random>

//NR library static link includes
#include "/Users/alexandregaillard/Documents/Compiler/lib/library/nr/nr.h"
#include "/Users/alexandregaillard/Documents/Compiler/lib/library/nr/ran1.cpp"


// Useful math function
#define max(a,b) ((a)>(b))?(a):(b)
#define min(a,b) ((a)<(b))?(a):(b)
#define interpol(x,y,z) (y+(x-floor(x))*(z-y))
#define inter1d(x1,y1,y2) ((1.0-(x1))*(y1)+(x1)*(y2))
#define inter2d(x1,x2,y11,y21,y12,y22) ((1.0-(x2))*((1.0-(x1))*(y11)+(x1)*(y21))+(x2)*((1.0-(x1))*(y12)+(x1)*(y22)))


//////////////////
//// OMP PARA ////
//////////////////

#if OMP == 1 || OMP == 2
#include <omp.h>
#endif

////////////////
/// GRID DEF ///
////////////////
#define maxgrid 300
#define maxprod 3
#define maxfirmtype 2
#define ifulldim (maxgrid*maxprod)
#define ifulldimE (maxgrid*maxprod*maxfirmtype)

#if BRC == 1
#define ifulldimKK (maxgrid*maxfirmtype*maxprod)
#define inxKK(igridindex,ygridindex,firmtype) (((ygridindex)*(maxgrid)*(maxfirmtype))+(igridindex)*(maxfirmtype)+firmtype)
#endif

#if BRC == 2 
#define ifulldimKK (maxgrid*maxfirmtype*maxprod)
#define inxKK(igridindex,ygridindex,firmtype) (((ygridindex)*(maxgrid)*(maxfirmtype))+(igridindex)*(maxfirmtype)+firmtype)
#endif

//Macro switch to uniform grid
//Uncomment to use an uniform grid
//#define iUniformGrid 1
#define Gridmin 0.0
#define Gridmax 5000.0

#ifdef iUniformGrid
#define phi(x) ((((Gridmax*1.0)-(Gridmin*1.0))/(maxgrid-1))*(x))
#define phiinv(x) ((x)*((maxagrid-1)/(Gridmax-Gridmin)))
#else
const double Echelle1=1.0;
const double grmin=(Gridmin/Echelle1)-1.0;
const double Exponen=log((Gridmax/Echelle1)-grmin)/(maxgrid-1);
#define phi(x) ( Echelle1*(exp(Exponen*(x))+grmin) )
#define phiinv(x) (log((x)/Echelle1-grmin)/Exponen)
#endif

double grid[maxgrid];

//access index
#define inx(igridindex,ygridindex) (((ygridindex)*(maxgrid))+(igridindex))
#define inxE(igridindex,ygridindex,firmtype) (((ygridindex)*(maxgrid)*(maxfirmtype))+(igridindex)*(maxfirmtype)+firmtype)

// state here should be entrepreneur / worker / unemployed ST / LT (1,2,3,4)




///////////////////////
/// CALIBRATION DEF ///
///////////////////////

// equilibrium variables
double rstar;
double wstar;
double lstar;
double ustar; // possibly useless.
double taxrate;

double rstars = 0.07;
double lstars = 0.8919;
double taxrates = 0.1417;

// Estimated Parameters
double betapar = 0.869;      //  discount factor
double phiparW = 1.853593;        //  Effort elasticity parameter for job search
double phiparE = 2.064395;         //  Effort elasticity parameter for Entrepreneurship search
double kappaE = 0.56;
double kappaEU = 0.56;
double kappaW = 6.488173;
double costpar = 0.0;         //  Entry cost
double fpar = 0.847430;           //  Borrowing constraint parameter (could be fixed to 0.75)
double xipar = 1.0;         //  Probability to get an idea
double zetapar =0.0;       //  Probability to lose the idea
double mupar = 0.025;          //  bankrupcy rate.
double nupar =0.86041;           //  Return to scale in entrepreneurial sector
double pibad = 0.18776;

// REGRESSION FOR UNEMPLOYEMENT PROBABILITY, G(), AND EMPLOYEMENT PROBABILITY (kappa) //
double etapar[maxprod] = {0.05, 0.031, 0.021}; //  Layoff rate
double kappaWU[maxprod];  //  Job matching rate

double alphakappaw1 = 1.2;
double alphakappaw2 = 0.2;
double g1 = 0.481418;
double g2 = 0.526987;
double g3 = 0.714549;

double stateg[maxfirmtype] = {0.415344, 1.0};

double pig1 = pibad;
double pig2 = pibad;
double pig3 = pibad;

double pig11 = pibad;
double pig22 = pibad;
double pig33 = pibad;

//double alphamu1 = 0.018; // Cannot be fixed, because many entrepreneur shut down to be worker !!
//double alphamu2 = -0.002; // Cannot be fixed, because many entrepreneur shut down to be worker !!
double gtrans[maxprod];

double mgtrans[maxprod][maxfirmtype] = {
    {pig1, 1-pig1},
    {pig2, 1-pig2},
    {pig3, 1-pig3}
    };

double mggtrans[maxprod][maxfirmtype] = {
    {pig11, 1-pig11},
    {pig22, 1-pig22},
    {pig33, 1-pig33}
    };

double prod[maxprod] = {0.25, 1.0, 1.75};   //  productivity worker.

// COULD USE GINI HERE...
double mprod[maxprod][maxprod] = {
    {0.780000,   0.220000,   0.000000},
    {0.110000,   0.780000,   0.110000},
    {0.000000,   0.220000,   0.780000}
};

// Fixed Parameters
const double minpar = 0.005;         //  minimum subsistence level
const double sigmapar =1.39;        //  Intertemp. Elasticity of subsitution (CRRA)
const double rhostar =0.4;         //  Unemployment insurance rate
const double TFPpar =1.0;           //  TFP
const double deltapar =0.06;        //  depreciation rate
const double pLTpar =1.0;          //  Probability to switch to Long-run Unemployment
const double PBnoUIVE = pLTpar;          //  Probability to switch to no UI in Entrepreneurship.
const double alphapar =0.36;        //  return to scale parameter
const double corptax =0.0;          //  tax paid by entrepreneurs
//
//double prod[maxprod] ={0.24675246, 0.44731655, 0.76540823, 1.3096984, 2.3742408};
//double mprod[maxprod][maxprod] = {
//    {0.737571957533869,   0.247281329173195,   0.014989047949883,   0.000157529629473,   0.000000135713580},
//    {0.194679850961718,   0.555490252744122,   0.232791251149989,   0.016914625083558,   0.000124020060613},
//    {0.011257411090059,   0.222075921776607,   0.533333334266667,   0.222075921776607,   0.011257411090059},
//    {0.000124020060613,   0.016914625083558,   0.232791251149989,   0.555490252744122,   0.194679850961718},
//    {0.000000135713580,   0.000157529629473,   0.014989047949883,   0.247281329173195,   0.737571957533869}
//};


// Utility function
#define utilc(x) pow((x),(1.0-sigmapar*1.0))/(1.0-sigmapar*1.0)

//disutility of effort
#define disutilityW(x) (-pow((x),phiparW))
#define disutilityE(x) (-pow((x),phiparE))

//marginal disutility of effort
#define dsearchdesutilityW(x) (phiparW*(pow((x),(phiparW-1.0))))
#define dsearchdesutilityE(x) (phiparE*(pow((x),(phiparE-1.0))))




////////////////////////////////
/// NUMERICAL DEF AND VALUES ///
////////////////////////////////

//maximum search effort
#define EffortMaxW 100.0
#define EffortMaxE 100.0

//maximum endogenous borrowing constraint in the economy
#define amaxendoK (Gridmax/fpar)

//value function cnvg criterion
#define epsilonValue 0.000001
const int maxiterVF = 200;

//endogenous borrowing constraint cnvg criterion
#define epsilonendoK 0.000001
const int maxiterendoK = 300;

//stationary distribution cnvg criterion
#define epsilonDist 0.00000001

//aggregate capital cnvg criterion
#define epsilonprice 0.00001
#define itermax 200

//aggregate capital cnvg criterion
#define epsilonparams 0.000001

//endogenous borrowing constraint relaxation
#define relaxendoK 0.3

//aggregate capital relaxation
#define relaxK 0.05

//aggregate capital relaxation
#define relaxL 0.2

//aggregate tax relaxation
#define relaxTax 0.2

// TOLERANCE MNBRACK
const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;

// TOLERANCE MNGOLDEN
const double TOL=1.0e-6;
const double R=0.61803399,C=1.0-R;

//seed de ran1
int idum;
int cseed;



////////////////////////////////////////
/// DEFINE TEMP FILE AND OUTPUT FILE ///
////////////////////////////////////////

char tempfileout[80];//="results\tempfile_20141202@1632.out";
const char valuefnout[]="valuefn.out";
const char endokout[]="valuefnVE.out";
const char saveout[]="save.out";
const char searchout[]="search.out";
const char searchout2[]="search2D.out";
const char distout[]="dist.out";
const char distVEfile[]="distVE.out";
const char distasset2[]="distasset.out";
const char disttransitionWWtoVE[]="disttransitionWWtoVE.out";
const char disttransitionVEtoWW[]="disttransitionVEtoWW.out";
const char distnecessity[]="distnecessity.out";
const char moment[]="moment.out";
const char BRC1ENDO[]="BRC1ENDO.out";

const char endoKfile[]="project_method/endoK.out";
const char pricepathfile[]="project_method/pricepath.out";

const char distVE_0file[]="project_method/distVE_0.out";
const char distWW1_0file[]="project_method/distWW1_0.out";
const char distWW0_0file[]="project_method/distWW0_0.out";
const char distUS0_0file[]="project_method/distUS0_0.out";
const char distUS1_0file[]="project_method/distUS1_0.out";
const char distUL0_0file[]="project_method/distUL0_0.out";
const char distUL1_0file[]="project_method/distUL1_0.out";

const char distVE_Tfile[]="project_method/distVE_T.out";
const char distWW1_Tfile[]="project_method/distWW1_T.out";
const char distWW0_Tfile[]="project_method/distWW0_T.out";
const char distUS0_Tfile[]="project_method/distUS0_T.out";
const char distUS1_Tfile[]="project_method/distUS1_T.out";
const char distUL0_Tfile[]="project_method/distUL0_T.out";
const char distUL1_Tfile[]="project_method/distUL1_T.out";

const char valueVETfile[]="project_method/valueVET.out";
const char valueVEcTfile[]="project_method/valueVEcT.out";
const char valueWW1Tfile[]="project_method/valueWW1T.out";
const char valueWW0Tfile[]="project_method/valueWW0T.out";
const char valueUS0Tfile[]="project_method/valueUS0T.out";
const char valueUS1Tfile[]="project_method/valueUS1T.out";
const char valueUL0Tfile[]="project_method/valueUL0T.out";
const char valueUL1Tfile[]="project_method/valueUL1T.out";
const char valueVEwUIcTfile[]="project_method/valueVEwUIcTfile.out";
const char valueVEwUITfile[]="project_method/valueVEwUIcTfile.out";

//////////////////////////////////////////
///// CRS SEARCH FOR OPTIMAL PARAMS //////
//////////////////////////////////////////
#if CRS == 1
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/sobol.cpp" // very time compiling
#endif

#define	nbsobol 44
#define nbpara 14
#define nbmoments 16

// FILE OUT //
const char PARA[]="PARA.out";
const char SOBOL[]="SOBOL.out";

// OBSERVED MOMENTS ////
double obsmoments[nbmoments];

    
// COVARIANCE MATRIX //
double covmat[nbmoments][nbmoments]; // identity matrix.

