// Entrepreneurship, Unemployment and Paid Worker //
// Sumudu Kankanamge, Alexandre Gaillard ///////////

/********************************** OBJECTIVE **************************************/
/* our objective is to create a model which explain the observed transition between SE -> WW and WW -> SE.
What is new in this literature is the fact that we account for different level of wealth and ability which
are crucial to explain the observed transitions. Especially, rich people are more likely to be able to set up
a business. 
The second objective is to assess the Unemployment Insurance in a general equilibirum model with heterogeneity and
labor market frictions. We want to assess the impact of the UI insurance in the switching probability to be entrepreneur
or worker. There is a potential tradeoff. When increasing the UI insurance, the unemployed people could have some disencentive 
effect to find a new job, since searching induce cost to him. On the other side, increasing the UI insurance allows the 
unemployed people to get more wealth and can potentially increase entrepreneurship.
We can also assess the effect of subsidy entrepreneurship. Subsidy should increase the incentive to search for an entrepreneurship and
decrease the searching effort for a job. This should increase the labor demand from entrepreneur and therefore decrease unemployment and
rises wages. Which induce an higher searching effort for a job. The effect is thus ambiguous in equilibrium.
--> The unemployment insurance should have disproportional effect among classes. Different ability level should respond differently 
to a rise in insurance or subsidy. Subsidy should create a larger fraction of "necessity entrepreneurs".
************************************************************************************/



//////////////////////////////////////////////

// HEADER FILE //
#include "EUP.h"


///////////////////////////////////////////////////////////
// Useful function, structure and mathematical algorithms//
///////////////////////////////////////////////////////////
namespace {
	inline void shft2(double &a, double &b, const double c)
	{
		a=b;
		b=c;
	}
    
	inline void shft3(double &a, double &b, double &c, const double d)
	{
		a=b;
		b=c;
		c=d;
	}
}

using namespace std;


// PARAMETER STRUCTURE FOR VALUE FUNCTION MNGOLDEN AND ZBRENT
// PARAMS STRUCTURE FOR MAIN (ASSET LEVEL NEXT PERIOD) //
struct paramsgolden //mngolden
{
    int igridST, ygridST, egridST, typeU; // typeU define whether the unemployed is long (=1) or short (=0)
    // next time value function.
    double *valueUL1nextST,*valueUS1nextST, *valueUL0nextST,*valueUS0nextST, *valueVEnextST, *valueVEcnextST, *valueWW1nextST,*valueWW0nextST, *valueVEwUInextST, *valueVEwUIcnextST;
    // Search function
    double searchUWoutST, searchUEoutST, searchoutST;
    // Endogenous capital level para:
    double endoKST;
};

// FOR INDIVIDUAL SEARCHING EFFORT (only one) (WORKER (SE) OR UNEMPLOYED (SW) OR ENTREPRENEUR (SW)) //
struct searchEorWstruct //zbrent to look for the optimal search level
{
    double diffValueNextfoc;
    int ygridparXX;
};

// FOR ENTREPRENEURS // (two ways: CD type or Buera (2011) type)
#if BRC == 1
struct searchfendoKYY //zbrent to look for the capital level invested
{
    int egridKYY, igridKYY, ygridKYY;
};
#endif

// CD type //
#if BRC == 2
struct searchfendoKYY //zbrent to look for the capital level invested
{
    double *valueULKYY, *valueVEKYY;
    int lgridKYY, egridKYY, ygridKYY;
};
#endif

// FOR INDIVIDUAL WITH two searching effort (UNEMPLOYED) //
struct searchfocEforW //zbrent to look for the optimal search level
{
    double Psw;
};
struct searchfocWwhenE //zbrent to look for the optimal search level
{
    double tempWW1val, tempWW1VEval, tempWW0val, tempUS1val, tempUS0val, tempUL1val, tempUL0val,tempUL1VEval, tempUS1VEval, tempWW1VEwUIval, searchUsEoutST;
    int ygridparXX, typeUU;
};

////////


void bascule(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorOUT[i]=VectorIN[i];}
}


///// MY OPTIMIZATION FUNCTIONS //////
#include "/Users/alexandregaillard/Documents/Compiler/lib/optimfun/entunemployed/mngolden1search.cpp"      /// to find minimum ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/optimfun/entunemployed/mymngolden.cpp"      /// to find minimum ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/optimfun/entunemployed/mngolden2search.cpp"      /// to find minimum ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/optimfun/entunemployed/zbrentNEW.cpp"   /// to find minimum sequentially ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/optimfun/entunemployed/myCRSwithoutparams.cpp"           /// To find best params ///

///// DISTRIBUTION FUNCTIONS //////
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/gini.cpp"   /// to compute gini///
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/toppercent.cpp"   /// to see top percent wealth ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/medianworth.cpp"   /// to get median worth ///
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/SMMfun.cpp"   /// to compute gini///

// OTHER FUNCTIONS //
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/readinput.cpp"
//////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
// END Useful function, structure and mathematical algorithms//
///////////////////////////////////////////////////////////////




//********************************************//
//**********VALUE FUNCTION ITERATION**********//
//********************************************//

/////////////////////////////////
////  PROBABILITY FUNCTIONS  ////
/////////////////////////////////

double piW(int ygridfoc, double searchval) {
    return (1.0-exp(-kappaW*(searchval)));
}
double piE(double searchval) {
    return (1.0-exp(-kappaE*(searchval)));
}
double piWU(int ygridfoc, double searchval) {
    return (1.0-exp(-kappaWU[ygridfoc]*(searchval)));
}
double piEU(double searchval) {
    return (1.0-exp(-kappaEU*(searchval)));
}


//////////////////////////////////
///////// FONCTION SEARCH ////////
//////////////////////////////////

double searchEwithWgivenfun(const double seval, void * params)
{
    struct searchfocEforW *focparamsWforE= (struct searchfocEforW *) params;
    double Pswval=focparamsWforE->Psw;
    
    double residu;
    
    residu =  exp(-kappaEU*seval) - dsearchdesutilityE(seval)/(betapar*kappaEU*(1-zetapar)*Pswval);
    
    return residu;
}


double searchWfun(const double xval, void * params)
{
    struct searchEorWstruct *focparamsW= (struct searchEorWstruct *) params;
    double diffValueNext=focparamsW->diffValueNextfoc;
    int ygridfoc=focparamsW->ygridparXX;
    
    double residu;
    
    //residu = kappaW*exp(-kappaW*xval)*diffValueNext-dsearchdesutilityW(xval); // becareful you change dsearchdesutility
    // if want to make diff prob in function of ability
    residu = kappaW*exp(-kappaW*xval)*diffValueNext-dsearchdesutilityW(xval);
    
    return residu;
}


double searchWfunU(const double xval, void * params)
{
    struct searchEorWstruct *focparamsW= (struct searchEorWstruct *) params;
    double diffValueNext=focparamsW->diffValueNextfoc;
    int ygridfoc=focparamsW->ygridparXX;
    
    double residu;
    
    //residu = kappaW*exp(-kappaW*xval)*diffValueNext-dsearchdesutilityW(xval); // becareful you change dsearchdesutility
    // if want to make diff prob in function of ability
    residu = kappaWU[ygridfoc]*exp(-kappaWU[ygridfoc]*xval)*diffValueNext-dsearchdesutilityW(xval);
    
    return residu;
}


double searchEfun(const double xval, void * params)
{
    struct searchEorWstruct *focparamsE= (struct searchEorWstruct *) params;
    double diffValueNext=focparamsE->diffValueNextfoc;
    int ygridfoc=focparamsE->ygridparXX;
    
    double residu;
    
    residu = kappaE*exp(-kappaE*xval)*diffValueNext-dsearchdesutilityE(xval); // becareful you change dsearchdesutility
    
    return residu;
}


double searchWwhenEknownfun(const double swval, void * params)
{
    struct searchfocWwhenE *focparamsWwhenE= (struct searchfocWwhenE *) params;
    double tempWW0=focparamsWwhenE->tempWW0val;
    double tempWW1=focparamsWwhenE->tempWW1val;
    double tempWW1VE=focparamsWwhenE->tempWW1VEval;
    double tempUS0=focparamsWwhenE->tempUS0val;
    double tempUS1=focparamsWwhenE->tempUS1val;
    double tempUL0=focparamsWwhenE->tempUL0val;
    double tempUL1=focparamsWwhenE->tempUL1val;
    double tempUS1VE=focparamsWwhenE->tempUS1VEval;
    double tempUL1VE=focparamsWwhenE->tempUL1VEval;
    #if SEA == 1
    double tempWW1VEwUI=focparamsWwhenE->tempWW1VEwUIval;
    #endif
    int ygridfoc=focparamsWwhenE->ygridparXX;
    int typeUval=focparamsWwhenE->typeUU;
    
    double diffValueNext, residu, testval1, testval2, testval3, testval4, searcheffVE;
//    double VERIF1, VERIF2;
    
    int iverif;
    
    struct searchfocEforW paramsXXX;
    
    if(typeUval == 0) {
    #if SEA == 1
    paramsXXX.Psw = max(0.000000001,piW(ygridfoc, swval)*(PBnoUIVE*tempWW1VE + (1 - PBnoUIVE)*tempWW1VEwUI  - tempWW1) + (1 - piW(ygridfoc, swval))*(((1-pLTpar)*tempUS1VE + pLTpar*tempUL1VE) - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1))); // I think this is OK.
    #endif
    #if SEA == 0
    paramsXXX.Psw = max(0.000000001,piW(ygridfoc, swval)*(tempWW1VE - tempWW1) + (1 - piW(ygridfoc, swval))*(((1-pLTpar)*tempUS1VE + pLTpar*tempUL1VE) - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1)));
    #endif
    }
    if(typeUval == 1) {
    paramsXXX.Psw = max(0.000000001,piW(ygridfoc, swval)*(tempWW1VE - tempWW1) + (1 - piW(ygridfoc, swval))*(tempUL1VE - tempUL1));
    }
    
    if(paramsXXX.Psw <= -0.0000001) {
    printf("There is a bug in US1");
    }
    
    testval1 = searchEwithWgivenfun(0.00001, &paramsXXX);
    testval2 = searchEwithWgivenfun(0.000011, &paramsXXX);
    testval3 = searchEwithWgivenfun(EffortMaxE, &paramsXXX);
    testval4 = searchEwithWgivenfun(EffortMaxE-0.00001, &paramsXXX);
    
    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffVE = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffVE = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffVE = EffortMaxE;
        printf("There is a bug:: EffortmaxE 3 riched %f %f %f\n", testval3, testval4, paramsXXX.Psw); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffVE = EffortMaxE;
        printf("There is a bug:: EffortmaxE 4 riched %f %f %f\n", testval3, testval4, paramsXXX.Psw); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffVE = zbrentNEW(searchEwithWgivenfun,0.0,EffortMaxE,(1.0e-10),&paramsXXX);
        }
    
//    if(paramsXXX.Psw >0.000001) {
//      printf("%f %f %f %f %f \n", testval1, testval2, swval, paramsXXX.Psw, searcheffVE); getchar();
//    }

    focparamsWwhenE->searchUsEoutST=searcheffVE;
    
    // Given the best SEval, find the best SWval function //
    if(typeUval == 0) {
    #if SEA == 1
    diffValueNext = max(0.000000000001,(betapar*((1-zetapar)*((1 - piEU(searcheffVE))*(tempWW1 - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1))+piEU(searcheffVE)*(pLTpar*tempWW1VE + (1-pLTpar)*tempWW1VEwUI - ((1-pLTpar)*tempUS1VE + pLTpar*tempUL1VE))) + zetapar*(tempWW0 - ((1-pLTpar)*tempUS0 + pLTpar*tempUL0)))));
    #endif
    #if SEA == 0
    diffValueNext = max(0.000000000001,(betapar*((1-zetapar)*((1 - piEU(searcheffVE))*(tempWW1 - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1))+piEU(searcheffVE)*(tempWW1VE - ((1-pLTpar)*tempUS1VE + pLTpar*tempUL1VE))) + zetapar*(tempWW0 - ((1-pLTpar)*tempUS0 + pLTpar*tempUL0)))));
    #endif
    }
    if(typeUval == 1) {
    diffValueNext = max(0.000000000001,(betapar*((1-zetapar)*((1 - piEU(searcheffVE))*(tempWW1 - tempUL1)+piEU(searcheffVE)*(tempWW1VE - tempUL1VE)) + zetapar*(tempWW0 - tempUL0))));
    }
    
//    
//    if(VERIF1 < -0.000001 & VERIF2 < -0.000001 & iverif == 0) {
//        printf("Verif1: %20.15f, Verif2: %20.15f, diffValueNext: %20.15f, For 0.0: Pw %20.15f Type: %d, SearchVE: %20.15f",VERIF1, VERIF2, diffValueNext, paramsXXX.Psw, typeUval, searcheffVE); getchar();
//    }
    
    residu = kappaWU[ygridfoc]*exp(-kappaWU[ygridfoc]*swval)*diffValueNext-dsearchdesutilityW(swval);
    
    return residu;
}

////////////////////////////////////////////////
///////// FONCTION BORROWING CONSTRAINT ////////
////////////////////////////////////////////////

#if BRC == 1 // Buera (2011) type //
double searchendoK(const double kval, void * params)
{
    struct searchfendoKYY *focparamsK = (struct searchfendoKYY *) params;
    int igridfoc=focparamsK->igridKYY;
    int egridfoc=focparamsK->egridKYY;
    int ygridfoc=focparamsK->ygridKYY;
    
    double residu;
    
    residu =  (1 - corptax)*(TFPpar*gtrans[ygridfoc]*stateg[egridfoc]*pow(kval, nupar) + (1 - deltapar)*kval - (1 + rstar)*(kval - grid[igridfoc])) - fpar*(1 - corptax)*(TFPpar*gtrans[ygridfoc]*stateg[egridfoc]*pow(kval, nupar) + (1 - deltapar)*kval);
    
    return residu;
}
#endif

#if BRC == 2
double searchendoK(const double kval, void * params)
{
    struct searchfendoKYY *searchendoKparamsXX= (struct searchfendoKYY *) params;
    double *valueULendoK=searchendoKparamsXX->valueULKYY;
    double *valueVEendoK=searchendoKparamsXX->valueVEKYY;
    int ygridendoK=searchendoKparamsXX->ygridKYY;
    int egridendoK=searchendoKparamsXX->egridKYY;
    int igridendoK=searchendoKparamsXX->lgridKYY;
    
	double khatnextix, khatnext, dkhatnext, ValULcheck, functionBRC;
    int khatnextixongrid;
    
	khatnextix = phiinv(fpar*kval);
	 
	khatnextixongrid = min((maxgrid - 1), (int)(floor(khatnextix)));
	khatnextixongrid = max(khatnextixongrid, 0);
	
	if(khatnextixongrid >= (maxgrid - 1)) {
		khatnextix = (maxgrid - 1) - 0.000000001;
		khatnextixongrid = min((maxgrid - 1), (int)(floor(khatnextix)));
		khatnextixongrid = max(khatnextixongrid, 0);
	}
	
	if(khatnextixongrid <= 0) {
		khatnextix = 0.000000001;
		khatnextixongrid = min((maxgrid - 1), (int)(floor(khatnextix)));
		khatnextixongrid = max(khatnextixongrid, 0);
	}
	
	khatnext = phi(khatnextix);
	
	dkhatnext = (khatnext - grid[khatnextixongrid])/(grid[khatnextixongrid + 1] - grid[khatnextixongrid]);
	
	
    // NOTE THAT RUN AWAY IS CONTEMPORANEROUS (SO YOU DON'T HAVE TO PUT ANY STOCHASTICITY HERE
	ValULcheck = inter1d(dkhatnext, valueULendoK[inx(khatnextixongrid,ygridendoK)], valueULendoK[inx(khatnextixongrid + 1,ygridendoK)]);
	
	functionBRC = ValULcheck - valueVEendoK[inxE(igridendoK,ygridendoK,egridendoK)];
	
	return functionBRC;
}
#endif



///////////////////////////////////////////////////
// FUNCTION CONSUMPTION (OR MAX WEALTH) FUNCTION //
///////////////////////////////////////////////////

double findconsWW(const double xval, const int igridfoc, const int ygridfoc)
{
	return ((1.0+rstar)*grid[igridfoc]+prod[ygridfoc]*wstar*(1.0-taxrate)-xval);
}

double findconsUS(const double xval, const int igridfoc, const int ygridfoc)
{
	return ((1.0+rstar)*grid[igridfoc]+prod[ygridfoc]*rhostar*wstar*(1.0-taxrate)-xval);
}

double findconsUL(const double xval, const int igridfoc)
{
	return ((1.0+rstar)*grid[igridfoc] + minpar*(1.0-taxrate) - xval);
}

double findconsVE(const double xval, const int igridfoc, const int egridfoc, const int ygridfoc, const double endoKval)
{
	return ((1.0 - corptax)*(TFPpar*gtrans[ygridfoc]*stateg[egridfoc]*pow(endoKval, nupar) - deltapar*endoKval - rstar*endoKval) + (1.0+rstar)*grid[igridfoc] - xval);
}

double findconsVEwUI(const double xval, const int igridfoc, const int egridfoc, const int ygridfoc, const double endoKval)
{
    double wealth, UI, entrepinc;
    
    UI = ((1.0+rstar)*grid[igridfoc]+prod[ygridfoc]*rhostar*wstar*(1.0-taxrate)-xval);
    entrepinc = ((1.0 - corptax)*(TFPpar*gtrans[ygridfoc]*stateg[egridfoc]*pow(endoKval, nupar) - deltapar*endoKval - rstar*endoKval) + (1.0+rstar)*grid[igridfoc] - xval);
    
    if(UI > entrepinc) {
        if(entrepinc > 0.0) {
            wealth = UI;
        }
    }
    
    if(UI <= entrepinc) {
        wealth = entrepinc;
    }
    
	return (wealth);
}
///////////////////////////////////////
// END CONS (OR MAX WEALTH) FUNCTION //
///////////////////////////////////////




/////////////////////
// VALUE FUNCTIONS //
/////////////////////

// WORKER //

// WITHOUT IDEA //
double focWWnoID(const double xval, void * params)
{
    struct paramsgolden *focparamsWW= (struct paramsgolden *) params;
    double *valueWW1nextfocWW=focparamsWW->valueWW1nextST;
    double *valueWW0nextfocWW=focparamsWW->valueWW0nextST;
    double *valueUS1nextfocWW=focparamsWW->valueUS1nextST;
    double *valueUS0nextfocWW=focparamsWW->valueUS0nextST;
    int igridfoc=focparamsWW->igridST;
    int ygridfoc=focparamsWW->ygridST;
    
    double cons,xgrid,dxgrid,valfoc,tempWW0,tempUS0,tempWW1,tempUS1;
    
    int ixgrid, index1, index2;

    cons = findconsWW(xval, igridfoc, ygridfoc);
    
    
    // ON THE GRID REDEFINITION //
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    // VERIFICATIONS //
    if ((xval>(Gridmax))){printf("focWW: (xval>(Gridmax)) %20.15f\n",xval);getchar();}
    if ((ixgrid>=(maxgrid-1))){printf("focWW: ixgrid>=(maxgrid-1) %d\n",ixgrid);getchar();}
    if ((cons<=0.0)){printf("focWW: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);getchar();}
    if ((xval<(0.0))){printf("focWW: (xval<(0.0)) %20.15f\n",xval);getchar();}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);
    
    tempUS0 = 0.0;
    tempWW0 = 0.0;
    tempUS1 = 0.0;
    tempWW1 = 0.0;
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempWW1 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueWW1nextfocWW[index1]),(valueWW1nextfocWW[index2]));
        tempUS1 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS1nextfocWW[index1]),(valueUS1nextfocWW[index2]));
        tempWW0 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueWW0nextfocWW[index1]),(valueWW0nextfocWW[index2]));
        tempUS0 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS0nextfocWW[index1]),(valueUS0nextfocWW[index2]));
   }

    valfoc=utilc(cons)+betapar*(xipar*((1.0-etapar[ygridfoc])*(tempWW1)+etapar[ygridfoc]*tempUS1)+(1-xipar)*(etapar[ygridfoc]*tempUS0 + (1.0 - etapar[ygridfoc])*tempWW0));
    
    valfoc=-valfoc;
    
    return valfoc;
}


// WITH AN IDEA //
double focWWwID(const double xval, void * params)
{
    struct paramsgolden *focparamsWW= (struct paramsgolden *) params;
    double *valueWW0nextfocWW=focparamsWW->valueWW0nextST;
    double *valueWW1nextfocWW=focparamsWW->valueWW1nextST;
    double *valueUS0nextfocWW=focparamsWW->valueUS0nextST;
    double *valueUS1nextfocWW=focparamsWW->valueUS1nextST;
    double *valueVEcnextfocWW=focparamsWW->valueVEcnextST;
    double *valueVEwUIcnextfocWW=focparamsWW->valueVEwUIcnextST;
    int igridfoc=focparamsWW->igridST;
    int ygridfoc=focparamsWW->ygridST;
    
    double cons,xgrid,dxgrid,valfoc,tempWW0,tempUS0,tempWW1,tempUS1, tempWW1VE, tempUS1VE, tempVEwUI, tempVE, searcheffVE, testval1, testval2, testval3, testval4;
    
    int ixgrid, index1, index2, index3, index4, iverif;

    cons = findconsWW(xval, igridfoc, ygridfoc);
    
    
    // ON THE GRID REDEFINITION //
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    // VERIFICATIONS //
    if ((xval>(Gridmax))){printf("focWW: (xval>(Gridmax)) %20.15f\n",xval);getchar();}
    if ((ixgrid>=(maxgrid-1))){printf("focWW: ixgrid>=(maxgrid-1) %d\n",ixgrid);getchar();}
    if ((cons<=0.0)){printf("focWW: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);getchar();}
    if ((xval<(0.0))){printf("focWW: (xval<(0.0)) %20.15f\n",xval);getchar();}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);
    
    tempUS0 = 0.0;
    tempWW0 = 0.0;
    tempUS1 = 0.0;
    tempWW1 = 0.0;
    tempVE = 0.0;
    tempWW1VE = 0.0;
    tempUS1VE = 0.0;
    
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempWW1 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueWW1nextfocWW[index1]),(valueWW1nextfocWW[index2]));
        tempUS1 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS1nextfocWW[index1]),(valueUS1nextfocWW[index2]));
        tempWW0 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueWW0nextfocWW[index1]),(valueWW0nextfocWW[index2]));
        tempUS0 += mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS0nextfocWW[index1]),(valueUS0nextfocWW[index2]));
        
        // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
        tempVE = 0.0;
        for(int e=0; e<maxfirmtype; e++) {
            index3 = inxE(ixgrid, y, e);
            index4 = inxE((ixgrid+1), y, e);
            
            tempVE += mgtrans[y][e]*inter1d(dxgrid,(valueVEcnextfocWW[index3]),(valueVEcnextfocWW[index4]));
        }
        
        tempVEwUI = 0.0;
        for(int e=0; e<maxfirmtype; e++) {
            index3 = inxE(ixgrid, y, e);
            index4 = inxE((ixgrid+1), y, e);
            
            tempVEwUI += mgtrans[y][e]*inter1d(dxgrid,(valueVEwUIcnextfocWW[index3]),(valueVEwUIcnextfocWW[index4]));
        }
        
        tempWW1VE += max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueWW1nextfocWW[index1]),(valueWW1nextfocWW[index2])),mprod[ygridfoc][y]*tempVE);
        
        #if SEA == 1
        tempUS1VE += max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS1nextfocWW[index1]),(valueUS1nextfocWW[index2])),mprod[ygridfoc][y]*tempVEwUI);
        #endif
        
        #if SEA == 0
        tempUS1VE += max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueUS1nextfocWW[index1]),(valueUS1nextfocWW[index2])),mprod[ygridfoc][y]*tempVE);
        #endif
   }
   
   struct searchEorWstruct paramsXX;
    
    // FIND THE OPTIMAL EFFORT LEVEL
    paramsXX.diffValueNextfoc = max(0.000000000001,(betapar*((1-zetapar)*((1-etapar[ygridfoc])*(tempWW1VE - tempWW1) + etapar[ygridfoc]*(tempUS1VE - tempUS1)))));
    paramsXX.ygridparXX = ygridfoc;
    
    if((betapar*((1-zetapar)*((1-etapar[ygridfoc])*(tempWW1VE - tempWW1) + etapar[ygridfoc]*(tempUS1VE - tempUS1)))) <= -0.0000001) {
    printf("There is a bug in WW"); getchar();
    }
    
    // Bracket the minimum (not sure to keep it) ---->
    testval1 = searchEfun(0.0, &paramsXX);
    testval2 = searchEfun(0.0001, &paramsXX);
    testval3 = searchEfun(EffortMaxE, &paramsXX);
    testval4 = searchEfun(EffortMaxE-0.0001, &paramsXX);
    
    iverif = 0;

    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffVE = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffVE = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffVE = EffortMaxE;
        printf("There is a bug:: EffortmaxE reached \n"); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffVE = EffortMaxE;
        printf("There is a bug:: EffortmaxE reached \n"); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffVE=zbrentNEW(searchEfun,0.0,EffortMaxE,(1.0e-10),&paramsXX);
        }
    
    focparamsWW->searchoutST=searcheffVE;
    
    // COMPUTE THE CONTINUATION VALUE
    
    valfoc=utilc(cons)+disutilityE(searcheffVE)+betapar*((1-zetapar)*(piE(searcheffVE)*(1 - etapar[ygridfoc])*tempWW1VE + etapar[ygridfoc]*(1 - piE(searcheffVE))*tempUS1 + etapar[ygridfoc]*piE(searcheffVE)*tempUS1VE + (1 - etapar[ygridfoc])*(1 - piE(searcheffVE))*tempWW1) + zetapar*(etapar[ygridfoc]*tempUS0 + (1.0 - etapar[ygridfoc])*tempWW0));
    
    valfoc=-valfoc;
    
    return valfoc;
}


// UNEMPLOYED (type is already included) //

// WITHOUT IDEA //
double focUnoID(const double xval, void * params)
{
    struct paramsgolden *focparamsU= (struct paramsgolden *) params;
    double *valueWW1nextfocU=focparamsU->valueWW1nextST;
    double *valueUS1nextfocU=focparamsU->valueUS1nextST;
    double *valueUL1nextfocU=focparamsU->valueUL1nextST;
    double *valueWW0nextfocU=focparamsU->valueWW0nextST;
    double *valueUS0nextfocU=focparamsU->valueUS0nextST;
    double *valueUL0nextfocU=focparamsU->valueUL0nextST;
    int igridfoc=focparamsU->igridST;
    int ygridfoc=focparamsU->ygridST;
    int typeUfoc=focparamsU->typeU;
    
    double cons,xgrid,dxgrid,tempWW1, tempUS1, tempWW0, tempUS0, tempUL1, tempUL0, valfoc,searcheffWW, testval1, testval2, testval3, testval4;
    
    int ixgrid, index1, index2, iverif;
    
    if(typeUfoc == 0) {
        cons = findconsUS(xval, igridfoc, ygridfoc); // short-run unemployed
    }
    
    if(typeUfoc == 1) {
        cons = findconsUL(xval, igridfoc); // long-run unemployed
    }
    
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if ((xval>(Gridmax))){printf("focU: (xval>(Gridmax)) %20.15f\n",xval);getchar();}
    if ((ixgrid>=(maxgrid-1))){printf("focU: ixgrid>=(maxgrid-1) %d\n",ixgrid);getchar();}
    if ((cons<0.0)){printf("focU: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);getchar();}
    if ((xval<(0.0))){printf("focU: (xval<(0.0)) %20.15f\n",xval);getchar();}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);
    
    
    tempWW1 = 0.0;
    tempWW0 = 0.0;
    tempUS1 = 0.0;
    tempUS0 = 0.0;
    tempUL1 = 0.0;
    tempUL0 = 0.0;
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempWW1 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2])));
        tempWW0 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW0nextfocU[index1]),(valueWW0nextfocU[index2])));
        tempUS1 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUS1nextfocU[index1]),(valueUS1nextfocU[index2])));
        tempUS0 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUS0nextfocU[index1]),(valueUS0nextfocU[index2])));
        tempUL1 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL1nextfocU[index1]),(valueUL1nextfocU[index2])));
        tempUL0 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL0nextfocU[index1]),(valueUL0nextfocU[index2])));
    }
    
    
    struct searchEorWstruct paramsXX;
    
    // FIND THE OPTIMAL EFFORT LEVEL
    if(typeUfoc == 0) {
    paramsXX.diffValueNextfoc = max(0.000000000001,(betapar*(xipar*(tempWW1 - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1)) + (1 - xipar)*(tempWW0 - ((1-pLTpar)*tempUS0 + pLTpar*tempUL0)))));
    }
    if(typeUfoc == 1) {
    paramsXX.diffValueNextfoc = max(0.000000000001,(betapar*(xipar*(tempWW1 - tempUL1) + (1 - xipar)*(tempWW0 - tempUL0))));
    }
    paramsXX.ygridparXX = ygridfoc;
    
    if((betapar*(xipar*(tempWW1 - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1)) + (1 - xipar)*(tempWW0 - ((1-pLTpar)*tempUS0 + pLTpar*tempUL0)))) <= -0.0000001) {
    printf("There is a bug in US"); getchar();
    }
    
    // Bracket the minimum (not sure to keep it) ---->
    testval1 = searchWfunU(0.0, &paramsXX);
    testval2 = searchWfunU(0.0001, &paramsXX);
    testval3 = searchWfunU(EffortMaxW, &paramsXX);
    testval4 = searchWfunU(EffortMaxW-0.0001, &paramsXX);
    
    iverif = 0;

    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffWW=zbrentNEW(searchWfunU,0.0,EffortMaxW,(1.0e-10),&paramsXX);
        }
    
    focparamsU->searchoutST=searcheffWW;
    
    
    // COMPUTE THE CONTINUATION VALUE
    if(typeUfoc == 0) {
    valfoc=utilc(cons)+disutilityW(searcheffWW)+betapar*(xipar*((1-piWU(ygridfoc, searcheffWW))*((1-pLTpar)*tempUS1 + pLTpar*tempUL1) + piWU(ygridfoc, searcheffWW)*tempWW1) + (1 - xipar)*(1 - piWU(ygridfoc, searcheffWW))*((1-pLTpar)*tempUS0 + pLTpar*tempUL0)+piWU(ygridfoc, searcheffWW)*tempWW0);
    }
    if(typeUfoc == 1) {
    valfoc=utilc(cons)+disutilityW(searcheffWW)+betapar*(xipar*((1-piWU(ygridfoc, searcheffWW))*tempUL1 + piWU(ygridfoc, searcheffWW)*tempWW1) + (1 - xipar)*(1 - piWU(ygridfoc, searcheffWW))*tempUL0+piWU(ygridfoc, searcheffWW)*tempWW0);
    }
    
    valfoc=-valfoc;
    
    //printf("focUS %d\t%20.15f\t%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igridfoc,xval,ixgrid,dxgrid,cons,utilc(cons),(betapar*(exitrate(searcheffUS))*(tempW)+betapar*(1.0-exitrate(searcheffUS))*tempU),valfoc,searcheffUS);
    
    
    return valfoc;

}


// WITH AN IDEA //
double focUwID(const double xval, void * params)
{
    struct paramsgolden *focparamsU= (struct paramsgolden *) params;
    double *valueWW1nextfocU=focparamsU->valueWW1nextST;
    double *valueUS1nextfocU=focparamsU->valueUS1nextST;
    double *valueUL1nextfocU=focparamsU->valueUL1nextST;
    double *valueWW0nextfocU=focparamsU->valueWW0nextST;
    double *valueUS0nextfocU=focparamsU->valueUS0nextST;
    double *valueUL0nextfocU=focparamsU->valueUL0nextST;
    double *valueVEcnextfocU=focparamsU->valueVEcnextST;
    #if SEA == 1
    double *valueVEwUIcnextfocU=focparamsU->valueVEwUIcnextST;
    #endif
    int igridfoc=focparamsU->igridST;
    int ygridfoc=focparamsU->ygridST;
    int typeUfoc=focparamsU->typeU;
    
    double cons,xgrid,dxgrid,tempWW1, tempWW0, tempWW1VE, tempUS1, tempUS0, tempUL1, tempUL0, tempVE, tempUS1VE, tempUL1VE, tempVEwUI, tempWW1VEwUI, valfoc, searcheffWW, searcheffVE, testval1, testval2, testval3, testval4;
    
    int ixgrid, index1, index2, index3, index4, iverif;
    
    if(typeUfoc == 0) {
        cons = findconsUS(xval, igridfoc, ygridfoc); // short-run unemployed
    }
    
    if(typeUfoc == 1) {
        cons = findconsUL(xval, igridfoc); // long-run unemployed
    }
    
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if ((xval>(Gridmax))){printf("focU: (xval>(Gridmax)) %20.15f\n",xval);getchar();}
    if ((ixgrid>=(maxgrid-1))){printf("focU: ixgrid>=(maxgrid-1) %d\n",ixgrid);getchar();}
    if ((cons<0.0)){printf("focU: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);getchar();}
    if ((xval<(0.0))){printf("focU: (xval<(0.0)) %20.15f\n",xval);getchar();}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);
    
    
    tempWW1 = 0.0;
    tempWW0 = 0.0;
    tempWW1VE = 0.0;
    tempWW1VEwUI = 0.0;
    tempUS1 = 0.0;
    tempUS0 = 0.0;
    tempUL1 = 0.0;
    tempUL0 = 0.0;
    tempUS1VE = 0.0;
    tempUL1VE = 0.0;
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempWW1 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2])));
        tempWW0 += mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW0nextfocU[index1]),(valueWW0nextfocU[index2])));

        tempUS1 += mprod[ygridfoc][y]*(inter1d(dxgrid,valueUS1nextfocU[index1],valueUS1nextfocU[index2]));
        tempUS0 += mprod[ygridfoc][y]*(inter1d(dxgrid,valueUS0nextfocU[index1],valueUS0nextfocU[index2]));
        tempUL1 += mprod[ygridfoc][y]*(inter1d(dxgrid,valueUL1nextfocU[index1],valueUL1nextfocU[index2]));
        tempUL0 += mprod[ygridfoc][y]*(inter1d(dxgrid,valueUL0nextfocU[index1],valueUL0nextfocU[index2]));
        
        tempVE = 0.0;
        for(int e=0; e<maxfirmtype; e++) {
            index3 = inxE(ixgrid, y, e);
            index4 = inxE((ixgrid+1), y, e);
            
            tempVE += mgtrans[y][e]*inter1d(dxgrid,(valueVEcnextfocU[index3]),(valueVEcnextfocU[index4]));
        }
        
        tempUL1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,valueUL1nextfocU[index1],valueUL1nextfocU[index2])),mprod[ygridfoc][y]*tempVE);
        
        // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
        if(typeUfoc == 0) {
        
            #if SEA == 0
            tempWW1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2]))),mprod[ygridfoc][y]*tempVE);
            tempUS1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,valueUS1nextfocU[index1],valueUS1nextfocU[index2])),mprod[ygridfoc][y]*tempVE);
            #endif
            
            #if SEA == 1
            tempVEwUI = 0.0;
            for(int e=0; e<maxfirmtype; e++) {
                index3 = inxE(ixgrid, y, e);
                index4 = inxE((ixgrid+1), y, e);
                
                tempVEwUI += mgtrans[y][e]*inter1d(dxgrid,(valueVEwUIcnextfocU[index3]),(valueVEwUIcnextfocU[index4]));
            }
            
            tempWW1VEwUI += max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2]))),mprod[ygridfoc][y]*tempVEwUI); //++= here there is something to correct
            // because you should have a probability PBnoUIVE to switch to tempWW1VE and (1 - PBnoUIVE)*tempWW1VEwUI. You don't respect it here.
            tempWW1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2]))),mprod[ygridfoc][y]*tempVE);
            tempUS1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,valueUS1nextfocU[index1],valueUS1nextfocU[index2])),mprod[ygridfoc][y]*tempVEwUI);
            #endif
            
        }
        
        if(typeUfoc == 1) {
            tempWW1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocU[index1]),(valueWW1nextfocU[index2]))),mprod[ygridfoc][y]*tempVE);
            tempUS1VE += max(mprod[ygridfoc][y]*(inter1d(dxgrid,valueUS1nextfocU[index1],valueUS1nextfocU[index2])),mprod[ygridfoc][y]*tempVE); // useless, cannot come back to US1 but still...
        }
        

    }
    
    
    // STARTING 2D searching effort
    
    
    struct searchfocWwhenE paramsXX;
    paramsXX.tempWW1val = tempWW1;
    paramsXX.tempWW0val = tempWW0;
    paramsXX.tempWW1VEval = tempWW1VE;
    #if SEA == 1
    paramsXX.tempWW1VEwUIval = tempWW1VEwUI;
    #endif
    paramsXX.tempUS1val = tempUS1;
    paramsXX.tempUS0val = tempUS0;
    paramsXX.tempUL1val = tempUL1;
    paramsXX.tempUL0val = tempUL0;
    paramsXX.tempUS1VEval = tempUS1VE;
    paramsXX.tempUL1VEval = tempUL1VE;
    paramsXX.ygridparXX = ygridfoc;
    paramsXX.typeUU = typeUfoc;
    
    testval1 = searchWwhenEknownfun(0.0, &paramsXX);
    testval2 = searchWwhenEknownfun(0.0001, &paramsXX);
    testval3 = searchWwhenEknownfun(EffortMaxW, &paramsXX);
    testval4 = searchWwhenEknownfun(EffortMaxW-0.0001, &paramsXX);
    
    iverif = 0;

    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffWW = zbrentNEW(searchWwhenEknownfun,0.0,EffortMaxW,(1.0e-10),&paramsXX);
        searcheffVE = paramsXX.searchUsEoutST;
        }
    
    focparamsU->searchUWoutST=searcheffWW;
    focparamsU->searchUEoutST=searcheffVE;
    
 
    // COMPUTE THE CONTINUATION VALUE
    if(typeUfoc == 0) {
    #if SEA == 1
    valfoc=utilc(cons)+disutilityW(searcheffWW)+disutilityE(searcheffVE)+betapar*((1-zetapar)*(piWU(ygridfoc,searcheffWW)*(1 - piEU(searcheffVE))*tempWW1+piEU(searcheffVE)*(1-piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS1VE + pLTpar*tempUL1VE) + piWU(ygridfoc, searcheffWW)*piEU(searcheffVE)*(PBnoUIVE*tempWW1VE + (1-PBnoUIVE)*tempWW1VEwUI) + (1 - piEU(searcheffVE))*(1 - piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS1 + pLTpar*tempUL1)) + zetapar*((1 - piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS0 + pLTpar*tempUL0) + piWU(ygridfoc, searcheffWW)*tempWW0));
    #endif
    
    #if SEA == 0
    valfoc=utilc(cons)+disutilityW(searcheffWW)+disutilityE(searcheffVE)+betapar*((1-zetapar)*(piWU(ygridfoc,searcheffWW)*(1 - piEU(searcheffVE))*tempWW1+piEU(searcheffVE)*(1-piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS1VE + pLTpar*tempUL1VE) + piWU(ygridfoc, searcheffWW)*piEU(searcheffVE)*(tempWW1VE) + (1 - piEU(searcheffVE))*(1 - piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS1 + pLTpar*tempUL1)) + zetapar*((1 - piWU(ygridfoc, searcheffWW))*((1 - pLTpar)*tempUS0 + pLTpar*tempUL0) + piWU(ygridfoc, searcheffWW)*tempWW0));
    #endif
    }
    if(typeUfoc == 1) {
    valfoc=utilc(cons)+disutilityW(searcheffWW)+disutilityE(searcheffVE)+betapar*((1-zetapar)*(piWU(ygridfoc,searcheffWW)*(1 - piEU(searcheffVE))*tempWW1+piEU(searcheffVE)*(1-piWU(ygridfoc, searcheffWW))*tempUL1VE + piWU(ygridfoc, searcheffWW)*piEU(searcheffVE)*tempWW1VE + (1 - piEU(searcheffVE))*(1 - piWU(ygridfoc, searcheffWW))*tempUL1) + zetapar*((1 - piWU(ygridfoc, searcheffWW))*tempUL0 + piWU(ygridfoc, searcheffWW)*tempWW0));
    }

    valfoc=-valfoc;
    
    //printf("focUS %d\t%20.15f\t%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igridfoc,xval,ixgrid,dxgrid,cons,utilc(cons),(betapar*(exitrate(searcheffUS))*(tempW)+betapar*(1.0-exitrate(searcheffUS))*tempU),valfoc,searcheffUS);
    
    
    return valfoc;

}


// ENTREPRENEUR

// value function of an entrepreneur (with level of investment $endoKfoc).
double focVE(const double xval, void * params)
{
    struct paramsgolden *focparamsVE= (struct paramsgolden *) params;
    double *valueVEnextfocVE=focparamsVE->valueVEnextST; // no cost when already employed.
    double *valueWW0nextfocVE=focparamsVE->valueWW0nextST;
    double *valueUL0nextfocVE=focparamsVE->valueUL0nextST;
    double *valueWW1nextfocVE=focparamsVE->valueWW1nextST;
    double *valueUL1nextfocVE=focparamsVE->valueUL1nextST;
    double endoKfoc=focparamsVE->endoKST;
    int ygridfoc=focparamsVE->ygridST;
    int igridfoc=focparamsVE->igridST;
    int egridfoc=focparamsVE->egridST;
    
    double cons,xgrid,dxgrid,valfoc,tempWW0, tempUL0, tempWW1VE, tempUL1VE, testval1, testval2, testval3, testval4, searcheffWW, tempVE;
    int iverif, ixgrid, ibug, index1, index2, index3, index4;
    
    ibug=0;
    
    cons = findconsVE(xval, igridfoc, egridfoc, ygridfoc, endoKfoc);

    //printf("%f %f %d %d %f", cons, xval, igridfoc, egridfoc, endoKfoc);
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if ((xval>(Gridmax))){printf("focVE: (xval>(Gridmax)) %20.15f\n",xval);ibug=1;}
    if ((ixgrid>=(maxgrid-1))){printf("focVE: ixgrid>=(maxgrid-1) %d\n",ixgrid);ibug=1;}
    if ((cons<0.0)){printf("focVE: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);ibug=1;}
    if ((xval<(0.0))){printf("focVE: (xval<(0.0)) %20.15f\n",xval);ibug=1;}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);

    tempWW1VE = 0.0;
    tempWW0 = 0.0;
    tempUL1VE = 0.0;
    tempUL0 = 0.0;
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempUL0+=mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL0nextfocVE[index1]),(valueUL0nextfocVE[index2])));
        tempWW0+=mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW0nextfocVE[index1]),(valueWW0nextfocVE[index2])));

        // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
        if(egridfoc == 1) {
            tempVE = 0.0;
            for(int e=0; e<maxfirmtype; e++) {
                index3 = inxE(ixgrid, y, e);
                index4 = inxE((ixgrid+1), y, e);
                
                tempWW1VE+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])));
                tempUL1VE+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])), mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
            }
        }
        
        if(egridfoc == 0) {
            index3 = inxE(ixgrid, y, 0);
            index4 = inxE((ixgrid+1), y, 0);
            
            tempWW1VE+=max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])));
            tempUL1VE+=max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])), mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
        }
        
    }
    
    struct searchEorWstruct paramsXX;
    
    paramsXX.diffValueNextfoc=max(0.000000000001, (betapar*((1-mupar)*(tempWW1VE-tempUL1VE) + mupar*(tempWW0 - tempUL0))));
    paramsXX.ygridparXX = ygridfoc;
    
    if((betapar*((1-mupar)*(tempWW1VE-tempUL1VE) + mupar*(tempWW0 - tempUL0))) <= -0.0000001) {
    printf("There is a bug in VE"); getchar();
    }
    
    // Bracket the minimum (not sure to keep it) ---->
    testval1 = searchWfun(0.0, &paramsXX);
    testval2 = searchWfun(0.0001, &paramsXX);
    testval3 = searchWfun(EffortMaxW, &paramsXX);
    testval4 = searchWfun(EffortMaxW-0.0001, &paramsXX);
    
    iverif = 0;

    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffWW=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsXX);
        }
    
    focparamsVE->searchoutST=searcheffWW;
    
    
    valfoc=utilc(cons)+disutilityW(searcheffWW)+betapar*(piW(ygridfoc,searcheffWW)*(1.0 - mupar)*tempWW1VE + piW(ygridfoc, searcheffWW)*mupar*tempWW0 + mupar*(1 - piW(ygridfoc, searcheffWW))*tempUL0 + (1 - mupar)*(1-piW(ygridfoc, searcheffWW))*tempUL1VE);
    
    valfoc=-valfoc;
    
    if (ibug==1) {
        printf("focVE %20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","igridfoc","xval","ixgrid","dxgrid","cons","utilc(cons)","endoKfoc","vfRHS","valfoc");
        printf("focVE %20d\t%20.15f\t%20d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igridfoc,xval,ixgrid,dxgrid,cons,utilc(cons),endoKfoc,betapar*(piW(ygridfoc,searcheffWW)*(1.0 - mupar)*tempWW1VE + piW(ygridfoc, searcheffWW)*mupar*tempWW0 + mupar*(1 - piW(ygridfoc, searcheffWW))*tempUL0 + (1 - mupar)*(1-piW(ygridfoc, searcheffWW))*tempUL1VE),valfoc);
        getchar();
    }
    
    
    return valfoc;
}





// value function of an entrepreneur (with level of investment $endoKfoc).
double focVEwUI(const double xval, void * params)
{
    struct paramsgolden *focparamsVE= (struct paramsgolden *) params;
    double *valueVEwUInextfocVE=focparamsVE->valueVEwUInextST; // no cost when already employed.
    double *valueVEnextfocVE=focparamsVE->valueVEnextST; // no cost when already employed.
    double *valueWW0nextfocVE=focparamsVE->valueWW0nextST;
    double *valueUS0nextfocVE=focparamsVE->valueUS0nextST;
    double *valueUL0nextfocVE=focparamsVE->valueUL0nextST;
    double *valueWW1nextfocVE=focparamsVE->valueWW1nextST;
    double *valueUL1nextfocVE=focparamsVE->valueUL1nextST;
    double endoKfoc=focparamsVE->endoKST;
    int ygridfoc=focparamsVE->ygridST;
    int igridfoc=focparamsVE->igridST;
    int egridfoc=focparamsVE->egridST;
    
    double cons,xgrid,dxgrid,valfoc,tempWW0, tempUS0, tempUL0, tempWW1VEwUI, tempUL1VEwUI, tempWW1VE, tempUL1VE, testval1, testval2, testval3, testval4, searcheffWW, tempVE;
    int iverif, ixgrid, ibug, index1, index2, index3, index4;
    
    ibug=0;
    
    // here you have an additional insurance if the level of entrepreneurial acitivty is too low.
    // The additional effect comes from the insurance if bankruptcy.
    cons = findconsVEwUI(xval, igridfoc, egridfoc, ygridfoc, endoKfoc);

    //printf("%f %f %d %d %f", cons, xval, igridfoc, egridfoc, endoKfoc);
    xgrid=phiinv((xval));
    ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
    ixgrid=max(0,ixgrid);
    
    if (ixgrid>=(maxgrid-1))
	{
        xgrid=(maxgrid-1)-0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if (ixgrid<=0)
	{
        xgrid=0.000000001;
        ixgrid=min((maxgrid-1),(int)(floor(xgrid)));
        ixgrid=max(0,ixgrid);
	}
    
    if ((xval>(Gridmax))){printf("focVE: (xval>(Gridmax)) %20.15f\n",xval);ibug=1;}
    if ((ixgrid>=(maxgrid-1))){printf("focVE: ixgrid>=(maxgrid-1) %d\n",ixgrid);ibug=1;}
    if ((cons<0.0)){printf("focVE: consoFOC<0.0 %20.15f\t%20.15f\n",xval,cons);ibug=1;}
    if ((xval<(0.0))){printf("focVE: (xval<(0.0)) %20.15f\n",xval);ibug=1;}
    
    
    dxgrid=(xval-grid[ixgrid])/(grid[(ixgrid+1)]-grid[ixgrid]);

    tempWW1VE = 0.0;
    tempWW0 = 0.0;
    tempUL1VE = 0.0;
    tempUS0 = 0.0;
    tempUL0 = 0.0;
    tempWW1VEwUI = 0.0;
    tempUL1VEwUI = 0.0;
    
    for(int y=0;y<maxprod;y++) {
    	index1 = inx(ixgrid, y);
    	index2 = inx((ixgrid+1), y);
        
        tempUS0+=mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUS0nextfocVE[index1]),(valueUS0nextfocVE[index2])));
        tempUL0+=mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL0nextfocVE[index1]),(valueUL0nextfocVE[index2])));
        tempWW0+=mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW0nextfocVE[index1]),(valueWW0nextfocVE[index2])));

        // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
        if(egridfoc == 1) {
            tempVE = 0.0;
            for(int e=0; e<maxfirmtype; e++) {
                index3 = inxE(ixgrid, y, e);
                index4 = inxE((ixgrid+1), y, e);
                
                tempWW1VEwUI+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEwUInextfocVE[index3]),(valueVEwUInextfocVE[index4])));
                tempUL1VEwUI+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEwUInextfocVE[index3]),(valueVEwUInextfocVE[index4])),mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
                tempUL1VE+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])),mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
                tempWW1VE+=max(mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*mggtrans[ygridfoc][e]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])));
                
            }
        }
        
        if(egridfoc == 0) {
            index3 = inxE(ixgrid, y, 0);
            index4 = inxE((ixgrid+1), y, 0);
            
            tempWW1VEwUI+=max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEwUInextfocVE[index3]),(valueVEwUInextfocVE[index4])));
            tempUL1VEwUI+=max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEwUInextfocVE[index3]),(valueVEwUInextfocVE[index4])), mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
            tempUL1VE+=max(mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])), mprod[ygridfoc][y]*(inter1d(dxgrid,(valueUL1nextfocVE[index1]),(valueUL1nextfocVE[index2]))));
            tempWW1VE+=max(mprod[ygridfoc][y]*(inter1d(dxgrid,(valueWW1nextfocVE[index1]),(valueWW1nextfocVE[index2]))),mprod[ygridfoc][y]*inter1d(dxgrid,(valueVEnextfocVE[index3]),(valueVEnextfocVE[index4])));
        }
        
    }
    
    struct searchEorWstruct paramsXX;
    
    paramsXX.diffValueNextfoc=max(0.000000000001, (betapar*((1-mupar)*(PBnoUIVE*(tempWW1VE) + (1 - PBnoUIVE)*(tempWW1VEwUI) - PBnoUIVE*(tempUL1VE) - (1 - PBnoUIVE)*tempUL1VEwUI) + mupar*(tempWW0 - tempUS0)))); // if voluntary choose to shut down --> become unemployed without insurance.
    paramsXX.ygridparXX = ygridfoc;
    
    if((betapar*((1-mupar)*(PBnoUIVE*(tempWW1VE) + (1 - PBnoUIVE)*(tempWW1VEwUI) - PBnoUIVE*(tempUL1VE) - (1 - PBnoUIVE)*tempUL1VEwUI) + mupar*(tempWW0 - tempUS0))) <= -0.0000001) {
        printf("There is a bug in VEwUI"); getchar();
    }
    
    // Bracket the minimum (not sure to keep it) ---->
    testval1 = searchWfun(0.0, &paramsXX);
    testval2 = searchWfun(0.0001, &paramsXX);
    testval3 = searchWfun(EffortMaxW, &paramsXX);
    testval4 = searchWfun(EffortMaxW-0.0001, &paramsXX);
    
    iverif = 0;

    if(testval1 > testval2 & testval1 < 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    if(testval1 < testval2 & testval1 > 0.0) {
        searcheffWW = 0.0;
        iverif = 100;
        }
    
    if(testval3 > testval4 & testval3 < 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    if(testval3 < testval4 & testval3 > 0.0) {
        searcheffWW = EffortMaxW;
        printf("There is a bug:: EffortmaxW reached \n"); getchar();
        iverif = 100;
        }
    
    if(iverif == 0) {
        searcheffWW=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsXX);
        }
    
    focparamsVE->searchoutST=searcheffWW;
    
    
    valfoc=utilc(cons)+disutilityW(searcheffWW)+betapar*(piW(ygridfoc,searcheffWW)*(1.0 - mupar)*(PBnoUIVE*tempWW1VE + (1 - PBnoUIVE)*tempWW1VEwUI) + piW(ygridfoc, searcheffWW)*mupar*tempWW0 + mupar*(1 - piW(ygridfoc, searcheffWW))*(tempUS0) + (1 - mupar)*(1-piW(ygridfoc, searcheffWW))*((1-PBnoUIVE)*tempUL1VEwUI + PBnoUIVE*tempUL1VE));
    
    valfoc=-valfoc;
    
    if (ibug==1) {
        printf("focVE %20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","igridfoc","xval","ixgrid","dxgrid","cons","utilc(cons)","endoKfoc","vfRHS","valfoc");
        printf("focVE %20d\t%20.15f\t%20d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igridfoc,xval,ixgrid,dxgrid,cons,utilc(cons),endoKfoc,betapar*(piW(ygridfoc,searcheffWW)*(1.0 - mupar)*(PBnoUIVE*tempWW1VE + (1 - PBnoUIVE)*tempWW1VEwUI) + piW(ygridfoc, searcheffWW)*mupar*tempWW0 + mupar*(1 - piW(ygridfoc, searcheffWW))*(tempUS0) + (1 - mupar)*(1-piW(ygridfoc, searcheffWW))*((1-PBnoUIVE)*tempUL1VEwUI + PBnoUIVE*tempUL1VE)),valfoc);
        getchar();
    }
    
    
    return valfoc;
}


/////////////////////////
// END VALUE FUNCTIONS //
/////////////////////////





//////////////////////////////
// VALUE FUNCTION ITERATION //
//////////////////////////////

// find value functions (Value function Iteration).
void VFI(double *valueWW0, double *valueWW1, double *valueUS0, double *valueUS1, double *valueUL0, double *valueUL1, double *valueVE, double *valueVEc, double *valueVEwUI, double *valueVEwUIc, double *endoK,double *saveWW0, double *saveWW1, double *saveUS0, double *saveUS1, double *saveUL0, double *saveUL1, double *saveVE, double *saveVEwUI, double *searchEffortWWUS0, double *searchEffortWWUL0, double *searchEffortWWUS1, double *searchEffortVEUS1, double *searchEffortWWUL1, double *searchEffortVEUL1, double *searchEffortWWVE,double *searchEffortVEWW, double *searchEffortWWVEwUI)
{
    int ygrid, igrid, egrid, e, y, iter, iterK, jcount, ivaltype;
    double critereVF,critereendoK, valmax, valmaxold;
    
    // For Worker
    double *valueWW0next, valfnmaxWW0, testvalueWW0, wealthWW0, savemaxWW0, amaxsaveWW0;
    int isavemaxWW0, icaseWW0;
    double *valueWW1next, valfnmaxWW1, testvalueWW1, wealthWW1, savemaxWW1, amaxsaveWW1, searchvalVEWW;
    int isavemaxWW1, icaseWW1;
    
    // For entrepreneur
    double *valueVEnext, *endoKnext, valfnmaxVE, savemaxVE, amaxsaveVE,testvalueVE, wealthVE,searchvalWWVE;
    double opticap, constcap, investmax,investmin,investgrid,dinvestgrid,BRCtempUmin,BRCtempUmax; // BRC (for CD type)
    double xonentry, xentry, dxentry, *valueVEcnext;
    int ixentry, ixinvgrid, isavemaxVE, icaseVE;
    double *valueVEwUInext, valfnmaxVEwUI, savemaxVEwUI, amaxsaveVEwUI,testvalueVEwUI, wealthVEwUI,searchvalWWVEwUI;
    double investmaxVEwUI,investminVEwUI,investgridVEwUI,dinvestgridVEwUI,BRCtempUminVEwUI,BRCtempUmaxVEwUI; // BRC (for CD type)
    double xonentryVEwUI, xentryVEwUI, dxentryVEwUI, *valueVEwUIcnext;
    int ixentryVEwUI, ixinvgridVEwUI, isavemaxVEwUI, icaseVEwUI;
    
    // For unemployed
    double *valueUS0next, *valueUL0next, valfnmaxUS0,savemaxUS0,amaxsaveUS0,searchvalWWUS0, testvalueUL0, wealthUL0,testvalueUS0, wealthUS0,valfnmaxUL0,savemaxUL0,amaxsaveUL0,searchvalWWUL0;
    int isavemaxUS0,isavemaxUL0, icaseUS0, icaseUL0;
    double *valueUS1next, *valueUL1next, valfnmaxUS1,savemaxUS1,amaxsaveUS1,searchvalWWUS1,searchvalVEUS1, testvalueUL1, wealthUL1,testvalueUS1, wealthUS1,valfnmaxUL1,savemaxUL1,amaxsaveUL1,searchvalWWUL1,searchvalVEUL1;
    int isavemaxUS1,isavemaxUL1, icaseUS1, icaseUL1;
    
    // temp val
    double tempU,tempB, tempUS0, tempVE, tempVE2, tempUL0, tempWW0, tempUL1VE, tempWW1, tempUL1, tempUS1, tempWW1VE, tempUS1VE, tempWW1VE2, tempUL1VE2, tempVEwUI, tempUS1VEwUI, tempWW1VEwUI, tempWW1VEwUI2, tempUL1VEwUI2;

    double testval1, testval2, testval3, testval4;
    int iverif;

    // initialize files
    FILE *valuefnoutfile,*saveoutfile,*searchoutfile,*endokoutfile;
    
    // initialize pointers
    valueVEnext = (double *) calloc((ifulldimE), sizeof(double));
    valueVEwUInext = (double *) calloc((ifulldimE), sizeof(double));
    valueWW0next = (double *) calloc((ifulldim), sizeof(double));
    valueUS0next = (double *) calloc((ifulldim), sizeof(double));
    valueUL0next = (double *) calloc((ifulldim), sizeof(double));
    valueWW1next = (double *) calloc((ifulldim), sizeof(double));
    valueUS1next = (double *) calloc((ifulldim), sizeof(double));
    valueUL1next = (double *) calloc((ifulldim), sizeof(double));
    valueVEcnext = (double *) calloc((ifulldimE), sizeof(double));
    valueVEwUIcnext = (double *) calloc((ifulldimE), sizeof(double));
    endoKnext = (double *) calloc((ifulldimKK), sizeof(double));

    
    struct paramsgolden params;
    struct searchEorWstruct paramsYY;
    struct searchfocWwhenE paramsUU;


#if BRC == 2
   struct searchfendoKYY paramsKK;
    // initialize max capital (CD TYPE)
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
       for(int egrid=0;egrid<maxfirmtype;egrid++)
        {
            for(igrid=0;igrid<maxgrid;igrid++)
            {
                opticap = pow(((1-corptax)*TFPpar*gtrans[egrid]*nupar)/(rstar + deltapar), 1/(1-nupar));
                constcap = amaxendoK;
                endoK[inxKK(igrid, ygrid, egrid)] = min(min(constcap,opticap),amaxendoK);
            }
        }
    }
    
    iterK=0;
    critereendoK=10.0;
    

    //boucle endo investment.
    while ((critereendoK > epsilonendoK) && iterK < maxiterendoK)
    {
#endif
        bascule(endoK,endoKnext,ifulldimKK);

        critereVF=1.0;
        
        iter=0;

#if indexPM == 0
        // boucle VFI
        while ((critereVF > epsilonValue)  && (iter < maxiterVF))
        {
#endif
            
            bascule(valueVE,valueVEnext,ifulldimE);
            bascule(valueVEc,valueVEcnext,ifulldimE); // with entry cost
            bascule(valueVEwUI,valueVEwUInext,ifulldimE);
            bascule(valueVEwUIc,valueVEwUIcnext,ifulldimE); // with entry cost
            bascule(valueWW0,valueWW0next,ifulldim);
            bascule(valueUS0,valueUS0next,ifulldim);
            bascule(valueUL0,valueUL0next,ifulldim);
            bascule(valueWW1,valueWW1next,ifulldim);
            bascule(valueUS1,valueUS1next,ifulldim);
            bascule(valueUL1,valueUL1next,ifulldim);
            
            #if OMP == 1
            omp_set_dynamic(0);     // Explicitly disable dynamic teams
            omp_set_num_threads(nbthread); // Use 4 threads for all consecutive parallel regions
           
            #pragma omp parallel
            {
            #pragma omp for private(params, paramsYY, paramsUU, ygrid, igrid, egrid, e, y, valfnmaxWW0, testvalueWW0, wealthWW0, savemaxWW0, amaxsaveWW0, isavemaxWW0, icaseWW0, valfnmaxWW1, testvalueWW1, wealthWW1, savemaxWW1, amaxsaveWW1, searchvalVEWW, isavemaxWW1, icaseWW1, valfnmaxVE, savemaxVE, amaxsaveVE,testvalueVE, wealthVE,searchvalWWVE, opticap, constcap, investmax,investmin,investgrid,dinvestgrid,BRCtempUmin,BRCtempUmax, xonentry, xentry, dxentry, ixentry, ixinvgrid, isavemaxVE, icaseVE, valfnmaxUS0,savemaxUS0,amaxsaveUS0,searchvalWWUS0, testvalueUL0, wealthUL0,testvalueUS0, wealthUS0,valfnmaxUL0,savemaxUL0,amaxsaveUL0,searchvalWWUL0, isavemaxUS0,isavemaxUL0, icaseUS0, icaseUL0, valfnmaxUS1,savemaxUS1,amaxsaveUS1,searchvalWWUS1,searchvalVEUS1, testvalueUL1, wealthUL1,testvalueUS1, wealthUS1,valfnmaxUL1,savemaxUL1,amaxsaveUL1,searchvalWWUL1,searchvalVEUL1, isavemaxUS1,isavemaxUL1, icaseUS1, icaseUL1, tempU,tempB, tempUS0, tempVE, tempVE2, tempUL0, tempWW0, tempUL1VE, tempWW1, tempUL1, tempUS1, tempWW1VE, tempUS1VE, tempWW1VE2, tempUL1VE2, tempVEwUI, tempUS1VEwUI, tempWW1VEwUI, tempWW1VEwUI2, tempUL1VEwUI2, testval1, testval2, testval3, testval4, iverif, valfnmaxVEwUI, savemaxVEwUI, amaxsaveVEwUI,testvalueVEwUI, wealthVEwUI,searchvalWWVEwUI, investmaxVEwUI,investminVEwUI,investgridVEwUI,dinvestgridVEwUI,BRCtempUminVEwUI,BRCtempUmaxVEwUI, xonentryVEwUI, xentryVEwUI, dxentryVEwUI, ixentryVEwUI, ixinvgridVEwUI, isavemaxVEwUI, icaseVEwUI)
            #endif
            
        
            
           for(ygrid=0;ygrid<maxprod;ygrid++)
            {
            
                    for(igrid=0;igrid<maxgrid;igrid++)
                    {
    
                        
                        tempUS0 = 0.0;
                        tempUL0 = 0.0;
                        tempWW0 = 0.0;
                        tempUS1 = 0.0;
                        tempUL1 = 0.0;
                        tempWW1 = 0.0;
                        tempWW1VE = 0.0;
                        tempUL1VE = 0.0;
                        tempUS1VE = 0.0;
                        tempUS1VEwUI = 0.0;
                        tempWW1VEwUI = 0.0;
                        
                        for(y=0; y<maxprod; y++) {
                        
                            // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
                            tempVE = 0.0;
                            for(e=0; e<maxfirmtype; e++) {
                                tempVE += mgtrans[y][e]*valueVEcnext[inxE(0, y, e)];
                            }
                            
                            tempUL1VE += max(mprod[ygrid][y]*valueUL1next[inx(0, y)], mprod[ygrid][y]*tempVE);
                            tempWW1VE += max(mprod[ygrid][y]*valueWW1next[inx(0, y)], mprod[ygrid][y]*tempVE);
                            
                            #if SEA == 1
                            tempVEwUI = 0.0;
                            for(e=0; e<maxfirmtype; e++) {
                                tempVEwUI += mgtrans[y][e]*valueVEwUIcnext[inxE(0, y, e)];
                            }
                            tempUS1VEwUI += max(mprod[ygrid][y]*valueUS1next[inx(0, y)], mprod[ygrid][y]*tempVEwUI);
                            tempWW1VEwUI += max(mprod[ygrid][y]*valueWW1next[inx(0, y)], mprod[ygrid][y]*tempVEwUI);
                            tempUS1VE += max(mprod[ygrid][y]*valueUS1next[inx(0, y)], mprod[ygrid][y]*tempVEwUI);
                            #endif
                            
                            #if SEA == 0
                            tempUS1VE += max(mprod[ygrid][y]*valueUS1next[inx(0, y)], mprod[ygrid][y]*tempVE);
                            #endif
                            
                            
                            tempWW1 += mprod[ygrid][y]*valueWW1next[inx(0, y)];
                            tempWW0 += mprod[ygrid][y]*valueWW0next[inx(0, y)];
                            tempUS1 += mprod[ygrid][y]*valueUS1next[inx(0, y)];
                            tempUS0 += mprod[ygrid][y]*valueUS0next[inx(0, y)];
                            tempUL0 += mprod[ygrid][y]*valueUL0next[inx(0, y)];
                            tempUL1 += mprod[ygrid][y]*valueUL1next[inx(0, y)];

                        }
                
                
                        paramsUU.tempWW1val = tempWW1;
                        paramsUU.tempWW0val = tempWW0;
                        paramsUU.tempUS1val = tempUS1;
                        paramsUU.tempUS0val = tempUS0;
                        paramsUU.tempUL1val = tempUL1;
                        paramsUU.tempUL0val = tempUL0;
                        paramsUU.tempUL1VEval = tempUL1VE;
                        
                        #if SEA == 0
                        paramsUU.tempWW1VEval = tempWW1VE;
                        paramsUU.tempUS1VEval = tempUS1VE;
                        #endif
                        
                        #if SEA == 1
                        paramsUU.tempWW1VEwUIval = tempWW1VEwUI;
                        paramsUU.tempWW1VEval = tempWW1VE;
                        paramsUU.tempUS1VEval = tempUS1VEwUI;
                        #endif
                        
                        params.igridST=igrid;
                        params.ygridST=ygrid;
                        params.valueVEnextST=valueVEnext;
                        params.valueVEcnextST=valueVEcnext;
                        params.valueVEwUInextST=valueVEwUInext;
                        params.valueVEwUIcnextST=valueVEwUIcnext;
                        params.valueWW0nextST=valueWW0next;
                        params.valueUS0nextST=valueUS0next;
                        params.valueUL0nextST=valueUL0next;
                        params.valueWW1nextST=valueWW1next;
                        params.valueUS1nextST=valueUS1next;
                        params.valueUL1nextST=valueUL1next;
                        
                        
                        // printf("I AM HERE BEFORE WW \n"); getchar();
        
                        ///////////////////////
                        //Solving the WW case//
                        ///////////////////////
                        
                        // WITHOUT AN IDEA //
              
                        icaseWW0 = 100;
                        
                        testvalueWW0=phi(0)+0.000000001;
                        wealthWW0=findconsWW(0.0,igrid,ygrid);

                        if (wealthWW0<=testvalueWW0 || wealthWW0 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                
                                testvalueWW0=wealthWW0-0.000000001;
                                icaseWW0=101;
                                
                                if (testvalueWW0<=0.0 || wealthWW0 <= 0.0)
                                {//wealth is so small that consumption is almost zero

                                    valfnmaxWW0=utilc(0.0000000001)+betapar*(xipar*((1.0-etapar[ygrid])*(tempWW1)+etapar[ygrid]*tempUS1)+(1-xipar)*(etapar[ygrid]*tempUS0 + (1.0 - etapar[ygrid])*tempWW0));
                                    
                                    valfnmaxWW0=-valfnmaxWW0;
                                    amaxsaveWW0=Gridmin;
                                    icaseWW0=2;
                                    
                                }
                            }
                        
                        
                        if (icaseWW0>=100)
                        {
                            if (focWWnoID(Gridmin,&params)<focWWnoID(testvalueWW0,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxWW0=focWWnoID(Gridmin,&params);
                                amaxsaveWW0=Gridmin;
                            }
                            else
                            {//cas ou le max est a droite de la borne inf
                                
                                //printf("here\n");
                                
                                // define the highest level of xval, which means that consume 0.000000001 and kept the rest for saving.
                                savemaxWW0=min(Gridmax,findconsWW((0.0),igrid, ygrid)-0.000000001);
                                
                                if (savemaxWW0>=(Gridmax))
                                {
                                    savemaxWW0=Gridmax-0.000000001;
                                }
                                
                                isavemaxWW0=(int)floor(phiinv(savemaxWW0));
                                
                                if ((isavemaxWW0<=0)){printf("isavemaxWW0<=0 %d\n",isavemaxWW0);getchar();}
                                
                                
                                // Si utility(today avec un peu moins) < utility(today avec un peu plus) alors intérêt encore à épargner (monotonicity).
                                // Recall that focWW give the opposite of the valuefunction (find min)
                                if (focWWnoID(phi((isavemaxWW0)),&params)>focWWnoID(savemaxWW0,&params))
                                { // corner solution at the top.
                                    //printf("savemax\n");
                                    valfnmaxWW0=focWWnoID(savemaxWW0, &params);
                                    amaxsaveWW0=savemaxWW0;
                                    
                                }
                                else
                                {// interior solution
                                    valfnmaxWW0=mygolden(Gridmin,Gridmin+0.000001,savemaxWW0,TOL,amaxsaveWW0,&params,focWWnoID);
                                }
                                
                                
                            }
                        }
                        
                        valueWW0[inx(igrid, ygrid)]  =   -valfnmaxWW0;
                        saveWW0[inx(igrid, ygrid)]   =   amaxsaveWW0;



                        // WITH AN IDEA
                        
                        icaseWW1 = 100;
                        
                        testvalueWW1=phi(0)+0.000000001;
                        wealthWW1=findconsWW(0.0,igrid,ygrid);

                        if (wealthWW1<=testvalueWW1 || wealthWW1 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                testvalueWW1=wealthWW1-0.000000001;
                                icaseWW1=101;
                                
                                if (testvalueWW1<=0.0 || wealthWW1 <= 0.0)
                                {//wealth is so small that consumption is almost zero
                                    
                                    paramsYY.diffValueNextfoc = max(0.000000000001,(betapar*((1-zetapar)*((1-etapar[ygrid])*(tempWW1VE - tempWW1) + etapar[ygrid]*(tempUS1VE - tempUS1)))));
                                    paramsYY.ygridparXX = ygrid;
                                    
                                    // Bracket the minimum (not sure to keep it) ---->
                                    testval1 = searchEfun(0.0, &paramsYY);
                                    testval2 = searchEfun(0.0001, &paramsYY);
                                    testval3 = searchEfun(EffortMaxE, &paramsYY);
                                    testval4 = searchEfun(EffortMaxE-0.0001, &paramsYY);
                                    
                                    iverif = 0;

                                    if(testval1 > testval2 & testval1 < 0.0) {
                                        searchvalVEWW = 0.0;
                                        iverif = 100;
                                        }
                                    if(testval1 < testval2 & testval1 > 0.0) {
                                        searchvalVEWW = 0.0;
                                        iverif = 100;
                                        }
                                    
                                    if(testval3 > testval4 & testval3 < 0.0) {
                                        searchvalVEWW = EffortMaxE;
                                        printf("There is a bug:: EffortmaxE reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    if(testval3 < testval4 & testval3 > 0.0) {
                                        searchvalVEWW = EffortMaxE;
                                        printf("There is a bug:: EffortmaxE reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    
                                    if(iverif == 0) {
                                        searchvalVEWW=zbrentNEW(searchEfun,0.0,EffortMaxE,(1.0e-10),&paramsYY);
                                        }
                                    
                                    valfnmaxWW1=utilc(0.0000000001)+disutilityE(searchvalVEWW)+betapar*((1-zetapar)*(piE(searchvalVEWW)*(1 - etapar[ygrid])*tempWW1VE + etapar[ygrid]*(1 - piE(searchvalVEWW))*tempUS1 + etapar[ygrid]*piE(searchvalVEWW)*tempUS1VE + (1 - etapar[ygrid])*(1 - piE(searchvalVEWW))*tempWW1) + zetapar*(etapar[ygrid]*tempUS0 + (1.0 - etapar[ygrid])*tempWW0));
                                    
                                    valfnmaxWW1=-valfnmaxWW1;
                                    amaxsaveWW1=Gridmin;
                                    icaseWW1=2;
                                }
                            }
                    
                        
                        if (icaseWW1>=100)
                        {
                            if (focWWwID(Gridmin,&params)<focWWwID(testvalueWW1,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxWW1=focWWwID(Gridmin,&params);
                                searchvalVEWW=params.searchoutST;
                                amaxsaveWW1=Gridmin;
                            }
                            else
                            {
                                //printf("here\n");
                                
                                savemaxWW1=min(Gridmax,findconsWW((0.0),igrid, ygrid)-0.000000001);
                                
                                if (savemaxWW1>=(Gridmax))
                                {
                                    savemaxWW1=Gridmax-0.000000001;
                                }
                                
                                
                                isavemaxWW1=(int)floor(phiinv(savemaxWW1));
                                
                                if ((isavemaxWW1<=0)){printf("isavemaxWW1<=0 %d\n",isavemaxWW1);getchar();}
                                
                                
                                if (focWWwID(phi((isavemaxWW1)),&params)>focWWwID(savemaxWW1,&params))
                                {//corner solution at the top
                                    //printf("savemax\n");
                                    valfnmaxWW1=focWWwID(savemaxWW1, &params);
                                    searchvalVEWW=params.searchoutST;
                                    amaxsaveWW1=savemaxWW1;
                                    
                                }
                                else
                                {//cas standard
                                    
                                    //printf("stdcase\n");

                                    valfnmaxWW1=mygolden1search(Gridmin,Gridmin+0.000001,savemaxWW1,TOL,amaxsaveWW1,searchvalVEWW,&params,focWWwID);
                                    
                                    //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,amaxsaveUS,searchvalUS);
                                }
                                
                                
                            }
                        }
                        
                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,0.0,searchvalUS);
                        valueWW1[inx(igrid, ygrid)]=-valfnmaxWW1;
                        saveWW1[inx(igrid, ygrid)]=amaxsaveWW1;
                        searchEffortVEWW[inx(igrid, ygrid)]=searchvalVEWW;
                        
                        
                        //printf("I AM HERE AFTER WW %f \n", searchvalVEWW); getchar();
                        



                        //////////////////////////////
                        //Solving the ST unemployed //
                        //////////////////////////////
                        icaseUS0 = 100;
                        
                        params.typeU = 0;
                        
                        // WITHOUT AN IDEA
                        
                        testvalueUS0=phi(0)+0.000000001;
                        wealthUS0=findconsUS(0.0,igrid,ygrid);

                        if (wealthUS0<=testvalueUS0 || wealthUS0 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                testvalueUS0=wealthUS0-0.000000001;
                                icaseUS0=101;
                                
                                if (testvalueUS0<=0.0 || wealthUS0 <= 0.0)
                                {//wealth is so small that consumption is almost zero
   
                                    // FIND THE OPTIMAL EFFORT LEVEL
                                    paramsYY.diffValueNextfoc = max(0.000000000001,(betapar*(xipar*(tempWW1 - ((1-pLTpar)*tempUS1 + pLTpar*tempUL1)) + (1 - xipar)*(tempWW0 - ((1-pLTpar)*tempUS0 + pLTpar*tempUL0)))));
                                    paramsYY.ygridparXX = ygrid;
                                    
                                    // Bracket the minimum (not sure to keep it) ---->
                                    testval1 = searchWfunU(0.0, &paramsYY);
                                    testval2 = searchWfunU(0.0001, &paramsYY);
                                    testval3 = searchWfunU(EffortMaxW, &paramsYY);
                                    testval4 = searchWfunU(EffortMaxW-0.0001, &paramsYY);
                                    
                                    iverif = 0;

                                    if(testval1 > testval2 & testval1 < 0.0) {
                                        searchvalWWUS0 = 0.0;
                                        iverif = 100;
                                        }
                                    if(testval1 < testval2 & testval1 > 0.0) {
                                        searchvalWWUS0 = 0.0;
                                        iverif = 100;
                                        }
                                    
                                    if(testval3 > testval4 & testval3 < 0.0) {
                                        searchvalWWUS0 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    if(testval3 < testval4 & testval3 > 0.0) {
                                        searchvalWWUS0 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    
                                    if(iverif == 0) {
                                        searchvalWWUS0=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsYY);
                                        }
                                    
                                
                                    valfnmaxUS0=utilc(0.0000000001)+disutilityW(searchvalWWUS0)+betapar*(xipar*((1-piWU(ygrid, searchvalWWUS0))*((1-pLTpar)*tempUS1 + pLTpar*tempUL1) + piWU(ygrid, searchvalWWUS0)*tempWW1) + (1 - xipar)*(1 - piWU(ygrid, searchvalWWUS0))*((1-pLTpar)*tempUS0 + pLTpar*tempUL0)+piWU(ygrid, searchvalWWUS0)*tempWW0);
                                    
                                    valfnmaxUS0=-valfnmaxUS0;
                                    amaxsaveUS0=Gridmin;
                                    icaseUS0=2;
                                }
                            }
                        
     
                        
                        if (icaseUS0>=100)
                        {
                            if (focUnoID(Gridmin,&params)<focUnoID(testvalueUS0,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxUS0=focUnoID(Gridmin,&params);
                                searchvalWWUS0=params.searchoutST;
                                amaxsaveUS0=Gridmin;
                            }
                            else
                            {
                                savemaxUS0=min(Gridmax,findconsUS((0.0),igrid, ygrid)-0.000000001);
                                
                                if (savemaxUS0>=(Gridmax))
                                {
                                    savemaxUS0=Gridmax-0.000000001;
                                }
                                
                                
                                isavemaxUS0=(int)floor(phiinv(savemaxUS0));
                                
                                if ((isavemaxUS0<0)){printf("isavemaxUS<=0 %d\n",isavemaxUS0);getchar();}
                                
                                
                                if (focUnoID(phi((isavemaxUS0)),&params)>focUnoID(savemaxUS0,&params))
                                {//corner solution at the top
                                    //printf("savemax\n");
                                    valfnmaxUS0=focUnoID(savemaxUS0, &params);
                                    searchvalWWUS0=params.searchoutST;
                                    amaxsaveUS0=savemaxUS0;
                                    
                                }
                                else
                                {//cas standard
                                    
                                    //printf("stdcase\n");

                                    valfnmaxUS0=mygolden1search(Gridmin,Gridmin+0.000001,savemaxUS0,TOL,amaxsaveUS0,searchvalWWUS0,&params,focUnoID);
                                    
                                    //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,amaxsaveUS,searchvalUS);
                                }
                                
                                
                            }
                        }
                        
                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,0.0,searchvalUS);
                        valueUS0[inx(igrid, ygrid)]=-valfnmaxUS0;
                        saveUS0[inx(igrid, ygrid)]=amaxsaveUS0;
                        searchEffortWWUS0[inx(igrid, ygrid)]=searchvalWWUS0;
                        
                        
                        // WITH AN IDEA
                        
                        icaseUS1 = 100;
                        
                        params.typeU = 0;
                        
                        testvalueUS1=phi(0)+0.000000001;
                        wealthUS1=findconsUS(0.0,igrid,ygrid);

                        if (wealthUS1<=testvalueUS1 || wealthUS1 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                testvalueUS1=wealthUS1-0.000000001;
                                icaseUS1=101;
                                
                                if (testvalueUS1<=0.0 || wealthUS1 <= 0.0)
                                {//wealth is so small that consumption is almost zero
                
                                    // STARTING 2D searching effort
                                    paramsUU.ygridparXX = ygrid;
                                    paramsUU.typeUU = params.typeU;
                                    
                                    testval1 = searchWwhenEknownfun(0.0, &paramsUU);
                                    testval2 = searchWwhenEknownfun(0.00001, &paramsUU);
                                    testval3 = searchWwhenEknownfun(EffortMaxW, &paramsUU);
                                    testval4 = searchWwhenEknownfun(EffortMaxW-0.00001, &paramsUU);
                                    
                                    iverif = 0;

                                    if(testval1 > testval2 & testval1 < 0.0) {
                                        searchvalWWUS1 = 0.0;
                                        iverif = 100;
                                        }
                                    if(testval1 < testval2 & testval1 > 0.0) {
                                        searchvalWWUS1 = 0.0;
                                        iverif = 100;
                                        }
                                    
                                    if(testval3 > testval4 & testval3 < 0.0) {
                                        searchvalWWUS1 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    if(testval3 < testval4 & testval3 > 0.0) {
                                        searchvalWWUS1 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    
                                    if(iverif == 0) {
                                        searchvalWWUS1 = zbrentNEW(searchWwhenEknownfun,0.0,EffortMaxW,(1.0e-10),&paramsUU);
                                        searchvalVEUS1 = paramsUU.searchUsEoutST;
                                        }
                                    
                                 
                                    // COMPUTE THE CONTINUATION VALUE
                                    #if SEA == 0
                                    valfnmaxUS1=utilc(0.0000000001)+disutilityW(searchvalWWUS1)+disutilityE(searchvalVEUS1)+betapar*((1-zetapar)*(piWU(ygrid,searchvalWWUS1)*(1 - piEU(searchvalVEUS1))*tempWW1+piEU(searchvalVEUS1)*(1-piWU(ygrid, searchvalWWUS1))*((1 - pLTpar)*tempUS1VE + pLTpar*tempUL1VE) + piWU(ygrid, searchvalWWUS1)*piEU(searchvalVEUS1)*tempWW1VE + (1 - piEU(searchvalVEUS1))*(1 - piWU(ygrid, searchvalWWUS1))*((1-pLTpar)*tempUS1 + pLTpar*tempUL1)) + zetapar*((1 - piWU(ygrid, searchvalWWUS1))*((1-pLTpar)*tempUS0 + pLTpar*tempUL0) + piWU(ygrid, searchvalWWUS1)*tempWW0));
                                    #endif
                                    
                                    #if SEA == 1
                                    valfnmaxUS1=utilc(0.0000000001)+disutilityW(searchvalWWUS1)+disutilityE(searchvalVEUS1)+betapar*((1-zetapar)*(piWU(ygrid,searchvalWWUS1)*(1 - piEU(searchvalVEUS1))*tempWW1+piEU(searchvalVEUS1)*(1-piWU(ygrid, searchvalWWUS1))*((1 - pLTpar)*tempUS1VEwUI + pLTpar*tempUL1VE) + piWU(ygrid, searchvalWWUS1)*piEU(searchvalVEUS1)*((1 - PBnoUIVE)*tempWW1VEwUI + PBnoUIVE*tempWW1VE) + (1 - piEU(searchvalVEUS1))*(1 - piWU(ygrid, searchvalWWUS1))*((1-pLTpar)*tempUS1 + pLTpar*tempUL1)) + zetapar*((1 - piWU(ygrid, searchvalWWUS1))*((1-pLTpar)*tempUS0 + pLTpar*tempUL0) + piWU(ygrid, searchvalWWUS1)*tempWW0));
                                    #endif
                                
                                    
                                    valfnmaxUS1=-valfnmaxUS1;
                                    amaxsaveUS1=Gridmin;
                                    icaseUS1=2;
                                }
                            }
                        
     
                        
                        if (icaseUS1>=100)
                        {
                            if (focUwID(Gridmin,&params)<focUwID(testvalueUS1,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxUS1=focUwID(Gridmin,&params);
                                searchvalWWUS1=params.searchUWoutST;
                                searchvalVEUS1=params.searchUEoutST;
                                amaxsaveUS1=Gridmin;
                            }
                            else
                            {
                                savemaxUS1=min(Gridmax,findconsUS((0.0),igrid, ygrid)-0.000000001);
                                
                                if (savemaxUS1>=(Gridmax))
                                {
                                    savemaxUS1=Gridmax-0.000000001;
                                }
                                
                                
                                isavemaxUS1=(int)floor(phiinv(savemaxUS1));
                                
                                if ((isavemaxUS1<0)){printf("isavemaxUS<=0 %d\n",isavemaxUS1);getchar();}
                                
                                
                                if (focUwID(phi((isavemaxUS1)),&params)>focUwID(savemaxUS1,&params))
                                {//corner solution at the top
                                    //printf("savemax\n");
                                    valfnmaxUS1=focUwID(savemaxUS1, &params);
                                    searchvalWWUS1=params.searchUWoutST;
                                    searchvalVEUS1=params.searchUEoutST;
                                    amaxsaveUS1=savemaxUS1;
                                    
                                }
                                else
                                {//cas standard
                                    
                                    //printf("stdcase\n");

                                    valfnmaxUS1=mygolden2search(Gridmin,Gridmin+0.000001,savemaxUS1,TOL,amaxsaveUS1,searchvalVEUS1, searchvalWWUS1,&params,focUwID);
                                    
                                    //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,amaxsaveUS,searchvalUS);
                                }
                                
                                
                            }
                        }
                        
                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,0.0,searchvalUS);
                        valueUS1[inx(igrid, ygrid)]=-valfnmaxUS1;
                        saveUS1[inx(igrid, ygrid)]=amaxsaveUS1;
                        searchEffortWWUS1[inx(igrid, ygrid)]=searchvalWWUS1;
                        searchEffortVEUS1[inx(igrid, ygrid)]=searchvalVEUS1;
                        
                        //printf("I AM HERE BEFORE VE %f %f \n", searchvalWWUS1, searchvalVEUS1); getchar();
                        
                        
                        ///////////////////////////////
                        // Solving the LT unemployed //
                        ///////////////////////////////
                        icaseUL0 = 100;
                        
                        params.typeU = 1;
                        
                        // WITHOUT AN IDEA
                        
                        testvalueUL0=phi(0)+0.000000001;
                        wealthUL0=findconsUL(0.0,igrid);

                        if (wealthUL0<=testvalueUL0 || wealthUL0 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                testvalueUL0=wealthUL0-0.000000001;
                                icaseUL0=101;
                                
                                if (testvalueUL0<=0.0 || wealthUL0 <= 0.0)
                                {//wealth is so small that consumption is almost zero
   
                                    // FIND THE OPTIMAL EFFORT LEVEL
                                    paramsYY.diffValueNextfoc = max(0.000000000001,(betapar*(xipar*(tempWW1 - tempUL1) + (1 - xipar)*(tempWW0 - tempUL0))));
                                    paramsYY.ygridparXX = ygrid;
                                    
                                    // Bracket the minimum (not sure to keep it) ---->
                                    testval1 = searchWfunU(0.0, &paramsYY);
                                    testval2 = searchWfunU(0.0001, &paramsYY);
                                    testval3 = searchWfunU(EffortMaxW, &paramsYY);
                                    testval4 = searchWfunU(EffortMaxW-0.0001, &paramsYY);
                                    
                                    iverif = 0;

                                    if(testval1 > testval2 & testval1 < 0.0) {
                                        searchvalWWUL0 = 0.0;
                                        iverif = 100;
                                        }
                                    if(testval1 < testval2 & testval1 > 0.0) {
                                        searchvalWWUL0 = 0.0;
                                        iverif = 100;
                                        }
                                    
                                    if(testval3 > testval4 & testval3 < 0.0) {
                                        searchvalWWUL0 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    if(testval3 < testval4 & testval3 > 0.0) {
                                        searchvalWWUL0 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    
                                    if(iverif == 0) {
                                        searchvalWWUL0=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsYY);
                                        }
                                    
                                    valfnmaxUL0=utilc(0.0000000001)+disutilityW(searchvalWWUL0)+betapar*(xipar*((1-piWU(ygrid, searchvalWWUL0))*tempUL1 + piWU(ygrid, searchvalWWUL0)*tempWW1) + (1 - xipar)*(1 - piWU(ygrid, searchvalWWUL0))*tempUL0+piWU(ygrid, searchvalWWUL0)*tempWW0);
                                    
                                    valfnmaxUL0=-valfnmaxUL0;
                                    amaxsaveUL0=Gridmin;
                                    icaseUL0=2;
                                }
                            }
                        
                        if (icaseUL0>=100)
                        {
                            if (focUnoID(Gridmin,&params)<focUnoID(testvalueUL0,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxUL0=focUnoID(Gridmin,&params);
                                searchvalWWUL0=params.searchoutST;
                                amaxsaveUL0=Gridmin;
                            }
                            else
                            {
                                savemaxUL0=min(Gridmax,findconsUL((0.0),igrid)-0.000000001);
                                
                                if (savemaxUL0>=(Gridmax))
                                {
                                    savemaxUL0=Gridmax-0.000000001;
                                }
                                
                                
                                isavemaxUL0=(int)floor(phiinv(savemaxUL0));
                                
                                if ((isavemaxUL0<0)){printf("isavemaxUL<=0 %d\n",isavemaxUL0);getchar();}
                                
                                
                                if (focUnoID(phi((isavemaxUL0)),&params)>focUnoID(savemaxUL0,&params))
                                {//corner solution at the top
                                    //printf("savemax\n");
                                    valfnmaxUL0=focUnoID(savemaxUL0, &params);
                                    searchvalWWUL0=params.searchoutST;
                                    amaxsaveUL0=savemaxUL0;
                                    
                                }
                                else
                                {//cas standard
                                    
                                    //printf("stdcase\n");

                                    valfnmaxUL0=mygolden1search(Gridmin,Gridmin+0.000001,savemaxUL0,TOL,amaxsaveUL0,searchvalWWUL0,&params,focUnoID);
                                    
                                    //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,amaxsaveUS,searchvalUS);
                                }
                                
                                
                            }
                        }
                        
                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,0.0,searchvalUS);
                        valueUL0[inx(igrid, ygrid)]=-valfnmaxUL0;
                        saveUL0[inx(igrid, ygrid)]=amaxsaveUL0;
                        searchEffortWWUL0[inx(igrid, ygrid)]=searchvalWWUL0;
                        
                        // WITH AN IDEA
                        
                        icaseUL1 = 100;
                        
                        params.typeU = 1;
                        
                        testvalueUL1=phi(0)+0.000000001;
                        wealthUL1=findconsUL(0.0,igrid);

                        if (wealthUL1<=testvalueUL1 || wealthUL1 <= 0.0) // actual saving is higher than wealth (not possible)
                            {
                                testvalueUL1=wealthUL1-0.000000001;
                                icaseUL1=101;
                                
                                if (testvalueUL1<=0.0 || wealthUL1 <= 0.0)
                                {//wealth is so small that consumption is almost zero
                
                                
                                    // STARTING 2D searching effort
                                    paramsUU.ygridparXX = ygrid;
                                    paramsUU.typeUU = params.typeU;
                                
                                    
                                    testval1 = searchWwhenEknownfun(0.0, &paramsUU);
                                    testval2 = searchWwhenEknownfun(0.00001, &paramsUU);
                                    testval3 = searchWwhenEknownfun(EffortMaxW, &paramsUU);
                                    testval4 = searchWwhenEknownfun(EffortMaxW-0.00001, &paramsUU);
                                    
                                    iverif = 0;

                                    if(testval1 > testval2 & testval1 < 0.0) {
                                        searchvalWWUL1 = 0.0;
                                        iverif = 100;
                                        }
                                    if(testval1 < testval2 & testval1 > 0.0) {
                                        searchvalWWUL1 = 0.0;
                                        iverif = 100;
                                        }
                                    
                                    if(testval3 > testval4 & testval3 < 0.0) {
                                        searchvalWWUL1 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    if(testval3 < testval4 & testval3 > 0.0) {
                                        searchvalWWUL1 = EffortMaxW;
                                        printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                        iverif = 100;
                                        }
                                    
                                    if(iverif == 0) {
                                        searchvalWWUL1 = zbrentNEW(searchWwhenEknownfun,0.0,EffortMaxW,(1.0e-10),&paramsUU);
                                        searchvalVEUL1 = paramsUU.searchUsEoutST;
                                        }
                                    
                                 
                                    // COMPUTE THE CONTINUATION VALUE
                                    valfnmaxUL1=utilc(0.0000000001)+disutilityW(searchvalWWUL1)+disutilityE(searchvalVEUL1)+betapar*((1-zetapar)*(piWU(ygrid,searchvalWWUL1)*(1 - piEU(searchvalVEUL1))*tempWW1+piEU(searchvalVEUL1)*(1-piWU(ygrid, searchvalWWUL1))*tempUL1VE + piWU(ygrid, searchvalWWUL1)*piEU(searchvalVEUL1)*tempWW1VE + (1 - piEU(searchvalVEUL1))*(1 - piWU(ygrid, searchvalWWUL1))*tempUL1) + zetapar*((1 - piWU(ygrid, searchvalWWUL1))*tempUL0 + piWU(ygrid, searchvalWWUL1)*tempWW0));
                                    
                                    valfnmaxUL1=-valfnmaxUL1;
                                    amaxsaveUL1=Gridmin;
                                    icaseUL1=2;
                                }
                            }
                        
     
                        
                        if (icaseUL1>=100)
                        {
                            if (focUwID(Gridmin,&params)<focUwID(testvalueUL1,&params))
                                
                            {// corner solution at the bottom
                                valfnmaxUL1=focUwID(Gridmin,&params);
                                searchvalWWUL1=params.searchUWoutST;
                                searchvalVEUL1=params.searchUEoutST;
                                amaxsaveUL1=Gridmin;
                            }
                            else
                            {
                                savemaxUL1=min(Gridmax,findconsUL((0.0),igrid)-0.000000001);
                                
                                if (savemaxUL1>=(Gridmax))
                                {
                                    savemaxUL1=Gridmax-0.000000001;
                                }
                                
                                
                                isavemaxUL1=(int)floor(phiinv(savemaxUL1));
                                
                                if ((isavemaxUL1<0)){printf("isavemaxUL<=0 %d\n",isavemaxUL1);getchar();}
                                
                                
                                if (focUwID(phi((isavemaxUL1)),&params)>focUwID(savemaxUL1,&params))
                                {//corner solution at the top
                                    //printf("savemax\n");
                                    valfnmaxUL1=focUwID(savemaxUL1, &params);
                                    searchvalWWUL1=params.searchUWoutST;
                                    searchvalVEUL1=params.searchUEoutST;
                                    amaxsaveUL1=savemaxUL1;
                                    
                                }
                                else
                                {//cas standard
                                    
                                    //printf("stdcase\n");

                                    valfnmaxUL1=mygolden2search(Gridmin,Gridmin+0.000001,savemaxUL1,TOL,amaxsaveUL1,searchvalVEUL1, searchvalWWUL1,&params,focUwID);
                                    
                                    //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,amaxsaveUS,searchvalUS);
                                }
                                
                                
                            }
                        }
                        
                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\n",igrid,valfnmaxUS,0.0,searchvalUS);
                        valueUL1[inx(igrid, ygrid)]=-valfnmaxUL1;
                        saveUL1[inx(igrid, ygrid)]=amaxsaveUL1;
                        searchEffortWWUL1[inx(igrid, ygrid)]=searchvalWWUL1;
                        searchEffortVEUL1[inx(igrid, ygrid)]=searchvalVEUL1;
                      
                        
                        //printf("I AM HERE BEFORE VE %f %f \n", searchvalWWUL1, searchvalVEUL1); getchar();

 
 
                        ///////////////////////////////
                        // Solving the ENTREPRENEUR  //
                        ///////////////////////////////
                        for(egrid=0;egrid<maxfirmtype;egrid++) {
                        
                            params.egridST = egrid;
                            
                            #if BRC == 1
                            params.endoKST=endoKnext[inxKK(igrid, ygrid, egrid)];
                            #endif
                            
                            #if BRC == 2
                            params.endoKST=endoKnext[inxKK(igrid, ygrid, egrid)];
                            #endif
                            
                            
                            // ENTREPRENEUR WITHOUT AN ACCESS TO UI ///
                        
                            tempWW1VE2 = 0.0;
                            tempUL1VE2 = 0.0;
                            
                            for(y = 0;y<maxprod;y++) {
                                // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
                                //tempVE2 = 0.0;
                                for(e=0; e<maxfirmtype; e++) {
                                if(egrid == 1) {
                                    tempWW1VE2 += max(mprod[ygrid][y]*mggtrans[ygrid][e]*valueWW1next[inx(0, y)], mprod[ygrid][y]*mggtrans[ygrid][e]*valueVEnext[inxE(0, y, e)]); // don't have to repay the cost
                                    tempUL1VE2 += max(mprod[ygrid][y]*mggtrans[ygrid][e]*valueUL1next[inx(0, y)], mprod[ygrid][y]*mggtrans[ygrid][e]*valueVEnext[inxE(0, y, e)]); // don't have to repay the cost
                                }
                                if(egrid == 0) {
                                    tempWW1VE2 += max(mprod[ygrid][y]*valueWW1next[inx(0, y)], mprod[ygrid][y]*valueVEnext[inxE(0, y, 0)]); // don't have to repay the cost
                                    tempUL1VE2 += max(mprod[ygrid][y]*valueUL1next[inx(0, y)], mprod[ygrid][y]*valueVEnext[inxE(0, y, 0)]); // don't have to repay the cost
                                }
                                }
                            }
        
                        
                            // first test if phi(1)+0.00000001 exhibit or not positive consumption (in case when we test for corner solution).
                            icaseVE = 100;
                            
                            testvalueVE=phi(0)+0.0000001;
                            
                            //printf("I AM HERE %f", egrid); getchar();
                            
                            #if BRC == 1
                                wealthVE=findconsVE(0.0, igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)]);
                            #endif
                            
                            #if BRC == 2
                                wealthVE=findconsVE(0.0, igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)]);
                            #endif
                        
                            if(igrid == 0) {
                                valueVE[inxE(igrid, ygrid, egrid)]=-1000;
                                saveVE[inxE(igrid, ygrid, egrid)]=0.0;
                                searchEffortWWVE[inxE(igrid, ygrid, egrid)]=0.0;
                                } else {
                            if (wealthVE<=testvalueVE || wealthVE <= 0.0) // actual saving is higher than wealth (not possible)
                                {
                                    testvalueVE=wealthVE-0.000000001;
                                    icaseVE=101;
                                    
                                    if (testvalueVE<=0.0 || wealthVE <= 0.0)
                                    {//wealth is so small that consumption is almost zero
                                        
                                        paramsYY.diffValueNextfoc = max(0.000000000001, (betapar*((1-mupar)*(tempWW1VE2-tempUL1VE2) + mupar*(tempWW0 - tempUL0))));
                                        paramsYY.ygridparXX = ygrid;
                                        
                                        // Bracket the minimum (not sure to keep it) ---->
                                        testval1 = searchWfun(0.0, &paramsYY);
                                        testval2 = searchWfun(0.0001, &paramsYY);
                                        testval3 = searchWfun(EffortMaxW, &paramsYY);
                                        testval4 = searchWfun(EffortMaxW-0.0001, &paramsYY);
                                        
                                        iverif = 0;

                                        if(testval1 > testval2 & testval1 < 0.0) {
                                            searchvalWWVE = 0.0;
                                            iverif = 100;
                                            }
                                        if(testval1 < testval2 & testval1 > 0.0) {
                                            searchvalWWVE = 0.0;
                                            iverif = 100;
                                            }
                                        
                                        if(testval3 > testval4 & testval3 < 0.0) {
                                            searchvalWWVE = EffortMaxW;
                                            printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                            iverif = 100;
                                            }
                                        if(testval3 < testval4 & testval3 > 0.0) {
                                            searchvalWWVE = EffortMaxW;
                                            printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                            iverif = 100;
                                            }
                                        
                                        if(iverif == 0) {
                                            searchvalWWVE=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsYY);
                                            }
            
                                        valfnmaxVE=utilc(0.0000000001)+disutilityW(searchvalWWVE)+betapar*(piW(ygrid,searchvalWWVE)*(1.0 - mupar)*tempWW1VE2 + piW(ygrid, searchvalWWVE)*mupar*tempWW0 + mupar*(1 - piW(ygrid, searchvalWWVE))*tempUL0 + (1 - mupar)*(1-piW(ygrid, searchvalWWVE))*tempUL1VE2);
        
                                        valfnmaxVE=-valfnmaxVE;
                                        amaxsaveVE=Gridmin;
                                        icaseVE=2;
                                        //printf("icaseVE=%d\t%d\t%20.15f\t%20.15f\t%20.15f\n",icaseVE,igrid,wealthVE,testvalueVE,(phi(0)+0.000001));
                                    }
                                }
                            
                            
                            if (icaseVE>=100)
                            {
                                if (focVE(Gridmin,&params)<focVE(testvalueVE,&params))
                                {// corner solution at the bottom
                                    valfnmaxVE=focVE(Gridmin,&params);
                                    searchvalWWVE=params.searchoutST;
                                    amaxsaveVE=Gridmin;
                                    //printf("%20.15f %20.15f %d %f", focVE(Gridmin,&params), focVE(testvalueVE,&params), egrid, params.endoKST); getchar();
                                }
                                else
                                {//cas ou le max est a droite de la borne inf
                                    
                                    //printf("here\n");
                                    
                                    #if BRC == 1
                                    savemaxVE=min(Gridmax,findconsVE((0.0),igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)])-0.000000001);
                                    #endif 
                                    
                                    #if BRC == 2
                                    savemaxVE=min(Gridmax,findconsVE((0.0),igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)])-0.000000001);
                                    #endif
                                    
                                    if (savemaxVE>=(Gridmax))
                                    {
                                        savemaxVE=Gridmax-0.000000001;
                                    }
                                    
                                    
                                    isavemaxVE=(int)floor(phiinv(savemaxVE));
                                    
                                    if ((isavemaxVE<=0)){printf("isavemaxVE<=0 %d\n",isavemaxVE);getchar();}
                                    
                                    
                                    if (focVE(savemaxVE-0.0000001,&params)>focVE(savemaxVE,&params))
                                    {// corner solution at the top
                                        //printf("savemax\n");
                                        valfnmaxVE=focVE(savemaxVE, &params);
                                        searchvalWWVE=params.searchoutST;
                                        amaxsaveVE=savemaxVE;
                                        
                                    }
                                    else
                                    {//cas standard
                                        
                                        valfnmaxVE=mygolden1search(Gridmin,Gridmin+0.000001,savemaxVE,TOL,amaxsaveVE,searchvalWWVE,&params,focVE);
                                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,endoKnext[inxKK(igrid, ygrid, egrid)], savemaxVE, valfnmaxVE,amaxsaveVE,searchvalVE); getchar();
                                    }
                                    
                                     //if ((isavemaxVE<0)){printf("isavemaxVE<0 %d\t%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,isavemaxVE,savemaxVE,focVE(Gridmin,&params),focVE(testvalueVE,&params),focVE(phi((isavemaxVE)),&params),focVE(savemaxVE,&params));getchar();}
                                    
                                }
                            }
                            
                            //printf("%f", -valfnmaxVE);getchar();
                            valueVE[inxE(igrid, ygrid, egrid)]=-valfnmaxVE;
                            saveVE[inxE(igrid, ygrid, egrid)]=amaxsaveVE;
                            searchEffortWWVE[inxE(igrid, ygrid, egrid)]=searchvalWWVE;
                            }
                            
                            
                            //printf("I AM HERE AFTER VE %f  \n", searchvalWWVE); getchar();
                            
                            
                            // Compute the value function after paying the fixed cost to enter entrepreneurship
                            xonentry = max(grid[igrid] - costpar, 0.000000001);
                            
                            if(grid[igrid] - costpar <= 0.0){
                                valueVEc[inxE(igrid, ygrid, egrid)] = -1000;
                            }
                            else
                            {
                                xentry=phiinv(xonentry);
                                
                                ixentry=min((maxgrid-1),(int)(floor(xentry)));
                                ixentry=max(0,ixentry);
                                
                                if (ixentry>=(maxgrid-1))
                                {
                                    xentry=(maxgrid-1)-0.000000001;
                                    ixentry=min((maxgrid-1),(int)(floor(xentry)));
                                    ixentry=max(0,ixentry);
                                }
                                if (ixentry<=0)
                                {
                                    xentry=0.000000001;
                                    ixentry=min((maxgrid-1),(int)(floor(xentry)));
                                    ixentry=max(0,ixentry);
                                }
                                
                                // get the weight for linear interpolation.
                                dxentry	= (xonentry-grid[ixentry])/(grid[(ixentry+1)]-grid[ixentry]);
                                
                                if ((xonentry>(Gridmax+0.0000001))){printf("VFI: (xongridentry>(Gridmax)) %20.15f\n",xonentry);getchar();}
                                if ((ixentry>=(maxgrid-1))){printf("VFI: ixgridentry>=(maxgrid-1) %d\n",ixentry);getchar();}
                                if ((xonentry<(0.0))){printf("VFI: (xongridentry<(0.0)) %20.15f\n",xonentry);getchar();}
            
                                
                                valueVEc[inxE(igrid, ygrid, egrid)] = inter1d(dxentry, valueVE[inxE(ixentry,ygrid,egrid)], valueVE[inxE(ixentry + 1,ygrid,egrid)]);
                            }
                            
                            
                            
                            
                            #if SEA == 1
                            // ENTREPRENEUR WITH ACCESS TO UI ///
                            
                            tempWW1VEwUI2 = 0.0;
                            tempUL1VEwUI2 = 0.0;
                            tempUL1VE2 = 0.0;
                            
                            for(y = 0;y<maxprod;y++) {
                                // COMPUTE FUTURE STOCHASTIC RETURN OF BEING ENTREPRENEUR
                                //tempVE2 = 0.0;
                                for(e=0; e<maxfirmtype; e++) {
                                if(egrid == 1) {
                                    tempWW1VEwUI2 += max(mprod[ygrid][y]*mggtrans[ygrid][e]*valueWW1next[inx(0, y)], mprod[ygrid][y]*mggtrans[ygrid][e]*valueVEwUInext[inxE(0, y, e)]); // don't have to repay the cost
                                    tempUL1VEwUI2 += max(mprod[ygrid][y]*mggtrans[ygrid][e]*valueUL1next[inx(0, y)], mprod[ygrid][y]*mggtrans[ygrid][e]*valueVEwUInext[inxE(0, y, e)]); // don't have to repay the cost
                                }
                                if(egrid == 0) {
                                    tempWW1VEwUI2 += max(mprod[ygrid][y]*valueWW1next[inx(0, y)], mprod[ygrid][y]*valueVEwUInext[inxE(0, y, 0)]); // don't have to repay the cost
                                    tempUL1VEwUI2 += max(mprod[ygrid][y]*valueUL1next[inx(0, y)], mprod[ygrid][y]*valueVEwUInext[inxE(0, y, 0)]); // don't have to repay the cost
                                }
                                }
                            }
        
    
                        
                            // first test if phi(1)+0.00000001 exhibit or not positive consumption (in case when we test for corner solution).
                            icaseVEwUI = 100;
                            
                            testvalueVEwUI=phi(0)+0.0000001;
                            
                            //printf("I AM HERE %f", egrid); getchar();
                            
                            #if BRC == 1
                                wealthVEwUI=findconsVEwUI(0.0, igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)]);
                            #endif
                            
                            #if BRC == 2
                                wealthVEwUI=findconsVEwUI(0.0, igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)]);
                            #endif
                        
                            if(igrid == 0) {
                                valueVEwUI[inxE(igrid, ygrid, egrid)]=-1000;
                                saveVEwUI[inxE(igrid, ygrid, egrid)]=0.0;
                                searchEffortWWVEwUI[inxE(igrid, ygrid, egrid)]=0.0;
                                } else {
                            if (wealthVEwUI<=testvalueVEwUI || wealthVEwUI <= 0.0) // actual saving is higher than wealth (not possible)
                                {
                                    testvalueVEwUI=wealthVEwUI-0.000000001;
                                    icaseVEwUI=101;
                                    
                                    if (testvalueVEwUI<=0.0 || wealthVEwUI <= 0.0)
                                    {//wealth is so small that consumption is almost zero
                                        
    
                                        paramsYY.diffValueNextfoc = max(0.000000000001, (betapar*((1-mupar)*(PBnoUIVE*(tempWW1VE2) + (1 - PBnoUIVE)*(tempWW1VEwUI2) - PBnoUIVE*(tempUL1VE2) - (1 - PBnoUIVE)*tempUL1VEwUI2) + mupar*(tempWW0 - tempUS0))));
                                        paramsYY.ygridparXX = ygrid;
                                        
                                        // Bracket the minimum (not sure to keep it) ---->
                                        testval1 = searchWfun(0.0, &paramsYY);
                                        testval2 = searchWfun(0.0001, &paramsYY);
                                        testval3 = searchWfun(EffortMaxW, &paramsYY);
                                        testval4 = searchWfun(EffortMaxW-0.0001, &paramsYY);
                                        
                                        iverif = 0;

                                        if(testval1 > testval2 & testval1 < 0.0) {
                                            searchvalWWVEwUI = 0.0;
                                            iverif = 100;
                                            }
                                        if(testval1 < testval2 & testval1 > 0.0) {
                                            searchvalWWVEwUI = 0.0;
                                            iverif = 100;
                                            }
                                        
                                        if(testval3 > testval4 & testval3 < 0.0) {
                                            searchvalWWVEwUI = EffortMaxW;
                                            printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                            iverif = 100;
                                            }
                                        if(testval3 < testval4 & testval3 > 0.0) {
                                            searchvalWWVEwUI = EffortMaxW;
                                            printf("There is a bug:: EffortmaxW reached \n"); getchar();
                                            iverif = 100;
                                            }
                                        
                                        if(iverif == 0) {
                                            searchvalWWVEwUI=zbrentNEW(searchWfun,0.0,EffortMaxW,(1.0e-10),&paramsYY);
                                            }

                                        valfnmaxVEwUI=utilc(0.0000000001)+disutilityW(searchvalWWVEwUI)+betapar*(piW(ygrid,searchvalWWVEwUI)*(1.0 - mupar)*(PBnoUIVE*tempWW1VE2 + (1 - PBnoUIVE)*tempWW1VEwUI2) + piW(ygrid, searchvalWWVEwUI)*mupar*tempWW0 + mupar*(1 - piW(ygrid, searchvalWWVEwUI))*tempUS0 + (1 - mupar)*(1-piW(ygrid, searchvalWWVEwUI))*((1-PBnoUIVE)*tempUL1VEwUI2 + PBnoUIVE*tempUL1VE2));
        
                                        valfnmaxVEwUI=-valfnmaxVEwUI;
                                        amaxsaveVEwUI=Gridmin;
                                        icaseVEwUI=2;
                                        //printf("icaseVE=%d\t%d\t%20.15f\t%20.15f\t%20.15f\n",icaseVE,igrid,wealthVE,testvalueVE,(phi(0)+0.000001));
                                    }
                                }
                            
                             
                            if (icaseVEwUI>=100)
                            {
                                if (focVEwUI(Gridmin,&params)<focVEwUI(testvalueVEwUI,&params))
                                {// corner solution at the bottom
                                    valfnmaxVEwUI=focVEwUI(Gridmin,&params);
                                    searchvalWWVEwUI=params.searchoutST;
                                    amaxsaveVEwUI=Gridmin;
                                    //printf("%20.15f %20.15f %d %f", focVE(Gridmin,&params), focVE(testvalueVE,&params), egrid, params.endoKST); getchar();
                                }
                                else
                                {//cas ou le max est a droite de la borne inf
                                    
                                    //printf("here\n");
                                    
                                    #if BRC == 1
                                    savemaxVEwUI=min(Gridmax,findconsVEwUI((0.0),igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)])-0.000000001);
                                    #endif 
                                    
                                    #if BRC == 2
                                    savemaxVEwUI=min(Gridmax,findconsVEwUI((0.0),igrid, egrid, ygrid, endoKnext[inxKK(igrid, ygrid, egrid)])-0.000000001);
                                    #endif
                                    
                                    if (savemaxVEwUI>=(Gridmax))
                                    {
                                        savemaxVEwUI=Gridmax-0.000000001;
                                    }
                                    
                                    
                                    isavemaxVEwUI=(int)floor(phiinv(savemaxVEwUI));
                                    
                                    if ((isavemaxVEwUI<=0)){printf("isavemaxVEwUI<=0 %d\n",isavemaxVEwUI);getchar();}
                                    
                                    
                                    if (focVEwUI(savemaxVEwUI-0.0000001,&params)>focVEwUI(savemaxVEwUI,&params))
                                    {// corner solution at the top
                                        //printf("savemax\n");
                                        valfnmaxVEwUI=focVE(savemaxVEwUI, &params);
                                        searchvalWWVEwUI=params.searchoutST;
                                        amaxsaveVEwUI=savemaxVEwUI;
                                        
                                    }
                                    else
                                    {//cas standard
                                        
                                        valfnmaxVEwUI=mygolden1search(Gridmin,Gridmin+0.000001,savemaxVEwUI,TOL,amaxsaveVEwUI,searchvalWWVEwUI,&params,focVEwUI);
                                        //printf("mygolden %d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,endoKnext[inxKK(igrid, ygrid, egrid)], savemaxVE, valfnmaxVE,amaxsaveVE,searchvalVE); getchar();
                                    }
                                    
                                     //if ((isavemaxVE<0)){printf("isavemaxVE<0 %d\t%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,isavemaxVE,savemaxVE,focVE(Gridmin,&params),focVE(testvalueVE,&params),focVE(phi((isavemaxVE)),&params),focVE(savemaxVE,&params));getchar();}
                                    
                                }
                            }
                            
                            //printf("%f", -valfnmaxVE);getchar();
                            valueVEwUI[inxE(igrid, ygrid, egrid)]=-valfnmaxVEwUI;
                            saveVEwUI[inxE(igrid, ygrid, egrid)]=amaxsaveVEwUI;
                            searchEffortWWVEwUI[inxE(igrid, ygrid, egrid)]=searchvalWWVEwUI;
                            }
                            
                            
                            //printf("I AM HERE AFTER VE %f  \n", searchvalWWVE); getchar();
                            
                            
                            // Compute the value function after paying the fixed cost to enter entrepreneurship
                            xonentryVEwUI = max(grid[igrid] - costpar, 0.000000001);
                            
                            if(grid[igrid] - costpar <= 0.0){
                                valueVEwUIc[inxE(igrid, ygrid, egrid)] = -1000;
                            }
                            else
                            {
                                xentryVEwUI=phiinv(xonentryVEwUI);
                                
                                ixentryVEwUI=min((maxgrid-1),(int)(floor(xentryVEwUI)));
                                ixentryVEwUI=max(0,ixentryVEwUI);
                                
                                if (ixentryVEwUI>=(maxgrid-1))
                                {
                                    xentryVEwUI=(maxgrid-1)-0.000000001;
                                    ixentryVEwUI=min((maxgrid-1),(int)(floor(xentryVEwUI)));
                                    ixentryVEwUI=max(0,ixentryVEwUI);
                                }
                                if (ixentryVEwUI<=0)
                                {
                                    xentryVEwUI=0.000000001;
                                    ixentryVEwUI=min((maxgrid-1),(int)(floor(xentryVEwUI)));
                                    ixentryVEwUI=max(0,ixentryVEwUI);
                                }
                                
                                // get the weight for linear interpolation.
                                dxentryVEwUI	= (xonentryVEwUI-grid[ixentryVEwUI])/(grid[(ixentryVEwUI+1)]-grid[ixentryVEwUI]);
                                
                                if ((xonentryVEwUI>(Gridmax+0.0000001))){printf("VFI: (xongridentryVEwUI>(Gridmax)) %20.15f\n",xonentryVEwUI);getchar();}
                                if ((ixentryVEwUI>=(maxgrid-1))){printf("VFI: ixgridentryVEwUI>=(maxgrid-1) %d\n",ixentryVEwUI);getchar();}
                                if ((xonentryVEwUI<(0.0))){printf("VFI: (xongridentryVEwUI<(0.0)) %20.15f\n",xonentryVEwUI);getchar();}
            
                                
                                valueVEwUIc[inxE(igrid, ygrid, egrid)] = inter1d(dxentryVEwUI, valueVEwUI[inxE(ixentryVEwUI,ygrid,egrid)], valueVEwUI[inxE(ixentryVEwUI + 1,ygrid,egrid)]);
                            }
                            #endif

                        
                    
                            
                            
                        } // end egrid
                        
                        
                        
                        //printf("index: %d \n", inx(igrid, ygrid));

               }//for(igrid=0;igrid<maxgrid;igrid++)
            }//for(ygrid=0;ygrid<maxprod;ygrid++)
        #if OMP == 1
        } // end pragma
        #endif
        
#if indexPM == 0
            
        //convergence check
        critereVF=0.0;
        jcount=0;
        valmax=0.0;
        valmaxold=0.0;
        ivaltype=-1;

        for(ygrid=0; ygrid<maxprod;ygrid++)
        {
            for(igrid=0;igrid<maxgrid;igrid++)
            {
                for(egrid=0;egrid<maxfirmtype;egrid++) {
                    if ((fabs(valueVE[inxE(igrid, ygrid, egrid)]-valueVEnext[inxE(igrid, ygrid, egrid)]))>=critereVF){jcount=inxE(igrid, ygrid, egrid);ivaltype=1;valmax=valueVE[inxE(igrid, ygrid, egrid)];valmaxold=valueVEnext[inxE(igrid, ygrid, egrid)];}
                    critereVF=(max(fabs(valueVE[inxE(igrid, ygrid, egrid)]-valueVEnext[inxE(igrid, ygrid, egrid)]), critereVF));
                    if ((fabs(valueVEwUI[inxE(igrid, ygrid, egrid)]-valueVEwUInext[inxE(igrid, ygrid, egrid)]))>=critereVF){jcount=inxE(igrid, ygrid, egrid);ivaltype=7;valmax=valueVEwUI[inxE(igrid, ygrid, egrid)];valmaxold=valueVEwUInext[inxE(igrid, ygrid, egrid)];}
                    critereVF=(max(fabs(valueVEwUI[inxE(igrid, ygrid, egrid)]-valueVEwUInext[inxE(igrid, ygrid, egrid)]), critereVF));
                }
                if ((fabs(valueWW0[inx(igrid, ygrid)]-valueWW0next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=0;valmax=valueWW0[inx(igrid, ygrid)];valmaxold=valueWW0next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueWW0[inx(igrid, ygrid)]-valueWW0next[inx(igrid, ygrid)]), critereVF));
                if ((fabs(valueUS0[inx(igrid, ygrid)]-valueUS0next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=2;valmax=valueUS0[inx(igrid, ygrid)];valmaxold=valueUS0next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueUS0[inx(igrid, ygrid)]-valueUS0next[inx(igrid, ygrid)]), critereVF));
                if ((fabs(valueUL0[inx(igrid, ygrid)]-valueUL0next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=3;valmax=valueUL0[inx(igrid, ygrid)];valmaxold=valueUL0next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueUL0[inx(igrid, ygrid)]-valueUL0next[inx(igrid, ygrid)]), critereVF));
                if ((fabs(valueUS1[inx(igrid, ygrid)]-valueUS1next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=4;valmax=valueUS1[inx(igrid, ygrid)];valmaxold=valueUS1next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueUS1[inx(igrid, ygrid)]-valueUS1next[inx(igrid, ygrid)]), critereVF));
                if ((fabs(valueUL1[inx(igrid, ygrid)]-valueUL1next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=5;valmax=valueUL1[inx(igrid, ygrid)];valmaxold=valueUL1next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueUL1[inx(igrid, ygrid)]-valueUL1next[inx(igrid, ygrid)]), critereVF));
                if ((fabs(valueWW1[inx(igrid, ygrid)]-valueWW1next[inx(igrid, ygrid)]))>=critereVF){jcount=inx(igrid, ygrid);ivaltype=6;valmax=valueWW1[inx(igrid, ygrid)];valmaxold=valueWW1next[inx(igrid, ygrid)];}
                critereVF=(max(fabs(valueWW1[inx(igrid, ygrid)]-valueWW1next[inx(igrid, ygrid)]), critereVF));
            }
        }
            
            
        iter++;
  
        #if BRC == 1
        printf("Rules CNVG %d\t%20.15f\t%d\t%d\t%20.15f\t%20.15f\n",iter,critereVF,jcount,ivaltype,valmax,valmaxold);
        #endif
        #if BRC == 2
        printf("Rules CNVG %d\t%d\t%20.15f\t%d\t%d\t%20.15f\t%20.15f\n",iterK,iter,critereVF,jcount,ivaltype,valmax,valmaxold);
        #endif
            

        }//while (critereVF > epsilonValue)
#endif
        
        
        /////////////////////////////////////////
        //// ENDOGENOUS BORROWING CONSTRAINT ////
        /////////////////////////////////////////

#if BRC == 2
        critereendoK=0.0;
        

        for(ygrid=0; ygrid<maxprod;ygrid++)
        {
            for(igrid=0;igrid<maxgrid;igrid++)
            {
                // CASE 1: WHERE ENTREPRENEUR ALWAYS CHOOSE TO SAVE MAX AMOUNT (SAVE BRCtempUmin for test)
                investmax=Gridmax;
                investgrid=phiinv((investmax));
                
                ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                ixinvgrid=max(0,ixinvgrid);
                if (ixinvgrid>=(maxgrid-1))
                {
                    investgrid=(maxgrid-1)-0.000000001;
                    ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                    ixinvgrid=max(0,ixinvgrid);
                }
                if (ixinvgrid<=0)
                {
                    investgrid=0.000000001;
                    ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                    ixinvgrid=max(0,ixinvgrid);
                }
                
                if ((investmax>(Gridmax))){printf("vfirun: (yval>(Gridmax)) %20.15f\n",investmax);getchar();}
                if ((ixinvgrid>=(maxgrid-1))){printf("vfirun: iygrid>=(maxygrid-1) %d\n",ixinvgrid);getchar();}
                if ((investmax<(0.0))){printf("vfirun: (yval<(0.0)) %20.15f\n",investmax);getchar();}
                
                dinvestgrid=(investmax-grid[ixinvgrid])/(grid[(ixinvgrid+1)]-grid[ixinvgrid]);
                
                // CRITERION I: COME BACK TO Long-run UNEMPLOYED IF RUN AWAY
                BRCtempUmax=inter1d(dinvestgrid,(valueUL0[inx(ixinvgrid, ygrid)]),(valueUL0[inx((ixinvgrid+1), ygrid)])); // contemporaneous decision
                
                
                
                //CASE 2: IF ALWAYS RUN AWAY (SAVE BRCtempUmin for test)
                investmin=Gridmin;
                investgrid=phiinv((investmin));
                
                ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                ixinvgrid=max(0,ixinvgrid);
                if (ixinvgrid>=(maxgrid-1))
                {
                    investgrid=(maxgrid-1)-0.000000001;
                    ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                    ixinvgrid=max(0,ixinvgrid);
                }
                if (ixinvgrid<=0)
                {
                    investgrid=0.000000001;
                    ixinvgrid=min((maxgrid-1),(int)(floor(investgrid)));
                    ixinvgrid=max(0,ixinvgrid);
                }
                
                if ((investmin>(Gridmax))){printf("vfirun: (yval>(Gridmax)) %20.15f\n",investmin);getchar();}
                if ((ixinvgrid>=(maxgrid-1))){printf("vfirun: iygrid>=(maxygrid-1) %d\n",ixinvgrid);getchar();}
                if ((investmin<(0.0))){printf("vfirun: (yval<(0.0)) %20.15f\n",investmin);getchar();}
                
                dinvestgrid=(investmin-grid[ixinvgrid])/(grid[(ixinvgrid+1)]-grid[ixinvgrid]);
                
                BRCtempUmin=inter1d(dinvestgrid,(valueUL0[inx(ixinvgrid, ygrid)]),(valueUL0[inx((ixinvgrid+1), ygrid)])); // contemporaneous decision
               
               
               
                // TEST IN WHICH REGION WE ARE (recall : zbrent doesn't love corner solutions)
                
                if ((valueVE[inxE(igrid, ygrid, egrid)]) >= BRCtempUmax)
                {//these entrepreneurs never escape with borrowed money: we let them borrom maximum allowable amount
                    endoK[inxKK(igrid,ygrid, egrid)]=Gridmax/fpar;
                    endoK[inxKK(igrid,ygrid, egrid)]=relaxendoK*endoKnext[inxKK(igrid, ygrid, egrid)]+(1.0-relaxendoK)*endoK[inxKK(igrid,ygrid, egrid)]; // relaxation
                }
                else if ((valueVE[inxE(igrid, ygrid, egrid)]) <= BRCtempUmin)
                {
                    endoK[inxKK(igrid,ygrid, egrid)]=Gridmin;
                    endoK[inxKK(igrid,ygrid, egrid)]=relaxendoK*endoKnext[inxKK(igrid, ygrid, egrid)]+(1.0-relaxendoK)*endoK[inxKK(igrid,ygrid, egrid)]; // relaxation
                }
                else
                {// CASE3 : we look for the value with which the entrepreneur escapes
                    paramsKK.lgridKYY=igrid;
                    paramsKK.egridKYY=ygrid;
                    paramsKK.valueULKYY=valueUL0;
                    paramsKK.valueVEKYY=valueVE;
                    
                    endoK[inxKK(igrid,ygrid, egrid)]=zbrentNEW(searchendoK,Gridmin,(Gridmax/fpar),(1.0e-10),&paramsKK);
                    endoK[inxKK(igrid,ygrid, egrid)]=relaxendoK*endoKnext[inxKK(igrid, ygrid, egrid)]+(1.0-relaxendoK)*endoK[inxKK(igrid,ygrid, egrid)]; // relaxation
                }
                
                
                critereendoK=(max(fabs(endoK[inxKK(igrid,ygrid, egrid)]-endoKnext[inxKK(igrid, ygrid, egrid)]), critereendoK));
                
                if(critereendoK == 0.0) {
                    printf("endoK %d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,endoKnext[inxKK(igrid, ygrid, egrid)],endoK[inxKK(igrid,ygrid, egrid)],fpar*endoKnext[inxKK(igrid, ygrid, egrid)],fpar*endoK[inxKK(igrid,ygrid, egrid)]);
                    printf("endoK %d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",igrid,valueUL0[inx(igrid, ygrid)],valueVE[inxE(igrid, ygrid, egrid)],BRCtempUmax,BRCtempUmin); getchar();
                    }
                
            }
        }

        iterK++;
        
        printf("-------------------------------------------- CNVG %d\t%d\t%20.15f\t%20.15f\n",iterK,iter,critereVF,critereendoK);
        
        //getchar();
        
    }//while (critereendoK > epsilonendoK)
#endif
    

    // SAVING DECISION
    saveoutfile=fopen(saveout, "w");
    setbuf ( saveoutfile , NULL );
    
    
    fprintf(saveoutfile,"%5s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","PROD","ASSET","saveWW0","saveUS0","saveUL0","saveWW1","saveUS1","saveUL1");

    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(saveoutfile,"%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,grid[igrid],saveWW0[inx(igrid, ygrid)],saveUS0[inx(igrid, ygrid)],saveUL0[inx(igrid, ygrid)],saveWW1[inx(igrid, ygrid)],saveUS1[inx(igrid, ygrid)],saveUL1[inx(igrid, ygrid)]);
        }
    }

    fclose(saveoutfile);
    
    
    // SEARCH EFFORT DECISION (1D)
    searchoutfile=fopen(searchout, "w");
    setbuf ( searchoutfile , NULL );
    
    fprintf(searchoutfile,"%5s\t%20s\t%20s\t%20s\t%20s\n","PROD","ASSET","searchEffortWWUS0","searchEffortWWUL0","searchEffortVEWW");
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(searchoutfile,"%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,grid[igrid],searchEffortWWUS0[inx(igrid, ygrid)],searchEffortWWUL0[inx(igrid, ygrid)],searchEffortVEWW[inx(igrid, ygrid)]);
    
        }
    }
    fclose(searchoutfile);
    
    
    // SEARCH EFFORT DECISION (2D)
    searchoutfile=fopen(searchout2, "w");
    setbuf ( searchoutfile , NULL );
    
    fprintf(searchoutfile,"%5s\t%20s\t%20s\t%20s\t%20s\t%20s\n","PROD","ASSET","searchEffortVEUS1","searchEffortWWUS1","searchEffortVEUL1","searchEffortWWUL1");
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(searchoutfile,"%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,grid[igrid],searchEffortVEUS1[inx(igrid, ygrid)],searchEffortWWUS1[inx(igrid, ygrid)],searchEffortVEUL1[inx(igrid, ygrid)],searchEffortWWUL1[inx(igrid, ygrid)]);
    
        }
    }
    fclose(searchoutfile);
    
    
    // FOR ENTREPRENEUR CAPITAL
    endokoutfile=fopen(endokout, "w");
    setbuf ( endokoutfile , NULL );
    fprintf(endokoutfile,"%5s\t%5s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","TYPE","PROD","ASSET","ValueVE","ValueVEc", "SaveVE","endoKVE", "SearchWWVE","ValueVEwUI","ValueVEwUIc", "SaveVEwUI","SearchWWVEwUI");

    for(egrid=0;egrid<maxfirmtype;egrid++) {
        for(int ygrid = 0; ygrid < maxprod ; ygrid++)
        {
            for(igrid=0;igrid<maxgrid;igrid++)
            {
                #if BRC == 1
                fprintf(endokoutfile,"%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,egrid, grid[igrid],valueVE[inxE(igrid,ygrid,egrid)],valueVEc[inxE(igrid,ygrid,egrid)],saveVE[inxE(igrid,ygrid,egrid)],endoK[inxKK(igrid,ygrid,egrid)],searchEffortWWVE[inxE(igrid,ygrid,egrid)],valueVEwUI[inxE(igrid,ygrid,egrid)],valueVEwUIc[inxE(igrid,ygrid,egrid)],saveVEwUI[inxE(igrid,ygrid,egrid)],searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]);
                #endif
                
                #if BRC == 2
                fprintf(endokoutfile,"%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,egrid, grid[igrid],valueVE[inxE(igrid,ygrid,egrid)],valueVEc[inxE(igrid,ygrid,egrid)],saveVE[inxE(igrid,ygrid,egrid)],endoK[inxKK(igrid,ygrid,egrid)],searchEffortWWVE[inxE(igrid,ygrid,egrid)],valueVEwUI[inxE(igrid,ygrid,egrid)],valueVEwUIc[inxE(igrid,ygrid,egrid)],saveVEwUI[inxE(igrid,ygrid,egrid)],searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]);
                #endif
            }
        }
    }
    
    fclose(endokoutfile);
    
    
    // VALUE FUNCTIONS
    valuefnoutfile=fopen(valuefnout, "w");
    setbuf ( valuefnoutfile , NULL );
    fprintf(valuefnoutfile,"%5s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","PROD","ASSET","valueWW0","valueUS0","valueUL0","valueWW1","valueUS1","valueUL1");

    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(valuefnoutfile,"%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,grid[igrid],valueWW0[inx(igrid, ygrid)],valueUS0[inx(igrid, ygrid)],valueUL0[inx(igrid, ygrid)], valueWW1[inx(igrid, ygrid)],valueUS1[inx(igrid, ygrid)],valueUL1[inx(igrid, ygrid)]);
        }
    }
    
    fclose(valuefnoutfile);
    
    
#if indexPMwrite == 1
    // VALUE for projected methods:
    valuefnoutfile=fopen(valueVETfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++) {
            fprintf(valuefnoutfile,"%20.15f\n",valueVE[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    valuefnoutfile=fopen(valueVEwUITfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++) {
            fprintf(valuefnoutfile,"%20.15f\n",valueVEwUI[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    fclose(valuefnoutfile);
    
    valuefnoutfile=fopen(valueVEcTfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++) {
            fprintf(valuefnoutfile,"%20.15f\n",valueVEc[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    valuefnoutfile=fopen(valueVEwUIcTfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++) {
            fprintf(valuefnoutfile,"%20.15f\n",valueVEc[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    fclose(valuefnoutfile);
    
    valuefnoutfile=fopen(valueWW1Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueWW1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);
    
    valuefnoutfile=fopen(valueWW0Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueWW0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);

    valuefnoutfile=fopen(valueUS1Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueUS1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);

    valuefnoutfile=fopen(valueUL1Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueUL1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);

    valuefnoutfile=fopen(valueUS0Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueUS0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);

    valuefnoutfile=fopen(valueUL0Tfile, "w");
    setbuf ( valuefnoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(valuefnoutfile,"%20.15f\n",valueUL0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(valuefnoutfile);
#endif

    
    free(valueVEnext);
    free(valueVEcnext);
    free(valueVEwUInext);
    free(valueVEwUIcnext);
    free(valueWW0next);
    free(valueUS0next);
    free(valueUL0next);
    free(valueWW1next);
    free(valueUS1next);
    free(valueUL1next);
    free(endoKnext);

    
}

//////////////////////////////////
// END VALUE FUNCTION ITERATION //
//////////////////////////////////






//******************************//
//**********SIMULATION**********//
//******************************//


// get resid //
void getresid(double amylevel,double *residout,int *ixgrid)
{
    assert(residout);
    assert(ixgrid);
    
    *ixgrid=min((maxgrid-1),(int)(floor(phiinv((amylevel)))));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>(maxgrid-1))
	{
        printf("getresid: (ixgrid>(maxgrid-1))");
        getchar();
	}
    
    if (*ixgrid<0)
	{
        printf("getresid: (ixgrid<0)");
        getchar();
	}
    
    if ((amylevel>(Gridmax))){printf("bascule: (xval>(Gridmax)) %20.15f\n",amylevel);getchar();}
    if ((amylevel<(0.0))){printf("bascule: (xval<(0.0)) %20.15f\n",amylevel);getchar();}
    
    
    *residout=(amylevel-grid[*ixgrid])/(grid[(*ixgrid+1)]-grid[*ixgrid]);
}


// get resid WITH COST //
void getresidwcost(double amylevel, int *ixgrid)
{
    double amywcost;
    
    assert(ixgrid);
    
    amywcost = max(0.000000001,amylevel - costpar);
    
    *ixgrid=min((maxgrid-1),(int)(floor(phiinv(amywcost))));
    *ixgrid=max(0,*ixgrid);

}




// Need to change function value
                             
#if indexPM == 1
void SIMULATION(double *valueWW0, double *valueWW1, double *valueUS0, double *valueUS1, double *valueUL0, double *valueUL1, double *valueVE, double *valueVEc, double *valueVEwUI, double *valueVEwUIc, double *endoK,double *saveWW0, double *saveWW1, double *saveUS0, double *saveUS1, double *saveUL0, double *saveUL1, double *saveVE, double *saveVEwUI, double *searchEffortWWUS0, double *searchEffortWWUL0, double *searchEffortWWUS1, double *searchEffortVEUS1, double *searchEffortWWUL1, double *searchEffortVEUL1, double *searchEffortWWVE,double *searchEffortVEWW,double *searchEffortWWVEwUI,double *totalassets, double *totallabor, double *optimaltax, double *generatedmoment, double *start_distVE, double *start_distWW1, double *start_distWW0, double *start_distUL0, double *start_distUL1, double *start_distUS1, double *start_distUS0, double *start_distVEwUI,  double *distVE, double *distWW1, double *distWW0, double *distUL0, double *distUL1, double *distUS1, double *distUS0, double *fracVEE, double *fracUU)
{
#endif

// Need to change function value
#if indexPM == 0
void SIMULATION(double *valueWW0, double *valueWW1, double *valueUS0, double *valueUS1, double *valueUL0, double *valueUL1, double *valueVE, double *valueVEc, double *valueVEwUI, double *valueVEwUIc, double *endoK,double *saveWW0, double *saveWW1, double *saveUS0, double *saveUS1, double *saveUL0, double *saveUL1, double *saveVE, double *saveVEwUI, double *searchEffortWWUS0, double *searchEffortWWUL0, double *searchEffortWWUS1, double *searchEffortVEUS1, double *searchEffortWWUL1, double *searchEffortVEUL1, double *searchEffortWWVE,double *searchEffortVEWW,double *searchEffortWWVEwUI, double *totalassets, double *totallabor, double *optimaltax, double *generatedmoment, double *start_distVE, double *start_distWW1, double *start_distWW0, double *start_distUL0, double *start_distUL1, double *start_distUS1, double *start_distUS0, double *start_distVEwUI)
{
#endif

#if indexPM == 0
double *distVEwUI, *distVE, *distWW1, *distWW0, *distUL0, *distUL1, *distUS1, *distUS0;
#endif

    // Entrepreneur distribution
    double *distVEwUIold, *distVEold, *distassetVEwUI, *saveresVEwUI, fracVEwUI, *distassetVE, *saveresVE, fracVEtoWW, fracVEtoUL, fracVE;
    int *isaveVE, *isaveVEwUI;
    double tempVEWW1_0, tempVEUS1_0, tempVEUL1_0, tempVEWW1_1, tempVEUS1_1, tempVEUL1_1, tempVEwUIUS1_0, tempVEwUIWW1_0, tempVEwUIUS1_1, tempVEwUIWW1_1;
    
    // Worker distribution
    double *distWW, *distWWold, *distassetWW, fracWWtoUS, fracWWtoVE, fracWW;
    double *distWW0old, *distassetWW0, *saveresWW0;
    double *distWW1old, *distassetWW1, *saveresWW1;
    int *isaveWW0;
    int *isaveWW1, *isaveWW1c;
    
    // Short run unemployed
    double *distUS, *distUSold, *distassetUS, fracUStoVE, fracUStoWW, fracUStoUL, fracUS;
    double *distUS0old, *distassetUS0, *saveresUS0;
    double *distUS1old, *distassetUS1, *saveresUS1;
    int *isaveUS0;
    int *isaveUS1, *isaveUS1c;
    
    // long run unemployed
    double *distUL, *distULold, *distassetUL, fracULtoVE, fracULtoWW, fracUL;
    double *distUL0old, *distassetUL0, *saveresUL0;
    double *distUL1old, *distassetUL1, *saveresUL1;
    int *isaveUL0;
    int *isaveUL1, *isaveUL1c;
    
    // Statistics
    double *distVEtoUL, *disttotal, *distasset, *distWWtoVE, *distVEtoWW, *fracVEpability, *fracWWpability, *aggWWtoVE, *aggVEtoWW, *necessityUS, *necessityUL, *necessityWW, necessityWWtoVE, necessityUStoVE, necessityULtoVE, top001p, top01p, top05p, top10p, top20p,bot20p,bot40p,bot60p,bot80p, mworthVE,mworthWW,mworthUS,mworthUL,capitalVE,productionVE,totalcapital,corproduction,totalproduction, KY, GINIval, distWWearning[maxprod], GINIearn, welfareUS, welfareUL, welfareWW, welfareVE, avgfirmsize, welfaretot, fracVEtoUS;
    
    // taxes
    double expenditure, revenubasis, corprevenu, debt;
    
    // Simulation need
    double critDist, distval, distvalold, verifdist;
    int ygrid, igrid, egrid, k, itr;
    
    FILE *tempfileoutfile;
    FILE *distoutfile;
    
    // INITIALIZE POINTERS
    distasset = (double *) calloc((maxgrid), sizeof(double));
    disttotal = (double *) calloc((ifulldim), sizeof(double));
    necessityUS = (double *) calloc((ifulldim), sizeof(double));
    necessityUL = (double *) calloc((ifulldim), sizeof(double));
    necessityWW = (double *) calloc((ifulldim), sizeof(double));
    distWWtoVE = (double *) calloc((ifulldim), sizeof(double));
    distVEtoWW = (double *) calloc((ifulldimE), sizeof(double));
    distVEtoUL = (double *) calloc((ifulldimE), sizeof(double));
    aggWWtoVE = (double *) calloc((maxprod), sizeof(double));
    aggVEtoWW = (double *) calloc((maxprod), sizeof(double));
    fracWWpability = (double *) calloc((maxprod), sizeof(double));
    fracVEpability = (double *) calloc((maxprod), sizeof(double));
    
#if indexPM == 0
    distUL0 = (double *) calloc((ifulldim), sizeof(double));
    distUS0 = (double *) calloc((ifulldim), sizeof(double));
    distWW0 = (double *) calloc((ifulldim), sizeof(double));
    distUL1 = (double *) calloc((ifulldim), sizeof(double));
    distUS1 = (double *) calloc((ifulldim), sizeof(double));
    distWW1 = (double *) calloc((ifulldim), sizeof(double));
    distVE = (double *) calloc((ifulldimE), sizeof(double));
    distVEwUI = (double *) calloc((ifulldimE), sizeof(double));
#endif

    distassetUS = (double *) calloc((maxgrid), sizeof(double));
    distassetUL = (double *) calloc((maxgrid), sizeof(double));
    distassetWW = (double *) calloc((maxgrid), sizeof(double));
    distassetVE = (double *) calloc((maxgrid), sizeof(double));
    distassetVEwUI = (double *) calloc((maxgrid), sizeof(double));
    
    distUL0old = (double *) calloc((ifulldim), sizeof(double));
    distUS0old = (double *) calloc((ifulldim), sizeof(double));
    distWW0old = (double *) calloc((ifulldim), sizeof(double));
    distUL1old = (double *) calloc((ifulldim), sizeof(double));
    distUS1old = (double *) calloc((ifulldim), sizeof(double));
    distWW1old = (double *) calloc((ifulldim), sizeof(double));
    distVEold = (double *) calloc((ifulldimE), sizeof(double));
    distVEwUIold = (double *) calloc((ifulldimE), sizeof(double));
    
    saveresUS0 = (double *) calloc((ifulldim), sizeof(double));
    saveresUL0 = (double *) calloc((ifulldim), sizeof(double));
    saveresWW0 = (double *) calloc((ifulldim), sizeof(double));
    saveresUS1 = (double *) calloc((ifulldim), sizeof(double));
    saveresUL1 = (double *) calloc((ifulldim), sizeof(double));
    saveresWW1 = (double *) calloc((ifulldim), sizeof(double));
    saveresVE = (double *) calloc((ifulldimE), sizeof(double));
    saveresVEwUI = (double *) calloc((ifulldimE), sizeof(double));
    
    isaveUS0 = (int *) calloc((ifulldim), sizeof(int));
    isaveUL0 = (int *) calloc((ifulldim), sizeof(int));
    isaveWW0 = (int *) calloc((ifulldim), sizeof(int));
    isaveUS1 = (int *) calloc((ifulldim), sizeof(int));
    isaveUL1 = (int *) calloc((ifulldim), sizeof(int));
    isaveWW1 = (int *) calloc((ifulldim), sizeof(int));
    isaveUS1c = (int *) calloc((ifulldim), sizeof(int));
    isaveUL1c = (int *) calloc((ifulldim), sizeof(int));
    isaveWW1c = (int *) calloc((ifulldim), sizeof(int));
    isaveVE = (int *) calloc((ifulldimE), sizeof(int));
    isaveVEwUI = (int *) calloc((ifulldimE), sizeof(int));

    
    // GUESS FIRST DISTRIBUTION
    bascule(start_distVE, distVE, ifulldimE);
    bascule(start_distVEwUI, distVEwUI, ifulldimE);
    bascule(start_distWW1, distWW1, ifulldim);
    bascule(start_distWW0, distWW0, ifulldim);
    bascule(start_distUS1, distUS1, ifulldim);
    bascule(start_distUS0, distUS0, ifulldim);
    bascule(start_distUL1, distUL1, ifulldim);
    bascule(start_distUL0, distUL0, ifulldim);
    
    // Compute Residuals from moving from one point to the next in the distribution.
    
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++)
            {
                getresid(saveVE[inxE(igrid, ygrid, egrid)],&saveresVE[inxE(igrid, ygrid, egrid)],&isaveVE[inxE(igrid, ygrid, egrid)]);
                getresid(saveVEwUI[inxE(igrid, ygrid, egrid)],&saveresVEwUI[inxE(igrid, ygrid, egrid)],&isaveVEwUI[inxE(igrid, ygrid, egrid)]);
            }
            
            getresid(saveUS0[inx(igrid, ygrid)],&saveresUS0[inx(igrid, ygrid)],&isaveUS0[inx(igrid, ygrid)]);
            getresid(saveUL0[inx(igrid, ygrid)],&saveresUL0[inx(igrid, ygrid)],&isaveUL0[inx(igrid, ygrid)]);
            getresid(saveWW0[inx(igrid, ygrid)],&saveresWW0[inx(igrid, ygrid)],&isaveWW0[inx(igrid, ygrid)]);
            getresid(saveUS1[inx(igrid, ygrid)],&saveresUS1[inx(igrid, ygrid)],&isaveUS1[inx(igrid, ygrid)]);
            getresid(saveUL1[inx(igrid, ygrid)],&saveresUL1[inx(igrid, ygrid)],&isaveUL1[inx(igrid, ygrid)]);
            getresid(saveWW1[inx(igrid, ygrid)],&saveresWW1[inx(igrid, ygrid)],&isaveWW1[inx(igrid, ygrid)]);
            
            // IF Has to pay entry cost then:
            getresidwcost(saveUS1[inx(igrid, ygrid)],&isaveUS1c[inx(igrid, ygrid)]);
            getresidwcost(saveUL1[inx(igrid, ygrid)],&isaveUL1c[inx(igrid, ygrid)]);
            getresidwcost(saveWW1[inx(igrid, ygrid)],&isaveWW1c[inx(igrid, ygrid)]);
            
            if(saveresWW1[inx(igrid, ygrid)] < -0.000000001 || saveresWW1[inx(igrid, ygrid)] > 1.000000001 || saveresWW0[inx(igrid, ygrid)] < -0.000000001 || saveresWW0[inx(igrid, ygrid)] > 1.000000001 || saveresUS1[inx(igrid, ygrid)] < -0.000000001 || saveresUS1[inx(igrid, ygrid)] > 1.000000001 ) {
            printf("MISTAKE"); getchar();
            }
    //        printf("%d\t %d\t%20.15f\t%20.15f\n",ygrid,igrid,saveUS[inx(igrid, ygrid)],((saveUS[inx(igrid, ygrid)]-grid[(int)floor(phiinv(saveUS[inx(igrid, ygrid)]))])/(grid[(int)floor(phiinv(saveUS[inx(igrid, ygrid)]))+1]-grid[(int)floor(phiinv(saveUS[inx(igrid, ygrid)]))])));
    
    //        printf("%d\t %d\t%20.15f\t%20.15f\n",ygrid,igrid,saveUSc[inx(igrid, ygrid)],((saveUSc[inx(igrid, ygrid)]-grid[(int)floor(phiinv(saveUSc[inx(igrid, ygrid)]))])/(grid[(int)floor(phiinv(saveUSc[inx(igrid, ygrid)]))+1]-grid[(int)floor(phiinv(saveUSc[inx(igrid, ygrid)]))])));
            
        }
    }
    //getchar();
    
    
    critDist=1.0;
    itr=0;
    
#if indexPM == 0
    while((critDist>epsilonDist))  // start loop on distribution.
    {
#endif
        itr++;
        
        bascule(distUS0,distUS0old,ifulldim);
        bascule(distUL0,distUL0old,ifulldim);
        bascule(distWW0,distWW0old,ifulldim);
        bascule(distUS1,distUS1old,ifulldim);
        bascule(distUL1,distUL1old,ifulldim);
        bascule(distWW1,distWW1old,ifulldim);
        bascule(distVE,distVEold,ifulldimE);
        bascule(distVEwUI,distVEwUIold,ifulldimE);


    
    // INITIALIZE ALL DISTRIBUTION TO ZERO
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++)
            {
                verifdist += distVE[inxE(igrid, ygrid, egrid)];
                verifdist += distVEwUI[inxE(igrid, ygrid, egrid)];
                distVE[inxE(igrid, ygrid, egrid)]=0.0;
                distVEtoWW[inxE(igrid, ygrid, egrid)]=0.0;
                distVEwUI[inxE(igrid, ygrid, egrid)]=0.0;
            }
            verifdist += distUS0[inx(igrid, ygrid)];
            verifdist += distUL0[inx(igrid, ygrid)];
            verifdist += distWW0[inx(igrid, ygrid)];
            verifdist += distUS1[inx(igrid, ygrid)];
            verifdist += distUL1[inx(igrid, ygrid)];
            verifdist += distWW1[inx(igrid, ygrid)];

            distUS0[inx(igrid, ygrid)]=0.0;
            distUL0[inx(igrid, ygrid)]=0.0;
            distWW0[inx(igrid, ygrid)]=0.0;
            distUS1[inx(igrid, ygrid)]=0.0;
            distUL1[inx(igrid, ygrid)]=0.0;
            distWW1[inx(igrid, ygrid)]=0.0;
            distWWtoVE[inx(igrid, ygrid)]=0.0;
            necessityUS[inx(igrid, ygrid)]=0.0;
            necessityUL[inx(igrid, ygrid)]=0.0;
            necessityWW[inx(igrid, ygrid)]=0.0;
        }
    }
 
    

    // INITIALIZE STATISTICS
     fracWWtoVE=0.0;
     fracWWtoUS=0.0;
     fracUStoWW=0.0;
     fracUStoVE=0.0;
     fracVEtoUL=0.0;
     fracVEtoWW=0.0;
     fracULtoWW=0.0;
     fracULtoVE=0.0;
     fracVEtoUS = 0.0;
     necessityUStoVE=0.0;
     necessityULtoVE=0.0;
     necessityWWtoVE=0.0;
        
    // some verification
    distval=0.0;
    distvalold=0.0;
        
        
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            //(by default, we assume that WW is prefered to Entrepreneur which is prefered to Unemployment if indifferent)
            
            for(k = 0; k < maxprod; k++)
            {
            
            tempVEWW1_0 = 0.0;
            tempVEUS1_0 = 0.0;
            tempVEUL1_0 = 0.0;
            tempVEWW1_1 = 0.0;
            tempVEUS1_1 = 0.0;
            tempVEUL1_1 = 0.0;
            
            // TO COMPUTE SEA NEXT //
            tempVEwUIUS1_0 = 0.0;
            tempVEwUIWW1_0 = 0.0;
            tempVEwUIUS1_1 = 0.0;
            tempVEwUIWW1_1 = 0.0;
            
            for(int e = 0; e < maxfirmtype; e++) {
                tempVEWW1_0 += mgtrans[k][e]*valueVEc[inxE(isaveWW1[inx(igrid, ygrid)], k, e)];
                tempVEUS1_0 += mgtrans[k][e]*valueVEc[inxE(isaveUS1[inx(igrid, ygrid)], k, e)];
                tempVEUL1_0 += mgtrans[k][e]*valueVEc[inxE(isaveUL1[inx(igrid, ygrid)], k, e)];
                tempVEWW1_1 += mgtrans[k][e]*valueVEc[inxE(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k, e)];
                tempVEUS1_1 += mgtrans[k][e]*valueVEc[inxE(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k, e)];
                tempVEUL1_1 += mgtrans[k][e]*valueVEc[inxE(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k, e)];
                
                #if SEA == 1
                tempVEwUIUS1_0 += mgtrans[k][e]*valueVEwUIc[inxE(isaveUS1[inx(igrid, ygrid)], k, e)];
                tempVEwUIWW1_0 += mgtrans[k][e]*valueVEwUIc[inxE(isaveWW1[inx(igrid, ygrid)], k, e)];
                tempVEwUIUS1_1 += mgtrans[k][e]*valueVEwUIc[inxE(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k, e)];
                tempVEwUIWW1_1 += mgtrans[k][e]*valueVEwUIc[inxE(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k, e)];
                #endif
            }
            
            ////////////////////////////
            // TRANSITION OF A WORKER //
            ////////////////////////////
            
            // WITHOUT IDEA
            
            // TRANSITION WORKER0 TO UNEMPLOYED0 //
            distUS0[inx(isaveWW0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1-xipar)*etapar[ygrid]*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            distUS0[inx(min((maxgrid-1),(isaveWW0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1-xipar)*etapar[ygrid]*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(1-xipar)*etapar[ygrid]*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(1-xipar)*etapar[ygrid]*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            
            // TRANSITION WORKER0 TO WORKER0 //
            distWW0[inx(isaveWW0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1-xipar)*(1.0 - etapar[ygrid])*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveWW0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1-xipar)*(1.0 - etapar[ygrid])*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];

            // TRANSITION WORKER0 TO UNEMPLOYED1 //
            distUS1[inx(isaveWW0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(xipar)*etapar[ygrid]*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            distUS1[inx(min((maxgrid-1),(isaveWW0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(xipar)*etapar[ygrid]*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(xipar)*etapar[ygrid]*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(xipar)*etapar[ygrid]*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            
            // TRANSITION WORKER0 TO WORKER1 //
            distWW1[inx(isaveWW0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(xipar)*(1.0-etapar[ygrid])*(1.0 - saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveWW0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(xipar)*(1.0-etapar[ygrid])*(saveresWW0[inx(igrid, ygrid)])*distWW0old[inx(igrid, ygrid)];
            
            // WITH AN IDEA
            
            // DOES NOT LOSE THEIR IDEA //
            
            // TRANSITION WORKER1 TO ENTREPRENEUR //
            for(int e = 0; e < maxfirmtype; e++) {
            
                // IF DOES NOT FALL IN UNEMPLOYMENT.
                distVE[inxE(isaveWW1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 > (valueWW1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveWW1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 > (valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracWWtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 > (valueWW1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracWWtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 > (valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                distWWtoVE[inx(isaveWW1c[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 > (valueWW1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distWWtoVE[inx(min((maxgrid-1),(isaveWW1c[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 > (valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                // IF FALL
                #if SEA == 0
                distVE[inxE(isaveWW1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 >= (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveWW1c[inx(igrid, ygrid)]+1)), k,e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracUStoVE +=  mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 >= (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE +=  mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveWW1[inx(igrid, ygrid)], k))] >= (tempVEWW1_0)) & (tempVEWW1_0 >= (valueUS1[(inx(isaveWW1[inx(igrid, ygrid)], k))])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                    necessityUS[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] >= (tempVEWW1_1)) & (tempVEWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                    necessityUS[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                }
                #endif
                
                #if SEA == 1
                distVEwUI[inxE(isaveWW1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_0 >= (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVEwUI[inxE(min((maxgrid-1),(isaveWW1c[inx(igrid, ygrid)]+1)), k,e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_0 >= (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveWW1[inx(igrid, ygrid)], k))] >= (tempVEwUIWW1_0)) & (tempVEwUIWW1_0 >= (valueUS1[(inx(isaveWW1[inx(igrid, ygrid)], k))])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                    necessityUS[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] >= (tempVEwUIWW1_1)) & (tempVEwUIWW1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                    necessityUS[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
                }
                #endif
                
            }
            
            // TRANSITION WORKER1 TO WORKER1 //
            distWW1[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - zetapar)*(1 - etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - zetapar)*(1 - etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            
            distWW1[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 <= (valueWW1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - zetapar)*(1 - etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            
            
            
            // TRANSITION WORKER1 TO UNEMPLOYED1 //
            distUS1[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            distUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*(1.0 - piE(searchEffortVEWW[inx(igrid, ygrid)]))*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            
            #if SEA == 0
            distUS1[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 < (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 < (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_0 < (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEWW1_1 < (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            #if SEA == 1
            distUS1[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_0 < (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_1 < (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_0 < (valueUS1[inx(isaveWW1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracWWtoUS += mprod[ygrid][k]*(1 - zetapar)*(etapar[ygrid])*piE(searchEffortVEWW[inx(igrid, ygrid)])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)]*(((tempVEwUIWW1_1 < (valueUS1[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            // LOSE ITS IDEA
            
            // TRANSITION WORKER1 TO WORKER0 //
            distWW0[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(zetapar)*(1 - etapar[ygrid])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(zetapar)*(1 - etapar[ygrid])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            
            // TRANSITION WORKER1 TO US0 //
            distUS0[inx(isaveWW1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(zetapar)*(etapar[ygrid])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            distUS0[inx(min((maxgrid-1),(isaveWW1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(zetapar)*(etapar[ygrid])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(zetapar)*(etapar[ygrid])*(1.0-saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            fracWWtoUS += mprod[ygrid][k]*(zetapar)*(etapar[ygrid])*(saveresWW1[inx(igrid, ygrid)])*distWW1old[inx(igrid, ygrid)];
            
            
            /////////////////////////////////
            // TRANSITION OF AN UNEMPLOYED //
            /////////////////////////////////
            
            // SHORT-RUN WITHOUT AN IDEA //
            
            // TRANSITION US0 TO WW0 //
            distWW0[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            // TRANSITION US0 TO US0 //
            distUS0[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distUS0[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            // TRANSITION US0 TO UL0 //
            distUL0[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(pLTpar)*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distUL0[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(pLTpar)*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            // TRANSITION US0 TO US1 //
            distUS1[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distUS1[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            // TRANSITION US0 TO UL1 //
            distUL1[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(pLTpar)*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distUL1[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(pLTpar)*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            // TRANSITION US0 TO WW1 //
            distWW1[inx(isaveUS0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveUS0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(1.0 - saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUS0[inx(igrid,ygrid)]))*(saveresUS0[inx(igrid, ygrid)])*distUS0old[inx(igrid, ygrid)];
            
            
            // LONG-RUN WITHOUT AN IDEA //
            
            // TRANSITION UL0 TO WW0 //
            distWW0[inx(isaveUL0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveUL0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(1-xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            
            // TRANSITION UL0 TO UL0 //
            distUL0[inx(isaveUL0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            distUL0[inx(min((maxgrid-1),(isaveUL0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1-xipar)*(1.0 - piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            
            // TRANSITION UL0 TO UL1 //
            distUL1[inx(isaveUL0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            distUL1[inx(min((maxgrid-1),(isaveUL0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(xipar)*(1.0 - piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            
            // TRANSITION UL0 TO WW1 //
            distWW1[inx(isaveUL0[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveUL0[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(1.0 - saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(xipar)*(piWU(ygrid, searchEffortWWUL0[inx(igrid,ygrid)]))*(saveresUL0[inx(igrid, ygrid)])*distUL0old[inx(igrid, ygrid)];
            
            
            
            // SHORT-RUN WITH AN IDEA //
            
            // LOSING ITS IDEA //
            
            // TRANSITION US1 TO WW0 //
            distWW0[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            
            // TRANSITION US1 TO US0 //
            distUS0[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distUS0[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            
            // TRANSITION US1 TO UL0 //
            distUL0[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(pLTpar)*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distUL0[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(pLTpar)*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            
            // DOESN'T LOSE THEIR IDEA //
            
            // TRANSITION US1 TO WW1 //
            #if SEA == 0
            distWW1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            #if SEA == 1
            distWW1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1 - PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1 - PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(1 - PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(1 - PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            
            distWW1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 <= (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(PBnoUIVE)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            
            distWW1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            fracUStoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            
            
            for(int e = 0; e < maxfirmtype; e++) {
            
                // TRANSITION US1 TO VE (but doesn't in LT) //
                #if SEA == 0
                distVE[inxE(isaveUS1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 > (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveUS1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0)); // this take into account ST and LT unemployed since VWW > VUS > VUL
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - pLTpar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 > (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - pLTpar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                distVE[inxE(isaveUS1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 >= (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveUS1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 >= (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveUS1[inx(igrid, ygrid)], k))] >= (tempVEUS1_0)) & (tempVEUS1_0 >= (valueUS1[(inx(isaveUS1[inx(igrid, ygrid)], k))])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUS[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] >= (tempVEUS1_1)) & (tempVEUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUS[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                #endif
                
                
                #if SEA == 1
                distVEwUI[inxE(isaveUS1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 > (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVEwUI[inxE(min((maxgrid-1),(isaveUS1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0)); // this take into account ST and LT unemployed since VWW > VUS > VUL
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - pLTpar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 > (valueWW1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - pLTpar)*(piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                distVEwUI[inxE(isaveUS1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 >= (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVEwUI[inxE(min((maxgrid-1),(isaveUS1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 >= (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveUS1[inx(igrid, ygrid)], k))] >= (tempVEwUIUS1_0)) & (tempVEwUIUS1_0 >= (valueUS1[(inx(isaveUS1[inx(igrid, ygrid)], k))])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUS[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] >= (tempVEwUIUS1_1)) & (tempVEwUIUS1_1 >= (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))
                {
                    necessityUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUS[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*mgtrans[k][e]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                #endif
                
                
                // TRANSITION US1 TO VE (but fall in LT) //
                distVE[inxE(isaveUS1c[inx(igrid, ygrid)], k,e)] += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 >= (valueUL1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveUS1c[inx(igrid, ygrid)]+1)), k,e)] += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 >= (valueUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 >= (valueUL1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracUStoVE += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 >= (valueUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveUS1[inx(igrid, ygrid)], k))] >= (tempVEUS1_0)) & (tempVEUS1_0 >= (valueUL1[(inx(isaveUS1[inx(igrid, ygrid)], k))])))
                {
                    necessityULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUL[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] >= (tempVEUS1_1)) & (tempVEUS1_1 >= (valueUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))
                {
                    necessityULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                    necessityUL[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
                }
                
                
            }
            
            
            // TRANSITION US1 TO US1 //
            #if SEA == 0
            distUS1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 < (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 < (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            #if SEA == 1
            distUS1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_0 < (valueUS1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEwUIUS1_1 < (valueUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            #endif
            
            distUS1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distUS1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1 - pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            
            // TRANSITION US1 TO UL1 //
            distUL1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_0 < (valueUL1[inx(isaveUS1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)]*(((tempVEUS1_1 < (valueUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            
            distUL1[inx(isaveUS1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(1.0 - saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            distUL1[inx(min((maxgrid-1),(isaveUS1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(pLTpar)*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUS1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUS1[inx(igrid,ygrid)]))*(saveresUS1[inx(igrid, ygrid)])*distUS1old[inx(igrid, ygrid)];
            

            
            // LONG-RUN WITH AN IDEA //
            
            // LOSING ITS IDEA //
            
            // TRANSITION UL1 TO WW0 //
            distWW0[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            distWW0[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            
            // TRANSITION UL1 TO UL0 //
            distUL0[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            distUL0[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];

            
            // DOESN'T LOSE THEIR IDEA //
            
            // TRANSITION UL1 TO WW1 //
            distWW1[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 <= (valueWW1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            fracULtoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 <= (valueWW1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            fracULtoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 <= (valueWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            
            distWW1[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            distWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            fracULtoWW += mprod[ygrid][k]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            
            // TRANSITION UL1 TO VE //
            for(int e = 0; e < maxfirmtype; e++) {
                distVE[inxE(isaveUL1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 > (valueWW1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveUL1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 > (valueWW1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 > (valueWW1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                
                distVE[inxE(isaveUL1c[inx(igrid, ygrid)], k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 >= (valueUL1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                distVE[inxE(min((maxgrid-1),(isaveUL1c[inx(igrid, ygrid)]+1)), k, e)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 >= (valueUL1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                fracULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 >= (valueUL1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
                fracULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 >= (valueUL1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
                
                if((valueWW1[(inx(isaveUL1[inx(igrid, ygrid)], k))] >= (tempVEUL1_0)) & (tempVEUL1_0 >= (valueUL1[(inx(isaveUL1[inx(igrid, ygrid)], k))])))
                {
                    necessityULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
                    necessityUL[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
                }
                
                if((valueWW1[(inx(isaveUL1[inx(igrid, ygrid)], k))] >= (tempVEUL1_1)) & (tempVEUL1_1 >= (valueUL1[(inx(isaveUL1[inx(igrid, ygrid)], k))])))
                {
                    necessityULtoVE += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
                    necessityUL[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*mgtrans[k][e]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
                }
            }
            
            
            // TRANSITION UL1 TO UL1 //
            distUL1[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_0 < (valueUL1[inx(isaveUL1[inx(igrid, ygrid)], k)])))?(1.0):(0.0));
            distUL1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)]*(((tempVEUL1_1 < (valueUL1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)])))?(1.0):(0.0));
            
            distUL1[inx(isaveUL1[inx(igrid, ygrid)], k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(1.0 - saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            distUL1[inx(min((maxgrid-1),(isaveUL1[inx(igrid, ygrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - zetapar)*(1.0 - piWU(ygrid, searchEffortWWUL1[inx(igrid,ygrid)]))*(1.0 - piEU(searchEffortVEUL1[inx(igrid,ygrid)]))*(saveresUL1[inx(igrid, ygrid)])*distUL1old[inx(igrid, ygrid)];
            
            
            
                
            ///////////////////////////////////
            // TRANSITION OF AN ENTREPRENEUR //
            ///////////////////////////////////
            for(egrid = 0; egrid < maxfirmtype; egrid++)
            {
            
            
            // ENTREPRENEUR WITHOUT THE RIGHT OF UI //
            
            
                if(egrid == 1) {
                for(int e = 0; e < maxfirmtype; e++) {
                    // TRANSITION ENTREPRENEUR TO WORKER1 //
                    distWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // TRANSITION ENTREPRENEUR TO WORKER0 //
                    distWW0[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distWW0[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                
                    // TRANSITION ENTREPRENEUR TO ENTREPRENEUR //
                
                    distVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] > (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] > (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] >= (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] >= (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                
                
                    // TRANSITION ENTREPRENEUR TO UNEMPLOYED (LT) //
                    distUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    
                    distUL0[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distUL0[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoUL[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    }
                }
                if(egrid == 0) {
                    // TRANSITION ENTREPRENEUR TO WORKER1 //
                    distWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // TRANSITION ENTREPRENEUR TO WORKER0 //
                    distWW0[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distWW0[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                
                    // TRANSITION ENTREPRENEUR TO ENTREPRENEUR //
                
                    distVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] > (valueWW1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] > (valueWW1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] >= (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] >= (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                
                
                    // TRANSITION ENTREPRENEUR TO UNEMPLOYED (LT) //
                    distUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    
                    distUL0[inx(isaveVE[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distUL0[inx(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUL += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUL += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoUL[inxE(isaveVE[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVE[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVE[inxE(igrid,ygrid,egrid)]))*(saveresVE[inxE(igrid, ygrid, egrid)])*distVEold[inxE(igrid, ygrid, egrid)];
            }
                
                
            // ENTREPRENEUR WITH THE RIGHT OF UI //
            // I NEED TO DISTINGUISH BETWEEN THOSE FALLING IN NO UI AND THE OTHER //
            #if SEA == 1
                if(egrid == 1) {
                for(int e = 0; e < maxfirmtype; e++) {
                    // TRANSITION ENTREPRENEUR TO WORKER1 //
                    // if does not fall in no UI
                    distWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // if fall
                    distWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // TRANSITION ENTREPRENEUR TO WORKER0 //
                    distWW0[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distWW0[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                
                    // TRANSITION ENTREPRENEUR TO ENTREPRENEUR //
                    // if does not lose UI
                    distVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] > (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] > (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] >= (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] >= (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    //if lose
                    distVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] > (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] > (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] >= (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] >= (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                
                
                    // TRANSITION ENTREPRENEUR TO UNEMPLOYED (LT) //
                    // if does not lose UI
                    distUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    //if lose UI
                    distUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, e)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, e)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    
                    distUS0[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distUS0[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUS += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUS += mprod[ygrid][k]*mggtrans[ygrid][e]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    
                    }
                }
                
                if(egrid == 0) {
                    // TRANSITION ENTREPRENEUR TO WORKER1 //
                    // if does not lose his UI //
                    distWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // if does not lose his UI //
                    distWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoWW += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] <= (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] <= (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    
                    // TRANSITION ENTREPRENEUR TO WORKER0 //
                    distWW0[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distWW0[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoWW += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distVEtoWW[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(mupar)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                
                
                    // TRANSITION ENTREPRENEUR TO ENTREPRENEUR //
                    // if does not fall without UI
                    distVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] > (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] > (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] >= (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] >= (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // if fall without UI.
                    distVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] > (valueWW1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] > (valueWW1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] >= (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] >= (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                
                
                    // TRANSITION ENTREPRENEUR TO UNEMPLOYED (LT) //
                    // if have UI
                    distUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(1 - PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVEwUI[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    // if lose UI
                    distUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    fracVEtoUL += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    distVEtoUL[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)])))?(1.0):(0.0));
                    distVEtoUL[inxE(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k, 0)] += mprod[ygrid][k]*(1.0 - mupar)*(PBnoUIVE)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)]*(((valueVE[inxE(isaveVEwUI[inxE(igrid, ygrid, egrid)], k, 0)] < (valueUL1[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)])))?(1.0):(0.0));
                    
                    
                    distUS0[inx(isaveVEwUI[inxE(igrid, ygrid, egrid)], k)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    distUS0[inx(min((maxgrid-1),(isaveVEwUI[inxE(igrid, ygrid, egrid)]+1)), k)] += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUS += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(1.0 - saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];
                    fracVEtoUS += mprod[ygrid][k]*(mupar)*(1.0 - piW(ygrid, searchEffortWWVEwUI[inxE(igrid,ygrid,egrid)]))*(saveresVEwUI[inxE(igrid, ygrid, egrid)])*distVEwUIold[inxE(igrid, ygrid, egrid)];

                }
            #endif
                
                
            } // end egrid
            
            }
              
        }
    } // end computation of ergodic transition
    
    
    
    critDist=0.0;
    distval=0.0;
        
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            distval+=distUS0[inx(igrid, ygrid)];
            distval+=distUL0[inx(igrid, ygrid)];
            distval+=distWW0[inx(igrid, ygrid)];
            distval+=distUS1[inx(igrid, ygrid)];
            distval+=distUL1[inx(igrid, ygrid)];
            distval+=distWW1[inx(igrid, ygrid)];
            
            distvalold+=distUS0old[inx(igrid, ygrid)];
            distvalold+=distUL0old[inx(igrid, ygrid)];
            distvalold+=distWW0old[inx(igrid, ygrid)];
            distvalold+=distUS1old[inx(igrid, ygrid)];
            distvalold+=distUL1old[inx(igrid, ygrid)];
            distvalold+=distWW1old[inx(igrid, ygrid)];
            
            critDist=(max(fabs(distUS0[inx(igrid, ygrid)]-distUS0old[inx(igrid, ygrid)]), critDist));
            critDist=(max(fabs(distUL0[inx(igrid, ygrid)]-distUL0old[inx(igrid, ygrid)]), critDist));
            critDist=(max(fabs(distWW0[inx(igrid, ygrid)]-distWW0old[inx(igrid, ygrid)]), critDist));
            critDist=(max(fabs(distUS1[inx(igrid, ygrid)]-distUS1old[inx(igrid, ygrid)]), critDist));
            critDist=(max(fabs(distUL1[inx(igrid, ygrid)]-distUL1old[inx(igrid, ygrid)]), critDist));
            critDist=(max(fabs(distWW1[inx(igrid, ygrid)]-distWW1old[inx(igrid, ygrid)]), critDist));
            
            for(egrid=0; egrid<maxfirmtype; egrid++) {
                distval+=distVE[inxE(igrid, ygrid, egrid)];
                distvalold+=distVEold[inxE(igrid, ygrid, egrid)];
                critDist=(max(fabs(distVE[inxE(igrid, ygrid, egrid)]-distVEold[inxE(igrid, ygrid, egrid)]), critDist));
                
                #if SEA == 1
                distval+=distVEwUI[inxE(igrid, ygrid, egrid)];
                distvalold+=distVEwUIold[inxE(igrid, ygrid, egrid)];
                critDist=(max(fabs(distVEwUI[inxE(igrid, ygrid, egrid)]-distVEwUIold[inxE(igrid, ygrid, egrid)]), critDist));
                #endif
            }
        }
    }
    
    // VERIF
//    double distoldWW0, distnewWW0, distoldVE, distnewVE, distoldUL0, distnewUL0, distoldUS0, distnewUS0, distoldWW1, distnewWW1, distoldUL1, distnewUL1, distoldUS1, distnewUS1;
//    
//        distoldWW0 = 0.0;
//        distoldUS0 = 0.0;
//        distoldUL0 = 0.0;
//        distoldWW1 = 0.0;
//        distoldUS1 = 0.0;
//        distoldUL1 = 0.0;
//        distoldVE = 0.0;
//        distnewWW0 = 0.0;
//        distnewUS0 = 0.0;
//        distnewUL0 = 0.0;
//        distnewWW1 = 0.0;
//        distnewUS1 = 0.0;
//        distnewUL1 = 0.0;
//        distnewVE = 0.0;
//    
//    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
//    {
//        for(igrid = 0.0; igrid < maxgrid; igrid++)
//        {
//        distoldWW0 += distWW0old[inx(igrid, ygrid)];
//        distoldVE += distVEold[inxE(igrid, ygrid, egrid)];
//        distoldUS0 += distUS0old[inx(igrid, ygrid)];
//        distoldUL0 += distUL0old[inx(igrid, ygrid)];
//        distnewWW0 += distWW0[inx(igrid, ygrid)];
//        distnewVE += distVE[inxE(igrid, ygrid, egrid)];
//        distnewUS0 += distUS0[inx(igrid, ygrid)];
//        distnewUL0 += distUL0[inx(igrid, ygrid)];
//        distoldWW1 += distWW1old[inx(igrid, ygrid)];
//        distoldUS1 += distUS1old[inx(igrid, ygrid)];
//        distoldUL1 += distUL1old[inx(igrid, ygrid)];
//        distnewWW1 += distWW1[inx(igrid, ygrid)];
//        distnewUS1 += distUS1[inx(igrid, ygrid)];
//        distnewUL1 += distUL1[inx(igrid, ygrid)];
//        }
//    }
//    
//    //if(distval < distvalold || distval > distvalold)
//    if(itr >= 1)
//    {
//        printf("WARNING distval %d\t%20.15f\t%20.15f\t%20.15f\n",itr, distval,distvalold,critDist);
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldWW0 - distnewWW0,distoldWW0,distnewWW0);
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldUS0 - distnewUS0,distoldUS0,distnewUS0);//
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldUL0 - distnewUL0,distoldUL0,distnewUL0);
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldWW1 - distnewWW1,distoldWW1,distnewWW1);//
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldUS1 - distnewUS1,distoldUS1,distnewUS1);//
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldUL1 - distnewUL1,distoldUL1,distnewUL1);//
//        printf("WARNING distval %20.15f\t%20.15f\t%20.15f\n",distoldVE - distnewVE,distoldVE,distnewVE);
//        printf("WARNING distval %20.15f\n",distnewVE+distnewUL1+distnewUS1+distnewWW1+distnewUL0+distnewUS0+distnewWW0);
//        getchar();
//    }
//    
//        
//    printf("Dist CNVG %d\t%20.15f\t%20.15f\n",itr,critDist,distval);
//    

#if indexPM == 0
    }//while((critDist>epsilonDist))


    printf("Dist CNVG %d\t%20.15f\t%20.15f\n",itr,critDist,distval);
#endif
    
    ////////////////////////
    // COMPUTE STATISTICS //
    ////////////////////////
    
    // print all distribution
    distoutfile=fopen(distout, "w");
    setbuf ( distoutfile , NULL );
    
    
    fprintf(distoutfile,"%5s\t%5s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","PROD","ASSET","WW0","US0","UL0","WW1","US1","UL1");
    
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
        fprintf(distoutfile,"%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",ygrid,igrid, distWW0[inx(igrid, ygrid)],distUS0[inx(igrid, ygrid)],distUL0[inx(igrid, ygrid)],distWW1[inx(igrid, ygrid)],distUS1[inx(igrid, ygrid)],distUL1[inx(igrid, ygrid)]);
        }
    }
    
    fclose(distoutfile);
    
    
    distoutfile=fopen(distVEfile, "w");
    setbuf ( distoutfile , NULL );
    
    
    fprintf(distoutfile,"%5s\t%5s\t%5s\t%20s\t%20s\n","TYPE","PROD","ASSET","VE", "VEwUI");
    
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0; egrid<maxfirmtype; egrid++)
            {
                fprintf(distoutfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\n",egrid, ygrid,igrid, distVE[inxE(igrid, ygrid, egrid)], distVEwUI[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    fclose(distoutfile);
    

    
    // compute distribution of wealth + compute total investment
    *totalassets= 0.0;
    *totallabor = 0.0;
    capitalVE=0.0;
    productionVE=0.0;
    fracUS = 0.0;
    fracUL = 0.0;
    fracWW = 0.0;
    fracVE = 0.0;
    expenditure = 0.0;
    revenubasis = 0.0;
    corprevenu  = 0.0;
    welfareUS = 0.0;
    welfareUL = 0.0;
    welfareWW = 0.0;
    welfareVE = 0.0;
    avgfirmsize = 0.0;
    welfaretot = 0.0;
    
    
    for(ygrid=0; ygrid<maxprod; ygrid++) {
        aggVEtoWW[ygrid] = 0.0;
        aggWWtoVE[ygrid] = 0.0;
        fracVEpability[ygrid] = 0.0;
        fracWWpability[ygrid] = 0.0;
        distWWearning[ygrid] = 0.0;
        }
    
    for(igrid=0;igrid<maxgrid;igrid++)
    {
    distasset[igrid]+=0.0;
    distassetVE[igrid]+=0.0;
    distassetWW[igrid]+=0.0;
    distassetUS[igrid]+=0.0;
    distassetUL[igrid]+=0.0;

        
        for(ygrid=0;ygrid<maxprod;ygrid++)
        {
        disttotal[inx(igrid, ygrid)]=0.0;
        
        for(egrid=0;egrid<maxfirmtype; egrid++) {
            disttotal[inx(igrid, ygrid)]+=distVE[inxE(igrid, ygrid, egrid)];
            fracVE+=distVE[inxE(igrid, ygrid, egrid)];
            distassetVE[igrid] += distVE[inxE(igrid, ygrid, egrid)];
            fracVEpability[ygrid] += distVE[inxE(igrid, ygrid, egrid)];
            aggVEtoWW[ygrid] += distVEtoWW[inxE(igrid, ygrid, egrid)];
            welfareVE += distVE[inxE(igrid, ygrid, egrid)]*valueVE[inxE(igrid, ygrid, egrid)];
            avgfirmsize += distVE[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid,ygrid, egrid)];
            
            #if BRC == 1
            *totalassets += distVE[inxE(igrid, ygrid, egrid)]*(grid[igrid]-endoK[inxKK(igrid,ygrid, egrid)]);
            capitalVE+=distVE[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid,ygrid, egrid)];
            productionVE+=distVE[inxE(igrid, ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid,ygrid, egrid)], nupar);
            corprevenu += distVE[inxE(igrid,ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid,ygrid, egrid)], nupar)*corptax;
            #endif
            
            #if BRC == 2
            *totalassets += distVE[inxE(igrid, ygrid, egrid)]*(grid[igrid]-endoK[inxKK(igrid, ygrid, egrid)]);
            capitalVE+=distVE[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid, ygrid, egrid)];
            productionVE+=distVE[inxE(igrid, ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid, ygrid, egrid)], nupar);
            corprevenu += distVE[inxE(igrid,ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid, ygrid, egrid)], nupar)*corptax;
            #endif
            
            #if SEA == 1
            disttotal[inx(igrid, ygrid)]+=distVEwUI[inxE(igrid, ygrid, egrid)];
            fracVE+=distVEwUI[inxE(igrid, ygrid, egrid)];
            distassetVE[igrid] += distVEwUI[inxE(igrid, ygrid, egrid)];
            fracVEpability[ygrid] += distVEwUI[inxE(igrid, ygrid, egrid)];
            welfareVE += distVEwUI[inxE(igrid, ygrid, egrid)]*valueVEwUI[inxE(igrid, ygrid, egrid)];
            avgfirmsize += distVEwUI[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid,ygrid, egrid)];
            
            #if BRC == 1
            *totalassets += distVEwUI[inxE(igrid, ygrid, egrid)]*(grid[igrid]-endoK[inxKK(igrid,ygrid, egrid)]);
            capitalVE+=distVEwUI[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid,ygrid, egrid)];
            productionVE+=distVEwUI[inxE(igrid, ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid,ygrid, egrid)], nupar);
            corprevenu += distVEwUI[inxE(igrid,ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid,ygrid, egrid)], nupar)*corptax;
            #endif
            
            #if BRC == 2
            *totalassets += distVEwUI[inxE(igrid, ygrid, egrid)]*(grid[igrid]-endoK[inxKK(igrid, ygrid, egrid)]);
            capitalVE+=distVEwUI[inxE(igrid, ygrid, egrid)]*endoK[inxKK(igrid, ygrid, egrid)];
            productionVE+=distVEwUI[inxE(igrid, ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid, ygrid, egrid)], nupar);
            corprevenu += distVEwUI[inxE(igrid,ygrid, egrid)]*TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid, ygrid, egrid)], nupar)*corptax;
            #endif
            #endif
        }
        
        disttotal[inx(igrid, ygrid)]+=distUS0[inx(igrid, ygrid)];
        disttotal[inx(igrid, ygrid)]+=distWW0[inx(igrid, ygrid)];
        disttotal[inx(igrid, ygrid)]+=distUL0[inx(igrid, ygrid)];
        disttotal[inx(igrid, ygrid)]+=distUS1[inx(igrid, ygrid)];
        disttotal[inx(igrid, ygrid)]+=distWW1[inx(igrid, ygrid)];
        disttotal[inx(igrid, ygrid)]+=distUL1[inx(igrid, ygrid)];
        
        *totalassets += distUS0[inx(igrid, ygrid)]*grid[igrid];
        *totalassets += distWW0[inx(igrid, ygrid)]*grid[igrid];
        *totalassets += distUL0[inx(igrid, ygrid)]*grid[igrid];
        *totalassets += distUS1[inx(igrid, ygrid)]*grid[igrid];
        *totalassets += distWW1[inx(igrid, ygrid)]*grid[igrid];
        *totalassets += distUL1[inx(igrid, ygrid)]*grid[igrid];
        
        *totallabor+=distWW0[inx(igrid,ygrid)]*prod[ygrid];
        *totallabor+=distWW1[inx(igrid,ygrid)]*prod[ygrid];
        
        distasset[igrid]+=disttotal[inx(igrid, ygrid)];
        
        fracUS+=distUS0[inx(igrid, ygrid)];
        fracWW+=distWW0[inx(igrid, ygrid)];
        fracUL+=distUL0[inx(igrid, ygrid)];
        fracUS+=distUS1[inx(igrid, ygrid)];
        fracWW+=distWW1[inx(igrid, ygrid)];
        fracUL+=distUL1[inx(igrid, ygrid)];
        
        distassetWW[igrid] += distWW0[inx(igrid, ygrid)];
        distassetUS[igrid] += distUS0[inx(igrid, ygrid)];
        distassetUL[igrid] += distUL0[inx(igrid, ygrid)];
        distassetWW[igrid] += distWW1[inx(igrid, ygrid)];
        distassetUS[igrid] += distUS1[inx(igrid, ygrid)];
        distassetUL[igrid] += distUL1[inx(igrid, ygrid)];
        
        welfareWW += distWW0[inx(igrid, ygrid)]*valueWW0[inx(igrid,ygrid)] + distWW1[inx(igrid, ygrid)]*valueWW1[inx(igrid,ygrid)];
        welfareUS += distUS0[inx(igrid, ygrid)]*valueUS0[inx(igrid,ygrid)] + distUS1[inx(igrid, ygrid)]*valueUS1[inx(igrid,ygrid)];
        welfareUL += distUL0[inx(igrid, ygrid)]*valueUL0[inx(igrid,ygrid)] + distUL1[inx(igrid, ygrid)]*valueUL1[inx(igrid,ygrid)];
        
        aggWWtoVE[ygrid] += distWWtoVE[inx(igrid, ygrid)];
        fracWWpability[ygrid] += distWW1[inx(igrid, ygrid)] + distWW0[inx(igrid, ygrid)];
        
        revenubasis += distWW0[inx(igrid,ygrid)]*wstar*prod[ygrid];
        expenditure += distUS0[inx(igrid,ygrid)]*wstar*rhostar*prod[ygrid];
        expenditure += distUL0[inx(igrid,ygrid)]*minpar;
        revenubasis += distWW1[inx(igrid,ygrid)]*wstar*prod[ygrid];
        expenditure += distUS1[inx(igrid,ygrid)]*wstar*rhostar*prod[ygrid];
        expenditure += distUL1[inx(igrid,ygrid)]*minpar;
        
        distWWearning[ygrid] += distWW0[inx(igrid,ygrid)] + distWW1[inx(igrid,ygrid)];
    
        }
    }
    welfaretot = welfareVE + welfareUL + welfareUS + welfareWW;
    
    welfareVE = welfareVE/fracVE;
    welfareWW = welfareWW/fracWW;
    welfareUS = welfareUS/fracUS;
    welfareUL = welfareUL/fracUL;
    
    avgfirmsize = avgfirmsize/fracVE;
    
#if indexPM == 1
*fracVEE = fracVE;
*fracUU = fracUL + fracUS;
#endif

    // COMPUTE AGGREGATES
    totalcapital=*totalassets+capitalVE;
    corproduction=pow((*totalassets),(alphapar))*pow((*totallabor), (1.0-alphapar));
    totalproduction=corproduction+productionVE;
    KY=totalcapital/totalproduction;
    
    
    
    // COMPUTE POOL OF ENTREPRENEURS
    double ENT1, ENT2, ENT3, SUM1, SUM2, SUM3, AVGFIRM1,AVGFIRM2, AVGFIRM3,AVGASSET1,AVGASSET2,AVGASSET3,ENT1frac,ENT2frac,ENT3frac;
    ENT1 = 0; ENT2 = 0; ENT3 = 0;
    SUM1 = 0; SUM2 = 0; SUM3 = 0;
    AVGFIRM1 = 0; AVGFIRM2 = 0; AVGFIRM3 = 0;
    AVGASSET1 = 0; AVGASSET2 = 0; AVGASSET3 = 0;

    #if SEA == 1
    for(igrid = 0; igrid < maxgrid; igrid++) {
        for(egrid = 0; egrid < maxfirmtype; egrid++) {
            ENT1 += distVE[inxE(igrid,0,egrid)] + distVEwUI[inxE(igrid,0,egrid)];
            ENT2 += distVE[inxE(igrid,1,egrid)] + distVEwUI[inxE(igrid,1,egrid)];
            ENT3 += distVE[inxE(igrid,2,egrid)] + distVEwUI[inxE(igrid,2,egrid)];
            
            SUM1 += distVE[inxE(igrid,0,egrid)] + distVEwUI[inxE(igrid,0,egrid)];
            SUM2 += distVE[inxE(igrid,1,egrid)] + distVEwUI[inxE(igrid,1,egrid)];
            SUM3 += distVE[inxE(igrid,2,egrid)] + distVEwUI[inxE(igrid,2,egrid)];
        }
        SUM1 += distWW1[inx(igrid,0)] + distWW0[inx(igrid,0)] +distUS1[inx(igrid,0)] +distUS0[inx(igrid,0)] +distUL1[inx(igrid,0)] +distUL0[inx(igrid,0)] ;
        SUM2 += distWW1[inx(igrid,1)] + distWW0[inx(igrid,1)] +distUS1[inx(igrid,1)] +distUS0[inx(igrid,1)] +distUL1[inx(igrid,1)] +distUL0[inx(igrid,1)] ;
        SUM3 += distWW1[inx(igrid,2)] + distWW0[inx(igrid,2)] +distUS1[inx(igrid,2)] +distUS0[inx(igrid,2)] +distUL1[inx(igrid,2)] +distUL0[inx(igrid,2)] ;
    }
    
    ENT1frac = ENT1/SUM1;
    ENT2frac = ENT2/SUM2;
    ENT3frac = ENT3/SUM3;
    
    for(igrid = 0; igrid < maxgrid; igrid++) {
        for(egrid = 0; egrid < maxfirmtype; egrid++) {
            AVGFIRM1 += (distVEwUI[inxE(igrid,0,egrid)] + distVE[inxE(igrid,0,egrid)])*endoK[inxE(igrid,0,egrid)];
            AVGFIRM2 += (distVEwUI[inxE(igrid,1,egrid)] + distVE[inxE(igrid,1,egrid)])*endoK[inxE(igrid,1,egrid)];
            AVGFIRM3 += (distVEwUI[inxE(igrid,2,egrid)] + distVE[inxE(igrid,2,egrid)])*endoK[inxE(igrid,2,egrid)];
            AVGASSET1 += (distVEwUI[inxE(igrid,1,egrid)] + distVE[inxE(igrid,0,egrid)])*grid[igrid];
            AVGASSET2 += (distVEwUI[inxE(igrid,1,egrid)] + distVE[inxE(igrid,1,egrid)])*grid[igrid];
            AVGASSET3 += (distVEwUI[inxE(igrid,2,egrid)] + distVE[inxE(igrid,2,egrid)])*grid[igrid];
        }
    }
    
    AVGFIRM1 = AVGFIRM1/ENT1;
    AVGFIRM2 = AVGFIRM2/ENT2;
    AVGFIRM3 = AVGFIRM3/ENT3;
    
    AVGASSET1 = AVGASSET1/ENT1;
    AVGASSET2 = AVGASSET2/ENT2;
    AVGASSET3 = AVGASSET3/ENT3;
    #endif
    
    
    

    // PRINT DISTRIBUTIONS
    distoutfile=fopen(distasset2, "w");
    setbuf ( distoutfile , NULL );
    
    fprintf(distoutfile,"%5s\t%20s\t%20s\t%20s\t%20s\t%20s\n","ASSET", "all", "VE","WW","US","UL");
    
    for(igrid=0;igrid<maxgrid;igrid++)
    {
        fprintf(distoutfile,"%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",grid[igrid],distasset[igrid],distassetVE[igrid],distassetWW[igrid],distassetUS[igrid],distassetUL[igrid]);
    }
    
    fclose(distoutfile);
    
    
    // PRINT TRANSITION DISTRIBUTION VE TO WW AND WW TO VE
    distoutfile=fopen(disttransitionVEtoWW, "w");
    setbuf ( distoutfile , NULL );
    
    fprintf(distoutfile,"%5s\t%5s\t%20s\t%20s\n","TYPE", "PROD", "ASSET", "VEtoWW");
    
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++)
            {
                fprintf(distoutfile,"%5d\t%5d\t%20.15f\t%20.15f\n",egrid, ygrid, grid[igrid], distVEtoWW[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    fclose(distoutfile);
    
    distoutfile=fopen(disttransitionWWtoVE, "w");
    setbuf ( distoutfile , NULL );
    
    fprintf(distoutfile,"%5s\t%5s\t%20s\n","PROD", "ASSET", "WWtoVE");
    
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(distoutfile,"%5d\t%20.15f\t%20.15f\n",ygrid, grid[igrid], distWWtoVE[inx(igrid, ygrid)]);
        }
    }
    
    fclose(distoutfile);
    
    
    // PRINT NECESSITY SHARE DISTRIBUTION UL US
    distoutfile=fopen(distnecessity, "w");
    setbuf ( distoutfile , NULL );
    
    fprintf(distoutfile,"%5s\t%5s\t%20s\t%20s\n","PROD", "ASSET", "US", "UL");
    
    for(ygrid=0;ygrid<maxprod;ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            fprintf(distoutfile,"%5d\t%20.15f\t%20.15f\t%20.15f\n",ygrid, grid[igrid], necessityUS[inx(igrid, ygrid)], necessityUL[inx(igrid, ygrid)]);
        }
    }
    
    fclose(distoutfile);


    
    // COMPUTE DIFFERENT TOP BOTTOM OF THE DISTRIBUTION
    GINIval = GINI(grid, distasset, maxgrid);
    top001p = 1 - toppercent(grid, distasset, maxgrid, 0.999);
    top01p = 1 - toppercent(grid, distasset, maxgrid, 0.99);
    top05p = 1 - toppercent(grid, distasset, maxgrid, 0.95);
    top10p = 1 - toppercent(grid, distasset, maxgrid, 0.90);
    top20p = 1 - toppercent(grid, distasset, maxgrid, 0.80);
    bot20p = toppercent(grid, distasset, maxgrid, 0.20);
    bot40p = toppercent(grid, distasset, maxgrid, 0.40);
    bot60p = toppercent(grid, distasset, maxgrid, 0.60);
    bot80p = toppercent(grid, distasset, maxgrid, 0.80);

                                                                                                                     
    // COMPUTE MEDIAN WORTH
    if (fracVE>0.0){medianworth(distassetVE,fracVE, maxgrid, &mworthVE);}else{mworthVE=0.0;}
    if (fracWW>0.0){medianworth(distassetWW,fracWW, maxgrid, &mworthWW);}else{mworthWW=0.0;}
    if (fracUS>0.0){medianworth(distassetUS,fracUS, maxgrid, &mworthUS);}else{mworthUS=0.0;}
    if (fracUL>0.0){medianworth(distassetUL,fracUL, maxgrid, &mworthUL);}else{mworthUL=0.0;}
    
    
    //compute the distribution of earnings (gini around 0.36)
    GINIearn = GINI(prod, distWWearning, maxprod);
    

    
    // COMPUTE INCOME MATRIX + tax revenu
//    totrevenue = 0.0;
//    totincome = 0.0;
//    for(j=0;j<maxgrid;j++)
//    {
//        for(ai=0;ai<maxprod;ai++)
//        {
//            incomet[inx(igrid, ygrid)][0] = vEEability[ai]*pow(investYY[indexY],rend) -  depre*investYY[indexY] - dintrate*(investYY[indexY] - phit[j]);
//            incomet[inx(igrid, ygrid)][1] = distinYY[indexY];
//            
//            }
//        }
//    }
//		
//    qsort(incomet, incomegrid, 2*sizeof(double), comparefun2);
//    qsort(incometax, incomegrid, sizeof(double), comparefun2);
//    ginigeninc = GINI(incomelvl, incomedist, incomegrid);
    
    
    ///////////////////////////////////////
    // COMPUTE TAXES AND GOVERNMENT DEBT //
    ///////////////////////////////////////
    
    #if debtadj == 0
    *optimaltax = (expenditure-corprevenu)/revenubasis;
    //printf("TAXXX : %f %f %f %f \n", optimaltax, expenditure, revenubasis, corprevenu);
    #endif
    
    #if debtadj == 1
    *optimaltax = taxdebtfixed;
    debt = (expenditure-corprevenu) - taxdebtfixed*revenubasis;
    debt = (1+rstar)*debt;
    *totalassets = *totalassets + debt;
    #endif


    //////////////////////
    // PRINT STATISTICS //
    //////////////////////
    
    tempfileoutfile=fopen(tempfileout, "a");
    setbuf ( tempfileoutfile , NULL );
    
    fprintf(tempfileoutfile,"-------EQUILIBRIUM PRICES AND AGGREGATE--------\n");
    fprintf(tempfileoutfile,"r %20.15f \n", rstar);
    fprintf(tempfileoutfile,"w %20.15f \n", wstar);
    fprintf(tempfileoutfile,"u %20.15f \n", fracUS + fracUL);
    
    
    fprintf(tempfileoutfile,"-------MOMENTS GENERATED & STATISTICS--------\n");
    fprintf(tempfileoutfile,"Ratio KY [3.04] %20.15f\n",KY); // 1st moment (Capital output ratio)
    if (mworthWW>0.0){fprintf(tempfileoutfile,"Ratio mworthVE/mworthWW [8] %20.15f\n",(mworthVE/mworthWW));}else {fprintf(tempfileoutfile,"Ratio mworthVE/mworthWW %20s\n","nan");} // 2nd moment (net VE/WW)
    fprintf(tempfileoutfile,"Fraction E [7.55] %20.15f\n",fracVE); // 3rd moment (Fraction of entrepreneurs)
    fprintf(tempfileoutfile,"(UL: %20.15f, US: %20.15f) Unemploy. rate [4- 6] %20.15f\n",fracUL, fracUS, fracUS + fracUL); // 4th moment (Unemployment rate)
    if ((fracVE+fracWW+fracUS+fracUL)>0.0){fprintf(tempfileoutfile,"Fraction of nonE>E (rlv to whole pop) [0.39] %20.15f\n",((fracUStoVE + fracULtoVE + fracWWtoVE)/(fracVE+fracWW+fracUS+fracUL)));}else{fprintf(tempfileoutfile,"Fraction of nonE>E (rlv to whole pop) %20s\n","nan");} // 5th moment
    if ((fracVE+fracWW+fracUS+fracUL)>0.0){fprintf(tempfileoutfile,"Fraction of U>E (rlv to new entrant) [20] %20.15f\n",(fracUStoVE + fracULtoVE)/(fracUStoVE + fracULtoVE + fracWWtoVE));}else{fprintf(tempfileoutfile,"Fraction of U>E (rlv to new entrant) [20] %20s\n","nan");} // 6th moment
    fprintf(tempfileoutfile,"Worker exit rate (rlv to workers) [3-6.5] %20.15f\n",(fracWWtoVE + fracWWtoUS)/(fracWW));
    fprintf(tempfileoutfile,"Nshare (US) %20.15f, \t Nshare(UL): %20.15f, \t Nshare [0.1]: %20.15f\n",(necessityUStoVE)/(fracUStoVE + fracULtoVE + fracWWtoVE), necessityULtoVE/(fracUStoVE + fracULtoVE + fracWWtoVE), (necessityULtoVE + necessityUStoVE + necessityWWtoVE)/(fracUStoVE + fracULtoVE + fracWWtoVE)); // 11th moment
//    fprintf(tempfileoutfile,"Nshare (US/US>E) %20.15f, \t Nshare(UL/UL>E): %20.15f, \t Nshare (U/U>E) : %20.15f\n",(necessityUStoVE)/(fracUStoVE), necessityULtoVE/(fracULtoVE), (necessityULtoVE + necessityUStoVE)/(fracUStoVE + fracULtoVE)); // 11th moment
    fprintf(tempfileoutfile,"(U>E/pop>E) [0.2] %20.15f\n",(fracUStoVE + fracULtoVE)/(fracUStoVE + fracWWtoVE + fracULtoVE));
    fprintf(tempfileoutfile,"Wealth Gini [0.84] %20.15f\n",GINIval); // 7th moment
        fprintf(tempfileoutfile,"Earning Gini [0.36] %20.15f\n",GINIearn); // 7bisth moment
    fprintf(tempfileoutfile,"Top 1 [41.8] %20.15f\n",top01p); // 8th moment
    fprintf(tempfileoutfile,"Zero wealth fraction [7- 13] %20.15f\n",(distasset[0])); // 9th moment
    if ((fracUS+fracUL)>0.0){fprintf(tempfileoutfile,"Fraction U>nonU (rlv to unemployed) [40- 65] %20.15f\n",(fracUStoVE + fracUStoWW + fracULtoWW + fracULtoVE)/(fracUS+fracUL));}else{fprintf(tempfileoutfile,"Fraction U>nonU out %20s\n","nan");} // 10th moment
    fprintf(tempfileoutfile,"\n");
    
    
    fprintf(tempfileoutfile,"-------TRANSITION Worker to entrepreneur / ABILITY--------\n");
    fprintf(tempfileoutfile,"very low\t middle \t very high\n");
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        fprintf(tempfileoutfile,"%20.15f\t",aggWWtoVE[ygrid]/fracWWpability[ygrid]);
    }
    fprintf(tempfileoutfile,"\n");
    
    
    fprintf(tempfileoutfile,"-------TRANSITION Entrepreneur to Worker / ABILITY--------\n");
        fprintf(tempfileoutfile,"very low \t middle \t very high\n");
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        fprintf(tempfileoutfile,"%20.15f\t",aggVEtoWW[ygrid]/fracVEpability[ygrid]);
    }
        fprintf(tempfileoutfile,"\n");
        
        
    fprintf(tempfileoutfile,"--------- TRANSITION MATRIX BTW CLASSES ------\n");
    fprintf(tempfileoutfile, "  \t  U     \t W     \t E \n");
    fprintf(tempfileoutfile, "U \t %3.4f \t %3.4f \t %3.4f \n", (1- (fracUStoVE + fracUStoWW + fracULtoWW + fracULtoVE)/(fracUS+fracUL)), (fracUStoWW + fracULtoWW)/(fracUS+fracUL), (fracUStoVE + fracULtoVE)/(fracUS+fracUL));
    fprintf(tempfileoutfile, "W \t %3.4f \t %3.4f \t %3.4f \n", (fracWWtoUS)/fracWW, (1 - (fracWWtoVE + fracWWtoUS)/fracWW), (fracWWtoVE)/fracWW);
    fprintf(tempfileoutfile, "E \t %3.4f \t %3.4f \t %3.4f \n", (fracVEtoUL)/fracVE, fracVEtoWW/fracVE, (1-(fracVEtoWW + fracVEtoUL)/fracVE));
    fprintf(tempfileoutfile,"---------------\n");


    fprintf(tempfileoutfile,"--------- WEALTH DISTRIBUTION STATISTICS ------\n");
    fprintf(tempfileoutfile,"Top 0.1 [22] %20.15f\n",top001p);
    fprintf(tempfileoutfile,"Top 1 [41.8] %20.15f\n",top01p); // 8th moment
    fprintf(tempfileoutfile,"Top 5 [?] %20.15f\n",top05p);
    fprintf(tempfileoutfile,"Top 10 [77.2] %20.15f\n",top10p);
    fprintf(tempfileoutfile,"Top 20 [?] %20.15f\n",top20p);
    fprintf(tempfileoutfile,"Bot 20 [?] %20.15f\n",bot20p);
    fprintf(tempfileoutfile,"Bot 40 [?] %20.15f\n",bot40p);
    fprintf(tempfileoutfile,"Bot 60 [?] %20.15f\n",bot60p);
    fprintf(tempfileoutfile,"Bot 80 [?] %20.15f\n",bot80p);
    fprintf(tempfileoutfile,"mworthVE %20.15f\n",mworthVE);
    fprintf(tempfileoutfile,"mworthWW %20.15f\n",mworthWW);
    fprintf(tempfileoutfile,"mworthUS %20.15f\n",mworthUS);
    fprintf(tempfileoutfile,"mworthUL %20.15f\n",mworthUL);
    fprintf(tempfileoutfile,"---------------\n");

     
    fprintf(tempfileoutfile,"--------- PRODUCTION & CAPITAL -----\n");
    fprintf(tempfileoutfile,"CapitalVE %20.15f\n",capitalVE);
    fprintf(tempfileoutfile,"ProductionVE %20.15f\n",productionVE);
    fprintf(tempfileoutfile,"Totalcapital %20.15f\n",totalcapital);
    fprintf(tempfileoutfile,"Totallabor %20.15f\n",*totallabor);
    fprintf(tempfileoutfile,"Corproduction %20.15f\n",corproduction);
    fprintf(tempfileoutfile,"Corcapital %20.15f\n",*totalassets);
    fprintf(tempfileoutfile,"Totalproduction %20.15f\n",totalproduction);
    fprintf(tempfileoutfile,"CapitalVE/Totalcapital %20.15f\n",(capitalVE/totalcapital));
    fprintf(tempfileoutfile,"---------------\n");

     
    fprintf(tempfileoutfile,"--------- TAXATION -----\n");
    fprintf(tempfileoutfile,"Corporate revenu %20.15f\n",corprevenu);
    fprintf(tempfileoutfile,"Expenditure %20.15f\n",expenditure);
    fprintf(tempfileoutfile,"Revenue Basis %20.15f\n",revenubasis);
    fprintf(tempfileoutfile,"Optimal tax %20.15f\n",*optimaltax);
    
    
    fprintf(tempfileoutfile,"-------MASS--------\n");
    fprintf(tempfileoutfile,"WW %20.15f \n", fracWW);
    fprintf(tempfileoutfile,"U %20.15f \n", fracUL + fracUS);
    fprintf(tempfileoutfile,"VE %20.15f \n", fracVE);
    
    fprintf(tempfileoutfile,"-------WELFARE--------\n");
    fprintf(tempfileoutfile,"UL %20.15f  US %20.15f  VE %20.15f  WW %20.15f  Wtot%20.15f  Avg Firm: %20.15f \n", welfareUL,welfareUS,welfareVE, welfareWW, welfaretot, avgfirmsize);
    
    fprintf(tempfileoutfile,"-------POOL OF ENT--------\n");
    fprintf(tempfileoutfile,"THETA1:  %20.15f, %20.15f, %20.15f  \n", ENT1frac, AVGFIRM1, AVGASSET1);
    fprintf(tempfileoutfile,"THETA2:  %20.15f, %20.15f, %20.15f  \n", ENT2frac, AVGFIRM2, AVGASSET2);
    fprintf(tempfileoutfile,"THETA3:  %20.15f, %20.15f, %20.15f  \n", ENT3frac, AVGFIRM3, AVGASSET3);


    fclose(tempfileoutfile);
    
    
    
    ////////////////////////
    // PRINT MOMENTS USED //
    ////////////////////////
    
    distoutfile=fopen(moment, "w");
    setbuf ( distoutfile , NULL );
    
    fprintf(distoutfile,"\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","KY", "E", "U", "TOP1", "GINI", "netVEWW", "ZERO", "NECESSITY","EEXIT", "UEXIT", "UtoE","WWEXIT");
    
    fprintf(distoutfile,"%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",KY, fracVE, fracUS + fracUL, top01p, GINIval, (mworthVE/mworthWW), distasset[0], (necessityULtoVE + necessityUStoVE)/(fracUStoVE + fracULtoVE + fracWWtoVE), (fracVEtoUL + fracVEtoWW)/fracVE, (fracUStoWW + fracUStoVE + fracULtoWW + fracULtoVE)/(fracUL + fracUS), (fracUStoVE + fracULtoVE)/(fracUS + fracUL), (fracWWtoVE + fracWWtoUS)/(fracWW));

    fclose(distoutfile);
    
    
    
    ////////////////////////////
    // RETURN LIST OF MOMENTS //
    ////////////////////////////
    generatedmoment[0] = KY/10;
    generatedmoment[1] = fracVE*10;
    generatedmoment[2] = (fracUS + fracUL)*10;
    generatedmoment[3] = top01p;
    generatedmoment[4] = GINIval;
    generatedmoment[5] = (fracUStoVE + fracULtoVE)/(fracUStoVE + fracULtoVE + fracWWtoVE);
    generatedmoment[6] = (fracUStoWW + fracUStoVE + fracULtoWW + fracULtoVE)/((fracUL + fracUS));
    generatedmoment[7] = (fracVEtoWW)/(fracVE);
    generatedmoment[8] = (fracWWtoVE + fracWWtoUS)/(fracWW);
    generatedmoment[9] = aggVEtoWW[0]/fracVEpability[0];
    generatedmoment[10] = aggVEtoWW[int(floor(maxprod/2))]/fracVEpability[int(floor(maxprod/2))];
    generatedmoment[11] = aggVEtoWW[(maxprod - 1)]/fracVEpability[(maxprod - 1)];
    generatedmoment[12] = (aggWWtoVE[0]/fracWWpability[0])*10;
    generatedmoment[13] = (aggWWtoVE[int(floor(maxprod/2))]/fracWWpability[int(floor(maxprod/2))])*10;
    generatedmoment[14] = (aggWWtoVE[(maxprod - 1)]/fracWWpability[(maxprod - 1)])*10;
    generatedmoment[15] = GINIearn;
    
#if indexPM == 0
    printf("SIMULATION: Kcorp: %20.10f\t, KVE: %20.10f\t, Corprod: %20.10f\t, totalprod: %20.10f\t, KY: %20.10f\n",*totalassets, capitalVE, corproduction, totalproduction, KY);
    printf("SIMULATION: fracVE: %20.10f\t, fracU: %20.10f\t, fracWW: %20.10f\n",fracVE, fracUS+fracUL, fracWW);
    printf("------------------------\n");
    printf("-------TRANSITION Worker to entrepreneur / ABILITY--------\n");
    printf("\t very low \t middle \t very high\n");
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        printf("%20.15f\t",aggWWtoVE[ygrid]/fracWWpability[ygrid]);
    }
    printf("\n");
    
    
    printf("-------TRANSITION Entrepreneur to Worker / ABILITY--------\n");
    printf("\t very low \t middle \t very high\n");
    for(ygrid = 0.0; ygrid < maxprod; ygrid++)
    {
        printf("%20.15f\t",aggVEtoWW[ygrid]/fracVEpability[ygrid]);
    }
    printf("\n");
#endif
    
#if indexPMwrite == 1
    // DISTRIBUTION for projected methods:
    distoutfile=fopen(distVE_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(egrid=0;egrid<maxfirmtype;egrid++) {
                fprintf(distoutfile,"%20.15f\n",distVE[inxE(igrid, ygrid, egrid)]);
            }
        }
    }
    
    fclose(distoutfile);
    
    
    distoutfile=fopen(distWW1_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distWW1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);
    
    distoutfile=fopen(distWW0_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distWW0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);

    distoutfile=fopen(distUS1_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distUS1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);

    distoutfile=fopen(distUL1_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distUL1[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);

    distoutfile=fopen(distUS0_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distUS0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);

    distoutfile=fopen(distUL0_Tfile, "w");
    setbuf ( distoutfile , NULL );
    
    for(int ygrid = 0; ygrid < maxprod ; ygrid++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            
            fprintf(distoutfile,"%20.15f\n",distUL0[inx(igrid, ygrid)]);
            
        }
    }
    
    fclose(distoutfile);
#endif

    
    
    // Clear memory
    free(distasset);
    free(disttotal);
    free(necessityUS);
    free(necessityUL);
    free(aggVEtoWW);
    free(aggWWtoVE);
    free(distWWtoVE);
    free(distVEtoWW);
    
#if indexPM == 0
    free(distUS0);
    free(distUL0);
    free(distWW0);
    free(distUS1);
    free(distUL1);
    free(distWW1);
    free(distVE);
    free(distVEwUI);
#endif
    
    free(distassetUS);
    free(distassetUL);
    free(distassetWW);
    free(distassetVE);

    free(distUS0old);
    free(distUL0old);
    free(distWW0old);
    free(distUS1old);
    free(distUL1old);
    free(distWW1old);
    free(distVEold);
    free(distVEwUIold);
        
    free(saveresUS0);
    free(saveresUL0);
    free(saveresWW0);
    free(saveresUS1);
    free(saveresUL1);
    free(saveresWW1);
    free(saveresVE);
    free(saveresVEwUI);
    
    free(isaveWW1);
    free(isaveVE);
    free(isaveVEwUI);
    free(isaveUS1);
    free(isaveUL1);
    free(isaveWW0);
    free(isaveUS0);
    free(isaveUL0);
    free(isaveWW1c);
    free(isaveUS1c);
    free(isaveUL1c);
    
}






///////// PROJECTED METHOD ////////////
#define inxTP(indext, igridindex, ygridindex) ((maxgrid*maxprod*indext)+(ygridindex)*(maxgrid)+(igridindex))
#define inxTPE(indext, igridindex, ygridindex, typeindex) (((indext*maxgrid*maxfirmtype*maxprod)+((ygridindex)*(maxgrid)*(maxfirmtype))+(igridindex)*(maxfirmtype)+typeindex))
#define inxPrice(indext, indexprice) (7*indext+indexprice)

void bascule2E(double *VectorIN, double *VectorOUT, int TT)
{
    int igrid, ygrid, egrid;
    for(igrid=0; igrid<maxgrid; igrid++){
        for(ygrid=0; ygrid<maxprod; ygrid++){
            for(egrid=0; egrid<maxfirmtype; egrid++){
                VectorOUT[inxTPE(TT, igrid, ygrid, egrid)]=VectorIN[inxE(igrid, ygrid, egrid)];
            }
        }
    }
}
    
void bascule2(double *VectorIN, double *VectorOUT, int TT)
{
    int igrid, ygrid;
    for(igrid=0; igrid<maxgrid; igrid++){
        for(ygrid=0; ygrid<maxprod; ygrid++){
            VectorOUT[inxTP(TT, igrid, ygrid)]=VectorIN[inx(igrid, ygrid)];
        }
    }
}

void bascule3E(double *VectorIN, double *VectorOUT, int TT)
{
    int igrid, ygrid, egrid;
    for(igrid=0; igrid<maxgrid; igrid++){
        for(ygrid=0; ygrid<maxprod; ygrid++){
            for(egrid=0; egrid<maxfirmtype; egrid++){
                VectorOUT[inxE(igrid, ygrid, egrid)]=VectorIN[inxTPE(TT, igrid, ygrid, egrid)];
            }
        }
    }
}

void bascule3(double *VectorIN, double *VectorOUT, int TT)
{
    int igrid, ygrid;
    for(igrid=0; igrid<maxgrid; igrid++){
        for(ygrid=0; ygrid<maxprod; ygrid++){
            VectorOUT[inx(igrid, ygrid)]=VectorIN[inxTP(TT, igrid, ygrid)];
        }
    }
}


#if indexPM == 1
void Projected_TP(int T) {

// The goal is the following
/* We know the final distribution and the starting distribution, so we know interest rate, wage and so on...
Given an initial guess of price path, we compute the optimal decision given V(t+1) and prices starting
from T-1 to 1.
At period 1, all is given, we know price, we know the current distribution of wealth and so on. We use it to simulate
next period wealth level. That is, knowning that there is a distribution F(at) of asset at time 1, we take it a given
in period 2, we compute how individual transit according to this initial distribution, policy computed and price guessed.
We don't need to iterate over the simulation in this configuration (and we should not since we don't want ergodicity).
Given the new distribution of asset at period 2, we compute the new prices. If prices are different = update the guess
and restart interation of policy function and so on until convergence of TP.
--> At the end, compare the final distribution to the ergodic distribution that we should find, if to much different, 
increase the number of transition periods T.
*/

printf("STARTING...\n");

int i, j;
FILE *pricefile;
double epsilonTP = 0.000001;
double relaxTP = 0.8;
double *pricepath, *nowpricepath, GOODNESS, critTP, *genmoments, taxout, capitalout, laborout, fracVE, fracU;
genmoments = (double *) calloc((nbmoments), sizeof(double));
pricepath = (double *) calloc((T*7), sizeof(double));
nowpricepath = (double *) calloc((T*7), sizeof(double));

double *endprice, *startingprice;
endprice = (double *) calloc((7), sizeof(double));
startingprice = (double *) calloc((7), sizeof(double));

startingprice[0] = 0.057688258463640;        endprice[0] = 0.058865912247287;
startingprice[1] = 1.200363390425356;        endprice[1] = 1.193659293834749;
startingprice[2] = 0.010446298367941;        endprice[2] = 0.010446195624903;
startingprice[3] = 5.036197604700030;        endprice[3] = 4.960771969604878;
startingprice[4] = 0.877810214404596;        endprice[4] = 0.878221293800238;
startingprice[5] = 0.074337548852433;        endprice[5] = 0.071585727995411;
startingprice[6] = 0.052458167206997;        endprice[6] = 0.054604554469187;

// GUESS THE PRICE PATH //
for(i = 0; i < 7; i++) {
    pricepath[inxPrice(0,i)] = startingprice[i];
    pricepath[inxPrice((T-1),i)] = endprice[i];
}

for(j = 1; j < (T-1); j++) {
    for(i = 0; i < 7; i++) {
        pricepath[inxPrice(j,i)] = startingprice[i] + j*(endprice[i] - startingprice[i])/(T-1);
    }
}

for(j = 0; j < T; j++) {
    for(i = 0; i < 7; i++) {
        printf("%f\t", pricepath[inxPrice(j,i)]);
    }
    printf("\n");
}


/* IMPORTATION PART */
// LAST PERIOD VALUE FUNCTION AND DISTRIBUTION //
double *valueVET, *valueWW1T, *valueWW0T, *valueUS1T, *valueUL1T, *valueUS0T, *valueUL0T, *valueVEcT;
valueVET = (double *) calloc((ifulldimE), sizeof(double));
valueVEcT = (double *) calloc((ifulldimE), sizeof(double));
valueWW1T = (double *) calloc((ifulldim), sizeof(double));
valueWW0T = (double *) calloc((ifulldim), sizeof(double));
valueUL1T = (double *) calloc((ifulldim), sizeof(double));
valueUL0T = (double *) calloc((ifulldim), sizeof(double));
valueUS1T = (double *) calloc((ifulldim), sizeof(double));
valueUS0T = (double *) calloc((ifulldim), sizeof(double));

readinput(valueVET, ifulldimE, "project_method/valueVET.out");
readinput(valueVEcT, ifulldimE, "project_method/valueVET.out");
readinput(valueWW1T, ifulldim, "project_method/valueWW1T.out");
readinput(valueWW0T, ifulldim, "project_method/valueWW0T.out");
readinput(valueUL1T, ifulldim, "project_method/valueUL1T.out");
readinput(valueUL0T, ifulldim, "project_method/valueUL0T.out");
readinput(valueUS1T, ifulldim, "project_method/valueUS1T.out");
readinput(valueUS0T, ifulldim, "project_method/valueUS0T.out");


double *distVE_T, *distWW1_T, *distWW0_T, *distUL0_T, *distUL1_T, *distUS1_T, *distUS0_T;
distVE_T = (double *) calloc((ifulldimE), sizeof(double));
distWW1_T = (double *) calloc((ifulldim), sizeof(double));
distWW0_T = (double *) calloc((ifulldim), sizeof(double));
distUS1_T = (double *) calloc((ifulldim), sizeof(double));
distUS0_T = (double *) calloc((ifulldim), sizeof(double));
distUL1_T = (double *) calloc((ifulldim), sizeof(double));
distUL0_T = (double *) calloc((ifulldim), sizeof(double));

readinput(distVE_T, ifulldimE, "project_method/distVE_T.out");
readinput(distWW1_T, ifulldim, "project_method/distWW1_T.out");
readinput(distWW0_T, ifulldim, "project_method/distWW0_T.out");
readinput(distUS0_T, ifulldim, "project_method/distUS0_T.out");
readinput(distUS1_T, ifulldim, "project_method/distUS1_T.out");
readinput(distUL0_T, ifulldim, "project_method/distUL0_T.out");
readinput(distUL1_T, ifulldim, "project_method/distUL1_T.out");

// FIRST PERIOD WEALTH DISTRIBUTION //
double *distVE_0, *distWW1_0, *distWW0_0, *distUL0_0, *distUL1_0, *distUS1_0, *distUS0_0;
distVE_0 = (double *) calloc((ifulldimE), sizeof(double));
distWW1_0 = (double *) calloc((ifulldim), sizeof(double));
distWW0_0 = (double *) calloc((ifulldim), sizeof(double));
distUS1_0 = (double *) calloc((ifulldim), sizeof(double));
distUS0_0 = (double *) calloc((ifulldim), sizeof(double));
distUL1_0 = (double *) calloc((ifulldim), sizeof(double));
distUL0_0 = (double *) calloc((ifulldim), sizeof(double));

readinput(distVE_0, ifulldimE, "project_method/distVE_0.out");
readinput(distWW1_0, ifulldim, "project_method/distWW1_0.out");
readinput(distWW0_0, ifulldim, "project_method/distWW0_0.out");
readinput(distUS0_0, ifulldim, "project_method/distUS0_0.out");
readinput(distUS1_0, ifulldim, "project_method/distUS1_0.out");
readinput(distUL0_0, ifulldim, "project_method/distUL0_0.out");
readinput(distUL1_0, ifulldim, "project_method/distUL1_0.out");

// IMPORT ENDOK //
double *endoK;
endoK = (double *) calloc((ifulldimKK), sizeof(double));
readinput(endoK, ifulldimKK, "project_method/endoK.out");
//////////////////////


/* DEFINE TEMP POLICY FUNCTION AND VALUE FUNCTIONS */
double *saveVE, *saveWW1, *saveWW0, *saveUS1, *saveUL1, *saveUS0, *saveUL0, *searchWW1, *searchVE, *searchEffortVEWW, *searchEffortWWVE, *searchEffortWWUS1, *searchEffortVEUS1, *searchEffortWWUL1, *searchEffortVEUL1, *searchEffortWWUS0 , *searchEffortWWUL0;

saveVE = (double *) calloc((ifulldimE*T), sizeof(double));
saveWW1 = (double *) calloc((ifulldim*T), sizeof(double));
saveWW0 = (double *) calloc((ifulldim*T), sizeof(double));
saveUS1 = (double *) calloc((ifulldim*T), sizeof(double));
saveUL1 = (double *) calloc((ifulldim*T), sizeof(double));
saveUS0 = (double *) calloc((ifulldim*T), sizeof(double));
saveUL0 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortVEWW = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortWWVE = (double *) calloc((ifulldimE*T), sizeof(double));
searchEffortWWUS1 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortVEUS1 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortWWUL1 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortVEUL1 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortWWUS0 = (double *) calloc((ifulldim*T), sizeof(double));
searchEffortWWUL0 = (double *) calloc((ifulldim*T), sizeof(double));

double *temp_saveVE, *temp_saveWW1, *temp_saveWW0, *temp_saveUS1, *temp_saveUL1, *temp_saveUS0, *temp_saveUL0, *temp_searchEffortVEWW, *temp_searchEffortWWVE, *temp_searchEffortWWUS1, *temp_searchEffortVEUS1, *temp_searchEffortWWUL1 , *temp_searchEffortVEUL1 , *temp_searchEffortWWUS0, *temp_searchEffortWWUL0;

temp_saveVE = (double *) calloc((ifulldimE), sizeof(double));
temp_saveWW1 = (double *) calloc((ifulldim), sizeof(double));
temp_saveWW0 = (double *) calloc((ifulldim), sizeof(double));
temp_saveUS1 = (double *) calloc((ifulldim), sizeof(double));
temp_saveUL1 = (double *) calloc((ifulldim), sizeof(double));
temp_saveUS0 = (double *) calloc((ifulldim), sizeof(double));
temp_saveUL0 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortVEWW = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortWWVE = (double *) calloc((ifulldimE), sizeof(double));
temp_searchEffortWWUS1 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortVEUS1 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortWWUL1 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortVEUL1 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortWWUS0 = (double *) calloc((ifulldim), sizeof(double));
temp_searchEffortWWUL0 = (double *) calloc((ifulldim), sizeof(double));

double *valueVE, *valueWW1, *valueWW0, *valueUS1, *valueUL1, *valueUL0,  *valueUS0, *valueVEc;

valueVE = (double *) calloc((ifulldimE*T), sizeof(double));
valueVEc = (double *) calloc((ifulldimE*T), sizeof(double));
valueWW1 = (double *) calloc((ifulldim*T), sizeof(double));
valueWW0 = (double *) calloc((ifulldim*T), sizeof(double));
valueUL1 = (double *) calloc((ifulldim*T), sizeof(double));
valueUL0 = (double *) calloc((ifulldim*T), sizeof(double));
valueUS1 = (double *) calloc((ifulldim*T), sizeof(double));
valueUS0 = (double *) calloc((ifulldim*T), sizeof(double));

double *temp_valueVE, *temp_valueWW1, *temp_valueWW0, *temp_valueUS1, *temp_valueUL1, *temp_valueUS0, *temp_valueUL0, *temp_valueVEc;

temp_valueVE = (double *) calloc((ifulldimE), sizeof(double));
temp_valueVEc = (double *) calloc((ifulldimE), sizeof(double));
temp_valueWW1 = (double *) calloc((ifulldim), sizeof(double));
temp_valueWW0 = (double *) calloc((ifulldim), sizeof(double));
temp_valueUL1 = (double *) calloc((ifulldim), sizeof(double));
temp_valueUL0 = (double *) calloc((ifulldim), sizeof(double));
temp_valueUS1 = (double *) calloc((ifulldim), sizeof(double));
temp_valueUS0 = (double *) calloc((ifulldim), sizeof(double));

double *distVE, *distWW1, *distWW0, *distUL1, *distUL0, *distUS0, *distUS1;

distVE = (double *) calloc((ifulldimE*T), sizeof(double));
distWW1 = (double *) calloc((ifulldim*T), sizeof(double));
distWW0 = (double *) calloc((ifulldim*T), sizeof(double));
distUL1 = (double *) calloc((ifulldim*T), sizeof(double));
distUL0 = (double *) calloc((ifulldim*T), sizeof(double));
distUS1 = (double *) calloc((ifulldim*T), sizeof(double));
distUS0 = (double *) calloc((ifulldim*T), sizeof(double));

double *temp_distVE, *temp_distWW1, *temp_distWW0, *temp_distUL1, *temp_distUL0, *temp_distUS0, *temp_distUS1;

temp_distVE = (double *) calloc((ifulldimE), sizeof(double));
temp_distWW1 = (double *) calloc((ifulldim), sizeof(double));
temp_distWW0 = (double *) calloc((ifulldim), sizeof(double));
temp_distUL1 = (double *) calloc((ifulldim), sizeof(double));
temp_distUL0 = (double *) calloc((ifulldim), sizeof(double));
temp_distUS1 = (double *) calloc((ifulldim), sizeof(double));
temp_distUS0 = (double *) calloc((ifulldim), sizeof(double));

double *temp_distVEnext, *temp_distWW1next, *temp_distWW0next, *temp_distUL1next, *temp_distUL0next, *temp_distUS0next, *temp_distUS1next;

temp_distVEnext = (double *) calloc((ifulldimE), sizeof(double));
temp_distWW1next = (double *) calloc((ifulldim), sizeof(double));
temp_distWW0next = (double *) calloc((ifulldim), sizeof(double));
temp_distUL1next = (double *) calloc((ifulldim), sizeof(double));
temp_distUL0next = (double *) calloc((ifulldim), sizeof(double));
temp_distUS1next = (double *) calloc((ifulldim), sizeof(double));
temp_distUS0next = (double *) calloc((ifulldim), sizeof(double));
///////////////////////////////////


// SAVE THE LAST PERIOD VALUE FUNCTION //
bascule2E(valueVET,valueVE,(T-1));
bascule2E(valueVEcT, valueVEc,(T-1));
bascule2(valueWW1T, valueWW1,(T-1));
bascule2(valueWW0T, valueWW0,(T-1));
bascule2(valueUS1T, valueUS1,(T-1));
bascule2(valueUL1T, valueUL1,(T-1));
bascule2(valueUS0T, valueUS0,(T-1));
bascule2(valueUL0T, valueUL0,(T-1));

// SET FIRST PERIOD DISTRIBUTION //
bascule2E(distVE_0, distVE, 0);
bascule2E(distWW0_0, distWW0, 0);
bascule2(distWW1_0, distWW1, 0);
bascule2(distUS1_0, distUS1, 0);
bascule2(distUS0_0, distUS0, 0);
bascule2(distUL1_0, distUL1, 0);
bascule2(distUL0_0, distUL0, 0);
    
// STARTING LOOP OVER PRICES //

critTP = 1;

while(critTP > epsilonTP) {

// Plug value in T into the vector corresponding in temp_value //
bascule(valueVET, temp_valueVE, ifulldimE);
bascule(valueVEcT, temp_valueVEc, ifulldimE);
bascule(valueWW1T, temp_valueWW1, ifulldim);
bascule(valueWW0T, temp_valueWW0, ifulldim);
bascule(valueUS1T, temp_valueUS1, ifulldim);
bascule(valueUL1T, temp_valueUL1, ifulldim);
bascule(valueUS0T, temp_valueUS0, ifulldim);
bascule(valueUL0T, temp_valueUL0, ifulldim);

// SET FIRST PERIOD DISTRIBUTION //
bascule(distVE_0, temp_distVE, ifulldimE);
bascule(distWW0_0, temp_distWW0, ifulldim);
bascule(distWW1_0, temp_distWW1, ifulldim);
bascule(distUS1_0, temp_distUS1, ifulldim);
bascule(distUS0_0, temp_distUS0, ifulldim);
bascule(distUL1_0, temp_distUL1, ifulldim);
bascule(distUL0_0, temp_distUL0, ifulldim);

printf("POLICY VALUE STARTING\n");
    // COMPUTE OPTIMAL POLICY FUNCTION -- STORE IT
    // i = 0 correspond to i = T-1
    for(int ii = 1; ii < T;  ii++) {
    
        // Set prices right now here //
        rstar = pricepath[inxPrice((T-1-ii), 0)];
        wstar = pricepath[inxPrice((T-1-ii), 1)];
        taxrate = pricepath[inxPrice((T-1-ii), 2)];
        
        printf("PERIOD: %d, R: %f, W: %f, T: %f \n", (T-1-ii), rstar, wstar, taxrate);
        
        VFI(temp_valueWW0, temp_valueWW1, temp_valueUS0, temp_valueUS1, temp_valueUL0, temp_valueUL1, temp_valueVE, temp_valueVEc, endoK, temp_saveWW0, temp_saveWW1, temp_saveUS0, temp_saveUS1, temp_saveUL0, temp_saveUL1, temp_saveVE, temp_searchEffortWWUS0, temp_searchEffortWWUL0, temp_searchEffortWWUS1, temp_searchEffortVEUS1, temp_searchEffortWWUL1, temp_searchEffortVEUL1, temp_searchEffortWWVE, temp_searchEffortVEWW);
        
        /// HERE YOU WILL NEED TOW TYPES OF VFI DEPENDING ON THE TIMING AT WHICH THE MEASURE IS  IMPLEMENTED //
        /* I THINK U SHOULD USE VFIwUI and VFI */
        
        bascule2E(temp_valueVE,valueVE,(T-1-ii));
        bascule2E(temp_valueVEc,valueVEc,(T-1-ii));
        bascule2(temp_valueWW1,valueWW1,(T-1-ii));
        bascule2(temp_valueWW0,valueWW0,(T-1-ii));
        bascule2(temp_valueUS1,valueUS1,(T-1-ii));
        bascule2(temp_valueUL1,valueUL1,(T-1-ii));
        bascule2(temp_valueUS0,valueUS0,(T-1-ii));
        bascule2(temp_valueUL0,valueUL0,(T-1-ii));
        
        bascule2E(temp_searchEffortWWVE,searchEffortWWVE,(T-1-ii));
        bascule2(temp_searchEffortVEWW,searchEffortVEWW,(T-1-ii));
        bascule2(temp_searchEffortVEUS1,searchEffortVEUS1,(T-1-ii));
        bascule2(temp_searchEffortVEUL1,searchEffortVEUL1,(T-1-ii));
        bascule2(temp_searchEffortWWUS1,searchEffortWWUS1,(T-1-ii));
        bascule2(temp_searchEffortWWUL1,searchEffortWWUL1,(T-1-ii));
        bascule2(temp_searchEffortWWUS0,searchEffortWWUS0,(T-1-ii));
        bascule2(temp_searchEffortWWUL0,searchEffortWWUL0,(T-1-ii));
        
        bascule2E(temp_saveVE,saveVE,(T-1-ii));
        bascule2(temp_saveWW1,saveWW1,(T-1-ii));
        bascule2(temp_saveWW0,saveWW0,(T-1-ii));
        bascule2(temp_saveUS1,saveUS1,(T-1-ii));
        bascule2(temp_saveUL1,saveUL1,(T-1-ii));
        bascule2(temp_saveUS0,saveUS0,(T-1-ii));
        bascule2(temp_saveUL0,saveUL0,(T-1-ii));
        
//        printf("%f %f\n", temp_saveWW0[4], temp_saveWW1[4]);
//        printf("%f \n", temp_saveVE[inxE(205, 0, 0)]);
//        printf("%f %f\n", saveWW0[inxTP((T-1-ii), 4, 0)], saveWW1[inxTP((T-1-ii), 4, 0)]);
//        printf("%f \n", saveVE[inxTPE((T-1-ii), 205, 0, 0)]);
//        getchar();
//        printf("%d %d %f %f", T-1-ii, inxTP(T-1-ii, 4, 0), saveWW0[inxTP(T-1-ii, 4, 0)], saveWW1[inxTP(T-1-ii, 4, 0)]); getchar();
        
    }
    
printf("\nSTARTING DISTRIBUTION...\n");
    // START SIMULATING STARTING FROM PERIOD 0, where prices are given, but policy have changed.//
    // j = 0 --> moment where policy have changed, but starting distribution is given //
    for(int jj = 0; jj < (T); jj++) {
    
        if(jj == (T-1)) {
        bascule3E(valueVE,temp_valueVE,(jj+1));
        bascule3E(valueVEc,temp_valueVEc,(jj+1));
        bascule3(valueWW1,temp_valueWW1,(jj+1));
        bascule3(valueWW0,temp_valueWW0,(jj+1));
        bascule3(valueUS1,temp_valueUS1,(jj+1));
        bascule3(valueUL1,temp_valueUL1,(jj+1));
        bascule3(valueUS0,temp_valueUS0,(jj+1));
        bascule3(valueUL0,temp_valueUL0,(jj+1));
        } else {
        bascule3E(valueVE,temp_valueVE,(jj));
        bascule3E(valueVEc,temp_valueVEc,(jj));
        bascule3(valueWW1,temp_valueWW1,(jj));
        bascule3(valueWW0,temp_valueWW0,(jj));
        bascule3(valueUS1,temp_valueUS1,(jj));
        bascule3(valueUL1,temp_valueUL1,(jj));
        bascule3(valueUS0,temp_valueUS0,(jj));
        bascule3(valueUL0,temp_valueUL0,(jj));
        }
        
        bascule3E(searchEffortWWVE, temp_searchEffortWWVE, jj);
        bascule3(searchEffortVEWW, temp_searchEffortVEWW, jj);
        bascule3(searchEffortVEUS1, temp_searchEffortVEUS1, jj);
        bascule3(searchEffortVEUL1, temp_searchEffortVEUL1, jj);
        bascule3(searchEffortWWUS1, temp_searchEffortWWUS1, jj);
        bascule3(searchEffortWWUL1, temp_searchEffortWWUL1, jj);
        bascule3(searchEffortWWUS0, temp_searchEffortWWUS0, jj);
        bascule3(searchEffortWWUL0, temp_searchEffortWWUL0, jj);
        
        bascule3E(saveVE, temp_saveVE, jj);
        bascule3(saveWW1, temp_saveWW1, jj);
        bascule3(saveWW0, temp_saveWW0, jj);
        bascule3(saveUS1, temp_saveUS1, jj);
        bascule3(saveUL1, temp_saveUL1, jj);
        bascule3(saveUS0, temp_saveUS0, jj);
        bascule3(saveUL0, temp_saveUL0, jj);
        
//        printf("%f %f\n", temp_saveWW0[4], temp_saveWW1[4]); getchar();
//        printf("%f %f\n", temp_saveVE[405], temp_saveVE[405]); getchar();
        
//        printf("I am here");getchar();
        
        SIMULATION(temp_valueWW0, temp_valueWW1, temp_valueUS0, temp_valueUS1, temp_valueUL0, temp_valueUL1, temp_valueVE, temp_valueVEc, endoK, temp_saveWW0, temp_saveWW1, temp_saveUS0, temp_saveUS1, temp_saveUL0, temp_saveUL1, temp_saveVE, temp_searchEffortWWUS0, temp_searchEffortWWUL0, temp_searchEffortWWUS1, temp_searchEffortVEUS1, temp_searchEffortWWUL1, temp_searchEffortVEUL1, temp_searchEffortWWVE, temp_searchEffortVEWW, &capitalout, &laborout, &taxout, genmoments, temp_distVE, temp_distWW1, temp_distWW0, temp_distUL0, temp_distUL1, temp_distUS1, temp_distUS0, temp_distVEnext, temp_distWW1next, temp_distWW0next, temp_distUL0next, temp_distUL1next, temp_distUS1next, temp_distUS0next, &fracVE, &fracU);
        
        nowpricepath[inxPrice((jj+1), 0)] = alphapar*pow(capitalout,(alphapar-1.0))*pow(laborout,(1.0-alphapar))-deltapar;
        nowpricepath[inxPrice((jj+1), 1)] = (1.0-alphapar)*pow(((nowpricepath[inxPrice((jj+1), 0)]+deltapar)/alphapar),(alphapar/(alphapar-1.0)));
        nowpricepath[inxPrice((jj+1), 2)] = taxout;
        nowpricepath[inxPrice((jj+1), 3)] = capitalout;
        nowpricepath[inxPrice((jj+1), 4)] = laborout;
        nowpricepath[inxPrice((jj+1), 5)] = fracVE;
        nowpricepath[inxPrice((jj+1), 6)] = fracU;
        
        bascule(temp_distVEnext, temp_distVE, ifulldimE);
        bascule(temp_distWW0next, temp_distWW0, ifulldim);
        bascule(temp_distWW1next, temp_distWW1, ifulldim);
        bascule(temp_distUL0next, temp_distUL0, ifulldim);
        bascule(temp_distUL1next, temp_distUL1, ifulldim);
        bascule(temp_distUS1next, temp_distUS1, ifulldim);
        bascule(temp_distUS0next, temp_distUS0, ifulldim);
        
        printf("PERIOD: %d \n", jj);
    }

//    for(j = 1; j < T; j++) {
//        for(i = 0; i < 5; i++) {
//            printf("%f\t", nowpricepath[inxPrice(j,i)]);
//        }
//        printf("\n");
//    }
//    
//    getchar();
    
    // HERE UPDATE PRICES WITH A RELAXATION //
    critTP = 0.0;

    for(int ii = 1; ii < (T-1); ii++) {
        for(int jj = 0; jj < 7; jj++) {
            critTP = max(critTP, abs(pricepath[inxPrice(ii,jj)] - nowpricepath[inxPrice(ii,jj)])/pricepath[inxPrice(ii,jj)]);
//            printf("%f %f %f", relaxTP, pricepath[inxPrice(ii,jj)], nowpricepath[inxPrice(ii,jj)]);getchar();
            pricepath[inxPrice(ii,jj)] = relaxTP*pricepath[inxPrice(ii,jj)]+(1.0 - relaxTP)*nowpricepath[inxPrice(ii,jj)];
//            printf("%f %f %f", relaxTP, pricepath[inxPrice(ii,jj)], nowpricepath[inxPrice(ii,jj)]);getchar();
        }
    }
    
    // HERE COMPUTE THE GOODNESS OF FIT FOR F //
    GOODNESS = 0.0;
    for(j = 0; j < 7; j++) {
        GOODNESS = max(GOODNESS, abs(pricepath[inxPrice((T-1), j)] - nowpricepath[inxPrice((T-1), j)]));
    }
    printf("GOODNESS: %f \n", GOODNESS);



    // SAVING DECISION
    pricefile=fopen(pricepathfile, "w");
    setbuf ( pricefile , NULL );
    
    
    fprintf(pricefile,"%20s\t%20s\t%20s\t%20s\t%20s\n","R", "W", "T", "K", "L");

    for(i = 0; i < (T); i++) {
        for(j = 0; j < 7; j++) {
            fprintf(pricefile,"%20.15f\t",pricepath[inxPrice(i, j)]);
        }
        fprintf(pricefile,"\n");
    }

    fclose(pricefile);
    
    printf("Crit cvng: %20.15f", critTP);
    

} // end while


}
#endif





double ComputeEquilibrium(double *PARAM)
{
    // SET PARAMETERS
    betapar = PARAM[0];
    phiparW = PARAM[1];
    phiparE = PARAM[2];
    alphakappaw1 = PARAM[3];
    kappaE = PARAM[4];
    g1 = PARAM[5];
    g2 = PARAM[6];
    g3 = PARAM[7];
    fpar = PARAM[8];
    nupar = PARAM[9];
    stateg[0] = PARAM[10];
    pibad = PARAM[11];
    kappaW = PARAM[12];
    mupar = PARAM[13];


    // REGRESSION //
    gtrans[0] = g1;
    gtrans[1] = g2;
    gtrans[2] = g3;
    
    pig1 = pibad;
    pig2 = pibad;
    pig3 = pibad;

    pig11 = pibad;
    pig22 = pibad;
    pig33 = pibad;
    
    mgtrans[0][0] = pig1; mgtrans[0][1] = 1-pig1; mgtrans[1][0] = pig2; mgtrans[1][1] = 1-pig2; mgtrans[2][0] = pig3; mgtrans[2][1] = 1-pig3;
    // INITIALIZE REGRESSIONS
    for(int agridindex = 0; agridindex < maxprod; agridindex ++) {
        kappaWU[agridindex] = alphakappaw1 + alphakappaw2*(agridindex+1);
    }
    
    kappaEU = kappaE;
    
    int iter, igrid, ygrid, iterprice;
    
    // MOMENT POINTERS //
    int iteration, iteration2;
    double *genmoments, GMMmin;
    genmoments = (double *) calloc((nbmoments), sizeof(double));
    
    // Entrepreneur
    double *valueVE, *saveVE, *valueVEc, *searchEffortWWVE, *valueVEwUI, *saveVEwUI, *valueVEwUIc, *searchEffortWWVEwUI, *endoK;
    
    // Worker
    double *valueWW0, *saveWW0, *valueWW1, *saveWW1, *searchEffortVEWW;
    
    // Short-run Unemployed
    double *valueUS0, *saveUS0, *searchEffortWWUS0, *valueUS1, *saveUS1, *searchEffortWWUS1, *searchEffortVEUS1;
    
    // Long-run unemployed
    double *valueUL0, *saveUL0, *searchEffortWWUL0, *valueUL1, *saveUL1, *searchEffortWWUL1, *searchEffortVEUL1;
    
    // STATISTICS
    double capitalin, capitalout, laborin, laborout, taxin, taxout, critprice, critparams, taxrateout;
    
    // STARTING DISTRIBUTION
    double *start_distVE, *start_distVEwUI, *start_distUS0, *start_distUS1, *start_distWW0, *start_distWW1, *start_distUL0, *start_distUL1;


    // START COUNTING TIME.
    time_t rawtime,timeofstart,timeofend;
    struct tm * timeinfo;
    
    char buffer [80];
    
    FILE	*tempfileoutfile;
    
    //stop watch
    time ( &rawtime );
    time (&timeofstart);
    timeinfo = localtime ( &rawtime );
    
    strftime (buffer,80,"temp/tempfile_%Y%m%d@%H%M%S.out",timeinfo);
    strncpy(tempfileout, buffer, sizeof(tempfileout) - 1);
    tempfileout[sizeof(tempfileout) - 1] = 0;

    
    tempfileoutfile=fopen(tempfileout, "w");
    setbuf ( tempfileoutfile , NULL );
    fprintf(tempfileoutfile,"===================Program Starting: %s\n",asctime(timeinfo));
    fprintf(tempfileoutfile,"------ Estimated Model parameters------\n");
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","betapar",betapar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\t%20.15f\t%20.15f\n","Etapar",etapar[0],etapar[1],etapar[2]);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\t%20.15f\t%20.15f\n","Kappa",kappaWU[0], kappaWU[1],kappaWU[2]);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\t%20.15f\t%20.15f\n","PROD",prod[0], prod[1],prod[2]);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\t%20.15f\t%20.15f\n","GTRANS",gtrans[0], gtrans[1], gtrans[2]);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","phiparW",phiparW);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","phiparE",phiparE);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","costpar",costpar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","nupar",nupar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","fpar",fpar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","mupar",mupar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","xipar",xipar);
    fprintf(tempfileoutfile,"\n");
    fprintf(tempfileoutfile,"------ Fixed Model parameters------\n");
    fprintf(tempfileoutfile,"%20s\t\t%20d\n","maxgrid",maxgrid);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","Gridmin",Gridmin);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","Gridmax",Gridmax);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","corptax",corptax);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","alphapar",alphapar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","deltapar",deltapar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","rhostar",rhostar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","sigmapar",sigmapar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","TFPpar",TFPpar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","pLTpar",pLTpar);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","EffortMaxU",EffortMaxW);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","EffortMaxE",EffortMaxE);
    fprintf(tempfileoutfile,"\n");
    fprintf(tempfileoutfile,"------ Convergence criterions------\n");
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","epsilonValue",epsilonValue);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","epsilonendoK",epsilonendoK);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","maxiterVF",maxiterVF);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","maxiterendoK",maxiterendoK);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","epsilonDist",epsilonDist);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","epsilonprice",epsilonprice);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","epsilonparams",epsilonparams);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","relaxendoK",relaxendoK);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","relaxK",relaxK);
    fprintf(tempfileoutfile,"%20s\t\t%20.15f\n","relaxTax",relaxTax);
    fprintf(tempfileoutfile,"----------------------------\n");
    fclose(tempfileoutfile);
    
    //tempfileout="tempfile_20141202@1632.out";
    
    // SHOW PARAMETERS USED //
    for(int pa=0; pa < nbpara; pa++) {
        printf("%f \t", PARAM[pa]);
    }
    printf("\n");
    
    // INITIALIZE POINTERS //
    
    valueVE = (double *) calloc((ifulldimE), sizeof(double));
    valueVEc = (double *) calloc((ifulldimE), sizeof(double));
    valueVEwUI = (double *) calloc((ifulldimE), sizeof(double));
    valueVEwUIc = (double *) calloc((ifulldimE), sizeof(double));
    valueWW0 = (double *) calloc((ifulldim), sizeof(double));
    valueUS0 = (double *) calloc((ifulldim), sizeof(double));
    valueUL0 = (double *) calloc((ifulldim), sizeof(double));
    valueWW1 = (double *) calloc((ifulldim), sizeof(double));
    valueUS1 = (double *) calloc((ifulldim), sizeof(double));
    valueUL1 = (double *) calloc((ifulldim), sizeof(double));
    endoK = (double *) calloc((ifulldimKK), sizeof(double));
    
    saveVE = (double *) calloc((ifulldimE), sizeof(double));
    saveVEwUI = (double *) calloc((ifulldimE), sizeof(double));
    saveWW0 = (double *) calloc((ifulldim), sizeof(double));
    saveUS0 = (double *) calloc((ifulldim), sizeof(double));
    saveUL0 = (double *) calloc((ifulldim), sizeof(double));
    saveWW1 = (double *) calloc((ifulldim), sizeof(double));
    saveUS1 = (double *) calloc((ifulldim), sizeof(double));
    saveUL1 = (double *) calloc((ifulldim), sizeof(double));
    
    start_distVE = (double *) calloc((ifulldimE), sizeof(double));
    start_distVEwUI = (double *) calloc((ifulldimE), sizeof(double));
    start_distWW1 = (double *) calloc((ifulldim), sizeof(double));
    start_distWW0 = (double *) calloc((ifulldim), sizeof(double));
    start_distUS1 = (double *) calloc((ifulldim), sizeof(double));
    start_distUS0 = (double *) calloc((ifulldim), sizeof(double));
    start_distUL1 = (double *) calloc((ifulldim), sizeof(double));
    start_distUL0 = (double *) calloc((ifulldim), sizeof(double));
    
    searchEffortWWVE = (double *) calloc((ifulldimE), sizeof(double));
    searchEffortWWVEwUI = (double *) calloc((ifulldimE), sizeof(double));
    searchEffortVEWW = (double *) calloc((ifulldim), sizeof(double));
    searchEffortWWUL0 = (double *) calloc((ifulldim), sizeof(double));
    searchEffortWWUS0 = (double *) calloc((ifulldim), sizeof(double));
    searchEffortWWUL1 = (double *) calloc((ifulldim), sizeof(double));
    searchEffortWWUS1 = (double *) calloc((ifulldim), sizeof(double));
    searchEffortVEUL1 = (double *) calloc((ifulldim), sizeof(double));
    searchEffortVEUS1 = (double *) calloc((ifulldim), sizeof(double));
    
    
    start_distUL1[0]=1.0;
    
    
    // INITIALIZE starting value for pointers (Recall, should start by a point above) //
    for(igrid=0;igrid<maxgrid;igrid++)
    {
        grid[igrid] = phi(igrid);
        for(int ygrid = 0.0; ygrid < maxprod; ygrid ++)
        {
            for(int egrid = 0.0; egrid < maxfirmtype; egrid++) {
                valueVE[inxE(igrid,ygrid,egrid)]=0;
                valueVEc[inxE(igrid,ygrid,egrid)]=0;
                valueVEwUI[inxE(igrid,ygrid,egrid)]=0;
                valueVEwUIc[inxE(igrid,ygrid,egrid)]=0;
                }
            valueWW0[inx(igrid,ygrid)]=0;
            valueUS0[inx(igrid,ygrid)]=0;
            valueUL0[inx(igrid,ygrid)]=0;
            valueWW1[inx(igrid,ygrid)]=0;
            valueUS1[inx(igrid,ygrid)]=0;
            valueUL1[inx(igrid,ygrid)]=0;
        }
    }
    
    // PRICE (OR COUNTERPART) INITIALIZATION (no need to incorporate it in params since defined in header) //
    rstar = rstars;
    lstar = lstars;
    taxrate = taxrates;
    // ustar = 0.045; (should be useful only for agg fluctuation)
    // Compute corresponding aggregates
    capitalin=pow(((rstar+deltapar)/(alphapar)),(1.0/(alphapar-1.0)))*lstar;
    wstar=(1.0-alphapar)*pow(((rstar+deltapar)/alphapar),(alphapar/(alphapar-1.0)));
    
    
    #if indexPM == 1
    Projected_TP(200);

    printf("PROJECTED METHOD FINISHED"); getchar();
    #endif
    

    // INITIALIZE starting value for pointers (Recall, should start by a point above) //
    #if BRC == 1
    struct searchfendoKYY paramsKK;

    double testK1, testK2, testK3, testK4, Kmax;
    int iverifK;
    

    for(igrid=0;igrid<maxgrid;igrid++)
    {
        for(ygrid = 0; ygrid<maxprod; ygrid++) {
    
        for(int egrid = 0.0; egrid < maxfirmtype; egrid++)
        {
            paramsKK.igridKYY = igrid;
            paramsKK.egridKYY = egrid;
            paramsKK.ygridKYY = ygrid;
            
            Kmax = pow((nupar*gtrans[ygrid]*stateg[egrid])/(rstar + deltapar),1/(1-nupar));

            testK1 = searchendoK(0.0, &paramsKK);
            testK2 = searchendoK(0.01, &paramsKK);
            testK3 = searchendoK(Kmax, &paramsKK);
            testK4 = searchendoK(Kmax - 0.01, &paramsKK);

            iverifK = 0.0;
            
            //printf("%f %F %f %f %f", testK1,testK2,testK3,testK4, Kmax);

            if(testK1 > testK2 & testK1 < 0.0) {
                endoK[inxKK(igrid, ygrid, egrid)] = 0.0;
                iverifK = 100;
                //printf("I am here 1");getchar();
                }
            if(testK1 < testK2 & testK1 > 0.0) {
                endoK[inxKK(igrid, ygrid, egrid)] = 0.0;
                iverifK = 100;
                //printf("I am here 2");getchar();
                }
            if(testK3 > testK4 & testK3 < 0.0) {
                endoK[inxKK(igrid, ygrid, egrid)] = Kmax;
                iverifK = 100;
                //printf("I am here 3");getchar();
                }
            if(testK3 < testK4 & testK3 > 0.0) {
                endoK[inxKK(igrid, ygrid, egrid)] = Kmax;
                iverifK = 100;
                //printf("I am here 4 %f", Kmax);getchar();
                }
            
            if(iverifK == 0) {
                endoK[inxKK(igrid, ygrid, egrid)] = zbrentNEW(searchendoK,0.0,Kmax,(1.0e-10),&paramsKK);
            }

            if(endoK[inxKK(igrid, ygrid, egrid)] < 0.0) {
                endoK[inxKK(igrid, ygrid, egrid)] = 0.0;
                }
    
    
        }
        }
    
    }
    
    FILE *endoKK;
    
    endoKK=fopen(BRC1ENDO, "w");
    setbuf ( endoKK , NULL );
    
    fprintf(endoKK,"%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","ygrid","egrid","igrid","Index", "Grid", "Endok", "Profit");
    
    for(igrid=0;igrid<maxgrid;igrid++)
    {
        for(int egrid = 0.0; egrid < maxfirmtype; egrid ++)
        {
        for(int ygrid = 0.0; ygrid < maxprod; ygrid ++)
        {
            fprintf(endoKK,"%20d\t%20d\t%20d\t%20d\t%20.15f\t%20.15f\t%20.15f\n",ygrid, egrid, igrid, inxKK(igrid, ygrid, egrid), grid[igrid], endoK[inxKK(igrid,ygrid, egrid)], ((1.0 - corptax)*(TFPpar*gtrans[egrid]*pow(endoK[inxKK(igrid,ygrid, egrid)], nupar) - deltapar*endoK[inxKK(igrid,ygrid, egrid)] - rstar*endoK[inxKK(igrid,ygrid, egrid)])));
        }
        }
    }

    fclose(endoKK);
    
#if indexPMwrite == 1
    endoKK=fopen(endoKfile, "w");
    setbuf ( endoKK , NULL );
    
    // SHOULD BE IN THE SAME ORDER AS THE READING FILE (i.E. the index).
    for(int ygrid = 0.0; ygrid < maxprod; ygrid ++)
    {
        for(igrid=0;igrid<maxgrid;igrid++)
        {
            for(int egrid = 0.0; egrid < maxfirmtype; egrid ++)
            {
                fprintf(endoKK,"%20.15f\n",endoK[inxKK(igrid,ygrid, egrid)]);
            }
        }
    }

    fclose(endoKK);
#endif
    
    #endif
    
    
    printf("Initialization OK \n");
    
    // critparams = 1.0;
    // while(critparams > epsilonparams) {
    
    critprice = 1.0;
    iterprice = 0.0;
    while(critprice > epsilonprice) // loop over r and w (and possibly u)
    {
        
        printf("STARTING MODEL \n");
        printf("rstar %f, wstar %f, capitalin %f, lstar, %f, taxrate %f \n",rstar,wstar, capitalin, lstar, taxrate);
        
        // INSURE that you are in an incomplete market case //
		if (betapar*(1.0+rstar)>=1 || rstar < 0.015 || capitalin < 0.0 || iterprice > itermax)
		{
            printf("beta condition %20.15f\t%20.15f\n",rstar,betapar*(1.0+rstar));
            
            critprice = 0.0;
     
            for(int l = 0.0; l < nbmoments; l++)
            {
                genmoments[l] = -10000;
            }
		} else {
        
        
        tempfileoutfile=fopen(tempfileout, "a");
        fprintf(tempfileoutfile,"R: %20.15f,\t W: %20.15f,\t K: %20.15f\t, L: %20.15f, \t TAX: %20.15f\n",rstar,wstar,capitalin, lstar, taxrate);
        fprintf(tempfileoutfile, "\n");
        fclose(tempfileoutfile);
        
        
        // VFI Computation
        VFI(valueWW0, valueWW1, valueUS0, valueUS1, valueUL0, valueUL1, valueVE, valueVEc, valueVEwUI, valueVEwUIc, endoK, saveWW0, saveWW1, saveUS0, saveUS1, saveUL0, saveUL1, saveVE, saveVEwUI, searchEffortWWUS0, searchEffortWWUL0, searchEffortWWUS1, searchEffortVEUS1, searchEffortWWUL1, searchEffortVEUL1, searchEffortWWVE, searchEffortVEWW, searchEffortWWVEwUI);
            
        printf("VFI ok\n");
        
        // SIMULATION
        #if indexPM == 0
        SIMULATION(valueWW0, valueWW1, valueUS0, valueUS1, valueUL0, valueUL1, valueVE, valueVEc, valueVEwUI, valueVEwUIc, endoK, saveWW0, saveWW1, saveUS0, saveUS1, saveUL0, saveUL1, saveVE, saveVEwUI, searchEffortWWUS0, searchEffortWWUL0, searchEffortWWUS1, searchEffortVEUS1, searchEffortWWUL1, searchEffortVEUL1, searchEffortWWVE, searchEffortVEWW, searchEffortWWVEwUI, &capitalout,&laborout,&taxout,genmoments, start_distVE, start_distWW1, start_distWW0, start_distUL0, start_distUL1, start_distUS1, start_distUS0, start_distVEwUI);
        #endif
        
        printf("SIMULATION ok\n");   //getchar();
            
        // Compute convergence
        critprice=fabs((laborout-lstar)/lstar);
        critprice=max(fabs((capitalout-capitalin)/capitalin), critprice);
        critprice=max(fabs((taxout-taxrate)/taxrate), critprice);
        
        printf("TAX: simul: %20.10f\t, old: %20.10f\t, new: %20.10f\n",taxout, taxrate, (relaxTax*taxout+(1.0-relaxTax)*taxrate));
        printf("Capital: simul: %20.10f\t, old: %20.10f\t, new: %20.10f\n",capitalout, capitalin, (relaxK*capitalout+(1.0-relaxK)*capitalin));
        printf("agglab cnvg: crit: %20.10f\t %20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\n",critprice,lstar,laborout,capitalout,capitalin, taxout, taxrate);

        // Compute Next Aggregate
        lstar = (relaxL*laborout+(1.0-relaxL)*lstar);
        taxrate = (relaxTax*taxout+(1.0-relaxTax)*taxrate);
        capitalin = (relaxK*capitalout+(1.0-relaxK)*capitalin);
        rstar=alphapar*pow(capitalin,(alphapar-1.0))*pow(lstar,(1.0-alphapar))-deltapar;
        wstar=(1.0-alphapar)*pow(((rstar+deltapar)/alphapar),(alphapar/(alphapar-1.0)));
        
        } // end on BETA condition
        
        iterprice++;
        
    }//while(critprice > epsilonprice)
    

    // COMPUTE THE MOMENTS //
	GMMmin = SMM(obsmoments, genmoments, covmat);
    
    // CREATE A LOG FILE ///
	FILE *logfile;
	logfile = fopen(tempfileout, "a");
    fprintf(logfile,"\n GMM :: %f\n\n", GMMmin);
	fclose(logfile);
    
    logfile = fopen(PARA, "a");
    fprintf(logfile,"\n GMM :: %f \t", GMMmin);
    for(int pa=0; pa < nbpara; pa++) {
        fprintf(logfile,"%f \t", PARAM[pa]);
    }
    fprintf(logfile,"\n");
	fclose(logfile);
	
    
    // SHOW RESULTS ///
    printf("\n GMM :: %f \t", GMMmin);
    
    
    //stop watch
    time ( &rawtime );
    time (&timeofend);
    timeinfo = localtime ( &rawtime );
    tempfileoutfile=fopen(tempfileout, "a");
    setbuf ( tempfileoutfile , NULL );
    fprintf(tempfileoutfile,"===================Program ended: %s in %f seconds\n",asctime(timeinfo),difftime(timeofend,timeofstart));
    fclose(tempfileoutfile);
    
    printf("\a\a\a");
 
    return GMMmin;

}
    






int main(int argc, char* argv[])
{
  

#if CRS == 0
double *PARAMETERS;
PARAMETERS = (double *) calloc((nbpara), sizeof(double));

PARAMETERS[0] = betapar;
PARAMETERS[1] = phiparW;
PARAMETERS[2] = phiparE;
PARAMETERS[3] = alphakappaw1;
PARAMETERS[4] = kappaE;
PARAMETERS[5] = g1;
PARAMETERS[6] = g2;
PARAMETERS[7] = g3;
PARAMETERS[8] = fpar;
PARAMETERS[9] = nupar;
PARAMETERS[10] = stateg[0];
PARAMETERS[11] = pibad;
PARAMETERS[12] = kappaW;
PARAMETERS[13] = mupar;

ComputeEquilibrium(PARAMETERS);

printf("PROGRAM FINISHED"); getchar();
#endif
    

#if CRS == 1

// SET OBSERVED MOMENT //
obsmoments[0] = 0.304;      // Capital - output ratio (/10)
obsmoments[1] = 1.09;      // Fraction of Entrepreneurs (*10)
obsmoments[2] = 0.5;        // Unemployment rate (*10)
obsmoments[3] = 0.4;        // Top 1%
obsmoments[4] = 0.84;        // Gini Wealth
obsmoments[5] = 0.2;        // Necessity share (fraction of U to E)
obsmoments[6] = 0.525;     // U exit rate
obsmoments[7] = 0.226;      // Fraction E to WW
obsmoments[8] = 0.0475;     // Worker Exit rate
obsmoments[9] = 0.22;       // U-shape ent 0
obsmoments[10] = 0.24;      // U-shape ent 2
obsmoments[11] = 0.22;       // U-shape ent 4
obsmoments[12] = 0.2;      // U-shape WW 0 (*10)
obsmoments[13] = 0.12;      // U-shape WW 2 (*10)
obsmoments[14] = 0.2;      // U-shape WW 4 (*10)
obsmoments[15] = 0.36;      // Gini earnings

for(int l=0.0; l<nbmoments; l++){
    covmat[l][l] = 1.0;
}


/////////////////////////////
// STARTING CRS ALGORITHME //
/////////////////////////////

    FILE *soboloutfile;
    FILE *logfile;
    
	logfile = fopen(PARA, "w");
    fprintf(logfile,"\n");
	fclose(logfile);
    

    double *Best_PARA, *param_max, *param_min, ini_val2[nbsobol][nbpara];
    Best_PARA = (double *) calloc((nbpara), sizeof(double));
    // PARAMETERS MIN MAX ///
    param_max = (double *) calloc((nbpara), sizeof(double));
    param_min = (double *) calloc((nbpara), sizeof(double));


    param_min[0] = 0.865;     param_max[0] = 0.9; // BETA
    param_min[1] = 1.1;    param_max[1] = 2.0;  // phiparW
    param_min[2] = 1.1;    param_max[2] = 3.0;  // phiparE
    param_min[3] = 1.0;     param_max[3] = 3.0;  // alphakappaw1
    param_min[4] = 0.3;    param_max[4] = 1.0;  // kappaE
    param_min[5] = 0.3;     param_max[5] = 0.6;  // g1
    param_min[6] = 0.35;     param_max[6] = 0.65;  // g2
    param_min[7] = 0.7;     param_max[7] = 0.75;  // g3
    param_min[8] = 0.6;    param_max[8] = 0.9;  // f
    param_min[9] = 0.8;   param_max[9] = 0.9; // nupar
    param_min[10] = 0.2;   param_max[10] = 0.8; // statebad
    param_min[11] = 0.05;   param_max[11] = 0.25; // pibad
    param_min[12] = 3.0;   param_max[12] = 8.0; // kappaW
    param_min[13] = 0.03;   param_max[13] = 0.15; // mupar

    /// INITIALIZE SOBOL SEQUENCE //////
    //SOBOL GENERATION
//    float *sobolinitial;
//    sobolinitial = i4_sobol_generate(nbsobol, nbpara, 100);
//    
//    //SAVE Sobol
//    soboloutfile=fopen(SOBOL, "w");
//    setbuf ( soboloutfile , NULL );
//    for(int i = 0; i < nbsobol*nbpara; i++) {
//       fprintf(soboloutfile,"%20.15f\n", sobolinitial[i]);
//    }
//    fclose(soboloutfile);

    // OR load sobol
    double *sobolinitial;
    sobolinitial = (double *) calloc((nbsobol*nbpara), sizeof(double));
    readinput(sobolinitial, nbsobol*nbpara, "SOBOL.out");

    
    // b) Compute the corresponding points from the sequence
    for(int i = 0; i < nbsobol; i++) {
        for(int j = 0; j < nbpara; j++) {
            ini_val2[i][j] = param_min[j] + sobolinitial[nbpara*i + j] * (param_max[j] - param_min[j]);
        }
    }
    
   // printf("%f %f %f", ini_val2[0][1], ini_val2[0][2], ini_val2[0][3]);
    
    // WITHOUT OMP //
    myCRS(param_min, param_max, ini_val2, Best_PARA, ComputeEquilibrium);
    
    // WITH OMP //
    #if OMP == 2
    myCRSomp(param_min, param_max, ini_val2, Best_PARA, ComputeEquilibrium);
    #endif


    for(int i = 0; i < nbpara; i++) {
       printf("%20.15f\t", Best_PARA[i]);
    }
    
    printf("PROGRAM FINISHED"); getchar();
#endif

return 0;
}
