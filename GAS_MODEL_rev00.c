#include "udf.h"
#include "stdio.h"
#include "ctype.h"
#include "stdarg.h"



#define rgas UNIVERSAL_GAS_CONSTANT

//NASA coeffcient NO2:
#define a1
#define a2
#define a3
#define a4
#define a5
#define a6
#define a7
#define a8

static int (*usersMessage)(char *,...);
static void (*usersError)(char *,...);

DEFINE_ON_DEMAND(I_do_nothing)
{
   /* this is a dummy function to allow us */
   /* to use the compiled UDFs utility   */
}


/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_error                                      */
/*------------------------------------------------------------*/

void RKEOS_error(int err, char *f, char *msg)
{
 if (err)
     usersError("RKEOS_error (%d) from function: %s\n%s\n",err,f,msg);
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_Setup                                      */
/*------------------------------------------------------------*/

void RKEOS_Setup(Domain *domain, cxboolean vapor_phase, char *filename, int
                 (*messagefunc)(char *format,...),
                 void (*errorfunc)(char *format,...))
{

  usersMessage = messagefunc;
  usersError = errorfunc;
  usersMessage("\nLoading Ideal Gas Model: %s\n", filename);
}

/*------------------------------------------------------------*/
/* równanie stanu gazu                                 */
/*------------------------------------------------------------*/

double RKEOS_pressure(double temp, double density)
{
    double v = 1./density;
  
    return rgas*temp/v;
}

/*------------------------------------------------------------*/
/* równanie stanu gazu sprowadzone do równania wielomianiowego
w funkcji objetosci                                     */
/*------------------------------------------------------------*/

double RKEOS_spvol(double temp, double press)
{
    double a1,a2,a3;


    a1 = press;
    a2 = -rgas*temp;

  return vv=a1*v+a2;
  
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_density                                    */
/*      Returns density given T and P                         */
/*------------------------------------------------------------*/

double RKEOS_density(cxboolean vapor_phase, double temp, double press, double yi[])
{
    return 1./RKEOS_spvol(temp, press); /* (Kg/m3) */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_dvdp                                       */
/*      Returns dv/dp given T and rho                         */
/*------------------------------------------------------------*/

double RKEOS_dvdp(double temp, double density)
{
    press = RKEOS_pressure(temp, density);
    v = 1./density;

    return v;
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_dvdt                                       */
/*      Returns dv/dT given T and rho                         */
/*------------------------------------------------------------*/

double RKEOS_dvdt(double temp, double density)
{
    press = RKEOS_pressure(temp, density);
    v = 1./density;

    return -rgas;
}



/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_specific_heat                              */
/*      Returns specific heat given T and rho                  */
/*------------------------------------------------------------*/

double RKEOS_specific_heat(double temp)
{
    NASA POLYNOMIAL

    return RKEOS_Cp_ideal_gas(temp)+delta_Cp; /* (J/Kg-K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_enthalpy                                   */
/*      Returns specific enthalpy given T and rho             */
/*------------------------------------------------------------*/

double RKEOS_enthalpy(double temp)
{
    NASA POLYNOMIAL

    return H_REF+RKEOS_H_ideal_gas(temp)+delta_h; /* (J/Kg) */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_entropy                                    */
/*      Returns entropy given T and rho                       */
/*------------------------------------------------------------*/

double RKEOS_entropy(double temp, double density, double P, double yi[])
{
    NASA POLYNOMIAL

    return S_REF+cp_integral+delta_s; /* (J/Kg-K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_mw                                         */
/*      Returns molecular weight                              */
/*------------------------------------------------------------*/

double RKEOS_mw(double yi[])
{
    return MWT; /* (Kg/Kmol) */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_speed_of_sound                             */
/*      Returns s.o.s given T and rho                         */
/*------------------------------------------------------------*/

double RKEOS_speed_of_sound(double temp)
{
    double kappa=1.4;

    return sqrt(kappa*rgas*temp); /* m/s */
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_rho_t     drho/dT                          */
/*------------------------------------------------------------*/

double RKEOS_rho_t(double temp, double density, double P, double yi[])
{
    press = RKEOS_pressure(temp, density);
	return press/rgas*l;
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_rho_p      drho/dp                         */
/*------------------------------------------------------------*/

double RKEOS_rho_p(double temp, double density,double P, double yi[])
{
	
    return -density*density*RKEOS_dvdp(temp,  density);
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_enthalpy_t                                 */
/*------------------------------------------------------------*/

double RKEOS_enthalpy_t(double temp, double density, double P, double yi[])
{
     return -density*density*RKEOS_dvdt(temp,  density);
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_enthalpy_p                                 */
/*------------------------------------------------------------*/

double RKEOS_enthalpy_p(double temp, double density, double P, double yi[])
{
    double v = 1./density;
    double dvdt = RKEOS_dvdt(temp, density);

    return v-temp*dvdt;
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_viscosity                                  */
/*------------------------------------------------------------*/

double RKEOS_viscosity(double temp, double density, double P, double yi[])
{
  //jakaś ambitna funkcja lepkości , jednostki mikropascalosekudny przeliczyć

    return mu=0.046;
}

/*------------------------------------------------------------*/
/* FUNCTION: RKEOS_thermal_conductivity                       */
/*------------------------------------------------------------*/

double RKEOS_thermal_conductivity(double temp, double density, double P, 
								 double yi[])
{
    double cp, mu;

    cp = RKEOS_Cp_ideal_gas(temp);
    mu = RKEOS_viscosity(temp, density, P, yi);

    return (cp+1.25*rgas)*mu;
}

/* Export Real Gas Functions to Solver */

UDF_EXPORT RGAS_Functions RealGasFunctionList =
{
    RKEOS_Setup,                 /* initialize           */
    RKEOS_density,               /* density              */
    RKEOS_enthalpy,              /* enthalpy             */
    RKEOS_entropy,               /* entropy              */
    RKEOS_specific_heat,         /* specific_heat        */
    RKEOS_mw,                    /* molecular_weight     */
    RKEOS_speed_of_sound,        /* speed_of_sound       */
    RKEOS_viscosity,             /* viscosity            */
    RKEOS_thermal_conductivity,  /* thermal_conductivity */
    RKEOS_rho_t,                 /* drho/dT |const p     */
    RKEOS_rho_p,                 /* drho/dp |const T     */
    RKEOS_enthalpy_t,            /* dh/dT |const p       */
    RKEOS_enthalpy_p             /* dh/dp |const T       */
};