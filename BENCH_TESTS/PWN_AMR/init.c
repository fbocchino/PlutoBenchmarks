#include "pluto.h"


static double IsentropicPressure (double rho, double s);

/* ************************************************************** */
void Init (double *v, double x1, double x2, double x3)
/* 
 * Set initial conditions for the Crab Nebula simulations
 * in 2D Cylindrical or 3D Cartesian coordinates.
 * 
 * Authors:   A. Mignone
 *            B. Olmi
 *            N. Bucciantini
 *
 * Last modified:   06.12.2023
 *  
 **************************************************************** */
{
  double r, s;
  double f, fk, fm, sncth;
  double vr, g, th, ph; 
  double bc, signz, Lp;
  double dph, rho0, p0, b0, Bphi;
  double L, c, q;
  double r_ej, m_ej, E_sn, rho_ej, v_ej, p_ej;
  double rho_ism, v_ism, p_ism;    


  double r_wind   = g_inputParam[RPW];
  double sig      = g_inputParam[SIGMA_0];
  double glf      = g_inputParam[GAMMA_0];
  double Lpwn     = g_inputParam[WL_0];



  // -- constants and units
  L = 9.48e17;
  c = 3.e10;
  q = 1.67e-24;
  
  r_ej    = 5.0;
  m_ej    = 6.e33;
  E_sn    = 1.e51;
  rho_ej  = 3.*m_ej/(4.*CONST_PI*(r_ej*r_ej*r_ej)*(L*L*L)*q);
  v_ej    = sqrt(10./3.*E_sn/(m_ej*c*c));
  p_ej    = IsentropicPressure(rho_ej, 1.e-10);
  rho_ism = 1.0;
  v_ism   = 0.0;
  p_ism   = 1.e-9*rho_ism;    

   
  // definition of eos and adiabatic index
  #if EOS == IDEAL
    g_gamma = 4.0/3.0;
  #endif

  // -- GEOMETRY AND ROTATIONS
  //-- only valid for CARTESIAN GEOMETRY --  
  r    = sqrt(x1*x1 + x2*x2 + x3*x3);           /* spherical r */
  ph   = atan2(x2,x1);                           /* phi */               
  th   = acos(x3/r);    
  s    = sin(th);       /* -- sin(theta) -- */
  if (th <= CONST_PI/2.){
     signz=1.0;
   } else{
     signz=-1.0;
   }
   
  // ---------------------------------------
  
  bc     = 0.03;
  Lp     = CONST_PI*(8.0/3.0 + bc);
  rho0   = Lpwn/(Lp*L*L*glf*glf*c*c*c*q);  
  p0     = 1.e-2*rho0;
  b0     = sqrt(Lpwn/Lp * 4.0*CONST_PI/c)/L; 




  /* -- Define Fluxes ----------------*/
  f    = ( sin(th)*sin(th) + bc );
	fk   = f/(1.0 + sig);
	fm   = f*sig/(1.0 + sig);


	if (r <= r_wind){
	
	  v[RHO]   = rho0/(r*r)*fk;
    vr       = sqrt(1.0 - 1.0/(glf*glf));
    v[PRS]   =  MAX(v[PRS],(v[RHO]*vr*vr)*1.e-3);
    v[TRC]   = 1.0;
    Bphi     = signz * b0/r * sqrt(fm);

 	} else if (r <= r_ej) {
  /* --------------------------------------
      Ejecta
     -------------------------------------- */

    v[RHO]   = rho_ej;
    vr       = v_ej*r/r_ej;
    v[PRS]   = p_ej;
    Bphi     = 0.0;
    v[TRC]   = 0.0;
      
    
  } else {
  /* -----------------------------------
      ISM
     ----------------------------------- */
  
    v[RHO]   = rho_ism;
    vr       = v_ism;
    v[PRS]   = p_ism;
    v[TRC]   = 0.0;
    Bphi     = 0.0;

  }




#if PHYSICS == RMHD
  v[BX1] = v[BX2] = v[BX3] = 0.0;
#endif


 // -----velocity--------
  v[VX1]  = vr*s*cos(ph);
  v[VX2]  = vr*s*sin(ph);
  v[VX3]  = vr*x3/r;
  // ----magnetic field-------
  v[BX1]  = - Bphi*sin(ph);
  v[BX2]  =   Bphi*cos(ph);
  v[BX3]  =   0.0 ; 


/***************************************************/
/***************************************************/


  #ifdef GLM_MHD
    v[PSI_GLM] = 0.0;
  #endif  
     

  g_smallDensity  = 1.e-10;
  g_smallPressure = 1.e-12;


}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{

}
/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 * PURPOSE
 *  
 *   Perform some pre-processing data
 *
 * ARGUMENTS
 *
 *   d:      the PLUTO Data structure.
 *   grid:   pointer to array of GRID structures  
 *
 **************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double vin[256], r;
 
/* -----------------------------------------------------------
    Define size of reset region (internal boundary).
    Since coarser levels may not even see this small internal
    radius, we floor r_reset to at least two zones.
  ------------------------------------------------------------ */
  
  double r_reset = g_inputParam[RRES];


  r_reset = MAX(r_reset,2.0*grid->x[IDIR][IBEG]);
  r_reset = MAX(r_reset,2.0*grid->x[JDIR][JBEG]); 
  r_reset = MAX(r_reset,2.0*grid->x[KDIR][KBEG]);  


  //print("r_reset 2 (%f), (%f), (%f), (%f)\n ", r_reset, 2.0*grid[IDIR].dx[IBEG],2.0*grid[JDIR].dx[JBEG], 2.0*grid[KDIR].dx[KBEG]);

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  
  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], g_smallDensity);
      d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], g_smallPressure);
      d->Vc[TRC][k][j][i] = MAX(d->Vc[TRC][k][j][i], g_smallDensity);

#if GEOMETRY == SPHERICAL
      r = x1[i];
#elif GEOMETRY == CYLINDRICAL
      r = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
#elif GEOMETRY == CARTESIAN
      r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
#endif
      if (r <= r_reset){
        Init (vin, x1[i], x2[j], x3[k]);
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = vin[nv];
        d->flag[k][j ][i] |= FLAG_INTERNAL_BOUNDARY;
      }else if (r < -0.06) {
        Init (vin, x1[i], x2[j], x3[k]);
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][k][j][i]
             + (vin[nv] - d->Vc[nv][k][j][i])*(0.06-r)/0.01;
        }
      }
    }
  }
}

/* ********************************************************************* */
double IsentropicPressure (double rho, double s)
/*!
 *  Compute pressure from given density (rho) and entropy (s).
 *  (see PLUTO/Doc/Notes/rmhd.pdf)
 *********************************************************************** */
{
  double rho23, c,t,p;

#if EOS == IDEAL
  p = s*pow(rho,g_gamma);
#elif EOS == TAUB
  rho23 = pow(rho,2.0/3.0);
  c     = 1.5*s*rho23;
  t     = 2.0*c/(1.0 + sqrt(1.0 + 4.0*(1.0 + c)*c));
  p     = 2.0/3.0*rho*t/sqrt(1.0 - t*t);
#endif
  return p;
    
}




