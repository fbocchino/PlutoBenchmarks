#include <stdlib.h>
#include <math.h>
#include "pluto.h"

/*   USER DEFINITIONS  */

/* *********************************************************************
                             Function Prototypes
   ********************************************************************* */

double linear_interp(double x, double x1, double x0, double y1, double y0);

void ambient_field(double, double, double, double *);
void ambient_hydro(double, double, double, double *);
void csm_abundance(double, double, double, double *);


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{

  double Bfield[3];
  double hydroV[8];
  double abundance[19];

  double mu;
  double rho_zone, pres_zone;
  double bx_zone, by_zone, bz_zone;
  double time_sh, trac_ej, trac_trs, trac_clp;
  double lagr_x1, lagr_x2, lagr_x3;
  double vx_zone, vy_zone, vz_zone, vel_rad;

/* ------------------------------------------------------------------
         Set PLUTO parameters
   ------------------------------------------------------------------ */

  g_minCoolingTemp = 5.e6;
  g_maxCoolingRate = 0.9;
  g_smallDensity   = 1.e-20;
  g_smallPressure  = 4.e-09;

  g_time = g_inputParam[initime];

/* ------------------------------------------------------------------
         Set units
   ------------------------------------------------------------------ */

  mu = MU_CSM/2.;

/* Other units, for reference only, defined as in PLUTO UG:  */
  double unit_pressure = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  double unit_time     = UNIT_LENGTH / UNIT_VELOCITY; //approx default 1000 years
  double unit_b        = UNIT_VELOCITY * sqrt(4.0 * CONST_PI * UNIT_DENSITY);
  double unit_temp     = mu*CONST_mp*unit_pressure/(2.0*UNIT_DENSITY*CONST_kB);

/* ------------------------------------------------------------------
         Derive initial conditions for CSM
   ------------------------------------------------------------------ */

  /* Define magnetic field ------------------------------------------ */

  ambient_field(x1, x2, x3, Bfield);

  bx_zone = Bfield[0]/unit_b;
  by_zone = Bfield[1]/unit_b;
  bz_zone = Bfield[2]/unit_b;

  /* Define ambient velocities -------------------------------------- */

  vx_zone = 0.0;
  vy_zone = 0.0;
  vz_zone = 0.0;
  vel_rad = 0.0;

  /* Define ambient pressure, density, and tracers ------------------ */

  ambient_hydro(x1, x2, x3, hydroV);

  rho_zone  = hydroV[0];
  pres_zone = hydroV[1];
  trac_trs  = hydroV[2];
  trac_clp  = hydroV[3];
  trac_ej   = hydroV[4];
  lagr_x1   = hydroV[5];
  lagr_x2   = hydroV[6];
  lagr_x3   = hydroV[7];

  time_sh   = 0.0;

/* ------------------------------------------------------------------ */

  /* Define SNR pressure, density, velocity ------------------------- */

  v[RHO]   = rho_zone;
  v[PRS]   = pres_zone;
  v[VX1]   = vx_zone;
  v[VX2]   = vy_zone;
  v[VX3]   = vz_zone;
  v[TRC]   = time_sh;
  v[TRC+1] = trac_trs;
  v[TRC+2] = trac_clp;
  v[TRC+3] = trac_ej;
  v[TRC+4] = lagr_x1;
  v[TRC+5] = lagr_x2;
  v[TRC+6] = lagr_x3;
  v[TRC+7] = vel_rad;

  /* Define abundances of CSM --------------------------------------- */

  csm_abundance(x1, x2, x3, abundance);

  v[TRC+8]  = abundance[0];   /* Ar36 */
  v[TRC+9]  = abundance[1];   /* C12 */
  v[TRC+10] = abundance[2];   /* Ca40 */
  v[TRC+11] = abundance[3];   /* Cr48 */
  v[TRC+12] = abundance[4];   /* Fe52 */
  v[TRC+13] = abundance[5];   /* Fe54 */
  v[TRC+14] = abundance[6];   /* H1 */
  v[TRC+15] = abundance[7];   /* He3 */
  v[TRC+16] = abundance[8];   /* He4 */
  v[TRC+17] = abundance[9];   /* Mg24 */
  v[TRC+18] = abundance[10];  /* N14 */
  v[TRC+19] = abundance[11];  /* Ne20 */
  v[TRC+20] = abundance[12];  /* neut */
  v[TRC+21] = abundance[13];  /* Ni56 */
  v[TRC+22] = abundance[14];  /* O16 */
  v[TRC+23] = abundance[15];  /* prot */
  v[TRC+24] = abundance[16];  /* S32 */
  v[TRC+25] = abundance[17];  /* Si28 */
  v[TRC+26] = abundance[18];  /* Ti44 */

  /* Define circumstellar magnetic field ---------------------------- */

  #if PHYSICS == MHD || PHYSICS == RMHD
    v[BX1] = bx_zone;
    v[BX2] = by_zone;
    v[BX3] = bz_zone;
  
    #ifdef STAGGERED_MHD
      us[AX1] = 0.0;
      us[AX2] = 0.0;
      us[AX3] = 0.0;
    #endif
  #endif

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

  int i,j,k,id;

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

/* -------------------------------------------------------------------------- */
/*         READ THE CSM CALCULATED WITH A DEDICATED HD SIMULATION             */

/* ---------- Density ---------- */
  id = InputDataOpen ("CSM/rho.0000.flt","CSM/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
  }
  InputDataClose(id);

/* -------------------------------------------------------------------------- */

  double small_domain = g_inputParam[s_dom];

/* ---------- time when the cell was shocked ---------- */
    id = InputDataOpen ("EXT/tr1.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) ){
         d->Vc[TRC][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- Density ---------- */
  id = InputDataOpen ("EXT/rho.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[RHO][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- vel x1 ---------- */
  id = InputDataOpen ("EXT/vx1.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[VX1][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);
  
/* ---------- vel x2 ---------- */
  id = InputDataOpen ("EXT/vx2.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[VX2][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- vel x3 & lagr. vel ---------- */
  id = InputDataOpen ("EXT/vx3.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[VX3][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);

       d->Vc[TRC+7][k][j][i] = sqrt((d->Vc[VX1][k][j][i]) * (d->Vc[VX1][k][j][i]) +
                                    (d->Vc[VX2][k][j][i]) * (d->Vc[VX2][k][j][i]) +
                                    (d->Vc[VX3][k][j][i]) * (d->Vc[VX3][k][j][i]));
    }
  }
  InputDataClose(id);

/* ---------- pressure ---------- */
  id = InputDataOpen ("EXT/prs.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[PRS][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- tracer of ejecta & time_sh & lagr. coord. ---------- */
  id = InputDataOpen ("EXT/tr4.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
         (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+3][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
       d->Vc[TRC+4][k][j][i] = x1[i];
       d->Vc[TRC+5][k][j][i] = x2[j];
       d->Vc[TRC+6][k][j][i] = x3[k];
    }
  }
  InputDataClose(id);

/* ==================================== start IF =================================== */
  if (g_inputParam[EXT] > 0){  

    #if PHYSICS == MHD || PHYSICS == RMHD

/* ---------- B x1 ---------- */
      id = InputDataOpen ("EXT/Bx1.0000.flt","EXT/grid.out"," ", 0, CENTER);
      TOT_LOOP(k,j,i){
        if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
             (x2[j] > -small_domain) && (x2[j] < small_domain) &&
             (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
           d->Vc[BX1][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
        }
      }
      InputDataClose(id);
  
/* ---------- B x2 ---------- */
      id = InputDataOpen ("EXT/Bx2.0000.flt","EXT/grid.out"," ", 0, CENTER);
      TOT_LOOP(k,j,i){
        if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
             (x2[j] > -small_domain) && (x2[j] < small_domain) &&
             (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
           d->Vc[BX2][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
        }
      }
      InputDataClose(id);
  
/* ---------- B x3 ---------- */
      id = InputDataOpen ("EXT/Bx3.0000.flt","EXT/grid.out"," ", 0, CENTER);
      TOT_LOOP(k,j,i){
        if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
             (x2[j] > -small_domain) && (x2[j] < small_domain) &&
             (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
           d->Vc[BX3][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
        }
      }
      InputDataClose(id);

    #endif

/* ---------- tracer of HII region ---------- */
    id = InputDataOpen ("EXT/tr2.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+1][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- tracer of ring material ---------- */
    id = InputDataOpen ("EXT/tr3.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+2][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- lagrangian coordinate x1 ---------- */
    id = InputDataOpen ("EXT/tr5.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+4][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- lagrangian coordinate x2 ---------- */
    id = InputDataOpen ("EXT/tr6.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+5][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- lagrangian coordinate x3 ---------- */
    id = InputDataOpen ("EXT/tr7.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+6][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);

/* ---------- velocity when the cell was shocked ---------- */
    id = InputDataOpen ("EXT/tr8.0000.flt","EXT/grid.out"," ", 0, CENTER);
    TOT_LOOP(k,j,i){
      if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
           (x2[j] > -small_domain) && (x2[j] < small_domain) &&
           (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
         d->Vc[TRC+7][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
      }
    }
    InputDataClose(id);
  }
/* ====================================== end IF ===================================== */

/* ---------- Ar 36 ---------- */
  id = InputDataOpen ("EXT/tr9.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+8][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- C 12 ---------- */
  id = InputDataOpen ("EXT/tr10.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+9][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Ca 40 ---------- */
  id = InputDataOpen ("EXT/tr11.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+10][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Cr 48 ---------- */
  id = InputDataOpen ("EXT/tr12.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+11][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Fe 52 ---------- */
  id = InputDataOpen ("EXT/tr13.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+12][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Fe 54 ---------- */
  id = InputDataOpen ("EXT/tr14.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+13][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- H 1 ---------- */
  id = InputDataOpen ("EXT/tr15.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+14][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- He 3 ---------- */
  id = InputDataOpen ("EXT/tr16.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+15][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- He 4 ---------- */
  id = InputDataOpen ("EXT/tr17.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+16][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Mg 24 ---------- */
  id = InputDataOpen ("EXT/tr18.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+17][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- N 14 ---------- */
  id = InputDataOpen ("EXT/tr19.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+18][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Ne 20 ---------- */
  id = InputDataOpen ("EXT/tr20.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+19][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Neut ---------- */
  id = InputDataOpen ("EXT/tr21.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+20][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Ni 56 ---------- */
  id = InputDataOpen ("EXT/tr22.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+21][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- O 16 ---------- */
  id = InputDataOpen ("EXT/tr23.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+22][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- prot ---------- */
  id = InputDataOpen ("EXT/tr24.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+23][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- S 32 ---------- */
  id = InputDataOpen ("EXT/tr25.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+24][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Si 28 ---------- */
  id = InputDataOpen ("EXT/tr26.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+25][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

/* ---------- Ti 44 ---------- */
  id = InputDataOpen ("EXT/tr27.0000.flt","EXT/grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i){
    if ( (x1[i] > -small_domain) && (x1[i] < small_domain) &&
         (x2[j] > -small_domain) && (x2[j] < small_domain) &&
         (x3[k] > -small_domain) && (x3[k] < small_domain) &&
             (d->Vc[TRC][k][j][i] > 1.e-10) ){
       d->Vc[TRC+26][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

}


/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/*
 *
 ****************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox * box, int side, Grid *grid) 
/*
 * Sets inflow boundary condition at the top boundary (side == X2_END)
 * and the stellar wind region when side == 0.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double xx1, xx2, xx3;
  double bx_zone, by_zone, bz_zone;
  double rho_zone, pres_zone, trac_trs, trac_clp, trac_ej;

  double Bfield[3];
  double hydroV[8];
  double abundance[19];

  double unit_pressure = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  double unit_b        = UNIT_VELOCITY * sqrt(4.0 * CONST_PI * UNIT_DENSITY);

  x1 = grid->xgc[IDIR];
  x2 = grid->xgc[JDIR];
  x3 = grid->xgc[KDIR];

/* ---------------------------------------------------------------------- */

  if (side == X1_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X2FACE){  /* -- y staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX2s][k][j][i] = Bfield[1]/unit_b;
      }

    #endif
    } else if (box->vpos == X3FACE){  /* -- z staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX3s][k][j][i] = Bfield[2]/unit_b;
      }

    #endif
    }
  }

  if (side == X1_END){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X2FACE){  /* -- y staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX2s][k][j][i] = Bfield[1]/unit_b;
      }

    #endif
    }else if (box->vpos == X3FACE){  /* -- z staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX3s][k][j][i] = Bfield[2]/unit_b;
      }

    #endif
    }

  }

  if (side == X2_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X1FACE){  /* -- x staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX1s][k][j][i] = Bfield[0]/unit_b;
      }

    #endif
    }else if (box->vpos == X3FACE){  /* -- z staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX3s][k][j][i] = Bfield[2]/unit_b;
      }

    #endif
    }

  }

  if (side == X2_END){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X1FACE){  /* -- x staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX1s][k][j][i] = Bfield[0]/unit_b;
      }

    #endif
    }else if (box->vpos == X3FACE){  /* -- z staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX3s][k][j][i] = Bfield[2]/unit_b;
      }

    #endif
    }

  }

  if (side == X3_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X1FACE){  /* -- x staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX1s][k][j][i] = Bfield[0]/unit_b;
      }

    #endif
    }else if (box->vpos == X2FACE){  /* -- y staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX2s][k][j][i] = Bfield[1]/unit_b;
      }

    #endif
    }

  }

  if (side == X3_END){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];

        ambient_field(xx1, xx2, xx3, Bfield);

        bx_zone = Bfield[0]/unit_b;
        by_zone = Bfield[1]/unit_b;
        bz_zone = Bfield[2]/unit_b;

        ambient_hydro(xx1, xx2, xx3, hydroV);

        rho_zone  = hydroV[0];
        pres_zone = hydroV[1];
        trac_trs  = hydroV[2];
        trac_clp  = hydroV[3];
        trac_ej   = hydroV[4];

        d->Vc[RHO][k][j][i]   = rho_zone;
        d->Vc[PRS][k][j][i]   = pres_zone;
        d->Vc[TRC][k][j][i]   = 0.0;
        d->Vc[TRC+1][k][j][i] = trac_trs;
        d->Vc[TRC+2][k][j][i] = trac_clp;
        d->Vc[TRC+3][k][j][i] = trac_ej;
        d->Vc[TRC+4][k][j][i] = 0.0;
        d->Vc[TRC+5][k][j][i] = 0.0;
        d->Vc[TRC+6][k][j][i] = 0.0;
        d->Vc[TRC+7][k][j][i] = 0.0;

        csm_abundance(xx1, xx2, xx3, abundance);

        d->Vc[TRC+8][k][j][i]  = abundance[0];   /* Ar36 */
        d->Vc[TRC+9][k][j][i]  = abundance[1];   /* C12 */
        d->Vc[TRC+10][k][j][i] = abundance[2];   /* Ca40 */
        d->Vc[TRC+11][k][j][i] = abundance[3];   /* Cr48 */
        d->Vc[TRC+12][k][j][i] = abundance[4];   /* Fe52 */
        d->Vc[TRC+13][k][j][i] = abundance[5];   /* Fe54 */
        d->Vc[TRC+14][k][j][i] = abundance[6];   /* H1 */
        d->Vc[TRC+15][k][j][i] = abundance[7];   /* He3 */
        d->Vc[TRC+16][k][j][i] = abundance[8];   /* He4 */
        d->Vc[TRC+17][k][j][i] = abundance[9];   /* Mg24 */
        d->Vc[TRC+18][k][j][i] = abundance[10];  /* N14 */
        d->Vc[TRC+19][k][j][i] = abundance[11];  /* Ne20 */
        d->Vc[TRC+20][k][j][i] = abundance[12];  /* neut */
        d->Vc[TRC+21][k][j][i] = abundance[13];  /* Ni56 */
        d->Vc[TRC+22][k][j][i] = abundance[14];  /* O16 */
        d->Vc[TRC+23][k][j][i] = abundance[15];  /* prot */
        d->Vc[TRC+24][k][j][i] = abundance[16];  /* S32 */
        d->Vc[TRC+25][k][j][i] = abundance[17];  /* Si28 */
        d->Vc[TRC+26][k][j][i] = abundance[18];  /* Ti44 */

        DIM_EXPAND(d->Vc[VX1][k][j][i] = 0.0; ,
                   d->Vc[VX2][k][j][i] = 0.0; ,
                   d->Vc[VX3][k][j][i] = 0.0;)

  #if PHYSICS == MHD || PHYSICS == RMHD

        DIM_EXPAND(d->Vc[BX1][k][j][i] = bx_zone; ,
                   d->Vc[BX2][k][j][i] = by_zone; ,
                   d->Vc[BX3][k][j][i] = bz_zone;)

  #endif

      }
    } else if (box->vpos == X1FACE){  /* -- x staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX1s][k][j][i] = Bfield[0]/unit_b;
      }

    #endif
    }else if (box->vpos == X2FACE){  /* -- y staggered field -- */
    #ifdef STAGGERED_MHD

      BOX_LOOP(box,k,j,i){
        xx1 = x1[i];
        xx2 = x2[j];
        xx3 = x3[k];
        ambient_field(xx1, xx2, xx3, Bfield);
        d->Vs[BX2s][k][j][i] = Bfield[1]/unit_b;
      }

    #endif
    }

  }

}


/*! ***************************************
 * USER FUNCTIONS
 * ****************************************
* */

double linear_interp(double x, double x0, double x1, double y0, double y1){
  double y;

  y = y0+(x-x0)*(y1-y0)/(x1-x0);

  return y;
}

// Binary search implementation for double*
int* nearest_indexes_by(double* array, int size, double value) {
  static int idxs[2] = {0, 0};
  int kmid;

  idxs[1] = size-1;

  while (idxs[0] != (idxs[1] - 1)){
    kmid = (idxs[0] + idxs[1])/2;
    if (value <= array[kmid]){
      idxs[1] = kmid;
    }else if (value > array[kmid]){
      idxs[0] = kmid;
    }
  }

  return idxs;
}

// Count values (total, rows, columns)
// TODO: skiplines
size_t* countValues(char* filename) {
    float tmp;
    char ch;
    static size_t vals[3] = {0,0,0};

    FILE* fp = fopen(filename, "r");
    while((ch = fgetc(fp))!=EOF && ch!='\n') {
        vals[0]+=fscanf(fp, "%f", &tmp);
    }
    vals[2] = vals[0];
    while((ch = fgetc(fp))!=EOF) {
        vals[0]+=fscanf(fp, "%f", &tmp);
    }
    vals[1] = vals[0]/vals[2];
    fclose(fp);

    return vals;
}


void ambient_field(double x1, double x2, double x3, double *Bfield) {

/*
 *
 * PURPOSE
 *    Calculate the initial ambient magnetic field
 *
 * ARGUMENTS
 *
 **************************************************************** */


  double vwind, r_rs, bref;

  double r0, omega;
  double rr1, rr2;
  double theta_w, phi_w;
  double b_r, b_phi;

  double dist    = sqrt(x1*x1 + x2*x2 + x3*x3);
  double dist_dp = sqrt(x1*x1 + x2*x2);

  double bx_zone = 0.0;
  double by_zone = 0.0;
  double bz_zone = 0.0;

  vwind  = g_inputParam[vel_w];
  r_rs   = g_inputParam[r_star];
  bref   = g_inputParam[B_ref];

  if (dist > 0){

    r0    = r_rs*6.96e10;
    omega = 2.*CONST_PI/(200000.*27.*24.*3600.);
//    omega = 2.*CONST_PI/(40.*27.*24.*3600.);

    rr1 = dist*UNIT_LENGTH;

    theta_w = acos(x3/dist);
    phi_w   = acos(x1/dist_dp);

    if (x2 < 0){
      phi_w = -phi_w;
    }

    b_r   =  bref*(r0/rr1)*(r0/rr1);
    b_phi = -bref*omega*r0/vwind*(r0/rr1)*sin(theta_w);

    Bfield[0] = (b_r*sin(theta_w)*cos(phi_w)-b_phi*sin(phi_w));
    Bfield[1] = (b_r*sin(theta_w)*sin(phi_w)+b_phi*cos(phi_w));
    Bfield[2] = (b_r*cos(theta_w));

  } else {

    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;

  }

}


void ambient_hydro(double x1, double x2, double x3, double *hydroV) {

/*
 *
 * PURPOSE
 *    Calculate the initial ambient hydro variables
 *
 * ARGUMENTS
 *
 **************************************************************** */

  double T_ism, mu;
  double dist, dist_dp;
  int nnb;

  double rho_zone, rho_flat, pres_zone, trac_trs, trac_clp, trac_ej;
  double xClump, yClump, zClump, dist_cl, nu, func_sm;
  double height, rho_disk, theta;
  double lagr_x1, lagr_x2, lagr_x3;

  static double n_clumps;
  static double rad_clp[1000], rho_clp[1000], trs_clp[1000], the_clp[1000], phi_clp[1000];
  static double ref_dist, M_dot, v_wind, rho_wind, r_int, rho_dsk, pow_dsk;
  static double m_rsg, v_rsg, AA, BB, CC;
  static double smooth_par;
  static double year_sec = 3.1536e+7;

  double rad_clp0, rho_clp0, trs_clp0, the_clp0, rr_pow;

  int seed;
  double rand_num;
  static int first_call = 1;

  FILE *clumps;

/* ----------------------------------------------------------
         Set units
   ---------------------------------------------------------- */

  mu = MU_CSM/2.;

  double unit_time     = UNIT_LENGTH/UNIT_VELOCITY;
  double unit_pressure = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
  double unit_temp     = MU_CSM*CONST_mp*unit_pressure/(2.0*UNIT_DENSITY*CONST_kB);

/* ---------------------------------------------------------- */

  if (first_call==1){
    first_call = 0;

    seed = 3000;            /* choose a seed value */
    srand(seed);            /* initialize the randomizer */

/* ---------------------- bkg wind -------------------------- */

    ref_dist  = g_inputParam[ref_w]/UNIT_LENGTH;
    M_dot     = g_inputParam[mdot_w]*CONST_Msun/year_sec;
    v_wind    = g_inputParam[vel_w];
    rho_wind  = M_dot/(4.*CONST_PI*(ref_dist*UNIT_LENGTH)*(ref_dist*UNIT_LENGTH)*v_wind);
    rho_wind /= UNIT_DENSITY;

/* ----------------------- disk ----------------------------- */

    r_int    = g_inputParam[rad_d]/UNIT_LENGTH;
    m_rsg    = 7.e-4*1.99e33/3.16e7;
    v_rsg    = 10.e5;
    pow_dsk  = 2.4;
    AA       = 1.0;
    BB       = 5.0;
    CC       = 1.05e7;

/* ---------------------- clumps ---------------------------- */

    n_clumps   = g_inputParam[nclumps];
    rad_clp0   = 0.002;                    // radius of clumps
    rho_clp0   = 10.0;                     // density contrast of clumps
    trs_clp0   = r_int;                    // average distance of clumps from the star
    the_clp0   = 30./180.*CONST_PI;        // MAX abs(latitude) of clumps
    smooth_par = 3.;                       // smoothing of the clumps

    clumps = fopen("clumps.txt", "w");

    for (nnb = 0; nnb < n_clumps; nnb++){

      rand_num = (double)rand()/( (double)RAND_MAX + (double)(1) );
      rad_clp[nnb] = rad_clp0*(1.0 + 0.3*(2.*rand_num-1.0));

      rand_num = (double)rand()/( (double)RAND_MAX + (double)(1) );
      rho_clp[nnb] = rho_clp0*( pow(10., (2.*rand_num-1.0) ) );

      rand_num = (double)rand()/( (double)RAND_MAX + (double)(1) );
      trs_clp[nnb] = trs_clp0*(1. + 6.*rand_num);

      rand_num = (double)rand()/( (double)RAND_MAX + (double)(1) );
      the_clp[nnb] = the_clp0*(2.*rand_num-1.0);

      rand_num = (double)rand()/( (double)RAND_MAX + (double)(1) );
      phi_clp[nnb] = 2.*CONST_PI*rand_num;

      fprintf(clumps, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", rad_clp[nnb], rho_clp[nnb], trs_clp[nnb], the_clp[nnb], phi_clp[nnb]);
    }

    fclose(clumps);

  }

/* ---------------------------------------------------------- */

//  T_ism = 1.e3;                       // ISM Temperature
  T_ism = 1.e5;                       // ISM Temperature

  dist    = sqrt(x1*x1 + x2*x2 + x3*x3);
  dist_dp = sqrt(x1*x1 + x2*x2);

  rho_flat = 0.1*MU_CSM*CONST_mp/UNIT_DENSITY;
  rho_zone = rho_wind*(ref_dist/dist)*(ref_dist/dist) + rho_flat;

  pres_zone  = 2.0*(rho_zone*UNIT_DENSITY)/(MU_CSM*CONST_mp)*CONST_kB*T_ism;
  pres_zone /= unit_pressure;

/* ---------------------- Add the disk ---------------------- */

  trac_trs = 0.0;

  if (dist > r_int) {

    theta     = atan(x1/x3);
    rr_pow    = pow( (dist*UNIT_LENGTH), pow_dsk);
    rho_disk  = CC*m_rsg/(4.*v_rsg*CONST_PI*rr_pow)*
                (1.0-AA*(exp(-2.*BB*cos(theta)*cos(theta))-1.0)/(exp(-2.*BB)-1.0));
    rho_disk  = rho_disk/UNIT_DENSITY;

    rho_zone = MAX(rho_zone, rho_disk);
    trac_trs = 1.0;

  }

/* --------------------- Add the clumps --------------------- */

  trac_clp = 0.0;

  for (nnb = 0; nnb < n_clumps; nnb++){

    xClump  = trs_clp[nnb]*cos(the_clp[nnb])*cos(phi_clp[nnb]);
    yClump  = trs_clp[nnb]*cos(the_clp[nnb])*sin(phi_clp[nnb]);
    zClump  = trs_clp[nnb]*sin(the_clp[nnb]);

    dist_cl = sqrt( (x1-xClump)*(x1-xClump) +
                    (x2-yClump)*(x2-yClump) +
                    (x3-zClump)*(x3-zClump) );

    nu      = 1.0/rho_clp[nnb];

    if (dist_cl <= rad_clp[nnb]){
       func_sm  = (nu-(nu-1.0)/cosh(smooth_par*pow( (dist_cl/rad_clp[nnb]), smooth_par) ) );
       rho_zone = rho_zone*rho_clp[nnb]*func_sm;
       trac_clp = 1.0;
    }
  }

//  pres_zone = rho_wind*(T_ism/unit_temp);

  lagr_x1 = 0.0;
  lagr_x2 = 0.0;
  lagr_x3 = 0.0;

  trac_ej = 0.0;

/* --------------------------------------------- */

  hydroV[0] = rho_zone;
  hydroV[1] = pres_zone;
  hydroV[2] = trac_trs;
  hydroV[3] = trac_clp;
  hydroV[4] = trac_ej;
  hydroV[5] = lagr_x1;
  hydroV[6] = lagr_x2;
  hydroV[7] = lagr_x3;

}

void csm_abundance(double x1, double x2, double x3, double *abundance) {

/*
 *
 * PURPOSE
 * Calculate the abundance of CSM
 *
 * ARGUMENTS
 *
 ***************************************************************** */


  abundance[0]  = 0.00000E+00;   /* Ar36 */
  abundance[1]  = 0.00000E+00;   /* C12 */
  abundance[2]  = 0.00000E+00;   /* Ca40 */
  abundance[3]  = 0.00000E+00;   /* Cr48 */
  abundance[4]  = 0.00000E+00;   /* Fe52 */
  abundance[5]  = 0.00000E+00;   /* Fe54 */
  abundance[6]  = 0.00000E+00;   /* H1 */
  abundance[7]  = 0.00000E+00;   /* He3 */
  abundance[8]  = 0.00000E+00;   /* He4 */
  abundance[9]  = 0.00000E+00;   /* Mg24 */
  abundance[10] = 0.00000E+00;   /* N14 */
  abundance[11] = 0.00000E+00;   /* Ne20 */
  abundance[12] = 0.00000E+00;   /* neut */
  abundance[13] = 0.00000E+00;   /* Ni56 */
  abundance[14] = 0.00000E+00;   /* O16 */
  abundance[15] = 0.00000E+00;   /* prot */
  abundance[16] = 0.00000E+00;   /* S32 */
  abundance[17] = 0.00000E+00;   /* Si28 */
  abundance[18] = 0.00000E+00;   /* Ti44 */

}

