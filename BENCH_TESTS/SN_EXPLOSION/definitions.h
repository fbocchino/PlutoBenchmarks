#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        27
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            14

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 SELECTIVE
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  EXT                            0
#define  AMPL                           1
#define  s_dom                          2
#define  initime                        3
#define  ref_w                          4
#define  mdot_w                         5
#define  vel_w                          6
#define  rad_d                          7
#define  size_d                         8
#define  rho_d                          9
#define  B_ref                          10
#define  r_star                         11
#define  nclumps                        12
#define  BOUND_SH                       13

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  MULTIPLE                       LOG
#define  WARNING_MESSAGES               NO
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MC_LIM
#define  UNIT_LENGTH                    3.09e18
#define  UNIT_DENSITY                   2.15253e-16
#define  UNIT_VELOCITY                  1.000E+8
#define  MU_CSM                         1.2889425

/* [End] user-defined constants (do not change this line) */
