#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   DIV_CLEANING

/* -- user-defined parameters (labels) -- */

#define  GAMMA_0                        0
#define  SIGMA_0                        1
#define  WL_0                           2
#define  TLIM_SWITCH                    3
#define  RPW                            4
#define  RRES                           5

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOCK_FLATTENING               MULTID
#define  LIMITER                        MC_LIM
#define  RECONSTRUCT_4VEL               YES
#define  RMHD_FAST_EIGENVALUES          YES
#define  STATIC_REFINEMENT              NO

/* [End] user-defined constants (do not change this line) */
