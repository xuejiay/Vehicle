//
// Created by Xuejiao Yang on 4/10/19.
//

#ifndef TEST_PRIVATEDATA_H
#define TEST_PRIVATEDATA_H


/* Dimensions number*/
#define num_state_variable              7     // number of state variables x in original ODE systems
#define num_state_redundant             3
#define num_state_redundant_advance     5
#define num_reference_input             2
#define num_uncertain_parameter         0
#define samples                         4
#define num_loop                        2
#define num_uncertainties               0
#define num_piecewise                   1
/* Tolerance for the solver*/
#define RTOL     RCONST(1.0e-6)   /* scalar relative tolerance            */
#define ATOL     RCONST(1.0e-6)    // scalar absolute tolerance for original system
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */

#define File_Output              0                                      //whether output data of sampling
#define File_Output_SDI          1                                     //whether output data of SDI bd
#define File_Output_DI           0                                      //whether output data of DI with constraints bd
#define File_Output_DI_advance   0                                      //whether output data of DI with constraints and z bd



/* Problem Time output option */
#define T0    RCONST(0.0)      /* initial time           */
#define Tf    RCONST(0.005)      /* first output time      */
#define TMULT RCONST(0.1)     /* output time factor     */
#define NOUT  4//749//850//840//574//765              /* number of output times */




#endif //TEST_PRIVATEDATA_H
