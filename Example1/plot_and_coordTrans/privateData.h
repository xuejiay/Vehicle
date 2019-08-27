//
// Created by Xuejiao Yang on 4/10/19.
//

#ifndef TEST_PRIVATEDATA_H
#define TEST_PRIVATEDATA_H


/* Dimensions number*/


/* Tolerance for the solver*/
#define RTOL     RCONST(1.0e-6)   /* scalar relative tolerance            */
#define ATOL     RCONST(1.0e-6)    // scalar absolute tolerance for original system
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */

#define File_Output              1                                      //whether output data of sampling


/* Problem Time output option */
#define T0    RCONST(0.0)      /* initial time           */
#define Tf    RCONST(0.005)      /* first output time      */
#define TMULT RCONST(0.1)     /* output time factor     */
#define NOUT  801//749//850//840//574//765              /* number of output times */




#endif //TEST_PRIVATEDATA_H
