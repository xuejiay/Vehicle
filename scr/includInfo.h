//
// Created by Xuejiao Yang on 4/23/19.
//

#ifndef TEST_INCLUDINFO_H
#define TEST_INCLUDINFO_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include "privateData.h"
#include "interval.hpp"

typedef mc::Interval IA;


double rhs(double t, double *x_vector, double *p, double *ref, int i);

void refinement(IA *x_vector, IA *p);

void refinement_advance(IA *x_vector, IA *p);

IA rhsI(double t, IA *x_vector, IA *p, IA *ref, int i);


struct memberInfoSDI {
    double X0[2 * num_state_variable];
    double Uncertainty[2 * num_uncertain_parameter];
    double InputReference[num_reference_input * num_piecewise];
};

struct memberInfoDI_constraints {
    double X0[2 * num_state_redundant];
    double Uncertainty[2 * num_uncertain_parameter];
    double InputReference[num_reference_input * num_piecewise];
    double NaturBD[2 * num_state_redundant];
};

struct memberInfoDI_cst_advance {
    double X0[2 * num_state_redundant_advance];
    double Uncertainty[2 * num_uncertain_parameter];
    double InputReference[num_reference_input * num_piecewise];
    double NaturBD[2 * num_state_redundant_advance];
};

#endif //TEST_INCLUDINFO_H
