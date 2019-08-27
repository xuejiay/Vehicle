//
// Created by Xuejiao Yang on 4/25/19.
//

#ifndef EXAMPLE1_MAINHEADER_H
#define EXAMPLE1_MAINHEADER_H

#include "includInfo.h"
int SDI(double* X0, double* Uncertainty, memberInfoSDI Input);
int DI_constraints(double* X0, double* Uncertainty, memberInfoDI_constraints Input);
int DI_constraints_advance(double* X0, double* Uncertainty, memberInfoDI_cst_advance Input);
int trajectories(double* X0, double* Uncertainty, memberInfoSDI Input);

#endif //EXAMPLE1_MAINHEADER_H
