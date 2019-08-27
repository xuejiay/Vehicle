//
// Created by Xuejiao Yang on 4/23/19.
//


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
using namespace std;

static void extendIntersection(IA *origint, IA refint);

static IA SQRinterval(IA SQinterval);

double rhs(double t, double *x_vector, double *p, double *ref, int i) {
    double fx;

    double wref = ref[0], vref = ref[1], xe = x_vector[0], ye = x_vector[1], theta_e = x_vector[2], d1 = p[0],
            d2 = p[1], k1 = 10, k2 = 6.4 * 0.001, k3 = 0.161;
    double w = wref + vref * (k2 * ye + k3 * sin(theta_e)) + d1;
    switch (i) {
        case 0:
            fx = w * ye - k1 * xe - d2;
            break;
        case 1:
            fx = -w * xe + vref * sin(theta_e);
            break;
        case 2:
            fx = wref - w;
            break;
    }


    return fx;

}


IA rhsI(double t, IA *x_vector, IA *p, IA *ref, int i) {
    IA fx;
    IA wref = ref[0], vref = ref[1], xe = x_vector[0], ye = x_vector[1], theta_e = x_vector[2], d1 = p[0],
            d2 = p[1], k1 = 10, k2 = 6.4 * 0.001, k3 = 0.161;

    //IA inter, z1 = x_vector[3], z2 = x_vector[4];
    IA w = wref + vref * (k2 * ye + k3 * sin(theta_e)) + d1;
    switch (i) {
        case 0: {
            fx = w * ye - k1 * xe - d2;
            break;
        }
        case 1: {
            fx = -w * xe + vref * sin(theta_e);
            break;
        }
        case 2: {
            fx = -vref * (k2 * ye + k3 * sin(theta_e)) - d1;
            break;
        }

        case 3: {
            IA xeSQ = p[2], sinSQ = p[3];
            //fx=-k1*pow(xe,2)-xe*d2-(vref*k3*pow(sin(theta_e),2))/k2+sin(theta_e)/k2*d1;
            fx = -k1 * xeSQ - xe * d2 - (vref * k3 * sinSQ) / k2 + sin(theta_e) / k2 * d1;
            //fx = 2*g1;
            break;
        }
        case 4: {
            IA v = vref * cos(theta_e) + k1 * xe + d2;
            fx = v * cos(x_vector[6]);
            break;
        }

        case 5: {
            IA v = vref * cos(theta_e) + k1 * xe + d2;
            fx = v * sin(x_vector[6]);
            break;
        }

        case 6: {
            IA w = wref + vref * (k2 * ye + k3 * sin(theta_e)) + d1;
            fx = w;
            break;
        }
        case 7: {
            fx = vref * cos(x_vector[9]);
            break;
        }

        case 8: {
            fx = vref * sin(x_vector[9]);
            break;
        }


        case 9: {
            fx = wref;
            break;
        }
    }


    return fx;

}


void refinement(IA x_vector[num_state_redundant], IA *p) {

    for (int l = 0; l < num_loop; l++) {
        // V = 1/2*（xe^2+ye^2）+(1-cos(theta_e))/k2
        IA k2 = 6.4 * 0.001, invers_V, xeSQR, yeSQR, theta_eCOS, invers_xeSQR, invers_yeSQR, invers_thetaeCOS;

        xeSQR = pow(x_vector[0], 2);
        yeSQR = pow(x_vector[1], 2);
        theta_eCOS = cos(x_vector[2]);


        // refine lyapunov function V
        invers_V = (xeSQR + yeSQR) / 2. + (1. - theta_eCOS) / k2;
        // Taking intersection with the original interval
        extendIntersection(&x_vector[3], invers_V);


        // refine xe^2
        invers_xeSQR = 2 * (x_vector[3] - (1. - theta_eCOS) / k2 - yeSQR / 2.);
        extendIntersection(&xeSQR, invers_xeSQR);
        p[2] = xeSQR;
        // refine xe by taking squreRoot
        IA xe_inverse = SQRinterval(xeSQR);
        extendIntersection(&x_vector[0], xe_inverse);

        // refine ye^2
        invers_yeSQR = 2 * (x_vector[3] - (1. - theta_eCOS) / k2 - xeSQR / 2.);
        extendIntersection(&yeSQR, invers_yeSQR);
        // refine xe by taking squreRoot
        IA ye_inverse = SQRinterval(yeSQR);
        extendIntersection(&x_vector[1], ye_inverse);


        // refine cos(theta_e)
        invers_thetaeCOS = 1. - (x_vector[3] - xeSQR / 2. - yeSQR / 2.) * k2;
        extendIntersection(&theta_eCOS, invers_thetaeCOS);
        p[3] = 1 - pow(theta_eCOS, 2);
        // refine theta_e
        IA invers_thetae = hull(acos(theta_eCOS), -acos(theta_eCOS));
        extendIntersection(&x_vector[2], invers_thetae);

    }


}


void refinement_advance(IA x_vector[num_state_redundant_advance], IA *p) {
    for (int l = 0; l < num_loop; l++) {
        // V = 1/2*（xe^2+ye^2）+(1-cos(theta_e))/k2
        IA k2 = 6.4 * 0.001, invers_V, xeSQR, yeSQR, theta_eCOS, invers_xeSQR, invers_yeSQR, invers_thetaeCOS;

        xeSQR = pow(x_vector[0], 2);
        yeSQR = pow(x_vector[1], 2);
        theta_eCOS = cos(x_vector[2]);


        // refine lyapunov function V
        invers_V = (xeSQR + yeSQR) / 2. + (1. - theta_eCOS) / k2;
        // Taking intersection with the original interval
        extendIntersection(&x_vector[3], invers_V);


        // refine xe^2
        invers_xeSQR = 2 * (x_vector[3] - (1. - theta_eCOS) / k2 - yeSQR / 2.);
        extendIntersection(&xeSQR, invers_xeSQR);
        p[2] = xeSQR;
        // refine xe by taking squreRoot
        IA xe_inverse = SQRinterval(xeSQR);
        extendIntersection(&x_vector[0], xe_inverse);

        // refine ye^2
        invers_yeSQR = 2 * (x_vector[3] - (1. - theta_eCOS) / k2 - xeSQR / 2.);
        extendIntersection(&yeSQR, invers_yeSQR);
        // refine xe by taking squreRoot
        IA ye_inverse = SQRinterval(yeSQR);
        extendIntersection(&x_vector[1], ye_inverse);


        // refine cos(theta_e)
        invers_thetaeCOS = 1. - (x_vector[3] - xeSQR / 2. - yeSQR / 2.) * k2;
        extendIntersection(&theta_eCOS, invers_thetaeCOS);
        p[3] = 1 - pow(theta_eCOS, 2);
        // refine theta_e
        IA invers_thetae = hull(acos(theta_eCOS), -acos(theta_eCOS));
        extendIntersection(&x_vector[2], invers_thetae);


        //xe = Y[0] ye = Y[1] thetae = Y[2] yr=Y[8] xr = Y[7] thetar = Y[9] yc = Y[5] xc=Y[4] thetac=Y[6]
        // Use errors to tight the original positions

        IA theta_center = x_vector[9] - x_vector[2];
        extendIntersection(&x_vector[6], theta_center);

        IA y_center = x_vector[8] - sin(x_vector[6]) * x_vector[0] - cos(x_vector[6]) * x_vector[1];
        extendIntersection(&x_vector[5], y_center);

        IA x_center = x_vector[7] - cos(x_vector[6]) * x_vector[0] + sin(x_vector[6]) * x_vector[1];
        extendIntersection(&x_vector[4], x_center);


        // use the original positions to refine errors
        xe_inverse = cos(x_vector[6]) * (x_vector[7] - x_vector[4]) + sin(x_vector[6]) * (x_vector[8] - x_vector[5]);
        extendIntersection(&x_vector[0], xe_inverse);

        ye_inverse = sin(x_vector[6]) * (x_vector[7] - x_vector[4]) + cos(x_vector[6]) * (x_vector[8] - x_vector[5]);
        extendIntersection(&x_vector[1], ye_inverse);

        invers_thetae = x_vector[9] - x_vector[6];
        extendIntersection(&x_vector[2], invers_thetae);
    }


}


static void extendIntersection(IA *origint, IA refint) {
    origint->l(min(max(origint->l(), refint.l()), origint->u()));
    origint->u(max(min(origint->u(), refint.u()), origint->l()));

}


static IA SQRinterval(IA SQinterval) {

    IA SQR_invers;
    if (SQinterval.u() > 0.0001) {

        SQR_invers = hull(sqrt(SQinterval), -sqrt(SQinterval));
    } else {
        SQR_invers.u(50. * SQinterval.u() + (0.01 - 5. * pow(10, -3)));
        SQR_invers.l(-50. * SQinterval.u() + (-0.01 + 5. * pow(10, -3)));
    }
    return SQR_invers;
}
