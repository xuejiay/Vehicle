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

    double k1 = 2.0, k2 = 3.0, k3 = 1.0, k4 = 10.;

    double et = x_vector[0], en = x_vector[1], efi = x_vector[2],
            v = x_vector[3], kdelta = x_vector[4], vd = x_vector[5], wref = x_vector[6],
            aref = ref[0], alphad = ref[1], kd = wref / vd;


    switch (i) {
        case 0: {
            fx = v * cos(efi) - vd * (1 - kd * en);
            break;
        }
        case 1: {
            fx = v * sin(efi) - vd * kd * et;
            break;
        }

        case 2: {
            fx = v * kdelta - vd * kd;
            break;
        }
        case 3: {
            double ev = v - vd;
            fx = aref - k1 * et - k3 * ev + k2 * pow(efi, 2) - efi * kd;
            break;
        }
        case 4: {
            double epsilon = kd - k1 * (et * (cos(efi) - 1) / efi + en * (sin(efi) / efi)) - k2 * efi,
                    kd_dot = pow((alphad * vd - aref * wref) / vd, 2.), et_dot = v * cos(efi) - vd * (1 - kd * en),
                    en_dot = v * sin(efi) - vd * kd * et, efi_dot = v * kdelta - vd * kd,
                    epsilon_dot = kd_dot - k1 * (et_dot * (cos(efi) - 1) / efi +
                                                 et * efi_dot * (-sin(efi) / efi - (cos(efi) - 1) / pow(efi, 2)) +
                                                 en_dot * sin(efi) / efi +
                                                 en * efi_dot * (efi * cos(efi) - sin(efi)) / pow(efi, 2)) -
                                  k2 * efi_dot,
                    edelta = kdelta - epsilon;
            fx = -efi * v + epsilon_dot - k4 * edelta;
            break;
        }
        case 5: {
            fx = aref;
            break;
        }
        case 6: {
            fx = alphad;
            break;
        }

    }


    return fx;

}


IA rhsI(double t, IA *x_vector, IA *p, IA *ref, int i) {
    IA fx;
    IA k1 = 2.0, k2 = 3.0, k3 = 1.0, k4 = 10.;
    IA et = x_vector[0], en = x_vector[1], efi = x_vector[2],
            v = x_vector[3], kdelta = x_vector[4], vd = x_vector[5], wref = x_vector[6],
            aref = ref[0], alphad = ref[1], kd = wref / vd;


    switch (i) {
        case 0: {
            fx = v * cos(efi) - vd * (1 - kd * en);
            break;
        }

        case 1: {
            fx = v * sin(efi) - vd * kd * et;
            break;
        }

        case 2: {
            fx = v * kdelta - vd * kd;
            break;
        }

        case 3: {
            IA ev = v - vd;
            fx = aref - k1 * et - k3 * ev + k2 * pow(efi, 2) - efi * kd;
            break;
        }

        case 4: {
            IA kd_dot = pow((alphad * vd - aref * wref) / vd, 2.), et_dot = v * cos(efi) - vd * (1 - kd * en),
                    en_dot = v * sin(efi) - vd * kd * et, efi_dot = v * kdelta - vd * kd;
            double efiL = efi.l();
            double efiU = efi.u();
            IA IntFunc1, IntFunc4;
            IntFunc1.l((cos(efiU) - 1) / efiU);
            IntFunc1.u((cos(efiL) - 1) / efiL);
            IntFunc4.l(pow((efiU * cos(efiU) - sin(efiU)) / efiU, 2));
            IntFunc4.l(pow((efiL * cos(efiL) - sin(efiL)) / efiL, 2));

            double U = max(abs(efiU), abs(efiL)), L = min(abs(efiL), abs(efiL));
            IA IntFunc2, IntFunc3;

            if (efiU * efiL <= 0) {
                IntFunc2 = IA(sin(U) / U, 1),
                        IntFunc3 = IA(-0.5, (cos(U) - 1) / pow(U, 2));
            } else {
                IntFunc2 = IA(sin(U) / U, sin(L) / L),
                        IntFunc3 = IA((cos(L) - 1) / pow(L, 2), (cos(U) - 1) / pow(U, 2));
            }
            if (efiU == 0. && efiL == 0.) {
                IntFunc2 = IA(1.), IntFunc3 = IA(-0.5);
            }


            IA epsilon = kd - k1 * (et * IntFunc1 + en * IntFunc2) - k2 * efi,
                    epsilon_dot = kd_dot - k1 * (et_dot * IntFunc1 - et * efi_dot * (IntFunc2 + IntFunc3) +
                                                 en_dot * IntFunc2 + en * efi_dot * IntFunc4) - k2 * efi_dot,
                    edelta = kdelta - epsilon;
            fx = -efi * v + epsilon_dot - k4 * edelta;
            break;
        }
        case 5: {
            fx = aref;
            break;
        }
        case 6: {
            fx = alphad;
            break;
        }

    }


    return fx;

}


void refinement(IA x_vector[num_state_redundant], IA *p) {
    //natural bound here
    for (int l = 0; l < num_loop; l++) {
        IA //e=x_vector[0], theta_e = x_vector[1], V = x_vector[2],
                a = 2, g2 = pow(a,
                                2), invers_V, eSQR, theta_eSQR, invers_eSQR, invers_thetaeSQR, invers_e, invers_thetae;

        eSQR = pow(x_vector[0], 2);
        theta_eSQR = pow(x_vector[1], 2);


        // refine lyapunov function V
        invers_V = (eSQR + theta_eSQR / g2) / 2;
        // Taking intersection with the original interval
        extendIntersection(&x_vector[2], invers_V);


        // refine e^2
        invers_eSQR = 2 * x_vector[2] - theta_eSQR / g2;
        extendIntersection(&eSQR, invers_eSQR);
        // refine e
        if (eSQR.u() > 0.0001) {

            invers_e = hull(sqrt(eSQR), -sqrt(eSQR));
        } else {
            invers_e.u(50. * eSQR.u() + (0.01 - 5. * pow(10., -3)));
            invers_e.l(-50. * eSQR.u() + (-0.01 + 5. * pow(10., -3)));
        }
        extendIntersection(&x_vector[0], invers_e);


        // refine theta_e^2
        invers_thetaeSQR = (2 * x_vector[2] - eSQR) * g2;
        extendIntersection(&theta_eSQR, invers_thetaeSQR);
        p[1] = theta_eSQR;
        // refine e
        if (theta_eSQR.u() > 0.0001) {

            invers_thetae = hull(sqrt(theta_eSQR), -sqrt(theta_eSQR));
        } else {
            invers_thetae.u(50. * theta_eSQR.u() + (0.01 - 5. * pow(10, -3)));
            invers_thetae.l(-50. * theta_eSQR.u() + (-0.01 + 5. * pow(10, -3)));
        }
        extendIntersection(&x_vector[1], invers_thetae);

    }


}


void refinement_advance(IA x_vector[num_state_redundant_advance], IA *p) {
    //natural bound here
    for (int l = 0; l < num_loop; l++) {
        IA //e=x_vector[0], theta_e = x_vector[1], V = x_vector[2],
                a = 2, g2 = pow(a,
                                2), invers_V, eSQR, theta_eSQR, invers_eSQR, invers_thetaeSQR, invers_e, invers_thetae;

        eSQR = pow(x_vector[0], 2);
        theta_eSQR = pow(x_vector[1], 2);


        // refine lyapunov function V
        invers_V = (eSQR + theta_eSQR / g2) / 2;
        // Taking intersection with the original interval
        extendIntersection(&x_vector[2], invers_V);


        // refine e^2
        invers_eSQR = 2 * x_vector[2] - theta_eSQR / g2;
        extendIntersection(&eSQR, invers_eSQR);
        // refine e
        if (eSQR.u() > 0.0001) {

            invers_e = hull(sqrt(eSQR), -sqrt(eSQR));
        } else {
            invers_e.u(50. * eSQR.u() + (0.01 - 5. * pow(10., -3)));
            invers_e.l(-50. * eSQR.u() + (-0.01 + 5. * pow(10., -3)));
        }
        extendIntersection(&x_vector[0], invers_e);


        // refine theta_e^2
        invers_thetaeSQR = (2 * x_vector[2] - eSQR) * g2;
        extendIntersection(&theta_eSQR, invers_thetaeSQR);
        p[1] = theta_eSQR;
        // refine e
        if (theta_eSQR.u() > 0.0001) {

            invers_thetae = hull(sqrt(theta_eSQR), -sqrt(theta_eSQR));
        } else {
            invers_thetae.u(50. * theta_eSQR.u() + (0.01 - 5. * pow(10, -3)));
            invers_thetae.l(-50. * theta_eSQR.u() + (-0.01 + 5. * pow(10, -3)));
        }
        extendIntersection(&x_vector[1], invers_thetae);


        //z1 = cos(theta)
        //z2 = tan(theta)


        IA invers_z1, invers_z2, z1SQR, z1SQR_invers, z2SQR, z2SQR_invers;

        invers_z1 = cos(x_vector[1]);
        invers_z2 = tan(x_vector[1]);
        extendIntersection(&x_vector[3], invers_z1);
        extendIntersection(&x_vector[4], invers_z2);



        //z2SQR + 1= 1/z1SQR
        z1SQR = pow(x_vector[3], 2);
        z2SQR = pow(x_vector[4], 2);

        z1SQR_invers = 1. / (z2SQR + 1.);
        extendIntersection(&z1SQR, z1SQR_invers);
        p[2] = z1SQR;

        z2SQR_invers = 1. / z1SQR - 1.;
        extendIntersection(&z2SQR, z2SQR_invers);


        //square root of z1 and z2
        IA z1_invers = SQRinterval(z1SQR);
        extendIntersection(&x_vector[3], z1_invers);

        IA z2_invers = SQRinterval(z2SQR);
        extendIntersection(&x_vector[4], z2_invers);


        //use z1 and z2 to refine for theta_e
        //z1 = cos(theta)
        //z2 = tan(theta)
        // avoid beyond -1,1

        x_vector[3].l(min(max(-1., x_vector[3].l()), 1.));
        x_vector[3].u(max(min(1., x_vector[3].u()), -1.));

        //x_vector[4].l(min(max(-1.,x_vector[4].l()),1.));
        //x_vector[4].u(max(min(1.,x_vector[4].u()),-1.));

        invers_thetae = atan(x_vector[4]);
        extendIntersection(&x_vector[1], invers_thetae);


        invers_thetae = hull(acos(x_vector[3]), -acos(x_vector[3]));
        extendIntersection(&x_vector[1], invers_thetae);


        //use z1 and z2 to refine for theta_e
        // z2*z1 = sin(theta)
        IA Z1Z2 = sin(x_vector[1]);
        //IA Z1Z2_invers=x_vector[4]*x_vector[3];
        extendIntersection(&Z1Z2, x_vector[4] * x_vector[3]);

        invers_thetae = asin(Z1Z2);
        extendIntersection(&x_vector[1], invers_thetae);

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
