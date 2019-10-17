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
                    kd_dot = pow((alphad * vd - aref * wref) / vd, 2), et_dot = v * cos(efi) - vd * (1 - kd * en),
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
    IA k1 = 2.0, k2 = 3.0, k3 = 1.0, k4 = 10., l = 2.0;
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
            IA kd_dot = pow((alphad * vd - aref * wref) / vd, 2), et_dot = v * cos(efi) - vd * (1 - kd * en),
                    en_dot = v * sin(efi) - vd * kd * et, efi_dot = v * kdelta - vd * kd;
            double efiL = efi.l();
            double efiU = efi.u();
            IA IntFunc1, IntFunc4;
            IntFunc1.l((cos(efiU) - 1) / efiU);
            IntFunc1.u((cos(efiL) - 1) / efiL);
            IntFunc4.l((efiU * cos(efiU) - sin(efiU)) / pow(efiU, 2));
            IntFunc4.l((efiL * cos(efiL) - sin(efiL)) / pow(efiL, 2));

            double U = max(abs(efiU), abs(efiL)), L = min(abs(efiL), abs(efiL));
            IA IntFunc2, IntFunc3;

            if (efiU * efiL <= 0) {
                IntFunc2 = IA(sin(U) / U, 1);
                IntFunc3 = IA(-0.5, (cos(U) - 1) / pow(U, 2));
            } else {
                IntFunc2 = IA(sin(U) / U, sin(L) / L);
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
        case 7: {
            IA ev = v-vd;
            fx = -vd*k2*pow(efi,2)-k3* pow(ev,2);
            break;
        }
        case 8: {
            IA ev = v-vd;
            fx = -vd*k2*pow(efi,2)-k3*pow(ev,2)-k4* pow(x_vector[9],2);
            break;
        }
        case 9: {
            fx = -efi*v-k4*x_vector[9];
            break;
        }
        case 10:{
            fx = v * cos(x_vector[12]);
            break;
        }
        case 11: {
            fx = v * sin(x_vector[12]);
            break;
        }
        case 12:{

            fx = v * tan(x_vector[13])/l;
            break;
        }
        case 13:{
            IA  kd_dot = pow((alphad * vd - aref * wref) / vd, 2), et_dot = v * cos(efi) - vd * (1 - kd * en),
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
                IntFunc2 = IA(sin(U) / U, 1);
                IntFunc3 = IA(-0.5, (cos(U) - 1) / pow(U, 2));
            } else {
                IntFunc2 = IA(sin(U) / U, sin(L) / L);
                IntFunc3 = IA((cos(L) - 1) / pow(L, 2), (cos(U) - 1) / pow(U, 2));
            }
            if (efiU == 0. && efiL == 0.) {
                IntFunc2 = IA(1.), IntFunc3 = IA(-0.5);
            }

            IA epsilon = kd - k1 * (et * IntFunc1 + en * IntFunc2) - k2 * efi,
                    epsilon_dot = kd_dot - k1 * (et_dot * IntFunc1 - et * efi_dot * (IntFunc2 + IntFunc3) +
                                                 en_dot * IntFunc2 + en * efi_dot * IntFunc4) - k2 * efi_dot,
                    edelta = kdelta - epsilon;
            IA w1=-efi*v+epsilon_dot-k4*edelta;
            fx = w1/(1.0/l+l*pow(kdelta,2));
            break;
        }
        case 14:{
            fx = vd*cos(x_vector[16]);
            break;
        }
        case 15:{
            fx = vd*sin(x_vector[16]);
            break;
        }
        case 16:{
            fx = wref;
            break;
        }

    }


    return fx;

}


void refinement(IA x_vector[num_state_redundant], IA *p) {
    IA vd = x_vector[5], invers_V, invers_edeltaSq, edeltaSq, invers_edelta;
    IA k1 = 2.0, k2 = 3.0, invers_etSq, etSq, invers_et, invers_enSq, enSq, invers_en,
    invers_efiSq, efiSq, invers_efi, invers_evSq, evSq, invers_ev;
    edeltaSq = pow(x_vector[9],2);
    etSq=pow(x_vector[0],2);
    enSq=pow(x_vector[1],2);
    efiSq=pow(x_vector[2],2);
    evSq=pow((x_vector[3]-vd),2);
    for (int l = 0; l < num_loop; l++) {


        //refine V and square of edelta using Vc

        invers_V = x_vector[8]-0.5*edeltaSq;
        extendIntersection(&x_vector[7], invers_V);

        invers_edeltaSq = 2*(x_vector[8]-x_vector[7]);
        extendIntersection(&edeltaSq, invers_edeltaSq);

        invers_edelta = SQRinterval(edeltaSq);
        extendIntersection(&x_vector[9], invers_edelta);


        // Use V to refine the other errors
        //et
        invers_etSq = (x_vector[7]*2-k1*enSq-efiSq-evSq)/k1;
        extendIntersection(&etSq, invers_etSq);
        invers_et = SQRinterval(etSq);
        extendIntersection(&x_vector[0], invers_et);

        //en
        invers_enSq = (x_vector[7]*2-k1*etSq-efiSq-evSq)/k1;
        extendIntersection(&enSq, invers_enSq);
        invers_en = SQRinterval(enSq);
        extendIntersection(&x_vector[1], invers_en);

        //efi
        invers_efiSq = x_vector[7]*2-k1*etSq-k1*enSq-evSq;
        extendIntersection(&efiSq, invers_efiSq);
        invers_efi = SQRinterval(efiSq);
        extendIntersection(&x_vector[2], invers_efi);

        //ev
        invers_evSq = x_vector[7]*2-k1*etSq-k1*enSq-efiSq;
        extendIntersection(&evSq, invers_evSq);
        invers_ev = SQRinterval(evSq);
        extendIntersection(&x_vector[3], invers_ev+vd);

        //edelta = kdelta - epsilon
        IA epsilon, invers_kdelta, IntFunc1, IntFunc2;

        IntFunc1 = IA((cos(x_vector[2].u())-1)/x_vector[2].u(),
                             (cos(x_vector[2].l())-1)/x_vector[2].l());
        double U = max(abs(x_vector[2].u()), abs(x_vector[2].l()));
        double L = min(abs(x_vector[2].u()), abs(x_vector[2].l()));
        if (x_vector[2].u()* x_vector[2].l()<=0){
            IntFunc2 = IA(sin(U)/U, 1.);
        }
        else {
            IntFunc2 = IA(sin(U)/U, sin(L)/L);
        }
        if (x_vector[2].u() == 0. && x_vector[2].l() ==0){
            IntFunc2 = IA(1.);
        }

        epsilon = x_vector[6]/vd-k1*(x_vector[0]*IntFunc1+x_vector[1]*IntFunc2)-k2*x_vector[2];

        invers_kdelta = x_vector[9] + epsilon;
        extendIntersection(&x_vector[4], invers_kdelta);

        invers_edelta = x_vector[4] - epsilon;
        extendIntersection(&x_vector[9], invers_edelta);

    }


}


void refinement_advance(IA x_vector[num_state_redundant_advance], IA *p) {
    IA vd = x_vector[5], invers_V, invers_edeltaSq, edeltaSq, invers_edelta;
    IA k1 = 2.0, k2 = 3.0, invers_etSq, etSq, invers_et, invers_enSq, enSq, invers_en,
            invers_efiSq, efiSq, invers_efi, invers_evSq, evSq, invers_ev;
    edeltaSq = pow(x_vector[9],2);
    etSq=pow(x_vector[0],2);
    enSq=pow(x_vector[1],2);
    efiSq=pow(x_vector[2],2);
    evSq=pow((x_vector[3]-vd),2);
    for (int l = 0; l < num_loop; l++) {


        //refine V and square of edelta using Vc

        invers_V = x_vector[8]-0.5*edeltaSq;
        extendIntersection(&x_vector[7], invers_V);

        invers_edeltaSq = 2*(x_vector[8]-x_vector[7]);
        extendIntersection(&edeltaSq, invers_edeltaSq);

        invers_edelta = SQRinterval(edeltaSq);
        extendIntersection(&x_vector[9], invers_edelta);


        // Use V to refine the other errors
        //et
        invers_etSq = (x_vector[7]*2-k1*enSq-efiSq-evSq)/k1;
        extendIntersection(&etSq, invers_etSq);
        invers_et = SQRinterval(etSq);
        extendIntersection(&x_vector[0], invers_et);

        //en
        invers_enSq = (x_vector[7]*2-k1*etSq-efiSq-evSq)/k1;
        extendIntersection(&enSq, invers_enSq);
        invers_en = SQRinterval(enSq);
        extendIntersection(&x_vector[1], invers_en);

        //efi
        invers_efiSq = x_vector[7]*2-k1*etSq-k1*enSq-evSq;
        extendIntersection(&efiSq, invers_efiSq);
        invers_efi = SQRinterval(efiSq);
        extendIntersection(&x_vector[2], invers_efi);

        //ev
        invers_evSq = x_vector[7]*2-k1*etSq-k1*enSq-efiSq;
        extendIntersection(&evSq, invers_evSq);
        invers_ev = SQRinterval(evSq);
        extendIntersection(&x_vector[3], invers_ev+vd);

        //edelta = kdelta - epsilon
        IA epsilon, invers_kdelta, IntFunc1, IntFunc2;

        IntFunc1 = IA((cos(x_vector[2].u())-1)/x_vector[2].u(),
                      (cos(x_vector[2].l())-1)/x_vector[2].l());
        double U = max(abs(x_vector[2].u()), abs(x_vector[2].l()));
        double L = min(abs(x_vector[2].u()), abs(x_vector[2].l()));
        if (x_vector[2].u()* x_vector[2].l()<=0){
            IntFunc2 = IA(sin(U)/U, 1.);
        }
        else {
            IntFunc2 = IA(sin(U)/U, sin(L)/L);
        }
        if (x_vector[2].u() == 0. && x_vector[2].l() ==0){
            IntFunc2 = IA(1.);
        }

        epsilon = x_vector[6]/vd-k1*(x_vector[0]*IntFunc1+x_vector[1]*IntFunc2)-k2*x_vector[2];

        invers_kdelta = x_vector[9] + epsilon;
        extendIntersection(&x_vector[4], invers_kdelta);

        invers_edelta = x_vector[4] - epsilon;
        extendIntersection(&x_vector[9], invers_edelta);




        // Use coordinate transformation as constraints
        double length =2.;

        IA invers_x1 = cos(x_vector[16])*x_vector[0]-sin(x_vector[16])*x_vector[1]+x_vector[14];
        extendIntersection(&x_vector[10], invers_x1);

        IA invers_x2 = sin(x_vector[16])*x_vector[0]+cos(x_vector[16])*x_vector[1]+x_vector[15];
        extendIntersection(&x_vector[11], invers_x2);

        IA invers_phi = x_vector[2]+x_vector[16];
        extendIntersection(&x_vector[12], invers_phi);

        IA invers_delta = atan(length*x_vector[4]);
        extendIntersection(&x_vector[13], invers_delta);


        invers_et = cos(x_vector[16])*(x_vector[10]-x_vector[14]) + sin(x_vector[16])*(x_vector[11]-x_vector[15]);
        extendIntersection(&x_vector[0], invers_et);

        invers_en = -sin(x_vector[16])*(x_vector[10]-x_vector[14]) + cos(x_vector[16])*(x_vector[11]-x_vector[15]);
        extendIntersection(&x_vector[1], invers_en);

        invers_efi =  x_vector[12] - x_vector[16];
        extendIntersection(&x_vector[2], invers_efi);

        invers_kdelta = tan(x_vector[13])/length;
        extendIntersection(&x_vector[4], invers_kdelta);

//        std::cout<< "inside tan" << length*x_vector[4]<<std::endl;
//        std::cout<< "inverse is" << invers_delta<<std::endl;
//        std::cout<<"refined is "<<x_vector[13]<<std::endl;
    }


}


static void extendIntersection(IA *origint, IA refint) {
    origint->l(min(max(origint->l(), refint.l()), origint->u()));
    origint->u(max(min(origint->u(), refint.u()), origint->l()));

}


static IA SQRinterval(IA SQinterval) {

    IA SQR_invers;
    if (SQinterval.u() > 0.0001) {
//        if (SQinterval.l()<0){
//            std::cout<<SQinterval.l()<<std::endl;
//
//        }

        SQR_invers = hull(sqrt(SQinterval), -sqrt(SQinterval));
    } else {
        SQR_invers.u(50. * SQinterval.u() + (0.01 - 5. * pow(10, -3)));
        SQR_invers.l(-50. * SQinterval.u() + (-0.01 + 5. * pow(10, -3)));
    }
    return SQR_invers;
}
