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

static void  extendIntersection(IA* origint, IA refint);
static IA SQRinterval(IA SQinterval);

double rhs(double t, double* x_vector, double* p, double* ref, int i) {
    double fx;
    double e = x_vector[0], theta_e = x_vector[1], v = p[0], k = ref[0], a = 2., xi = 0.7, epsilon, g2, g1;
    g2 = pow(a,2);
    switch (i) {
        case 0:
            fx = (1-k*e)*tan(theta_e);
            break;
        case 1:
            g1 = 2.*xi*a*sqrt(pow(v,2)+epsilon);
            fx = (1.-k*e)*(-g1*theta_e/v/cos(theta_e)-g2*tan(theta_e)*e/theta_e);
            break;
    }





    return fx;

}



IA rhsI(double t, IA* x_vector, IA* p, IA* ref, int i) {
    IA fx;
    IA e = x_vector[0], theta_e = x_vector[1], v = p[0], k = ref[0], a = 2., xi = 0.7, epsilon=0.1, g2, g1;
    g2 = pow(a,2);
    IA inter, z1 = x_vector[3], z2 = x_vector[4];
    IA ab_the= abs(theta_e);
    switch (i) {
        case 0:
            fx = (1-k*e)*tan(theta_e);
            break;
        case 1:
            g1 = 2*xi*a*sqrt(pow(v,2)+epsilon);
            fx = (1-k*e)*(-2*xi*a*sqrt(1.+0.1/pow(v,2))*theta_e/cos(theta_e)-g2*tan(theta_e)*e/theta_e);
//            fx = (1-k*e)*(-2.8*sqrt(1.+0.1/pow(v,2))*theta_e/cos(theta_e)-4.*e*(1.+pow(theta_e,2)/3));
//            fx = (1-k*e)*(-5.7*sqrt(1.+0.1/pow(v,2))*theta_e/cos(theta_e)-g2*tan(theta_e)*e/theta_e);

            break;
        case 2:
            g1 = 2*xi*a*sqrt(pow(v,2)+epsilon);
            fx = -2*xi*a*sqrt(1.+0.1/pow(v,2))/g2*p[1]*(1-k*e)/cos(theta_e);
//            fx = -5.7*sqrt(1.+0.1/pow(v,2))/g2*p[1]*(1-k*e)/cos(theta_e);
            //fx = -g1/g2*pow(theta_e,2)*(1-k*e)/v/cos(theta_e);
            break;
        case 3:
            g1 = 2*xi*a*sqrt(pow(v,2)+epsilon);

            if (abs(theta_e.l())<pow(10,-8) && abs(theta_e.u())<pow(10,-8)){
                //if (theta_e.l()==0 && theta_e.u()==0){
                inter.l(1.);
                inter.u(1.);
            }
            else if(theta_e.l()* theta_e.u()<=pow(10,-8)){
                inter.l(1.);
                inter.u(tan(ab_the.u())/ab_the.u());
            }
            else {
                inter.l(tan(ab_the.l())/ab_the.l());
                inter.u(tan(ab_the.u())/ab_the.u());
            }

            fx=-z1*z2*(1-k*e)*(-g1*theta_e/v/z1-g2*inter*e);

            //fx = 2*g1;
            break;
        case 4:
            g1 = 2*xi*a*sqrt(pow(v,2)+epsilon);

            if  (abs(theta_e.l())<pow(10,-8) && abs(theta_e.u())<pow(10,-8)) {
                //if (theta_e.l()==0 && theta_e.u()==0){
                inter.l(1.);
                inter.u(1.);
            }
            else if (theta_e.l()* theta_e.u()<=pow(10,-8)){
                inter.l(1.);
                inter.u(tan(ab_the.u())/ab_the.u());
            }
            else {
                inter.l(tan(ab_the.l())/ab_the.l());
                inter.u(tan(ab_the.u())/ab_the.u());
            }

            fx = (1-k*e)*(-g1*theta_e/v/z1-g2*inter*e)/p[2];
            break;
    }






    return fx;

}





void refinement(IA x_vector[num_state_redundant], IA* p){
    //natural bound here
    for(int l=0; l<num_loop; l++){
        IA //e=x_vector[0], theta_e = x_vector[1], V = x_vector[2],
                a=2, g2=pow(a,2), invers_V, eSQR, theta_eSQR, invers_eSQR, invers_thetaeSQR, invers_e, invers_thetae;

        eSQR=pow(x_vector[0],2);
        theta_eSQR=pow(x_vector[1],2);


        // refine lyapunov function V
        invers_V = (eSQR+theta_eSQR/g2)/2;
        // Taking intersection with the original interval
        extendIntersection(& x_vector[2], invers_V);


        // refine e^2
        invers_eSQR = 2*x_vector[2]-theta_eSQR/g2;
        extendIntersection(& eSQR, invers_eSQR);
        // refine e
        if (eSQR.u() > 0.0001) {

            invers_e = hull(sqrt(eSQR), -sqrt(eSQR));
        }
        else{
            invers_e.u(50.*eSQR.u()+(0.01-5.*pow(10.,-3)));
            invers_e.l(-50.*eSQR.u()+(-0.01+5.*pow(10.,-3)));
        }
        extendIntersection(& x_vector[0], invers_e);


        // refine theta_e^2
        invers_thetaeSQR = (2*x_vector[2]-eSQR)*g2;
        extendIntersection(& theta_eSQR, invers_thetaeSQR);
        p[1] = theta_eSQR;
        // refine e
        if (theta_eSQR.u() > 0.0001) {

            invers_thetae = hull(sqrt(theta_eSQR), -sqrt(theta_eSQR));
        }
        else{
            invers_thetae.u(50.*theta_eSQR.u()+(0.01-5.*pow(10,-3)));
            invers_thetae.l(-50.*theta_eSQR.u()+(-0.01+5.*pow(10,-3)));
        }
        extendIntersection(& x_vector[1], invers_thetae);

    }


}



void refinement_advance(IA x_vector[num_state_redundant_advance], IA* p){
    //natural bound here
    for(int l=0; l<num_loop; l++){
        IA //e=x_vector[0], theta_e = x_vector[1], V = x_vector[2],
                a=2, g2=pow(a,2), invers_V, eSQR, theta_eSQR, invers_eSQR, invers_thetaeSQR, invers_e, invers_thetae;

        eSQR=pow(x_vector[0],2);
        theta_eSQR=pow(x_vector[1],2);


        // refine lyapunov function V
        invers_V = (eSQR+theta_eSQR/g2)/2;

        // Taking intersection with the original interval
        extendIntersection(& x_vector[2], invers_V);


        // refine e^2
        invers_eSQR = 2*x_vector[2]-theta_eSQR/g2;
        extendIntersection(& eSQR, invers_eSQR);
        // refine e
        if (eSQR.u() > 0.0001) {

            invers_e = hull(sqrt(eSQR), -sqrt(eSQR));
        }
        else{
            invers_e.u(50.*eSQR.u()+(0.01-5.*pow(10.,-3)));
            invers_e.l(-50.*eSQR.u()+(-0.01+5.*pow(10.,-3)));
        }
        extendIntersection(& x_vector[0], invers_e);


        // refine theta_e^2
        invers_thetaeSQR = (2*x_vector[2]-eSQR)*g2;
        extendIntersection(& theta_eSQR, invers_thetaeSQR);
        p[1] = theta_eSQR;
        // refine e
        if (theta_eSQR.u() > 0.0001) {

            invers_thetae = hull(sqrt(theta_eSQR), -sqrt(theta_eSQR));
        }
        else{
            invers_thetae.u(50.*theta_eSQR.u()+(0.01-5.*pow(10,-3)));
            invers_thetae.l(-50.*theta_eSQR.u()+(-0.01+5.*pow(10,-3)));
        }
        extendIntersection(& x_vector[1], invers_thetae);


        //z1 = cos(theta)
        //z2 = tan(theta)


        IA invers_z1, invers_z2, z1SQR, z1SQR_invers, z2SQR, z2SQR_invers;

        invers_z1 = cos(x_vector[1]);
        invers_z2 = tan(x_vector[1]);
        extendIntersection(& x_vector[3], invers_z1);
        extendIntersection(& x_vector[4], invers_z2);



        //z2SQR + 1= 1/z1SQR
        z1SQR=pow(x_vector[3],2);
        z2SQR=pow(x_vector[4],2);

        z1SQR_invers = 1./(z2SQR+ 1.);
        extendIntersection(& z1SQR, z1SQR_invers);
        p[2] = z1SQR;

        z2SQR_invers = 1./z1SQR - 1.;
        extendIntersection(& z2SQR, z2SQR_invers);


        //square root of z1 and z2
        IA z1_invers = SQRinterval(z1SQR);
        extendIntersection(& x_vector[3], z1_invers);

        IA z2_invers = SQRinterval(z2SQR);
        extendIntersection(& x_vector[4], z2_invers);


        //use z1 and z2 to refine for theta_e
        //z1 = cos(theta)
        //z2 = tan(theta)
        // avoid beyond -1,1

        x_vector[3].l(min(max(-1.,x_vector[3].l()),1.));
        x_vector[3].u(max(min(1.,x_vector[3].u()),-1.));

        //x_vector[4].l(min(max(-1.,x_vector[4].l()),1.));
        //x_vector[4].u(max(min(1.,x_vector[4].u()),-1.));

        invers_thetae = atan(x_vector[4]);
        extendIntersection(& x_vector[1], invers_thetae);


        invers_thetae = hull(acos(x_vector[3]),-acos(x_vector[3]));
        extendIntersection(& x_vector[1], invers_thetae);


        //use z1 and z2 to refine for theta_e
        // z2*z1 = sin(theta)
        IA Z1Z2 = sin(x_vector[1]);
        //IA Z1Z2_invers=x_vector[4]*x_vector[3];
        extendIntersection(& Z1Z2, x_vector[4]*x_vector[3]);

        invers_thetae = asin(Z1Z2);
        extendIntersection(& x_vector[1], invers_thetae);

    }


}








static void extendIntersection(IA* origint, IA refint){
    origint->l( min(max(origint->l(),refint.l()),origint->u()));
    origint->u( max(min(origint->u(),refint.u()),origint->l()));

}


static IA SQRinterval(IA SQinterval){

    IA SQR_invers;
    if (SQinterval.u() > 0.0001) {

        SQR_invers = hull(sqrt(SQinterval), -sqrt(SQinterval));
    }
    else{
        SQR_invers.u(50.*SQinterval.u()+(0.01-5.*pow(10,-3)));
        SQR_invers.l(-50.*SQinterval.u()+(-0.01+5.*pow(10,-3)));
    }
    return SQR_invers;
}
