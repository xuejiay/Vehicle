//
// Created by Xuejiao Yang on 4/24/19.
//

#ifndef EXAMPLE1_GET_INFO_H
#define EXAMPLE1_GET_INFO_H

#include "math.h"
#include "privateData.h"
#include "includInfo.h"

//memberInfoSDI InputDataSDI = {{0.5,M_PI/12, 1, M_PI/6},{5.,6.},{1./30}};
//
//memberInfoDI_constraints InputDataDI_cst = {{0.5,M_PI/12,1./2.*(pow(0.5,2)+1./4.*pow(M_PI/12,2)), 1., M_PI/6,1./2.*(1.+1./4.*pow(M_PI/6,2))}
//        ,{5.,6.},{1./30},{-100000,-100000,-100000,100000,100000,100000}};
//
//memberInfoDI_cst_advance  InputDataDI_cst_advance = {{0.5,M_PI/12,1./2.*(pow(0.5,2)+1./4.*pow(M_PI/12,2)),  cos(M_PI/6), tan(M_PI/12),
//                               1., M_PI/6,         1./2.*(1.+1./4.*pow(M_PI/6,2)), cos(M_PI/12), tan(M_PI/6)}
//        ,{5.,6.},{1./30},{-100000,-100000,-100000,-100000,-100000,100000,100000,100000,100000,100000}};
//
//

//memberInfoSDI InputDataSDI = {{0.8,M_PI/12, 1, M_PI/6},{5.,6.},{1./30,-1./30}};
//
//memberInfoDI_constraints InputDataDI_cst = {{0.8,M_PI/12,1./2.*(pow(0.8,2)+1./4.*pow(M_PI/12,2)), 1., M_PI/6,1./2.*(1.+1./4.*pow(M_PI/6,2))}
//        ,{5.,6.},{1./30,-1./30},{-100000,-100000,-100000,100000,100000,100000}};
//
//memberInfoDI_cst_advance  InputDataDI_cst_advance = {{0.8,M_PI/12,1./2.*(pow(0.8,2)+1./4.*pow(M_PI/12,2)),  cos(M_PI/6), tan(M_PI/12),
//                                                             1., M_PI/6,         1./2.*(1.+1./4.*pow(M_PI/6,2)), cos(M_PI/12), tan(M_PI/6)}
//        ,{5.,6.},{1./30,-1./30},{-100000,-100000,-100000,-100000,-100000,100000,100000,100000,100000,100000}};
//
//


memberInfoSDI InputDataSDI = {{0.8,M_PI/12, 1, M_PI/6},{5.,6.},{1./30}};

memberInfoDI_constraints InputDataDI_cst = {{0.8,M_PI/12,1./2.*(pow(0.8,2)+1./4.*pow(M_PI/12,2)), 1., M_PI/6,1./2.*(1.+1./4.*pow(M_PI/6,2))}
        ,{5.,6.},{1./30},{-100000.,-100000.,0.,100000.,100000.,100000.}};

memberInfoDI_cst_advance  InputDataDI_cst_advance = {{0.8,M_PI/12,1./2.*(pow(0.8,2)+1./4.*pow(M_PI/12,2)),  cos(M_PI/6), tan(M_PI/12),
                                                             1., M_PI/6,         1./2.*(1.+1./4.*pow(M_PI/6,2)), cos(M_PI/12), tan(M_PI/6)}
        ,{5.,6.},{1./30},{-100000.,-100000.,0.,-100000.,-100000.,100000.,100000.,100000.,100000.,100000.}};




/*
memberInfoSDI InputDataSDI = {{1,M_PI/6, 1, M_PI/6},{5.,6.},{1./30,-1./30}};

memberInfoDI_constraints InputDataDI_cst = {{1,M_PI/6,1./2.*(pow(1.,2)+1./4.*pow(M_PI/6,2)),
                                             1., M_PI/6,1./2.*(1.+1./4.*pow(M_PI/6,2))}
        ,{5.,6.},{1./30,-1./30},{-100000,-100000,-100000,100000,100000,100000}};

memberInfoDI_cst_advance  InputDataDI_cst_advance = {{1,M_PI/6,1./2.*(pow(1.,2)+1./4.*pow(M_PI/6,2)),  cos(M_PI/6), tan(M_PI/6),
                                                             1., M_PI/6,         1./2.*(1.+1./4.*pow(M_PI/6,2)), cos(M_PI/6), tan(M_PI/6)}
        ,{5.,6.},{1./30,-1./30},{-100000,-100000,-100000,-100000,-100000,100000,100000,100000,100000,100000}};


*/









#endif //EXAMPLE1_GET_INFO_H
