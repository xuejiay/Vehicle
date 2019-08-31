
#include "MainHeader.h"
#include "includInfo.h"
#include "Get_Info.h"

//#include "interval.hpp"

//void rhs(double fx[2], double x[2], double w[2]);
using namespace std;
double get_cpu_time();

int main() {
    double cpu0  = get_cpu_time();

//    trajectories(InputDataSDI.X0, InputDataSDI.Uncertainty,InputDataSDI);
    SDI(InputDataSDI.X0, InputDataSDI.Uncertainty,InputDataSDI);
    //cout << pow(M_PI/12,2)/4 << endl;
//    DI_constraints(InputDataDI_cst.X0, InputDataDI_cst.Uncertainty,InputDataDI_cst);
//    DI_constraints_advance(InputDataDI_cst_advance.X0, InputDataDI_cst_advance.Uncertainty, InputDataDI_cst_advance);
    double cpu1  = get_cpu_time();
    cout << "CPU Time  of trajectory simulation is  " << cpu1  - cpu0  << endl;
    return (0);
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
