//
// Created by Xuejiao Yang on 4/29/19.
//

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


using namespace std;


static int check_flag(void *flagvalue, const char *funcname, int opt);
static int freal (realtype t, N_Vector x, N_Vector xdot, void *user_data);
static double rhs(double t, double* x_vector,  double* ref, int i);



typedef struct {// define variables that are uncertain: initial conditions and/or parameters
    double referenceInput[1];
} *UserData;





int main(){

    double X0[3] = {0.,0.,0.}, ref[2] ={1./30,-1./30};

    int num_state_for_odesolver=3; //ODE state solution +   desired state solution from open loop ODE with given control input

    N_Vector x, abstol;
    realtype reltol, t, tout;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    int flag, iout;

    A=NULL;
    LS=NULL;
    abstol = NULL;
    cvode_mem = NULL;
    x = N_VNew_Serial(num_state_for_odesolver);


    // Create serial vector of length NEQ for I.C. and abstol
    if (check_flag((void *) x, "N_VNew_Serial", 0)) return (1);
    abstol = N_VNew_Serial(num_state_for_odesolver);
    if (check_flag((void *) abstol, "N_VNew_Serial", 0)) return (1);


    UserData data;
    data = (UserData) malloc(sizeof *data);
    if(check_flag((void *)data, "malloc", 2)) return(1);
    reltol = RTOL;
    for(int i=1;i<=num_state_for_odesolver;i++){
        Ith(abstol,i)=ATOL; }



    //Write results into csv files.
    ofstream outfile;
    if(File_Output){
        outfile.open ("ref.csv");

    }



    for (int i = 0; i < 3  ; i++) {
        Ith(x,i+1) = X0[i];
    }

    for(int pctime=0;pctime<2;pctime++){

        double tfinal = Tf+TMULT*NOUT;


        (*data).referenceInput[0]=ref[pctime];


        cvode_mem =CVodeCreate(CV_ADAMS,CV_FUNCTIONAL);
        if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) return (1);
        flag = CVodeInit(cvode_mem, freal, T0, x);
        if (check_flag(&flag, "CVodeInit", 1)) return (1);
        flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
        if (check_flag(&flag, "CVodeSVtolerances", 1)) return (1);
        flag = CVodeSetUserData(cvode_mem, data);
        if (check_flag(&flag, "CVodeSetUserData", 1)) return (1);

        // Create dense SUNMatrix for use in linear solves
        A = SUNDenseMatrix(num_state_for_odesolver, num_state_for_odesolver);
        if (check_flag((void *) A, "SUNDenseMatrix", 0)) return (1);

        // Create dense SUNLinearSolver object for use by CVode
        LS = SUNDenseLinearSolver(x, A);
        if (check_flag((void *) LS, "SUNDenseLinearSolver", 0)) return (1);

        // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
        flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
        if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return (1);

        iout = 0;
        tout = Tf;

        // create a 1*2n matrix to store the time and all bounds at each time step
        double time_step_Data[num_state_for_odesolver + 1];
        while (1) {
            flag = CVode(cvode_mem, tout, x, &t, CV_NORMAL);




            if(File_Output){
                time_step_Data[0] = t+tfinal*pctime;
                for (int i = 1; i <= num_state_for_odesolver; i++) {
                    time_step_Data[i] = Ith(x, i);
                }


                for (int count = 0; count < num_state_for_odesolver + 1; count++) {
                    outfile << time_step_Data[count] << ",";
                }
                outfile << endl;

            }



            if (check_flag(&flag, "CVode", 1)) break;
            if (flag == CV_SUCCESS) {
                iout++;
                tout += TMULT;
            }
            if (iout == NOUT) break;
        } //end of integration while loop







    }

    outfile.close();


    // Free y and abstol vectors
    N_VDestroy(x);
    N_VDestroy(abstol);

    // Free integrator memory
    CVodeFree(&cvode_mem);

    // Free the linear solver memory
    SUNLinSolFree(LS);

    // Free the matrix memory
    SUNMatDestroy(A);

    return (0);
}



static int freal (realtype t, N_Vector x, N_Vector xdot, void *user_data) {

    //double X[2*num_state_variable]; // state variables, desired state variables



    UserData data;
    data = (UserData) user_data;
    double X[3];
    for (int i = 0; i < 3 ; i++) {
        X[i] = NV_Ith_S(x,i);
    }

    double fx;

    //1) compute RHS of original state variables
    for (int i = 0; i < 3; i++) { //loop from fi=1 to fi=n
        //flat the upper bound to lower bound

        //flat the lower bound to upper bound
        fx=rhs(t,X,(*data).referenceInput,i);


        NV_Ith_S(xdot, i)=fx;


    } //End of computing bounds of original right hand side


    return(0);
}


int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    // Check if SUNDIALS function returned NULL pointer - no memory allocated
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

        // Check if flag < 0
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return(1); }}

        // Check if function returned NULL pointer - no memory allocated
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}

double rhs(double t, double* x_vector,  double* ref, int i) {
    double fx,k = ref[0],thetad = x_vector[0];
    switch (i) {
        case 0:
            fx = k;
            break;
        case 1:
            fx = cos(thetad);
            break;
        case 2:
            fx = sin(thetad);
            break;
    }





    return fx;

}


