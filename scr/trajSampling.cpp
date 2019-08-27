//
// Created by Xuejiao Yang on 4/12/19.
//
//#include "boost/numeric/interval.hpp"     // Interval arithmetic


#include "includInfo.h"
//#include "NonIntFunc.h"
//#include "includInfo.h"
//#include "rhs.h"
using namespace std;



//double rhs(double t, double* x_vector, double* p, int i);
static int check_flag(void *flagvalue, const char *funcname, int opt);
static int freal (realtype t, N_Vector x, N_Vector xdot, void *user_data);
//static double rhs(double t, double* x_vector, double* p, int i);
typedef struct {// define variables that are uncertain: initial conditions and/or parameters
    double state[num_state_variable], uncertainty[num_uncertainties], referenceInput[num_reference_input]; //x0,y0,v;
} *UserData;

int trajectories(double X0[2*num_state_variable], double Uncertainty[2*num_uncertain_parameter], memberInfoSDI Input){


    int num_state_for_odesolver=num_state_variable; //ODE state solution +   desired state solution from open loop ODE with given control input
    //  extern double X0;
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


    /* Create serial vector of length NEQ for I.C. and abstol */
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
        outfile.open ("States.csv");

    }


    //Generate random numbers
    int total_number_samples = (int)round(pow(samples,num_state_variable)) ; //int y = (int)round(x);
    int total_number_uncertainties = (int)round(pow(samples,num_uncertain_parameter));

    for(int i=0;i<total_number_samples; i++) {
        int ix = 0;
        double next_digit = i + 0.;
        for (int n=num_state_variable-1; n>=0; n--) {
            ix = (int) (next_digit  / pow(samples, n));
            (*data).state[n] = X0[n] + (X0[n + num_state_variable] - X0[n]) / (samples - 1) * ix;
            next_digit -= ix * pow(samples,n);


        }
        for (int j = 0; j < total_number_uncertainties; j++) {
            int ip = 0;
            double next_digit_uncertain = j + 0.;
            for (int m = num_uncertain_parameter - 1; m >= 0; m--) {
                ip = (int) ( next_digit_uncertain/ pow(samples, m));
                (*data).uncertainty[m] = Uncertainty[m] +
                                         (Uncertainty[m + num_uncertain_parameter] - Uncertainty[m]) / (samples - 1) *
                                         ip;

                next_digit_uncertain -= ip * pow(samples,m);
            }


            for (int n = 0; n < num_state_variable  ; n++) {
                NV_Ith_S(x, n) =(*data).state[n];
            }

            for(int pctime=0;pctime<num_piecewise;pctime++) {

                double tfinal = Tf + TMULT * NOUT;

                for (int j = 0; j < num_reference_input; j++) {
                    (*data).referenceInput[j] = Input.InputReference[pctime * num_reference_input + j];

                }


                cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
                if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) return (1);
                flag = CVodeInit(cvode_mem, freal, T0, x);
                if (check_flag(&flag, "CVodeInit", 1)) return (1);
                flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
                if (check_flag(&flag, "CVodeSVtolerances", 1)) return (1);
                flag = CVodeSetUserData(cvode_mem, data);
                if (check_flag(&flag, "CVodeSetUserData", 1)) return (1);

                /* Create dense SUNMatrix for use in linear solves */
                A = SUNDenseMatrix(num_state_for_odesolver, num_state_for_odesolver);
                if (check_flag((void *) A, "SUNDenseMatrix", 0)) return (1);

                /* Create dense SUNLinearSolver object for use by CVode */
                LS = SUNDenseLinearSolver(x, A);
                if (check_flag((void *) LS, "SUNDenseLinearSolver", 0)) return (1);

                /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
                flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
                if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return (1);

                iout = 0;
                tout = Tf;

                // create a 1*2n matrix to store the time and all bounds at each time stepdouble time_step_Data[num_state_for_odesolver + 1];
                double time_step_Data[num_state_for_odesolver + 1];
                while (1) {
                    flag = CVode(cvode_mem, tout, x, &t, CV_NORMAL);
                    if (File_Output) {
                        time_step_Data[0] = t+tfinal*pctime;
                        for (int iw = 1; iw <= num_state_for_odesolver; iw++) {
                            time_step_Data[iw] = Ith(x, iw);
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


                } //loop of the ode solver




            }





            } //loop of uncertain parameters


    } //loop of the initial conditions



    outfile.close();



    /* Free y and abstol vectors */
    N_VDestroy(x);
    N_VDestroy(abstol);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */
    SUNLinSolFree(LS);

    /* Free the matrix memory */
    SUNMatDestroy(A);

    return (0);
}





static int freal (realtype t, N_Vector x, N_Vector xdot, void *user_data) {

    double X[num_state_variable]; // state variables, desired state variables

    UserData data;
    data = (UserData) user_data;

    double P[num_uncertain_parameter];
    for (int i = 0; i < num_uncertain_parameter ; i++) {
        P[i] = (*data).uncertainty[i];
    }




    for (int i = 0; i < num_state_variable ; i++) {
        X[i] = NV_Ith_S(x,i);
    }


    double fx;

    //1) compute RHS of original state variables
    for (int i = 0; i < num_state_variable; i++) { //loop from fi=1 to fi=n

        fx=rhs(t,X,P,(*data).referenceInput,i);
        NV_Ith_S(xdot, i)=fx;

        //NV_Ith_S(xdot, i+num_state_variable)=fx;

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