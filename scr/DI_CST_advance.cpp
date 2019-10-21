//
// Created by Xuejiao Yang on 4/23/19.
//


#include "includInfo.h"


using namespace std;


static int check_flag(void *flagvalue, const char *funcname, int opt);

static int freal (realtype t, N_Vector x, N_Vector xdot, void *user_data);

typedef struct {// define variables that are uncertain: initial conditions and/or parameters
    IA uncertainties[num_uncertainties], stateI[num_state_redundant_advance],
    referenceInput[num_reference_input], natural_bound[num_state_redundant_advance];
} *UserData;



int DI_constraints_advance(double* X0, double* Uncertainty, memberInfoDI_cst_advance Input){


    int num_state_for_odesolver=2*num_state_redundant_advance; //ODE state solution +   desired state solution from open loop ODE with given control input
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
    if(File_Output_DI_advance){
        outfile.open ("StatesBD_DI_CST_advance.csv");

    }


    for(int j=0; j<num_uncertain_parameter; j++){
        (*data).uncertainties[j].l(Uncertainty[j]);
        (*data).uncertainties[j].u(Uncertainty[j+num_uncertain_parameter]);
    }
    for (int i = 0; i < 2*num_state_redundant_advance  ; i++) {
        Ith(x,i+1) = X0[i];
        (*data).natural_bound[i].l(Input.NaturBD[i]);
        (*data).natural_bound[i].u(Input.NaturBD[i+num_state_redundant_advance]);
    }

    for(int pctime=0;pctime<num_piecewise;pctime++) {

        double tfinal = Tf + TMULT * (NOUT-1);

        for (int j = 0; j < num_reference_input; j++) {
            (*data).referenceInput[j].l(Input.InputReference[pctime * num_reference_input + j]);
            (*data).referenceInput[j].u(Input.InputReference[pctime * num_reference_input + j]);
        }


//        std::cout << (*data).referenceInput[1] << std::endl;

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

            //x satisfies the constraint -- using refinement algorithm

            for(int j=0;j<num_state_redundant_advance;j++){
                (*data).stateI[j].l(NV_Ith_S(x,j));
                (*data).stateI[j].u(NV_Ith_S(x,j+num_state_redundant_advance));
            }
            refinement_advance((*data).stateI, (*data).uncertainties);
            //Assign value to x
            for(int j=0;j<num_state_redundant_advance;j++){
                NV_Ith_S(x,j) = (*data).stateI[j].l();
                NV_Ith_S(x,j+num_state_redundant_advance) = (*data).stateI[j].u();
            }

            if(File_Output_DI_advance){
                time_step_Data[0] = t+tfinal*pctime;

                for (int i = 1; i <= num_state_for_odesolver; i++) {
                    time_step_Data[i] = Ith(x, i);
                    //time_step_Data[i] = (*data).stateI[i-1].l();
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

        //for(int i = 0; i < 2*num_state_redundant_advance  ; i++) {
         //   X0[i] = Ith(x,i+1);
        //}
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

    //double X[2*num_state_redundant]; // state variables, desired state variables



    UserData data;
    data = (UserData) user_data;



    IA fxl, fxu;
    int flat_flag[2 * num_state_redundant_advance] = {1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,0};
//    int flat_flag[2 * num_state_redundant_advance] = {1,1,1,1,1,1,1,0,0,0};

    // 1) compute RHS of original state variables
    for (int i = 0; i < num_state_redundant_advance; i++) { //loop from fi=1 to fi=n

        for(int j=0;j<num_state_redundant_advance;j++){
            (*data).stateI[j].l(NV_Ith_S(x,j));
            (*data).stateI[j].u(NV_Ith_S(x,j+num_state_redundant_advance));
        }
        if(flat_flag[i] == 1){
            //flat the upper bound to lower bound
            (*data).stateI[i].u(NV_Ith_S(x, i));
            refinement_advance((*data).stateI, (*data).uncertainties);
        }
        fxl=rhsI(t,(*data).stateI, (*data).uncertainties,(*data).referenceInput,i);


        //flat the lower bound to upper bound
        for(int j=0;j<num_state_redundant_advance;j++){
            (*data).stateI[j].l(NV_Ith_S(x,j));
            (*data).stateI[j].u(NV_Ith_S(x,j+num_state_redundant_advance));
        }
        if(flat_flag[i] == 1) {
            (*data).stateI[i].l(NV_Ith_S(x, i + num_state_redundant_advance));
            refinement_advance((*data).stateI, (*data).uncertainties);
        }
        fxu=rhsI(t,(*data).stateI, (*data).uncertainties,(*data).referenceInput,i);


        NV_Ith_S(xdot, i)=fxl.l();
        NV_Ith_S(xdot, i+num_state_redundant_advance)=fxu.u();

    } // End of computing bounds of original right hand side

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



