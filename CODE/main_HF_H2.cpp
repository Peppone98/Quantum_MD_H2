#include "definitions.h"
#include "xc.h"

using namespace std;


int main (int argc, char *argv[]){

    /**** Create space for matrices and vectors ****/
	int n, n_iter = 20;
	R R_A, R_B;
    double Q[2*N][2*N][2*N][2*N];
    double X[CP_iter];
    double E_1s=0., lambda=0., E0 = 0., E_old = 0., dE0_dX = 0., T_N = 0., T_e = 0.0;
	gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *S_auxiliary = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *H = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *F = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *U = gsl_matrix_calloc(2*N, 2*N); 
    gsl_matrix *V = gsl_matrix_calloc(2*N, 2*N);
	gsl_vector *c = gsl_vector_alloc(2*N);
    gsl_vector *c_old = gsl_vector_alloc(2*N);
	gsl_vector_set_all(c, 1.);
    gsl_matrix *V_xc = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *dVxc_dX = gsl_matrix_alloc(2*N, 2*N);

    /**** Initial nuclei positions ****/
    R_A.x = 0., R_A.y = 0., R_A.z = 0.;
    R_B.x = 1.2, R_B.y = 0., R_B.z = 0.;

    /**** Fill S and compute V ****/
    create_S(S, R_A, R_B);
    gsl_matrix_memcpy(S_auxiliary, S);
	diag_S(S_auxiliary, V);

    /**** Fill H and Q ****/
    one_body_H(H, R_A, R_B);
    build_Q(Q, R_A, R_B);

    /**** Initial conditions ****/
    X[0] = sqrt(scalar_prod(R_A, R_B));
    X[1] = X[0];
    double norm = normalization(c, S);
    gsl_vector_memcpy(c_old, c);
	



    /**** List of options: ****/
    /**** SC_Hartree_Fock ****/
    /**** SC_DFT ****/
    /**** Evolve_coefficients ****/
    /**** CPMD_HF ****/
    /**** CG_min ****/
    /**** CG_CP ****/
    /**** CG_CP_superposition ****/
    /**** EX_CORR ****/
    /**** XC_DER ****/
    /**** CPMD_DFT ****/
    /**** ADAPTIVE ****/
    /**** XC_DER_ADAPTIVE ****/





    
    /**** SELF-CONSISTENT HARTREE-FOCK ****/
    if (string(argv[1]) == "SC_HF"){

        ofstream myfile;
	    myfile.open("Energy_profile.txt", ios :: out | ios :: trunc);

        for(int i=0; i<24; i++){
            
            /**** Set the interatomic distance ****/
            R_B.x = 0.8 + i*0.1;
            X[0] = sqrt(scalar_prod(R_A, R_B));
            std::cout << "Interatomic distance X = " << X[0] << endl;

            /**** Refill S and compute V with new value of R_B.x ****/
            create_S(S, R_A, R_B);
            gsl_matrix_memcpy(S_auxiliary, S);
            diag_S(S_auxiliary, V);

            /**** Refill H and Q ****/
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            std::cout << "   " << endl;
            std::cout << "HF energies: " << endl;

            E0 = 1.0, E_old = 0.0;
            while(fabs(E0 - E_old) > 1E-7){

                /**** Build Fock matrix with the added exchange part ****/
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);

                /**** tmp_F is needed beacuse solve_FC_eSC destroys lower part of F ****/
                gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
                gsl_matrix_memcpy(tmp_F, F);
                E_1s = solve_FC_eSC(tmp_F, V, U);

                /**** Compute the energy and print it ****/
                E_old = E0;
                E0 = compute_E0(F, H, c) + 1./X[0];
                std::cout << E0 << endl;

                /**** Get the c vector from eigenvector matrix U ****/
                gsl_matrix_get_col(c, U, 0);
                gsl_matrix_free(tmp_F);
            }

            /**** Print the coefficients ****/
            std::cout << "   " << endl;
            std::cout << "Vector of coefficients" << endl;
            for(n=0; n<2*N; n++){
                std::cout << gsl_vector_get(c, n) << endl;
            }
            
            /**** Print energy to screen and fill the .txt file ****/
            std::cout << "   " << endl;
            std::cout << "Total HF energy: " << E0 << " H" << endl;
            std::cout << "*******************************************************" << endl;
            std::cout << "   " << endl;
            myfile << X[0] << "       " << E0 << endl;
        }  
        myfile.close(); 
    }










    /**** SELF-CONSISTENT DFT ****/
    if (string(argv[1]) == "SC_DFT"){

        /**** Get the name of the functional ****/
        xc_func_type func;
        xc_func_init(&func, FUNCTIONAL_C, XC_UNPOLARIZED);
        std::cout << "************ XC FUNCTIONAL ************" << endl;
        std::cout << func.info->name << endl;
        std::cout << "  " << endl;
        ofstream myfile;
	    myfile.open("Energy_profile.txt", ios :: out | ios :: trunc);

        for(int i=0; i<24; i++){
            
            /**** Set the interatomic distance ****/
            R_B.x = 0.8 + i*0.1;
            X[0] = sqrt(scalar_prod(R_A, R_B));
            std::cout << "Interatomic distance X = " << X[0] << endl;

            /**** Refill S and compute V with new value of R_B.x ****/
            create_S(S, R_A, R_B);
            gsl_matrix_memcpy(S_auxiliary, S);
            diag_S(S_auxiliary, V);

            /**** Refill H and Q ****/
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            std::cout << "   " << endl;
            std::cout << "DFT energies: " << endl;

            E0 = 1.0, E_old = 0.0;
            while(fabs(E0 - E_old) > 1E-7){

                /**** Build Fock matrix with the added exchange part ****/
                two_body_F(Q, c, F);
                gsl_matrix_scale(F, 1. + a_x);
                gsl_matrix_add(F, H);

                /**** Adding the exchange and correlation ****/
                string s = "V_xc";
                Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[0], s);
                gsl_matrix_add(F, V_xc);

                /**** tmp_F is needed beacuse solve_FC_eSC destroys lower part of F ****/
                gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
                gsl_matrix_memcpy(tmp_F, F);
                E_1s = solve_FC_eSC(tmp_F, V, U);

                /**** Compute the energy and print it ****/
                E_old = E0;
                E0 = compute_E0(F, H, c) + 1./X[0];
                std::cout << E0 << endl;

                /**** Get the c vector from eigenvector matrix U ****/
                gsl_matrix_get_col(c, U, 0);
                gsl_matrix_free(tmp_F);
            }

            /**** Print the coefficients ****/
            std::cout << "   " << endl;
            std::cout << "Vector of coefficients" << endl;
            for(n=0; n<2*N; n++){
                std::cout << gsl_vector_get(c, n) << endl;
            }
            
            /**** Print energy to screen and fill the .txt file ****/
            std::cout << "   " << endl;
            std::cout << "Total DFT energy: " << E0 << " H" << endl;
            std::cout << "*******************************************************" << endl;
            std::cout << "   " << endl;
            myfile << X[0] << "       " << E0 << endl;
        }  
        myfile.close();

        std::cout << "Employed functional: " << func.info->name << endl;
        std::cout << "E_xc[n] = (1 - a_x)*E_HF_x[n] + a_x*E_x[n] + E_c[n]" << endl;
        std::cout << "Value of a_x: " << a_x << endl; 
    }









    /**** BORN OPPENHEIMER MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "BOMD_HF"){

        /**** Start with very short equilibration with CPMD to get the first lambda ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Get the first reference value for lambda_shake ****/
        double lambda_shake = lambda;
        std::cout << lambda_shake << endl;

        /**** Create the files for storing the energies and the coefficients ****/
        string X_en = "MD_BO_X_energies.txt";
        string coeff = "MD_BO_coeff.txt";
        ofstream X_en_file;
        ofstream coeff_file;
        X_en_file.open(X_en, ios :: out | ios :: trunc);
        coeff_file.open(coeff, ios :: out | ios :: trunc);

        for(n=1; n<400; n++){

            /**** Fill H and Q for electronic problem ****/
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            /**** Create S and the V matrix, needed for the SC procedure ****/
            create_S(S, R_A, R_B);
            gsl_matrix_memcpy(S_auxiliary, S);
            diag_S(S_auxiliary, V);

            /**** SC while loop: here we get the c coefficients ****/
            E0 = 1.0, E_old = 0.0;
            std::cout << "******** SC procedure *******" << endl;
            while(fabs(E0 - E_old) > 1E-7){

                /**** Build Fock matrix with the added exchange part (self-consistent part!) ****/
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);

                /**** tmp_F is needed beacuse solve_FC_eSC destroys lower part of F ****/
                gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
                gsl_matrix_memcpy(tmp_F, F);
                E_1s = solve_FC_eSC(tmp_F, V, U);

                /**** Compute the energy and print it ****/
                E_old = E0;
                E0 = compute_E0(F, H, c) + 1./X[n];
                std::cout << E0 << endl;

                /**** Get the c vector from eigenvector matrix U ****/
                gsl_matrix_get_col(c, U, 0);
                gsl_matrix_free(tmp_F);
            }

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill H, Q and F for nuclear problem ****/
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X using the newly computed lambda_shake ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            double denominator = dsigma_dX(c, R_A, R_B);
            std::cout << "Before: " << R_B.x << endl;
            partial_update_shake(X[n], X[n - 1], &R_B, dE0_dX, lambda_shake, denominator);
            std::cout << "After: " << R_B.x << endl;
            double new_denominator = denominator*dsigma_dX(c, R_A, R_B);
            double numerator = sigma(c, R_A, R_B);
            double new_lambda = lambda_shake + numerator/new_denominator;
            std::cout << "Correction to lambda: " << numerator/new_denominator << endl;
            X[n + 1] = 2.0*X[n] - X[n - 1] - h_N*h_N/M_N*(2.0*dE0_dX + new_lambda*denominator);
            R_B.x = X[n + 1];
            std::cout << "New X: " << R_B.x << endl;
            std::cout << "  " << endl;
            X_en_file << X[n] << "    " << E0 << endl;
            std::cout << "lambda_shake: " << new_lambda << endl;    
        }

        X_en_file.close();
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Born-Oppenheimer Molecular Dynamics has been executed for " << CP_iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }









    /**** EVOLUTION OF THE COEFFICIENTS C(t) AT FIXED X ****/
    if(string(argv[1]) == "Evolve_coefficients"){
        ofstream myfile;
	    myfile.open("Energies_C_evolution.txt", ios :: out | ios :: trunc);

        /**** MD cycle: compute new c vector and keep F updated ****/
        for(n=0; n<50; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            myfile << compute_E0(F, H, c) + 1./X[0] << endl;
            lambda = update_c(F, S, c, c_old);
        }
        myfile.close();

        /**** Print the coefficients ****/
        std::cout << "Vector of coefficients: " << endl;
        for(n=0; n<2*N; n++){
            std::cout << gsl_vector_get(c, n) << endl;
        }
    }













    /**** CAR PARRINELLO MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "CPMD_HF"){
        double v = 0.0, v_max = 0.03;

        /**** Choose the number of trajectories to generate ****/
        int N_traj;
        string response;
        std::cout << "Number of trajectories: ";
        std::cin >> N_traj;
        std::cout << "Do you want to save the coefficients c? (Y or N)";
        std::cin >> response;


        /**** Create the files for storing the energies and the coefficients ****/
        string coeff = "outputs/HF_traj/CPMD_HF_coeff.txt";
        ofstream coeff_file;
        coeff_file.open(coeff, ios :: out | ios :: trunc);

        /**** Seed for random positions ****/
        srand((unsigned) time(NULL));

        for(int i=0; i<N_traj; i++){

            /**** Define the output trajectory files ****/
            string X_en = "outputs/HF_traj/CPMD_HF_" + to_string(i) + ".txt";
            ofstream X_en_file;
            X_en_file.open(X_en, ios :: out | ios :: trunc);
            
            /**** Generate random initial positions ****/
            X[1] = 1.0 + 1.0*rand()/RAND_MAX;
            R_B.x = X[1];

            /**** Generate random initial velocities ****/
            v = -v_max + 2.0*v_max*rand()/RAND_MAX;
            X[0] = X[1] - v*h_N; 

            std::cout << X[1] << "  " << X[0] << endl;

            for(n=1; n<CP_iter-1; n++){

                /**** Fill S, H and Q for electronic problem ****/
                create_S(S, R_A, R_B);
                one_body_H(H, R_A, R_B);
                build_Q(Q, R_A, R_B);

                /**** Evolve the c vector of coefficients ****/
                double n_times = h_N/h;
                for(int k=0; k<n_times; k++){
                    two_body_F(Q, c, F);
                    gsl_matrix_add(F, H);
                    lambda = update_c(F, S, c, c_old);
                }

                if(response == "Y"){
                    /**** Print the coefficients to file ****/
                    for(int i=0; i<2*N; i++){
                        coeff_file << gsl_vector_get(c, i) <<  "    ";
                    }
                    coeff_file << endl;
                }

                /**** Compute the correct energy after the update of F ****/
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);
                E0 = compute_E0(F, H, c) + 1./X[n];

                /**** Fill S, H, Q and F for nuclear problem ****/
                create_dS_dX(S, R_A, R_B);
                one_body_dH_dX(H, R_A, R_B, X[n]);
                build_dQ_dX(Q, R_A, R_B, X[n]);
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);

                /**** Update X (also R_B.x is updated) ****/
                dE0_dX = compute_dE0_dX(F, H, c, X[n]);
                X[n + 1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[n], X[n-1]);

                /**** Save relevant quantities ****/
                T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
                T_e = Orbital_kinetic_en(R_A, R_B, c);
                std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << T_e << endl;
                X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << T_e << endl;
            }

            X_en_file.close();
            std::cout << "Trajectory " << to_string(i) << " has been written in HF_traj" << endl;
        }
        coeff_file.close();

        Print_density(c, X[CP_iter - 2]);
        

        /**** Print to screen the record of the simulation ****/
        std::cout << "Car-Parrinello Molecular Dynamics has been executed for " << CP_iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Electronic step h: " << h << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Electronic fictitious mass: " << m << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Electronic damping: " << gamma_el << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
    }





    





    /**** CONJUGATE GRADIENT - NUCLEI MOTION WITH CP ****/
    if(string(argv[1]) == "CG_CP"){

        /**** Start with little equilibration cycle with CPMD ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }
        
        /**** The first lambda_CP comes from the CP equilibration ****/
        double lambda_CP = lambda;

        /**** Threshold for conjugate gradient convergence ****/
        double eps = 0.001;

        /**** Open files to write the X, energies and coefficients ****/
        string X_en = "CG_CP_X_energies.txt";
        string coeff = "CG_CP_coeff.txt";
        ofstream X_en_file;
        ofstream coeff_file;
	    X_en_file.open(X_en, ios :: out | ios :: trunc);
        coeff_file.open(coeff, ios :: out | ios :: trunc);

        /**** Create Hessian, b, and increment of C ****/
        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);

        for(n=1; n<CP_iter-1; n++){

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Fill the Hessian and b (minus the gradient of E) ****/
            Get_Hessian_and_b(Hessian, b, Q, S, F, c);

            /**** Apply the CG routine to get the increment Delta_c ****/
            gsl_vector_set_all(Delta_c, 0.0);
            Conj_grad(Hessian, b, Delta_c, eps);

            /**** Copy c in c_old and then update c (adding Delta_c) ****/
            gsl_vector_memcpy(c_old, c);
            gsl_vector_add(c, Delta_c);

            /**** Update the c (in CP fashion) ****/
            lambda_CP = Get_lambda_CP(S, c, c_old, lambda_CP);
            norm = Get_norm_C_cg(S, c);

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X using the newly computed lambda_CP ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda_CP, dE0_dX, X[n], X[n-1]);
            std::cout << "  " << endl;
            X_en_file << X[n] << "    " << E0 << endl;
            std::cout << "lambda_CP: " << lambda_CP << endl;
            std::cout << "normalisation: " << norm << endl;

            /**** Prepare the S, H, Q and F for new electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
        }

        X_en_file.close();
        coeff_file.close();

        std::cout << "Car-Parrinello with conjugate gradient has been executed for " << CP_iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }








    /**** CG AND CP SUPERPOSITION ****/
    if(string(argv[1]) == "CG_CP_superposition"){
        int i, j;
        ofstream scal_prod_file;
        scal_prod_file.open("scal_prod.txt");

        /**** Open the file stream ****/ 
        string CG = "CG_CP_coeff.txt";
        string MD = "MD_CP_coeff.txt";
        ifstream CG_file;
        ifstream CP_file;
        CG_file.open(CG);
        CP_file.open(MD);

        /**** Check if opening a file failed ****/ 
        if (CG_file.fail() || CP_file.fail()) {
            std::cout << "Error opening the coefficients file" << endl;
            exit(1);
        }else{
            double c_CG[2*N], c_CP[2*N];
            double scal_prod = 0.;
            string line;

            /**** Read the lines of the two files ****/
            while(getline(CG_file, line)){
                for(i=0; i<2*N; i++){
                    CG_file >> c_CG[i];
                    CP_file >> c_CP[i];
                }

                /**** Compute the scalar product using S ****/
                scal_prod = 0.;
                for(i=0; i<2*N; i++){
                    for(j=0; j<2*N; j++){
                        scal_prod += gsl_matrix_get(S, i, j)*c_CP[i]*c_CG[j];
                    }
                }
                scal_prod_file << scal_prod << endl;
            }
            std::cout << "Coefficients files " << MD << " and " << CG << " successfully read" << endl;
            std::cout << "File scal_prod.txt successfully written" << endl;
        }
        CG_file.close();
        CP_file.close();
        scal_prod_file.close();
    }
    







    /**** EXCHANGE AND CORRELATION PART ****/
    if(string(argv[1]) == "EX_CORR"){
        /**** Little equilibration cycle with CPMD ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        
        /**** Create the XC matrix ****/
        create_Ex_Corr(V_xc, R_A, R_B, c, X[0]);

        std::cout << "Matrix V_xc: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(V_xc, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

        /**** Print density ****/
        Print_density(c, X[0]);

        /**** Print all the integrand functions defined over the (rho, z) plane ****/
        Print_integrand(0, 0, R_A, R_B, c, X[0]);
        Print_integrand(0, 1, R_A, R_B, c, X[0]);
        Print_integrand(0, 2, R_A, R_B, c, X[0]);
        Print_integrand(0, 3, R_A, R_B, c, X[0]);
        Print_integrand(1, 1, R_A, R_B, c, X[0]);
        Print_integrand(1, 2, R_A, R_B, c, X[0]);
        Print_integrand(1, 3, R_A, R_B, c, X[0]);
        Print_integrand(2, 2, R_A, R_B, c, X[0]);
        Print_integrand(2, 3, R_A, R_B, c, X[0]);
        Print_integrand(3, 3, R_A, R_B, c, X[0]);
    }






    /**** EXCHANGE AND CORRELATION DERIVATIVES PART ****/
    if(string(argv[1]) == "XC_DER"){
        /**** Little equilibration cycle with CPMD ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Print density ****/
        Print_density_derivative(c, X[0]);

        /**** Create the derivative of XC matrix ****/
        create_dVxc_dX(dVxc_dX, R_A, R_B, c, X[0]);

        std::cout << "Matrix dVxc_dX: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(dVxc_dX, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

    }










    /**** CPMD WITH DFT ****/
    if(string(argv[1]) == "CPMD_DFT"){

        /**** Choose the number of steps (must be less than CP_iter) ****/
        int N_steps;
        std::cout << "Enter the number of CPMD_DFT steps: ";
        std::cin >> N_steps;

        /**** If a previous trajectory has been produced, then append the output ****/
        string X_en = "MD_CP_DFT_X_energies.txt";
        string coeff = "MD_CP_DFT_coeff.txt";
        ifstream X_energies_file;
        ifstream coefficients_file;
        X_energies_file.open(X_en);
        coefficients_file.open(coeff);
        double read_coeff[2*N];
        double last_X;

        /**** Check if opening a file failed ****/ 
        if (coefficients_file.fail() || X_energies_file.fail()) {
            std::cout << "Coefficients file or the X_energies file NOT FOUND" << endl;
            std::cout << "*********** Starting a NEW simulation ***********" << endl;
        }else{
            string line_coeff, line_X;

            /**** Read the coefficient file until the last line ****/
            while(getline(coefficients_file, line_coeff)){
                for(int i=0; i<2*N; i++){
                    coefficients_file >> read_coeff[i];
                }
            }

            /**** Print the coefficients to screen ****/
            std::cout << "****** Warning: starting from a given set of coefficients c ****" << endl;
            std::cout << "Coefficients read from file: " << endl;
            for(int i=0; i<2*N; i++){    
                std::cout << read_coeff[i] << endl;
            }
            
            /**** Save the read coefficients in the gsl_vector ****/
            for(int i=0; i<2*N; i++){
                gsl_vector_set(c, i, read_coeff[i]);
            }

            /**** To restart the c evolution on the run ****/
            gsl_vector_memcpy(c_old, c);

            /**** Save the lines as strings in contents_X ****/
            vector<string> contents_X;
            while(!X_energies_file.eof()){
                getline(X_energies_file, line_X);
                contents_X.push_back(line_X);
            }

            /**** Save the last and the second last lines ****/
            int number_of_lines = contents_X.size();
            std::cout << "Second last line: " << contents_X.at(number_of_lines-3) << endl;
            std::cout << "Last line: " << contents_X.at(number_of_lines-2) << endl;
            /**** For the user: leave an empty line at the end of the file ****/

            /**** Create a string object with a stream ****/
            stringstream s1(contents_X.at(number_of_lines - 3));
            stringstream s2(contents_X.at(number_of_lines - 2));
            string record_1, record_2;
            s1 >> record_1;
            s2 >> record_2;
            X[0] = atof(record_1.c_str());
            X[1] = atof(record_2.c_str());
            R_B.x = X[1]; 

            /**** Remove the unnecessary last line from the X_en file ****/
            Remove_last_line(X_en);

        }/**** End of the else statement ****/

        /**** Create the files for storing the energies and the coefficients ****/
        ofstream X_en_file;
        ofstream coeff_file;

        /**** If the X_en_file already exists, then append the output ****/
        if(X_en_file.fail() || coeff_file.fail()){
            X_en_file.open(X_en, ios :: out | ios :: trunc);
            coeff_file.open(coeff, ios :: out | ios :: trunc);
        }else{
            X_en_file.open(X_en, ios_base::app);
            coeff_file.open(coeff, ios_base::app);
        }

        std::cout << "X" << "      " << "Energy (H)" << endl;

        for(n=1; n<=N_steps; n++){

            /**** Fill S, H and Q for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            string s = "V_xc";
            /**** Evolve the c vector of coefficients ****/
            double n_times = h_N/h;
            for(int k=0; k<n_times; k++){
                two_body_F(Q, c, F);
                gsl_matrix_scale(F, 1. + a_x);

                /**** The XC part has to be recomputed because it strictly depends on the vector c ****/
                Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[n], s);
                gsl_matrix_add(F, H);
                gsl_matrix_add(F, V_xc);
                lambda = update_c(F, S, c, c_old);
            }

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);

            /**** Adding the exchange and correlation ****/
            create_Ex_Corr(V_xc, R_A, R_B, c, X[n]);
            gsl_matrix_add(F, V_xc);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            create_dVxc_dX(dVxc_dX, R_A, R_B, c, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);
            gsl_matrix_add(F, dVxc_dX);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[n], X[n-1]);
            std::cout << "  " << endl;

            /**** Print to X_en file and to screen ****/
            X_en_file << X[n] << "    " << E0 << endl;
            std::cout << X[n] << "    " << E0 << endl;
        }

        /**** Print the X[n+1] to restart the simulation ****/
        X_en_file << X[N_steps + 1] << endl;  

        X_en_file.close();
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Car-Parrinello Molecular Dynamics has been executed for " << N_steps << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Electronic step h: " << h << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Electronic fictitious mass: " << m << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Electronic damping: " << gamma_el << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        xc_func_type func;
        xc_func_init(&func, FUNCTIONAL_C, XC_UNPOLARIZED);
        std::cout << "XC functional employed: " << func.info->name << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }









    /**** ADAPTIVE INTEGRATION ****/
    if(string(argv[1]) == "ADAPTIVE"){
        /**** Little equilibration cycle with CPMD ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Create the derivative of XC matrix ****/
        string s = "V_xc";
        Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[0], s);

        std::cout << "Matrix V_xc: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(V_xc, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

    }







    /**** XC DERIVATIVES (ADAPTIVE) ****/
    if(string(argv[1]) == "XC_DER_ADAPTIVE"){
        /**** Little equilibration cycle with CPMD ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Create the derivative of XC matrix ****/
        string s = "dVxc_dX";
        Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[0], s);

        std::cout << "Matrix dVxc_dX: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(dVxc_dX, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

    }


    /**** End of the main function ****/
}




