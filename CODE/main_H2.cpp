#include "definitions.h"
#include "xc.h"

using namespace std;


int main (int argc, char *argv[]){

    /**** Create space for matrices and vectors ****/
	int n;
	R R_A, R_B;
    double Q[2*N][2*N][2*N][2*N];
    double X[iter];
    double E_1s=0., lambda=0., E0 = 0., E_old = 0., dE0_dX = 0., T_N = 0., E_ee = 0.0, f_kin_en = 0.0, E_one_body = 0.0, E_xc = 0.0;
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
    gsl_matrix *dS_dX = gsl_matrix_alloc(2*N, 2*N);

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
    /**** BOMD_HF ****/
    /**** BOMD_DFT ****/
    /**** Evolve_coefficients (CP electronic energy evolution in single step) ****/ 
    /**** CPMD_HF ****/
    /**** CPMD_DFT ****/
    /**** CGMD_HF ****/
    /**** CGMD_DFT ****/
    /**** EX_CORR (prints the V_xc integrands) ****/ 
    /**** CPMD_HF_TRAJECTORIES (computes multiple trajectories from random initial conditions) ****/





    
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

        /**** Create the files for storing the energies and the coefficients ****/
        string X_en = "BOMD_HF_X_energies.txt";
        string coeff = "BOMD_HF_coeff.txt";
        ofstream X_en_file;
        ofstream coeff_file;
        X_en_file.open(X_en, ios :: out | ios :: trunc);
        coeff_file.open(coeff, ios :: out | ios :: trunc);

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "Force" << endl;

        for(n=1; n<iter - 1; n++){

            /**** Fill H and Q for electronic problem ****/
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            /**** Create S and the V matrix, needed for the SC procedure ****/
            create_S(S, R_A, R_B);
            gsl_matrix_memcpy(S_auxiliary, S);
            diag_S(S_auxiliary, V);

            /**** SC while loop: here we get the c coefficients ****/
            E0 = 1.0, E_old = 0.0;
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
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E0 = compute_E0(F, H, c) + 1./X[n];
            

            /**** Fill H, Q and F for nuclear problem ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Compute energy gradient & update the X with shake ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, E_1s, dE0_dX, X[n], X[n - 1], "BO_CG");

            /**** Print relevant quantities ****/
            T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << dE0_dX << endl;
        }

        X_en_file.close();
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Born-Oppenheimer Molecular Dynamics has been executed for " << iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }













    /**** BOMD WITH DFT ****/
    if(string(argv[1]) == "BOMD_DFT"){

        /**** Choose the number of steps (must be less than iter=2000) ****/
        int N_steps;
        std::cout << "Enter the number of BOMD_DFT steps: ";
        std::cin >> N_steps;

        /**** If a previous trajectory has been produced, then append the output ****/
        string X_en = "BOMD_DFT_X_energies.txt";
        string coeff = "BOMD_DFT_coeff.txt";
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

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "E_xc" << " || " << "Force" << endl;

        for(n=1; n<=N_steps; n++){

            /**** Fill S and V ****/
            create_S(S, R_A, R_B);
            gsl_matrix_memcpy(S_auxiliary, S);
            diag_S(S_auxiliary, V);

            /**** Fill H and Q for electronic problem ****/
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            string s = "V_xc";
            /**** Evolve the c vector of coefficients ****/
            E0 = 0.0, E_old = 1.0;
            while(fabs(E0 - E_old) > 1E-7){

                /**** Build Fock matrix with the added exchange part ****/
                two_body_F(Q, c, F);
                gsl_matrix_scale(F, 1. + a_x);
                gsl_matrix_add(F, H);

                /**** Adding the exchange and correlation ****/
                string s = "V_xc";
                Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[n], s);
                gsl_matrix_add(F, V_xc);

                /**** tmp_F is needed beacuse solve_FC_eSC destroys lower part of F ****/
                gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
                gsl_matrix_memcpy(tmp_F, F);
                E_1s = solve_FC_eSC(tmp_F, V, U);

                /**** Compute the energy ****/
                E_old = E0;
                E0 = compute_E0(F, H, c) + 1./X[0];

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
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);

            /**** Adding the exchange and correlation ****/
            create_Ex_Corr(V_xc, R_A, R_B, c, X[n]);
            gsl_matrix_add(F, V_xc);

            /**** Computing el-el, one-body, XC energy & DFT ground state ****/
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E_xc = XC_energy(c, V_xc);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(dS_dX, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            create_dVxc_dX(dVxc_dX, R_A, R_B, c, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);
            gsl_matrix_add(F, dVxc_dX);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, dS_dX, &R_B, E_1s, dE0_dX, X[n], X[n - 1], "BO_CG");

            /**** Print to X_en file and to screen ****/
            T_N = Nuclear_kinetic_en(X[n], X[n-1]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << dE0_dX << endl;
        }

        /**** Print the X[n+1] to restart the simulation ****/
        X_en_file << X[N_steps + 1] << endl;  

        X_en_file.close();
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Born-Oppenheimer Molecular Dynamics has been executed for " << N_steps << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        xc_func_type func;
        xc_func_init(&func, FUNCTIONAL_C, XC_UNPOLARIZED);
        std::cout << "XC functional employed: " << func.info->name << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }













    /**** EVOLUTION OF THE COEFFICIENTS C(t) AT FIXED X ****/
    if(string(argv[1]) == "Evolve_coefficients"){
        ofstream myfile;
	    myfile.open("Energies_C_evolution.txt", ios :: out | ios :: trunc);

        /**** MD cycle: compute new c vector and keep F updated ****/
        int steps = h_N/h;
        for(n=0; n<steps; n++){
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

        /**** Create the files for storing the energies and the coefficients ****/
        string coeff = "CPMD_HF_coeff.txt";
        string X_en = "CPMD_HF_X_energies.txt";
        ofstream coeff_file;
        ofstream X_en_file;
        coeff_file.open(coeff, ios :: out | ios :: trunc);
        X_en_file.open(X_en, ios :: out | ios :: trunc);

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "Fictitious kinetic energy"  << " || " << "Force" << endl;           

        for(n=1; n<iter-1; n++){
            /**** Keep c_old updated ****/
            gsl_vector_memcpy(c_old, c);

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

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(dS_dX, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, dS_dX, &R_B, lambda, dE0_dX, X[n], X[n-1], "CP");

            /**** Compute fictitious kinetic energy ****/
            f_kin_en = Fictitious_kin_energy(c, c_old, S, dS_dX, X[n + 1], X[n], R_A, R_B);

            /**** Save relevant quantities ****/
            T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << f_kin_en << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << f_kin_en << "    " << dE0_dX << endl;
        

        }
        coeff_file.close();
        X_en_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Car-Parrinello Molecular Dynamics has been executed for " << iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Electronic step h: " << h << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Electronic fictitious mass: " << m << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Electronic damping: " << gamma_el << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
    }












    /**** CPMD WITH DFT ****/
    if(string(argv[1]) == "CPMD_DFT"){

        /**** Choose the number of steps (must be less than iter=2000) ****/
        int N_steps;
        std::cout << "Enter the number of CPMD_DFT steps: ";
        std::cin >> N_steps;

        /**** If a previous trajectory has been produced, then append the output ****/
        string X_en = "CPMD_DFT_X_energies.txt";
        string coeff = "CPMD_DFT_coeff.txt";
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

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "E_xc" << " || " << "Fic. kin. energy" << " || " << "Force" << endl;

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
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E_xc = XC_energy(c, V_xc);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(dS_dX, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            create_dVxc_dX(dVxc_dX, R_A, R_B, c, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);
            gsl_matrix_add(F, dVxc_dX);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, dS_dX, &R_B, lambda, dE0_dX, X[n], X[n-1], "CP");

            /**** Compute fictitious kinetic energy ****/
            f_kin_en = Fictitious_kin_energy(c, c_old, S, dS_dX, X[n + 1], X[n], R_A, R_B);

            /**** Print to X_en file and to screen ****/
            T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << f_kin_en << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << f_kin_en << "    " << dE0_dX << endl;
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














    /**** CONJUGATE GRADIENT MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "CGMD_HF"){

        /**** Small equilibration with Car-Parrinello ****/
        create_S(S, R_A, R_B);
        one_body_H(H, R_A, R_B);
        build_Q(Q, R_A, R_B);
        double n_times = h_N/h;
        for(int k=0; k<n_times; k++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Threshold for conjugate gradient convergence ****/
        double eps = 1e-4;

        /**** Open files to write the X, energies and coefficients ****/
        string X_en = "CGMD_HF_X_energies.txt";
        string coeff = "CGMD_HF_coeff.txt";
        ofstream X_en_file;
        ofstream coeff_file;
	    X_en_file.open(X_en, ios :: out | ios :: trunc);
        coeff_file.open(coeff, ios :: out | ios :: trunc);

        /**** Create Hessian, b, and increment of C ****/
        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "Force" << endl;

        for(n=1; n<iter; n++){

            /**** Prepare the S, H, Q and F for new electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Fill the Hessian and b (minus the gradient of E) ****/
            Get_Hessian_and_b(Hessian, b, Q, S, F, c);

            /**** Apply the CG routine to get the increment Delta_c ****/
            gsl_vector_set_all(Delta_c, 0.0);
            Conj_grad(Hessian, b, Delta_c, eps);

            /**** Update c (adding Delta_c) ****/
            gsl_vector_add(c, Delta_c);

            /**** Update F and compute the correct energy after the update of c ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Compute lambda needed in the X evolution ****/
            lambda = Get_lambda_CG(F, c);

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Fill H, Q and F for nuclear problem (for the calculus of dE0/dX) ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X using shake ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[n], X[n - 1], "BO_CG");

            /**** Print relevant quantities ****/
            T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << dE0_dX << endl;
            
        }

        X_en_file.close();
        coeff_file.close();

        std::cout << "Conjugate Gradient Molecular Dynamics has been executed for " << iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }









    /**** CPMD WITH DFT ****/
    if(string(argv[1]) == "CGMD_DFT"){

        /**** Small equilibration with Car-Parrinello ****/
        create_S(S, R_A, R_B);
        one_body_H(H, R_A, R_B);
        build_Q(Q, R_A, R_B);
        double n_times = h_N/h;
        for(int k=0; k<n_times; k++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Choose the number of steps (must be less than iter=2000) ****/
        int N_steps;
        std::cout << "Enter the number of CGMD_DFT steps: ";
        std::cin >> N_steps;

        /**** If a previous trajectory has been produced, then append the output ****/
        string X_en = "CGMD_DFT_X_energies.txt";
        string coeff = "CGMD_DFT_coeff.txt";
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

        /**** Create Hessian, b, and increment of C ****/
        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "E_xc" << " || " << "Force" << endl;

        for(n=1; n<=N_steps; n++){

            /**** Fill S, H, Q and F for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Add the XC contribution ****/
            string s = "V_xc";
            Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[n], s);
            gsl_matrix_add(F, V_xc);
            
            /**** Fill the Hessian and b (minus the gradient of E) ****/
            Get_Hessian_and_b(Hessian, b, Q, S, F, c);

            /**** Apply the CG routine to get the increment Delta_c ****/
            gsl_vector_set_all(Delta_c, 0.0);
            Conj_grad(Hessian, b, Delta_c, 0.001);

            /**** Update c (adding Delta_c) ****/
            gsl_vector_add(c, Delta_c);

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

            /**** Computing el-el, one-body, XC energy & DFT ground state ****/
            E_ee = Electron_electron_en(c, F);
            E_one_body = One_body(c, H);
            E_xc = XC_energy(c, V_xc);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Compute lambda needed in the X evolution ****/
            lambda = Get_lambda_CG(F, c);

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(dS_dX, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            create_dVxc_dX(dVxc_dX, R_A, R_B, c, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_scale(F, 1. + a_x);
            gsl_matrix_add(F, H);
            gsl_matrix_add(F, dVxc_dX);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, dS_dX, &R_B, lambda, dE0_dX, X[n], X[n - 1], "BO_CG");

            /**** Print to X_en file and to screen ****/
            T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
            std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << dE0_dX << endl;
            X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << E_xc << "    " << dE0_dX << endl;
        }

        /**** Print the X[n+1] to restart the simulation ****/
        X_en_file << X[N_steps + 1] << endl;  

        X_en_file.close();
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Conjugate Gradient Molecular Dynamics has been executed for " << N_steps << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
        xc_func_type func;
        xc_func_init(&func, FUNCTIONAL_C, XC_UNPOLARIZED);
        std::cout << "XC functional employed: " << func.info->name << endl;
        std::cout << "  " << endl;
        std::cout << "The X coordinate and the energies have been written to " << X_en << endl;
        std::cout << "The C coefficients have been written to " << coeff << endl;
    }










    /**** SUPERPOSITION ****/
    if(string(argv[1]) == "Scal_prod"){
        int i, j;
        double prod1, prod2;
        string set1, set2, record_1, record_2, line_X1, line_X2;
        gsl_vector *c1 = gsl_vector_alloc(2*N);
        gsl_vector *c2 = gsl_vector_alloc(2*N);
        gsl_vector *Sc = gsl_vector_alloc(2*N);

        std::cout << "Select two sets of coefficients between CP, BO and CG" << endl;
        std::cout << "First set: ";
        std::cin >> set1;
        std::cout << "Second set: ";
        std::cin >> set2;

        /**** Name of the output file ****/
        ofstream scal_prod_file;
        scal_prod_file.open("Scal_prod_" + set1 + "_" + set2 + ".txt");

        /**** Open the file stream ****/ 
        string SET1 = set1 + "MD_HF_coeff.txt";
        string SET2 = set2 + "MD_HF_coeff.txt";
        string X_EN1 = set1 + "MD_HF_X_energies.txt";
        string X_EN2 = set2 + "MD_HF_X_energies.txt";

        /**** The idea is to get the whole set of X first and then do a cycle ****/
        /**** for each line of the of the coefficients file ****/
        ifstream SET1_file;
        ifstream SET2_file;
        ifstream X_EN1_file;
        ifstream X_EN2_file;
        SET1_file.open(SET1);
        SET2_file.open(SET2); 
        X_EN1_file.open(X_EN1);
        X_EN2_file.open(X_EN2);  

        /**** Make the X_EN1 file an input stream ****/
        vector<string> contents_X1;
        while(!X_EN1_file.eof()){
            getline(X_EN1_file, line_X1);

            /**** Get the lines until the end of the file ****/
            contents_X1.push_back(line_X1);
        }

        /**** Repeat for the second file ****/
        vector<string> contents_X2;
        while(!X_EN2_file.eof()){
            getline(X_EN2_file, line_X2);
            contents_X2.push_back(line_X2);
        }

        /**** Check if opening a coefficients file failed ****/ 
        if (SET1_file.fail() || SET2_file.fail()) {
            std::cout << "Error opening the coeff or X_energies files" << endl;
            exit(1);
        }else{
            double c_SET1[2*N], c_SET2[2*N];
            string line;
            int counts;

            /**** Read the lines of the two coefficients files ****/
            while(getline(SET1_file, line)){

                /**** Capture the lines *****/
                for(i=0; i<2*N; i++){
                    SET1_file >> c_SET1[i];
                    SET2_file >> c_SET2[i];
                }

                /**** Set the vectors ****/
                for(i=0; i<2*N; i++){
                    gsl_vector_set(c1, i, c_SET1[i]);
                    gsl_vector_set(c2, i, c_SET2[i]);
                }

                /**** Account for the number of the current line ****/
                counts++;

                /**** Set R_B for the current line ****/
                stringstream s1(contents_X1.at(counts - 1));
                s1 >> record_1;
                R_B.x = atof(record_1.c_str());

                /**** Fill the S matrix with the current value of X ****/
                create_S(S, R_A, R_B);

                /**** Compute the first scalar product C_2^T * S_1 * C_1 ****/
                gsl_blas_dgemv(CblasNoTrans, 1., S, c1, 0., Sc);
                gsl_blas_ddot(c2, Sc, &prod1);

                /**** Repeat everything ****/
                stringstream s2(contents_X2.at(counts - 1));
                s2 >> record_2;
                R_B.x = atof(record_2.c_str());
                create_S(S, R_A, R_B);
                gsl_blas_dgemv(CblasNoTrans, 1., S, c1, 0., Sc);
                gsl_blas_ddot(c2, Sc, &prod2);
                scal_prod_file << fabs(prod1) << "    " << fabs(prod2) << endl; 
            }

            std::cout << "Coefficients files " << SET1 << " and " << SET2 << " successfully read" << endl;
            std::cout << "File .txt with scalar products successfully written" << endl;
        }
        SET1_file.close();
        SET2_file.close();
        X_EN1_file.close();
        X_EN2_file.close();
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

        std::cout << "Matrix V_xc computed with STANDARD SIMPSON: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(V_xc, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

        std::cout << "  " << endl;
        std::cout << "*************************************************************************" << endl;
        std::cout << "  " << endl;
        /**** Create the XC matrix ****/
        Adaptive_Ex_Corr(V_xc, dVxc_dX, R_A, R_B, c, X[0], "V_xc");

        std::cout << "Matrix V_xc computed with ADAPTIVE SIMPSON: " << endl;
        for(int k=0; k<2*N; k++){
            for(int j=0; j<2*N; j++){
                std::cout << gsl_matrix_get(V_xc, k, j) << "   ";
            }
            std::cout << "  " << endl;
        }

        std::cout << "  " << endl;
        
        /**** Print density ****/
        Print_density(c, X[0]);
        Print_density_derivative(c, X[0]);
        Print_density_dz(c, X[0], R_A, R_B);

        /**** Print all the integrand functions defined over the (rho, z) plane ****/
        Print_integrand(0, 0, R_A, R_B, c, X[0]);
        Print_integrand(0, 1, R_A, R_B, c, X[0]);
        Print_integrand(0, 2, R_A, R_B, c, X[0]);
        Print_integrand(0, 3, R_A, R_B, c, X[0]);
        Print_integrand(1, 1, R_A, R_B, c, X[0]);
        Print_integrand(1, 2, R_A, R_B, c, X[0]);
        Print_integrand(1, 3, R_A, R_B, c, X[0]);
        Print_integrand(2, 2, R_A, R_B, c, X[0]);
        Print_integrand(2, 3, R_A, R_A, c, X[0]);
        Print_integrand(3, 3, R_A, R_B, c, X[0]);
        Print_integrand_dX(0, 0, R_A, R_B, c, X[0]);
        Print_integrand_dX(0, 0, R_A, R_B, c, X[0]);
        Print_integrand_dX(0, 1, R_A, R_B, c, X[0]);
        Print_integrand_dX(0, 2, R_A, R_B, c, X[0]);
        Print_integrand_dX(0, 3, R_A, R_B, c, X[0]);
        Print_integrand_dX(1, 1, R_A, R_B, c, X[0]);
        Print_integrand_dX(1, 2, R_A, R_B, c, X[0]);
        Print_integrand_dX(1, 3, R_A, R_B, c, X[0]);
        Print_integrand_dX(2, 2, R_A, R_B, c, X[0]);
        Print_integrand_dX(2, 3, R_A, R_A, c, X[0]);
        Print_integrand_dX(3, 3, R_A, R_B, c, X[0]);
    }

















    /**** EXTRA : CAR PARRINELLO MULTIPLE TRAJECTORIES ****/
    if(string(argv[1]) == "CPMD_HF_TRAJECTORIES"){

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

        std::cout << "X" << " || " << "Electronic energy" << " || " << "Nuclear kin. energy" << " || " << "E_ee energy" << " || " << "One-body energy" << " || " << "Fictitious kinetic energy"  << " || " << "Force" << endl;

        /**** Seed for random positions ****/
        srand((unsigned) time(NULL));

        for(int i=0; i<N_traj; i++){

            /**** Define the output trajectory files ****/
            string X_en = "outputs/HF_traj/CPMD_HF_" + to_string(i) + ".txt";
            ofstream X_en_file;
            X_en_file.open(X_en, ios :: out | ios :: trunc);
            
            /**** Generate random initial positions ****/
            X[0] = 1.2;
            X[1] = X[0];
            R_B.x = X[1];            

            for(n=1; n<iter-1; n++){
                /**** Keep c_old updated ****/
                gsl_vector_memcpy(c_old, c);

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
                E_ee = Electron_electron_en(c, F);
                E_one_body = One_body(c, H);
                E0 = compute_E0(F, H, c) + 1./X[n];

                /**** Fill S, H, Q and F for nuclear problem ****/
                create_dS_dX(dS_dX, R_A, R_B);
                one_body_dH_dX(H, R_A, R_B, X[n]);
                build_dQ_dX(Q, R_A, R_B, X[n]);
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);

                /**** Update X (also R_B.x is updated) ****/
                dE0_dX = compute_dE0_dX(F, H, c, X[n]);
                X[n + 1] = evolve_X(c, dS_dX, &R_B, lambda, dE0_dX, X[n], X[n-1], "CP");

                /**** Compute fictitious kinetic energy ****/
                f_kin_en = Fictitious_kin_energy(c, c_old, S, dS_dX, X[n + 1], X[n], R_A, R_B);

                /**** Save relevant quantities ****/
                T_N = Nuclear_kinetic_en(X[n + 1], X[n]);
                std::cout << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << f_kin_en << "    " << dE0_dX << endl;
                X_en_file << X[n] << "    " << E0 << "    " << T_N << "    " << E_ee << "    " << E_one_body << "    " << f_kin_en << "    " << dE0_dX << endl;
            }

            X_en_file.close();
            std::cout << "Trajectory " << to_string(i) << " has been written in HF_traj" << endl;
        }
        coeff_file.close();

        /**** Print to screen the record of the simulation ****/
        std::cout << "Car-Parrinello Molecular Dynamics has been executed for " << iter << " steps" << endl;
        std::cout << "  " << endl;
        std::cout << "Parameters of the simulation: " << endl;
        std::cout << "Electronic step h: " << h << endl;
        std::cout << "Nuclear step h_N: " << h_N << endl;
        std::cout << "Electronic fictitious mass: " << m << endl;
        std::cout << "Nuclear mass: " << M_N << endl;
        std::cout << "Electronic damping: " << gamma_el << endl;
        std::cout << "Nuclear damping: " << gamma_N << endl;
    }


    /**** End of the main function ****/
}




