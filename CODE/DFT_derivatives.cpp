/********** EXCHANGE AND CORRELATION DERIVATIVES FOR CPMD ********/

#include "definitions.h" 
#include "xc.h"

using namespace std;


double density_derivative(double rho, double z, gsl_vector *c, double X){
    double dn_dX = 0;
    int p, q;
    double c_p, c_q, K, exp_factor_p, exp_factor_q;

    /**** Return the value of dn/dX, which is needed in the nuclear equations of CPMD ****/
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            c_p = gsl_vector_get(c, p);
            c_q = gsl_vector_get(c, q);
            dn_dX += 2.0*(a[p] + a[q])*(z - X)*c_p*c_q*(exp(-(a[p] + a[q])*(rho*rho + (z - X)*(z - X))));

            /**** Note that K is a function of X. So here we have the derivative of a product ****/
            K = exp(-a[p]*a[q]*X*X/(a[p] + a[q]));
            exp_factor_p = exp(-(a[p] + a[q])*(rho*rho + (z - a[p]*X/(a[p] + a[q]))*(z - a[p]*X/(a[p] + a[q]))));
            exp_factor_q = exp(-(a[p] + a[q])*(rho*rho + (z - a[q]*X/(a[p] + a[q]))*(z - a[q]*X/(a[p] + a[q]))));
            dn_dX -= 2.0*a[p]*a[q]*X/(a[p] + a[q])*c_p*c_q*K*(exp_factor_p + exp_factor_q);
            dn_dX += 2.0*a[p]*(z - a[p]*X/(a[p] + a[q]))*c_p*c_q*K*exp_factor_p;
            dn_dX += 2.0*a[q]*(z - a[q]*X/(a[p] + a[q]))*c_p*c_q*K*exp_factor_q;
        }
    }

    return dn_dX;
}


double dchi_p_chi_q_dX(double rho, double z, double alpha, double beta, R R_A, R R_B, double X){
    double derivative = 0.0, exp_factor_p, exp_factor_q;

    /**** Compute the factor K (which depends implicitly on X) ****/
    double prefactor = K(alpha, beta, R_A, R_B);

    /**** Three cases ****/
    if(R_A.x == R_B.x && R_B.x == 0.0){
        derivative = 0.0;
    }
    if(R_A.x == R_B.x && R_B.x != 0.0){
        derivative = 2.0*prefactor*(alpha + beta)*(z - X)*exp(-(alpha + beta)*(rho*rho + (z - X)*(z - X)));
    }
    if(R_A.x != R_B.x && R_B.x == 0.0){
        exp_factor_p = exp(-(alpha + beta)*(rho*rho + (z - alpha*X/(alpha + beta))*(z - alpha*X/(alpha + beta))));
        derivative = 2.0*prefactor*alpha*(z - alpha*X/(alpha + beta))*exp_factor_p;
        derivative -= 2.0*prefactor*alpha*beta*X/(alpha + beta)*exp_factor_p;
    }
    if(R_A.x != R_B.x && R_B.x == 1.0){
        exp_factor_q = exp(-(alpha + beta)*(rho*rho + (z - beta*X/(alpha + beta))*(z - beta*X/(alpha + beta))));
        derivative = 2.0*prefactor*beta*(z - beta*X/(alpha + beta))*exp_factor_q;
        derivative -= 2.0*prefactor*alpha*beta*X/(alpha + beta)*exp_factor_q;
    }
    return derivative;
}




double Integrand_dX(double rho, double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Compute the density ****/
    double n = density(rho, z, c, X);

    /**** Compute the exchange part (namely, dv_x/dn) ****/
    double v_x = 0.0, v_c = 0.0, dvx_dn = 0.0, dvc_dn = 0.0;
    xc_func_type functional_x;
    xc_func_init(&functional_x, FUNCTIONAL_X, XC_UNPOLARIZED);
    xc_lda_vxc(&functional_x, 1, &n, &v_x);
    xc_lda_fxc(&functional_x, 1, &n, &dvx_dn);
    xc_func_end(&functional_x);

    /**** Compute the correlation part (namely, dv_c/dn) ****/
    xc_func_type functional_c;
    xc_func_init(&functional_c, FUNCTIONAL_C, XC_UNPOLARIZED);
    xc_lda_vxc(&functional_c, 1, &n, &v_c);
    xc_lda_fxc(&functional_c, 1, &n, &dvc_dn);
    xc_func_end(&functional_c);

    /**** These are used in the integral. Note the a_x ****/
    double v_xc = a_x*v_x + v_c;
    double dvxc_dn = a_x*dvx_dn + dvc_dn;

    /**** Compute prefactor and R_C ****/
    double prefactor = K(alpha, beta, R_A, R_B);
    R R_C = R_weighted(alpha, beta, R_A, R_B);

    /**** Sum the two contributions to the integral ****/
    double integrand = prefactor*exp(-(alpha + beta)*(rho*rho + (z - R_C.x)*(z - R_C.x)))*dvxc_dn*rho*density_derivative(rho, z, c, X);
    integrand += dchi_p_chi_q_dX(rho, z, alpha, beta, R_A, R_B, X)*v_xc*rho;

    return integrand;
}




double Simpson_rho_dX(double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Extremes of integration ****/
    double m = alpha + beta;
    double a = 0.0;
    double b = 2.0 + 1/m;
    
    /**** First and last points ****/
    double f_a = Integrand_dX(a, z, alpha, beta, R_A, R_B, c, X);
    double f_b = Integrand_dX(b, z, alpha, beta, R_A, R_B, c, X);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 100;
    double drho = (b-a)/N_mesh;

    /**** Simpson integration ****/
    double rho = 0.0, Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
        rho = i*drho;
        f_i = Integrand_dX(rho, z, alpha, beta, R_A, R_B, c, X);
        f_i_plus_1 = Integrand_dX(rho + drho, z, alpha, beta, R_A, R_B, c, X);
        Simpson_sum = Simpson_sum + (4.*f_i + 2.*f_i_plus_1);
    }

    Simpson_sum = (Simpson_sum + f_a + f_b)*(drho/3.);
    return Simpson_sum;
}




double Simpson_z_dX(double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Extremes of integration ****/
    double m = alpha + beta;
    double a = -1.5 - 1/m;
    double b = 2.5 + 1/m;
    
    /**** First and last points ****/
    double f_a = Simpson_rho_dX(a, alpha, beta, R_A, R_B, c, X);
    double f_b = Simpson_rho_dX(b, alpha, beta, R_A, R_B, c, X);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 100;
    double dz = (b-a)/N_mesh;

    /**** Simpson integration ****/
    double z = 0.0, Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
    	
    	/**** The z coordinate starts at a finite value ****/
        z = a + i*dz;
        f_i = Simpson_rho_dX(z, alpha, beta, R_A, R_B, c, X);
        f_i_plus_1 = Simpson_rho_dX(z + dz, alpha, beta, R_A, R_B, c, X);
        Simpson_sum = Simpson_sum + (4.*f_i + 2.*f_i_plus_1);
    }

    Simpson_sum = (Simpson_sum + f_a + f_b)*(dz/3.);
    return Simpson_sum;
}



void create_dVxc_dX(gsl_matrix *dVxc_dX, R R_A, R R_B, gsl_vector *c, double X){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){ 

            /**** Factor 2*pi due to the integration over theta ****/
			val1 = 2.0*pi*Simpson_z_dX(a[p], a[q], R_A, R_A, c, X);
	        gsl_matrix_set(dVxc_dX, p, q, val1);
            gsl_matrix_set(dVxc_dX, q, p, val1);
            gsl_matrix_set(dVxc_dX, p + N, q + N, val1);
            gsl_matrix_set(dVxc_dX, q + N, p + N, val1);

            val2 = 2.0*pi*Simpson_z_dX(a[p], a[q], R_A, R_B, c, X);
	        gsl_matrix_set(dVxc_dX, p, q + N, val2);
            gsl_matrix_set(dVxc_dX, q, p + N, val2);
            gsl_matrix_set(dVxc_dX, p + N, q, val2);
            gsl_matrix_set(dVxc_dX, q + N, p, val2);

		}
	}
}



void Print_density_derivative(gsl_vector *c, double X){
    ofstream myfile;
	myfile.open("Density_derivative.txt", ios :: out | ios :: trunc);

    int i, j, N_mesh = 100;
    double z, rho;

    /**** Define the mesh for rho ****/
    double rho_a = 0.0, rho_b = 2.0;
    double drho = (rho_b - rho_a)/N_mesh;

    /**** Define the mesh for z ****/
    double z_a = -2.0, z_b = 3.0;
    double dz = (z_b - z_a)/N_mesh;

    for(i=0; i<N_mesh; i++){
        rho = i*drho;
        for(j=0; j<N_mesh; j++){
            z = z_a + j*dz;
            myfile << density_derivative(rho, z, c, X) << "       ";
        }
        myfile << "  " << endl;
    }

    myfile.close();
    std::cout << "Density derivative written successfully in a txt file" << endl;
}



void Print_integrand_dX(int p, int q, R R_A, R R_B, gsl_vector *c, double X){
    string name = "Integrands/Integrand_dX";
    string underscore = "_";
    string txt = ".txt";
    string complete_name = name + to_string(p) + underscore + to_string(q) + txt;
    ofstream myfile;
	myfile.open(complete_name, ios :: out | ios :: trunc);

    int i, j, N_mesh = 100;
    double z, rho;

    /**** Define the mesh for rho ****/
    double m = a[p] + a[q];
    double rho_a = 0.0, rho_b = 2.0 + 1/m;
    double drho = (rho_b - rho_a)/N_mesh;

    /**** Define the mesh for z ****/
    double z_a = -1.5 - 1/m;
    double z_b = 2.5 + 1/m;
    double dz = (z_b - z_a)/N_mesh;

    for(i=0; i<N_mesh + 10; i++){
        rho = i*drho;
        for(j=0; j<N_mesh; j++){
            z = z_a + j*dz;
            myfile << Integrand_dX(rho, z, a[p], a[q], R_A, R_B, c, X) << "       ";
        }
        myfile << "  " << endl;
    }

    myfile.close();
}