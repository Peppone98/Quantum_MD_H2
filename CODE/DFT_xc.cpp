/********** EXCHANGE AND CORRELATION PART DFT ********/

#include "definitions.h" 
#include "xc.h"

using namespace std;



double density(double rho, double z, gsl_vector *c, double X){
    double n = 0;
    int p, q;
    double c_p, c_q, K, exp_factor_p, exp_factor_q;

    /**** Return the value of n(rho, z), which is needed for the v_xc(n) ****/
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            c_p = gsl_vector_get(c, p);
            c_q = gsl_vector_get(c, q);
            
            n += c_p*c_q*(exp(-(a[p] + a[q])*(rho*rho + z*z)));
            n += c_p*c_q*(exp(-(a[p] + a[q])*(rho*rho + (z - X)*(z - X))));

            /**** K factor, i.e., exp(-a_p*a_q*|R_A - R_B|^2/(a_p + a_q))****/
            K = exp(-a[p]*a[q]*X*X/(a[p] + a[q]));
            exp_factor_p = exp(-(a[p] + a[q])*(rho*rho + (z - a[p]*X/(a[p] + a[q]))*(z - a[p]*X/(a[p] + a[q]))));
            exp_factor_q = exp(-(a[p] + a[q])*(rho*rho + (z - a[q]*X/(a[p] + a[q]))*(z - a[q]*X/(a[p] + a[q]))));
            n += c_p*c_q*K*(exp_factor_p + exp_factor_q);
        }
    }

    return n;
}



double Integrand(double rho, double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Compute the density ****/
    double n = density(rho, z, c, X);

    /**** Compute the exchange part ****/
    double v_x = 0.0, v_c =0.0;
    xc_func_type functional_x;
    xc_func_init(&functional_x, FUNCTIONAL_X, XC_UNPOLARIZED);
    xc_lda_vxc(&functional_x, 1, &n, &v_x);
    xc_func_end(&functional_x);

    /**** Compute the correlation part ****/
    xc_func_type functional_c;
    xc_func_init(&functional_c, FUNCTIONAL_C, XC_UNPOLARIZED);
    xc_lda_vxc(&functional_c, 1, &n, &v_c);
    xc_func_end(&functional_c);

    /**** This dE_dn is used in the integral. Note the a_x ****/
    double dE_dn = a_x*v_x + v_c;

    /**** Compute prefactor and R_C ****/
    double prefactor = K(alpha, beta, R_A, R_B);
    R R_C = R_weighted(alpha, beta, R_A, R_B);

    /**** The motion is restricted to the x-axis, so we consider only R_C.x ****/
    return prefactor*exp(-(alpha + beta)*(rho*rho + (z - R_C.x)*(z - R_C.x)))*dE_dn*rho;
}



double Simpson_rho(double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Extremes of integration ****/
    double m = alpha + beta;
    double a = 0.0;
    double b = 2.0 + 1/m;
    
    /**** First and last points ****/
    double f_a = Integrand(a, z, alpha, beta, R_A, R_B, c, X);
    double f_b = Integrand(b, z, alpha, beta, R_A, R_B, c, X);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 100;
    double drho = (b-a)/N_mesh;

    /**** Simpson integration ****/
    double rho = 0.0, Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
        rho = i*drho;
        f_i = Integrand(rho, z, alpha, beta, R_A, R_B, c, X);
        f_i_plus_1 = Integrand(rho + drho, z, alpha, beta, R_A, R_B, c, X);
        Simpson_sum = Simpson_sum + (4.*f_i + 2.*f_i_plus_1);
    }

    Simpson_sum = (Simpson_sum + f_a + f_b)*(drho/3.);
    return Simpson_sum;
}




double Simpson_z(double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Extremes of integration ****/
    double m = alpha + beta;
    double a = -1.5 - 1/m;
    double b = 2.5 + 1/m;
    
    /**** First and last points ****/
    double f_a = Simpson_rho(a, alpha, beta, R_A, R_B, c, X);
    double f_b = Simpson_rho(b, alpha, beta, R_A, R_B, c, X);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 100;
    double dz = (b-a)/N_mesh;

    /**** Simpson integration ****/
    double z = 0.0, Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
    	
    	/**** The z coordinate starts at a finite value ****/
        z = a + i*dz;
        f_i = Simpson_rho(z, alpha, beta, R_A, R_B, c, X);
        f_i_plus_1 = Simpson_rho(z + dz, alpha, beta, R_A, R_B, c, X);
        Simpson_sum = Simpson_sum + (4.*f_i + 2.*f_i_plus_1);
    }

    Simpson_sum = (Simpson_sum + f_a + f_b)*(dz/3.);
    return Simpson_sum;
}




void create_Ex_Corr(gsl_matrix *V_xc, R R_A, R R_B, gsl_vector *c, double X){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){ 

            /**** Factor 2*pi due to the integration over theta ****/
			val1 = 2.0*pi*Simpson_z(a[p], a[q], R_A, R_A, c, X);
	        gsl_matrix_set(V_xc, p, q, val1);
            gsl_matrix_set(V_xc, q, p, val1);
            gsl_matrix_set(V_xc, p + N, q + N, val1);
            gsl_matrix_set(V_xc, q + N, p + N, val1);

            val2 = 2.0*pi*Simpson_z(a[p], a[q], R_A, R_B, c, X);
	        gsl_matrix_set(V_xc, p, q + N, val2);
            gsl_matrix_set(V_xc, q, p + N, val2);
            gsl_matrix_set(V_xc, p + N, q, val2);
            gsl_matrix_set(V_xc, q + N, p, val2);

		}
	}
}




void Print_density(gsl_vector *c, double X){
    ofstream myfile;
	myfile.open("Density.txt", ios :: out | ios :: trunc);

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
            myfile << density(rho, z, c, X) << "       ";
        }
        myfile << "  " << endl;
    }

    myfile.close();
}



void Print_integrand(int p, int q, R R_A, R R_B, gsl_vector *c, double X){
    string name = "Integrands/Integrand";
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
            myfile << Integrand(rho, z, a[p], a[q], R_A, R_B, c, X) << "       ";
        }
        myfile << "  " << endl;
    }

    myfile.close();
}

