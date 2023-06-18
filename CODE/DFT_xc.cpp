/********** EXCHANGE AND CORRELATION PART DFT ********/

#include "definitions.h" 
#include "xc.h"

using namespace std;



double density(double rho, double z, gsl_vector *c, double X){
    double n = 0;
    int p, q;
    double c_p, c_q;

    /**** Return the value of n(rho, z), which is needed for e_c(n) ****/
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            c_p = gsl_vector_get(c, p);
            c_q = gsl_vector_get(c, q);
            
            n += c_p*c_q*(exp(-(a[p] + a[q])*(rho*rho + z*z)));
            n += c_p*c_q*(exp(-(a[p] + a[q])*(rho*rho + (z - X)*(z - X))));
            n += 2.*c_p*c_q*(exp(-a[p]*a[q]*X*X/(a[p] + a[q]))*exp(-(a[p] + a[q])*(rho*rho + (z - a[q]*X/(a[p] + a[q]))*(z - a[q]*X/(a[p] + a[q])))));
        }
    }

    return n;
}



double Integrand(double rho, double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Compute the density ****/
    double n = density(rho, z, c, X);

    /**** Compute the exchange part ****/
    double e_x=0.0, e_c=0.0, dEx_dn=0.0, dEc_dn=0.0;
    xc_func_type functional_x;
    xc_func_init(&functional_x, XC_LDA_X, XC_UNPOLARIZED);
    xc_lda_exc_vxc(&functional_x, 1, &n, &e_x, &dEx_dn);

    /**** Compute the correlation part ****/
    xc_func_type functional_c;
    xc_func_init(&functional_c, XC_LDA_C_PZ, XC_UNPOLARIZED);
    xc_lda_exc_vxc(&functional_c, 1, &n, &e_c, &dEc_dn);

    /**** This quantity will be used in the integral ****/
    double dE_dn = dEx_dn + dEc_dn;

    /**** Compute correlation energy density ****/
    double prefactor = 2.*pi*K(alpha, beta, R_A, R_B);
    R R_C = R_weighted(alpha, beta, R_A, R_B);

    /**** The motion is restricted to the x-axis, so we consider only R_C.x ****/
    return prefactor*exp(-(alpha + beta)*(rho*rho + (z - R_C.x)*(z - R_C.x)))*dE_dn*rho;
}



double Simpson_rho(double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    /**** Extremes of integration ****/
    double a = 0.0;
    double b = 2.0;
    
    /**** First and last points ****/
    double f_a = Integrand(a, z, alpha, beta, R_A, R_B, c, X);
    double f_b = Integrand(b, z, alpha, beta, R_A, R_B, c, X);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 50;
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
    int i, N_mesh = 50;
    double dz = (b-a)/N_mesh;

    /**** Simpson integration ****/
    double z = 0.0, Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
        z = i*dz;
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

			val1 = Simpson_z(a[p], a[q], R_A, R_A, c, X);
	        gsl_matrix_set(V_xc, p, q, val1);
            gsl_matrix_set(V_xc, q, p, val1);
            gsl_matrix_set(V_xc, p + N, q + N, val1);
            gsl_matrix_set(V_xc, q + N, p + N, val1);

            val2 = Simpson_z(a[p], a[q], R_A, R_B, c, X);
	        gsl_matrix_set(V_xc, p, q + N, val2);
            gsl_matrix_set(V_xc, q, p + N, val2);
            gsl_matrix_set(V_xc, p + N, q, val2);
            gsl_matrix_set(V_xc, q + N, p, val2);

		}
	}
}



