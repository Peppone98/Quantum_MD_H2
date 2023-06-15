/********** EXCHANGE AND CORRELATION PART DFT ********/

#include "definitions.h"
#include "global.h"
#include "xc.h"
#include "gsl/gsl_integration.h"

using namespace std;


void copy_global_variables(R R_A, R R_B, gsl_vector *c){
    R_global_A.x = R_A.x;
    R_global_B.x = R_B.x;
    R_global_A.y = R_A.y;
    R_global_B.y = R_B.y;
    R_global_A.z = R_A.z;
    R_global_B.z = R_B.z;
    gsl_vector_memcpy(c_global, c);
}


double density(double rho, double x){
    double n = 0, c_p, c_q;
    int p, q;

    /**** Return the value of n(rho, x), which is needed for e_c(n) ****/
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            c_p = gsl_vector_get(c_global, p);
            c_q = gsl_vector_get(c_global, q);
            
            /**** Sum the contributions to the density ****/
            n += c_p*c_q*phi(rho, x, a[p], R_global_A)*phi(rho, x, a[q], R_global_A);
            n += c_p*c_q*phi(rho, x, a[p], R_global_B)*phi(rho, x, a[q], R_global_B);
            n += 2.*c_p*c_q*phi(rho, x, a[p], R_global_A)*phi(rho, x, a[q], R_global_B);
        }
    }
    return n;
}



double phi(double rho, double x, int alpha, R R){
    double phi = 0.0;
    return phi = exp(-alpha*(rho*rho + (x - R.x)*(x - R.x)));
}


double phi_product(double rho, double x, double alpha, double beta, double quadrant){
    double phi_mu=0.0, phi_nu=0.0;

    /**** Identify the right couple of basis functions ****/
    if(quadrant == 1.0){
        phi_mu = phi(rho, x, alpha, R_global_A);
        phi_nu = phi(rho, x, beta, R_global_A);
    }
    if(quadrant == 2.0){
        phi_mu = phi(rho, x, alpha, R_global_A);
        phi_nu = phi(rho, x, beta, R_global_B);
    }
    if(quadrant == 3.0){
        phi_mu = phi(rho, x, alpha, R_global_B);
        phi_nu = phi(rho, x, beta, R_global_A);
    }
    if(quadrant == 4.0){
        phi_mu = phi(rho, x, alpha, R_global_B);
        phi_nu = phi(rho, x, beta, R_global_B);
    }
    return phi_mu*phi_nu;
}


double func(double x, void * params){
    double phi_mu=0.0, phi_nu=0.0;

    /**** Setting all the parameters ****/
    double rho = (( double *) params)[0];
    double alpha = (( double *) params)[1];
    double beta = (( double *) params)[2];

    /**** The quadrant can be 1.0, 2.0, 3.0, 4.0 ****/
    double quadrant = (( double *) params)[3];

    /**** Compute the needed quantities ****/
    double n = density(rho, x);

    /**** Compute the e_xc and its derivative de/dn ****/
    double e_x=0.0, e_c=0.0, dex_dn=0.0, dec_dn=0.0;
    xc_func_type functional_x;
    xc_func_init(&functional_x, XC_LDA_X, XC_UNPOLARIZED);
    xc_lda_exc_vxc(&functional_x, 1, &n, &e_x, &dex_dn);

    xc_func_type functional_c;
    xc_func_init(&functional_c, XC_LDA_C_PZ, XC_UNPOLARIZED);
    xc_lda_exc_vxc(&functional_c, 1, &n, &e_c, &dec_dn);

    double e_xc = e_x + e_c;
    double de_dn = dex_dn + dec_dn;

    /**** Return the integrand function ****/
    double phi_mu_phi_nu = 0.0;
    phi_mu_phi_nu = phi_product(rho, x, alpha, beta, quadrant);
    return phi_mu_phi_nu*(e_xc + de_dn*n);
}




double outIntegr(double rho, void * params){

    /**** Fix the rho as a parameter ****/
    ((double *) params)[0] = rho;

    /**** Initialize the workspace for integration over x ****/
    int WRKSPC_SIZE = 50;
    gsl_integration_workspace *wrkspc = gsl_integration_workspace_alloc(WRKSPC_SIZE);
    double result = 0.0, error = 0.0;

    /**** Link to the function to integrate ****/
    gsl_function intgr_func; 
    intgr_func.function = &func;
    intgr_func.params = params;

    /**** Integration over x in [-2, 3] using quadrature ****/
    gsl_integration_qags(&intgr_func, -2, 3, 0, 1e-3, WRKSPC_SIZE, wrkspc, &result, &error);
    gsl_integration_workspace_free(wrkspc);
    return result;
}



double integrate(double alpha, double beta, double quadrant){

    /**** Create a set of parameters (each parameter in the same position as in "double func(...)") ****/
    double params[4];
    params[1] = alpha;
    params[2] = beta;
    params[3] = quadrant;

    /**** Allocate space for integration over rho ****/
    int WRKSPC_SIZE = 50;
    gsl_integration_workspace *wrkspc = gsl_integration_workspace_alloc(WRKSPC_SIZE);
    double result = 0.0, error = 0.0;

    /**** Link to the function that has been integrated over dx ****/
    gsl_function intgr_func;
    intgr_func.function = &outIntegr;
    intgr_func.params = params;

    /**** Integration over rho ****/
    gsl_integration_qags(&intgr_func, 0, 3, 0, 1e-3, WRKSPC_SIZE, wrkspc, &result, &error);
    gsl_integration_workspace_free(wrkspc);
    return result;
}



void create_Ex_Corr(gsl_matrix *V_xc, R R_A, R R_B, gsl_vector *c){
	int p, q;
	double val1=0., val2=0., val3=0., val4=0.;

	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){ 

			val1 = integrate(a[p], a[q], 1.0);
	        gsl_matrix_set(V_xc, p, q, val1);
            gsl_matrix_set(V_xc, q, p, val1);

            val2 = integrate(a[p], a[q], 2.0);
            gsl_matrix_set(V_xc, p + N, q, val2);
            gsl_matrix_set(V_xc, q + N, p, val2);

            val3 = integrate(a[p], a[q], 3.0);
            gsl_matrix_set(V_xc, p, q + N, val3);
            gsl_matrix_set(V_xc, q, p + N, val3);

            val4 = integrate(a[p], a[q], 4.0);
            gsl_matrix_set(V_xc, p + N, q + N, val4);
            gsl_matrix_set(V_xc, q + N, p + N, val4);
		}
	}
}