/********** EXCHANGE AND CORRELATION PART DFT ********/

#include "definitions.h" 
#include "gsl/gsl_integration.h"
#include "xc.h"

using namespace std;


double density(double rho, double z, gsl_vector *c, double X){
    double n = 0;
    int p, q;
    double c_p, c_q;

    /**** Return the value of n(rho, z), which is needed for the v_xc(n) ****/
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



double Integrand(double z, void *p){
    /**** Get the parameters contained in the struct ****/
    my_params *params = (my_params *)p;
    double rho = params->rho;
    double alpha = params->alpha;
    double beta = params->beta;
    R R_A = params->R_A;
    R R_B = params->R_B;
    gsl_vector *c = gsl_vector_alloc(2*N);
    gsl_vector_memcpy(c, params->c);
    double X = params->X;

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

    /**** This quantity takes part in the integrand function ****/
    double dE_dn = dEx_dn + dEc_dn;

    /**** Compute prefactors and R_c ****/
    double prefactor = K(alpha, beta, R_A, R_B);
    R R_C = R_weighted(alpha, beta, R_A, R_B);

    /**** The motion is restricted to the x-axis, so we consider only R_C.x ****/
    return prefactor*exp(-(alpha + beta)*(rho*rho + (z - R_C.x)*(z - R_C.x)))*dE_dn*rho;
}




double Integrate_over_z(my_params params){
    int WRKSPC_SIZE = 100;
    gsl_integration_workspace *wrkspc = gsl_integration_workspace_alloc(WRKSPC_SIZE);
    double result = 0.0, error = 0.0;
    gsl_function intgr_func;
    intgr_func.function = &Integrand;
    intgr_func.params = &params;
    gsl_integration_qagi(&intgr_func, 0, 1e-7, WRKSPC_SIZE, wrkspc, &result, &error);
    gsl_integration_workspace_free(wrkspc);
    return result;
}




double Simpson_over_rho(double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X){
    double rho;

    /**** Extremes of integration ****/
    double rho_0 = 0.0;
    double rho_N = 2.0;

    /**** Define the parameters ****/
    my_params par = {rho_0, alpha, beta, R_A, R_B, c, X};

    /**** First and last points ****/
    double f_0 = Integrate_over_z(par);
    par.rho = rho_N;
    double f_N = Integrate_over_z(par);

    /**** Set the mesh resolution ****/
    int i, N_mesh = 10000;
    double drho = (rho_N - rho_0)/N_mesh;

    /**** Simpson integration ****/
    double Simpson_sum = 0.0, f_i = 0.0, f_i_plus_1 = 0.0;
    for(i=1; i<=N_mesh-1; i=i+2){
        par.rho = rho_0 + i*drho;
        f_i = Integrate_over_z(par);
        par.rho = rho_0 + (i+1)*drho;
        f_i_plus_1 = Integrate_over_z(par);
        Simpson_sum = Simpson_sum + (4.*f_i + 2.*f_i_plus_1);
    }

    Simpson_sum = (Simpson_sum + f_0 + f_N)*(drho/3.);
    return Simpson_sum;
}





void create_Ex_Corr(gsl_matrix *V_xc, R R_A, R R_B, gsl_vector *c, double X){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){ 

            /**** Factor 2*pi due to the integration over theta ****/
			val1 = 2.0*pi*Simpson_over_rho(a[p], a[q], R_A, R_A, c, X);
	        gsl_matrix_set(V_xc, p, q, val1);
            gsl_matrix_set(V_xc, q, p, val1);
            gsl_matrix_set(V_xc, p + N, q + N, val1);
            gsl_matrix_set(V_xc, q + N, p + N, val1);

            val2 = 2.0*pi*Simpson_over_rho(a[p], a[q], R_A, R_B, c, X);
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

