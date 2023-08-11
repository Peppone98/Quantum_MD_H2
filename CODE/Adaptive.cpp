/********** SIMPSON ADAPTIVE INTEGRATION PART ********/

#include "definitions.h" 


double Get_Simpson_rho(double a, double b, double z, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X, string s){

    double m=0.0, fm=0.0, fa=0.0, fb=0.0;

    if(s == "V_xc"){
        /**** The extremese of integration are a and b ****/
        m = (a + b)/2.0;
        fm = Integrand(m, z, alpha, beta, R_A, R_B, c, X);
        fa = Integrand(a, z, alpha, beta, R_A, R_B, c, X);
        fb = Integrand(b, z, alpha, beta, R_A, R_B, c, X);
    }
    if(s == "dVxc_dX"){
        /**** The extremese of integration are a and b ****/
        m = (a + b)/2.0;
        fm = Integrand_dX(m, z, alpha, beta, R_A, R_B, c, X);
        fa = Integrand_dX(a, z, alpha, beta, R_A, R_B, c, X);
        fb = Integrand_dX(b, z, alpha, beta, R_A, R_B, c, X);
    }

    /**** Return the "whole" integral on the [a, b] interval ****/
    return (b - a)/6.0*(fa + 4.0*fm + fb);
}



/**** Adaptive Simpson's Rule, recursive core on the rho axis ****/
double Adaptive_Simpsons_rho(double a, double b, double z, double eps, double whole, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X, string s){

    /**** The coordinate z here is treated as another parameter ****/
    double m = (a + b)/2.0;
    double left  = Get_Simpson_rho(a, m, z, alpha, beta, R_A, R_B, c, X, s);
    double right = Get_Simpson_rho(m, b, z, alpha, beta, R_A, R_B, c, X, s);

    /**** Difference between the sum of the separated integrals and the total integral ****/
    double delta = left + right - whole;

    /**** Lyness 1969 + Richardson extrapolation, see article ****/
    if(fabs(delta) <= 15*eps){
        return left + right + delta/15;
    }else{
        double new_left = Adaptive_Simpsons_rho(a, m, z, eps/2.0, left, alpha, beta, R_A, R_B, c, X, s);
        double new_right = Adaptive_Simpsons_rho(m, b, z, eps/2.0, right, alpha, beta, R_A, R_B, c, X, s);
        return new_left + new_right;
    }
}



double Get_Simpson_z(double a, double b, double eps, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X, string s){

    /**** Extremes of integration for rho ****/
    double k = alpha + beta;
    double rho_a = 0.0;
    double rho_b = 2.0 + 1/k;

    /**** Midpoint of the z-grid ****/
    double m = (a + b)/2.0;

    /**** Compute the first "whole" integral over rho with z = a, b, m ****/
    double whole_a = Get_Simpson_rho(rho_a, rho_b, a, alpha, beta, R_A, R_B, c, X, s);
    double whole_b = Get_Simpson_rho(rho_a, rho_b, b, alpha, beta, R_A, R_B, c, X, s);
    double whole_m = Get_Simpson_rho(rho_a, rho_b, m, alpha, beta, R_A, R_B, c, X, s);

    /**** Start the first iterative cycle ****/
    double fa = Adaptive_Simpsons_rho(rho_a, rho_b, a, eps, whole_a, alpha, beta, R_A, R_B, c, X, s);
    double fb = Adaptive_Simpsons_rho(rho_a, rho_b, b, eps, whole_b, alpha, beta, R_A, R_B, c, X, s);
    double fm = Adaptive_Simpsons_rho(rho_a, rho_b, m, eps, whole_m, alpha, beta, R_A, R_B, c, X, s);

    /**** Now, using the integrals over rho, we do the analog of "Get_Simpson_rho", but on z ****/
    return (b - a)/6.0*(fa + 4.0*fm + fb);
}



/**** Adaptive Simpson's Rule, recursive core on the z axis ****/
double Adaptive_Simpsons_z(double a, double b, double eps, double whole, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X, string s){

    double m = (a + b)/2.0;
    double left  = Get_Simpson_z(a, m, eps, alpha, beta, R_A, R_B, c, X, s);
    double right = Get_Simpson_z(m, b, eps, alpha, beta, R_A, R_B, c, X, s);

    /**** Difference between the sum of the separated integrals and the total integral ****/
    double delta = left + right - whole;

    /**** Lyness 1969 + Richardson extrapolation, see article ****/
    if(fabs(delta) <= 15*eps){
        return left + right + delta/15;
    }else{
        double new_left = Adaptive_Simpsons_z(a, m, eps/2.0, left, alpha, beta, R_A, R_B, c, X, s);
        double new_right = Adaptive_Simpsons_z(m, b, eps/2.0, right, alpha, beta, R_A, R_B, c, X, s);
        return new_left + new_right;
    }
}



double Adaptive_integration(double eps, double alpha, double beta, R R_A, R R_B, gsl_vector *c, double X, string s){

    /**** Set the extremes for z-integration ****/
    double k = alpha + beta;
    double a = -1.5 - 1/k;
    double b = 2.5 + 1/k;

    /**** We need the initial whole integral to start the iterative cycle ****/
    double whole = Get_Simpson_z(a, b, eps, alpha, beta, R_A, R_B, c, X, s);
    double I = Adaptive_Simpsons_z(a, b, eps, whole, alpha, beta, R_A, R_B, c, X, s);
    return I;
}



void Adaptive_Ex_Corr(gsl_matrix *V_xc, gsl_matrix *dVxc_dX, R R_A, R R_B, gsl_vector *c, double X, string s){
	int p, q;
	double val1=0., val2=0.;

    /**** The idea is to set an higher precision for the integration of fastly decaying functions ****/
    double eps[4][4] = {
        {1E-6, 1E-6, 1E-6, 1E-6},
        {1E-6, 1E-5, 1E-5, 1E-4},
        {1E-6, 1E-5, 1E-4, 1E-4},
        {1E-6, 1E-4, 1E-4, 1E-3}
    };

    if(s == "V_xc"){
        for(p=0; p<N; p++){
            for(q=0; q<=p; q++){ 

                /**** Factor 2*pi due to the integration over theta ****/
                val1 = 2.0*pi*Adaptive_integration(eps[p][q], a[p], a[q], R_A, R_A, c, X, s);
                gsl_matrix_set(V_xc, p, q, val1);
                gsl_matrix_set(V_xc, q, p, val1);
                gsl_matrix_set(V_xc, p + N, q + N, val1);
                gsl_matrix_set(V_xc, q + N, p + N, val1);

                val2 = 2.0*pi*Adaptive_integration(eps[p][q], a[p], a[q], R_A, R_B, c, X, s);
                gsl_matrix_set(V_xc, p, q + N, val2);
                gsl_matrix_set(V_xc, q, p + N, val2);
                gsl_matrix_set(V_xc, p + N, q, val2);
                gsl_matrix_set(V_xc, q + N, p, val2);
            }
        }
    }
    if(s == "dVxc_dX"){
        for(p=0; p<N; p++){
            for(q=0; q<=p; q++){ 

                /**** Factor 2*pi due to the integration over theta ****/
                val1 = 2.0*pi*Adaptive_integration(eps[p][q], a[p], a[q], R_A, R_A, c, X, s);
                gsl_matrix_set(dVxc_dX, p, q, val1);
                gsl_matrix_set(dVxc_dX, q, p, val1);
                gsl_matrix_set(dVxc_dX, p + N, q + N, val1);
                gsl_matrix_set(dVxc_dX, q + N, p + N, val1);

                val2 = 2.0*pi*Adaptive_integration(eps[p][q], a[p], a[q], R_A, R_B, c, X, s);
                gsl_matrix_set(dVxc_dX, p, q + N, val2);
                gsl_matrix_set(dVxc_dX, q, p + N, val2);
                gsl_matrix_set(dVxc_dX, p + N, q, val2);
                gsl_matrix_set(dVxc_dX, q + N, p, val2);
            }
        }
    }
}