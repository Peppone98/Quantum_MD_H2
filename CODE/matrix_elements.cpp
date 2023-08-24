/******** ROUTINES FOR MATRIX ELEMENTS *******/

#include "definitions.h"

using namespace std; 

double K(double alpha, double beta, R R_A, R R_B){
    return exp(-alpha*beta*scalar_prod(R_A, R_B)/(alpha + beta));
}

double scalar_prod(R R_A, R R_B){
    double x = R_A.x - R_B.x;
    double y = R_A.y - R_B.y;
    double z = R_A.z - R_B.z;
    return x*x + y*y + z*z;
}

R R_weighted(double alpha, double beta, R R_A, R R_B){
    R R_w;
    R_A.x = alpha*R_A.x, R_A.y = alpha*R_A.y, R_A.z = alpha*R_A.z;
    R_B.x = beta*R_B.x, R_B.y = beta*R_B.y, R_B.z = beta*R_B.z;
    R_w.x = (R_A.x + R_B.x)/(alpha + beta);
    R_w.y = (R_A.y + R_B.y)/(alpha + beta);
    R_w.z = (R_A.z + R_B.z)/(alpha + beta);
    return R_w;
}

double overlap(double alpha, double beta, R R_A, R R_B){
	return pow(pi/(alpha + beta), 1.5)*K(alpha, beta, R_A, R_B);
}


double laplacian(double alpha, double beta, R R_A, R R_B){
    double tmp = (3. + 2.*log(K(alpha, beta, R_A, R_B)));
    return alpha*beta/(alpha + beta)*tmp*overlap(alpha, beta, R_A, R_B);
}


double F0(double x){
    if(x == 0.){
        return 1;
    }else{
        return sqrt(pi/(4.*x))*erf(sqrt(x));
    }
}


double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C){
    double tmp = -2.*pi*K(alpha, beta, R_A, R_B)/(alpha + beta);
    R R_P = R_weighted(alpha, beta, R_A, R_B);
    return tmp*F0((alpha + beta)*scalar_prod(R_P, R_C));
}


double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime){
    double tmp = 2.*pow(pi, 2.5)/((alpha + alpha_prime)*(beta + beta_prime)*sqrt(alpha + alpha_prime + beta + beta_prime));
    tmp = tmp*K(alpha, alpha_prime, R_A, R_A_prime)*K(beta, beta_prime, R_B, R_B_prime);
    R R_P, R_Q;
    R_P = R_weighted(alpha, alpha_prime, R_A, R_A_prime);
    R_Q = R_weighted(beta, beta_prime, R_B, R_B_prime);
    return tmp*F0(-log(K(alpha + alpha_prime, beta + beta_prime, R_P, R_Q)));
}
