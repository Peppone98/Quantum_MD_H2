
#include "definitions.h"

using namespace std; 

/** forse prima posso fare un test con l'oscillatore armonico **/

/** impose the boundary conditions y_0 = y_M = 0 **/

void solve_SE(){
    gsl_vector *eval = gsl_vector_alloc(M);
    gsl_matrix *evec = gsl_matrix_alloc(M, M); 
    gsl_matrix *A = gsl_matrix_alloc(M, M);
    gsl_matrix_set_zero(A);
    double *y, x;
    y = new double[M];
    int i; 

    /**** filling the matrix ****/
    for(i=0; i<M; i++){
        x = -L/2. + i*h;
        if(i<M){
            gsl_matrix_set(A, i, i, 1./(h*h) + v(x));
        }
        if(i<M-2){
            gsl_matrix_set(A, i, i+1, -0.5/(h*h));
        }
        if(i>0){
            gsl_matrix_set(A, i, i-1, -0.5/(h*h));
        }
    }

    cout << gsl_matrix_get(A, 1, 0) << endl;
    cout << gsl_matrix_get(A, M-1, M-2) << endl;


    /** la A sembra essere troppo semplice per dare risultati ragionevoli */
    create_eval_evec(A, evec, eval);

    cout << gsl_matrix_get(A, 2, 3) << endl;

    ofstream solution;
    solution.open("Numerov_sol.txt", ios :: out | ios :: trunc);
    for(i=0; i<M; i++){
        y[i] = gsl_matrix_get(evec, i, 0);
        solution << y[i] << endl;
    }
    solution.close();
    double E = gsl_vector_get(eval, 1);
    cout << E << endl;
}


void prepare_grid(R *r){
    int i;
    for(i=0; i<M; i++){
        r[i].x = -L/2. + i*h;
        r[i].y = -L/2. + i*h;
        r[i].z = -L/2. + i*h;
    }
}


double v(double x){
    return 0.5*x*x;
}


void prepare_density(double p[M][M][M]){
    int i, j, k;
    for(i=0; i<M/4; i++){
        for(j=0; j<M/4; j++){
            for(k=0; k<M/4; k++){
                p[i][j][k] = 0.;
                p[i + 3*M/4][j + 3*M/4][k + 3*M/4] = 0.;
            }
        }
    }

    for(i=M/4; i<3*M/4; i++){
        for(j=M/4; j<3*M/4; j++){
            for(k=M/4; k<3*M/4; k++){
                p[i][j][k] = 1.;
            }
        }
    }
}


double compute_v_KS(int I, int J, int K, R *r, double p[M][M][M]){
   int i, j, k;
   double integral=0., dx, dy, dz;
   for(i=0; i<I; i++){
        for(j=0; j<J; j++){
            for(k=0; k<K; k++){
                dx = r[i].x - r[I].x;
                dy = r[j].y - r[J].y;
                dz = r[k].z - r[K].z;
                integral += h*h*h*p[i][j][k]/sqrt(dx*dx + dy*dy + dz*dz);
            }
        }
    }
    return integral; 
}