

#include <armadillo>

using namespace arma;
using namespace std;

int main(){
    int N = 5;
    sp_mat X(N, N);
    X = speye(N, N);
    for(int i=0; i<N-1; i++){
        X(i, i+1) = -2.;
    }

    vec eigval;
    mat eigvec;
    eig_gen(eigval, eigvec, X);
    eigvec.print("Eigen Vector");
    eigval.print("Eigen Value");
}