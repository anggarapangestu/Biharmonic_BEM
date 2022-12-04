#include "BEM.hpp"

void calcBEM::solve_F(element& elm){
    // Initialize the matrix
    int N = elm.num;
    Eigen::MatrixXd A_Mat = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd B_Mat = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd dFdn = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(N);

    // Fill the matrix
    double _a, _k, L, _Aij, _Bij, _Cij, _Dij;
    for (int i = 0; i < N; i++){
        // The evaluated point
        double xi = elm.xm[i];
        double yi = elm.ym[i];
        for (int j = 0; j < N; j++){
            // The point take into account
            double xj = elm.xm[j];
            double yj = elm.ym[j];
            double xn = elm.xn[i];
            double yn = elm.yn[i];

            // Calculating the matrix element
            _a = this->calc_a(xi,yi,xj,yj,xn,yn);
            _k = this->calc_k(xi,yi,xj,yj,xn,yn);
            _Aij = this->calc_Aij(_a,_k,L);
            _Bij = this->calc_Bij(_a,_k,L);
            A_Mat(i, j) = _Aij;
            B_Mat(i, j) = _Bij;
        }
        // Update the matrix
        dFdn(i) = elm.dFdn[i];
        F(i) = elm.F[i];
    }

    // Solve the matrix by using matrix invertion
    Eigen::VectorXd bi = B_Mat * dFdn;
    Eigen::VectorXd Fi = A_Mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi);

    // Update the value F on the element
    for (int i = 0; i < N; i++){
        elm.F[i] = Fi(i);
    }
}

void calcBEM::solve_phi(element& elm){
    // calcBEM::
}

void calcBEM::calculate_internal_phi(intElement& intElm, const element& elm){
    // caclBEM::
}