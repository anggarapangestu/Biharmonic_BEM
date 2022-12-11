#include "BEM.hpp"

// Matrix Element calculation
double calcBEM::calc_Aij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Aij, _tan;

    // Calculate internal variable
    _tan = std::atan2(a_ij*L, ((a_ij * a_ij) + (k_ij * k_ij) - (L * L)/4.0));
    
    // Calculate A element
    Aij = 1.0/(2.0 * M_PI) * _tan;
    return Aij;
}
double calcBEM::calc_Bij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Bij, var1, var2, _tan;

    // Calculate internal variable
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow((k_ij + L/2.0), 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow((k_ij - L/2.0), 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    
    // Calculate B element
    Bij = 1/(2.0 * M_PI) * 
    ( 
          (k_ij + L/2.0) * var1
        - (k_ij - L/2.0) * var2
        + a_ij * _tan
        - L
    );
    return Bij;
}
double calcBEM::calc_Cij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Cij, var1, var2, _tan, par1, par2;

    // Calculate internal variable
    par1 = k_ij + L/2.0;
    par2 = k_ij - L/2.0;
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par1, 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par2, 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    
    // Calculate B element
    Cij = a_ij/(8.0 * M_PI) * 
    ( 
          2.0 * par1 * var1
        - 2.0 * par2 * var2
        + 2.0 * _tan 
        - 3.0 * L
    );
    return Cij;
}
double calcBEM::calc_Dij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Dij, var1, var2, _tan, par1, par2;

    // Calculate internal variable
    par1 = k_ij + L/2.0;
    par2 = k_ij - L/2.0;
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par1, 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par2, 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    
    // Calculate B element
    Dij = 1/(8.0 * M_PI) * 
    ( 
          (std::pow(a_ij, 2) * par1 + (1.0/3.0) * std::pow(par1, 3)) * var1
        - (std::pow(a_ij, 2) * par2 + (1.0/3.0) * std::pow(par2, 3)) * var2
        - (4.0/9.0) * std::pow(par1, 3)
        + (4.0/9.0) * std::pow(par2, 3)
        - (5.0/3.0) * std::pow(a_ij, 2) * L
        + (2.0/3.0) * std::pow(a_ij, 3) * _tan
    );
    return Dij;
}

// Local parameter calculation tools
double calcBEM::calc_a(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n){
    // Calculation of a\
       -> the projection length of distance x_0 to x_m\
       -> a = (x_m - x_0) dot hat(n)
    double a;
    a = (x_m - x_0) * x_n + (y_m - y_0) * y_n;
    return a;
}
double calcBEM::calc_k(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n){
    // Calculation of k\
       -> the length from nearest point intersection to x_m\
       -> k = (x_m - x_0) dot hat(ell)
    double k;
    k = -(x_m - x_0) * y_n + (y_m - y_0) * x_n;
    return k;
}

// Matrix operator
void calcBEM::swap_col(Eigen::MatrixXd& A, Eigen::MatrixXd& B, int col){
    // Perform the colomn swap
    double _temp;
    for (int i = 0; i < A.rows(); i++){
        _temp = -A(i,col);
        A(i,col) = -B(i,col);
        B(i,col) = _temp;
    }
    std::cout << "[LOG] Matrix Swapped\n";
}
