#include "BEM.hpp"

// Matrix Element calculation
// Representative to integral(dGdn * dL)
double calcBEM::calc_Aij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Aij, _tan;

    // Calculate internal variable
    _tan = std::atan2(a_ij*L, ((a_ij * a_ij) + (k_ij * k_ij) - (L * L)/4.0));
    if (a_ij == 0 && k_ij == 0){
        _tan = 0;
    }
    
    // Calculate A element
    Aij = (1.0/(2.0 * M_PI)) * _tan;
    return Aij;
}
// Representative to integral(G * dL)
double calcBEM::calc_Bij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Bij, var1, var2, _tan;

    // Calculate internal variable
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow((k_ij + L/2.0), 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow((k_ij - L/2.0), 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    if (a_ij == 0 && k_ij == 0){
        _tan = 0;
    }

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
// Representative to integral(dWdn * dL)
double calcBEM::calc_Cij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Cij, var1, var2, _tan, par1, par2;

    // Calculate internal variable
    par1 = k_ij + L/2.0;
    par2 = k_ij - L/2.0;
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par1, 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par2, 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    if (a_ij == 0 && k_ij == 0){
        _tan = 0;
    }
    
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
// Representative to integral(W * dL)
double calcBEM::calc_Dij(double a_ij, double k_ij, double L){
    // Internal variable declaration
    double Dij, var1, var2, _tan, par1, par2;

    // Calculate internal variable
    par1 = k_ij + L/2.0;
    par2 = k_ij - L/2.0;
    var1 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par1, 2) ));
    var2 = std::log(std::sqrt( std::pow(a_ij, 2) + std::pow(par2, 2) ));
    _tan = std::atan2(a_ij*L, (std::pow(a_ij,2) + std::pow(k_ij,2) - std::pow(L/2.0,2)));
    if (a_ij == 0 && k_ij == 0){
        _tan = 0;
    }
    
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

// Matrix Element calculation
// Representative to Dewangga Code: K1
double calcBEM::calc_G_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L){
    // Mathematical Expression of function G
    // G = 1/(2pi) * ln|x - x_0|
    // G_dL = G * L
    
    // Calculate G_dL
    double G_dL, _r;
    _r = std::sqrt( std::pow((x - x_0), 2) + std::pow((y - y_0), 2) );
    if (_r == 0){   // Singular case
        G_dL = (L / (2.0 * M_PI)) * (std::log(L/2.0) - 1);
    }else{          // Non singular case
        G_dL = (L / (2.0 * M_PI)) * std::log(_r);
    }

    return G_dL;
}
// Representative to Dewangga Code: K2
double calcBEM::calc_W_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L){
    // Mathematical Expression of function W
    // W = 1/(8pi) * |x - x_0|^2 * ((ln|x - x_0|) - 1)
    // W_dL = W * L
    
    // Calculate G_dL
    double W_dL, _r;
    _r = std::sqrt( std::pow((x - x_0), 2) + std::pow((y - y_0), 2) );
    if (_r == 0){   // Singular case
        W_dL = (1.0 / (8.0 * M_PI)) * std::pow((L/2.0), 3) * ((2.0/3.0) * std::log(L/2.0) - (8.0/9.0));
    }else{          // Non singular case
        W_dL = (L / (8.0 * M_PI)) * std::pow(_r, 2) * (std::log(_r) - 1);
    }

    return W_dL;
}
// Representative to Dewangga Code: G1
double calcBEM::calc_dGdn_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L){
    // Mathematical Expression of function DelG
    // DelG = 1/(2pi) * (x - x_0) / |x - x_0|^2
    // dG_dn_dL = DelG(dot)hat(n) * L
    
    // Calculate dG_dn_dL
    double dG_dn_dL, _r, x_dot_n;
    _r = std::sqrt( std::pow((x - x_0), 2) + std::pow((y - y_0), 2) );
    x_dot_n = (x - x_0) * x_n + (y - y_0) * y_n;

    if (_r == 0){   // Singular case
        dG_dn_dL = 0;
    }else{          // Non singular case
        dG_dn_dL = (L / (2.0 * M_PI)) * (x_dot_n / std::pow(_r, 2));
    }

    return dG_dn_dL;
}
// Representative to Dewangga Code: G2
double calcBEM::calc_dWdn_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L){
    // Mathematical Expression of function DelG
    // DelG = 1/(8pi) (ln|x - x_0|^2 - 1) * (x - x_0)
    // dG_dn_dL = DelG(dot)hat(n) * L
    
    // Calculate dG_dn_dL
    double dW_dn_dL, _r, x_dot_n;
    _r = std::sqrt( std::pow((x - x_0), 2) + std::pow((y - y_0), 2) );
    x_dot_n = (x - x_0) * x_n + (y - y_0) * y_n;

    if (_r == 0){   // Singular case
        dW_dn_dL = 0;
    }else{          // Non singular case
        dW_dn_dL = (L / (8.0 * M_PI)) * x_dot_n * (2 * std::log(_r) - 1);
    }

    return dW_dn_dL;
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
    // std::cout << "[LOG] Matrix Swapped\n";
}


// Calculate the other F boundary value
void calcBEM::TEST_BEM(element& elm, std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM testing on calc. F ...\n");
    clock_t _time = clock();

    // Initialize the matrix
    Eigen::MatrixXd A_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::MatrixXd B_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::VectorXd A_Vec = Eigen::VectorXd::Zero(this->N);       // Left : Value is calculated
    Eigen::VectorXd B_Vec = Eigen::VectorXd::Zero(this->N);       // Right: Value is given

    // ================= Fill the BEM matrix =================
    // *******************************************************
    double _a, _k, L, _Aij, _Bij;
    for (int i = 0; i < this->N; i++){
        double xi, yi;  // Current element evaluated
        double xj, yj;  // Iteration element
        double xn, yn;  // Normal vector

        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            xi = elm.xm[i];
            yi = elm.ym[i];
        }else{
            // Inner geometry
            xi = in_elm[elmGIN[i]].xm[elmID[i]];
            yi = in_elm[elmGIN[i]].ym[elmID[i]];
        }

        // Calculate all element
        for (int j = 0; j < this->N; j++){
            // Update the iterated element
            if (elmGIN[j] < 0){
                // Base geometry
                xj = elm.xm[j];
                yj = elm.ym[j];
                xn = elm.xn[j];
                yn = elm.yn[j];
                L  = elm.L[j];
            }else{
                // Inner geometry
                xj = in_elm[elmGIN[j]].xm[elmID[j]];
                yj = in_elm[elmGIN[j]].ym[elmID[j]];
                xj = in_elm[elmGIN[j]].xn[elmID[j]];
                yj = in_elm[elmGIN[j]].yn[elmID[j]];
                L  = in_elm[elmGIN[j]].L[elmID[j]];
            }

            // Calculating the matrix element
            if (Par::opt_BEM == 1){
                _a = this->calc_a(xi,yi,xj,yj,xn,yn);
                _k = this->calc_k(xi,yi,xj,yj,xn,yn);
                _Aij = this->calc_Aij(_a,_k,L);
                _Bij = this->calc_Bij(_a,_k,L);
            }else if (Par::opt_BEM == 2){
                _Aij = this->calc_dGdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Bij = this->calc_G_dL(xi,yi,xj,yj,xn,yn,L);
            }    
            
            // Matrix assignment
            if (i == j){
                A_Mat(i, j) = _Aij - 0.5;
                std::cout << "> A (" << i << ", " << j << ") : " << _Aij << std::endl;
            }else{
                A_Mat(i, j) = _Aij;
            }
            B_Mat(i, j) = _Bij;
        }
    }

    save.write_Matrix(A_Mat, "A_TEST");
    save.write_Matrix(B_Mat, "B_TEST");

    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Calculating F comp. time           [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

