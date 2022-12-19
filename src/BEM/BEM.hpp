#ifndef BOUNDARY_ELEMENT_CALCULATION
#define BOUNDARY_ELEMENT_CALCULATION

#ifndef INCLUDE_VARIABLE
#include "../../variable.hpp"
#endif

#ifndef SAVE_DATA
#include "../Saving/save_data.hpp"
#endif

#include "../Eigen/Dense"

class calcBEM
{
private:
    // The class private method belong here
    // ************************************
    #define M_PI 3.14159265358979323846

    // Internal object variable
    dataSaving save;            // For saving the BEM matrix

    // Internal variable
    int N;                      // The matrix row and col length
    std::vector<int> elmID;     // The element ID
    std::vector<int> elmGIN;    // The element geometry position\
                                    (-1:= base, 0:= inner_1, and forth ...)

    // Element BEM matrix calculation tools [TYPE 1]
    double calc_Aij(double a_ij, double k_ij, double L);
    double calc_Bij(double a_ij, double k_ij, double L);
    double calc_Cij(double a_ij, double k_ij, double L);
    double calc_Dij(double a_ij, double k_ij, double L);
    
    // Local parameter calculation tools [TYPE 1]
    double calc_a(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n);
    double calc_k(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n);

    // Element BEM matrix calculation tools [TYPE 2]
    double calc_G_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L);
    double calc_W_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L);
    double calc_dGdn_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L);
    double calc_dWdn_dL(double x_0, double y_0, double x, double y, double x_n, double y_n, double L);

    // Matrix operator
    void swap_col(Eigen::MatrixXd& A, Eigen::MatrixXd& B, int col);
    
public:
    // The class public method belong here
    // ***********************************

    // The BEM class internal variable initialization
    void Define_BEM(const element& elm, const std::vector<element>& in_elm);
    
    // Biharmonic Solver
    void solve_F(element& elm, std::vector<element>& in_elm);
    void solve_phi(element& elm, std::vector<element>& in_elm);
    void calculate_internal_phi(intElement& intElm, const element& elm, const std::vector<element>& in_elm);
    
    // Temperature Solver
    void solve_T(element& elm, std::vector<element>& in_elm);
    void calculate_internal_T(intElement& intElm, const element& elm, const std::vector<element>& in_elm);

    // For TESTING !!!
    void TEST_BEM(element& elm, std::vector<element>& in_elm);
    void CALC_theta();
};

#endif