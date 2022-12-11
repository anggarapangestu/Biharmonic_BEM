#ifndef BOUNDARY_ELEMENT_CALCULATION
#define BOUNDARY_ELEMENT_CALCULATION

#ifndef INCLUDE_VARIABLE
#include "../../variable.hpp"
#endif

#include "../Eigen/Dense"

class calcBEM
{
private:
    // The class private method belong here
    #define M_PI 3.14159265358979323846
    // Element BEM matrix calculation tools
    double calc_Aij(double a_ij, double k_ij, double L);
    double calc_Bij(double a_ij, double k_ij, double L);
    double calc_Cij(double a_ij, double k_ij, double L);
    double calc_Dij(double a_ij, double k_ij, double L);
    
    // Local parameter calculation tools
    double calc_a(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n);
    double calc_k(double x_0, double y_0, double x_m, double y_m, double x_n, double y_n);
public:
    // The class public method belong here
    void solve_F(element& elm);
    void solve_phi(element& elm);
    void calculate_internal_phi(intElement& intElm, const element& elm);
};

// - Solve the Laplace of F
// - Solve the Poisson of p (phi - Airy Stress Function)
#endif

/*

class className
{
private:
    // The class private method belong here
public:
    // The class public method belong here
};

*/
