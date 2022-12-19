#ifndef INCLUDE_VARIABLE
#define INCLUDE_VARIABLE

#ifndef SETTINGs
#include "setting.hpp"
#endif

// Boundary Element class: Store the property data of each boundary element
class element
{
private:
    // The private variable belongs here
public:
    // The public variable belongs here
    // Number of the node
    int num;

    // Midpoint position (x,y)
    std::vector<double> xm;     // Midpoint x coordinate
    std::vector<double> ym;     // Midpoint y coordinate
    std::vector<double> L;      // Panel length
    
    // Panel normal unit vector (nx, ny)
    std::vector<double> xn;
    std::vector<double> yn;

    // Raw Boundary condition
    std::vector<double> Tx;     // Traction in x direction
    std::vector<double> Ty;     // Traction in y direction
    
    // Boundary condition Dirichet and Neumann Biharmonic
    std::vector<double> F;      // Airy stress laplacian : F = Del^2(phi)
    std::vector<double> p;      // Airy stress function  : phi
    std::vector<double> dFdn;   // Airy stress laplacian normal plane derivative
    std::vector<double> dpdn;   // Airy stress function normal plane derivative
    std::vector<bool> F_type;   // The type of given boundary condition for F   (true: neumann, false: dirichlet)
    std::vector<bool> p_type;   // The type of given boundary condition for phi (true: neumann, false: dirichlet)

    // Boundary condition Dirichet and Neumann Temperature
    std::vector<double> T;      // Temperature property T
    std::vector<double> dTdn;   // Temperature normal plane derivative
    std::vector<bool> T_type;   // The type of given boundary condition for T   (true: neumann, false: dirichlet)
};

// Internal Element class: Store the property data of each node inside the domain
class intElement
{
private:
    // The private variable belongs here
public:
    // The public variable belongs here
    // Number of the node
    int num;

    // Midpoint position (x,y)
    std::vector<double> x;      // Midpoint x coordinate
    std::vector<double> y;      // Midpoint y coordinate
    std::vector<double> s;      // The element spacing size
    std::vector<double> R;      // Distance to nearest panel

    // The BEM basic function
    std::vector<double> phi;      // Airy stress function (phi)
    std::vector<double> phi_an;   // Airy stress function analytic solution (phi)
    std::vector<double> T;        // Temperature (T)
    std::vector<double> T_an;     // Temperature analytic solution (T)
    
    // ============= STRESS =============
    // Stress properties (cartesian)
    std::vector<double> s_xx;   // Normal stress in x-direction
    std::vector<double> s_yy;   // Normal stress in y-direction
    std::vector<double> s_zz;   // Normal stress in z-direction
    std::vector<double> t_xy;   // Shear stress in xy-plane

    // Stress properties (cylindrical)
    std::vector<double> s_rr;   // Normal stress in r-direction
    std::vector<double> s_tt;   // Normal stress in θ-direction
    std::vector<double> t_rt;   // Shear stress in rθ-plane
    
    // ============= STRAIN =============
    // Strain properties (cartesian)
    std::vector<double> e_xx;   // Normal strain in x-direction
    std::vector<double> e_yy;   // Normal strain in y-direction
    std::vector<double> e_zz;   // Normal strain in z-direction
    std::vector<double> e_xy;   // Shear strain in xy-plane

    // Strain properties (cylindrical)
    std::vector<double> e_rr;   // Normal strain in r-direction
    std::vector<double> e_tt;   // Normal strain in θ-direction
    std::vector<double> e_rt;   // Shear strain in rθ-plane

    // ============= Displacement =============
    // Displacement properties
    std::vector<double> u;      // Displacement in x direction
    std::vector<double> v;      // Displacement in y direction
};

#endif

/* DATA STORAGE
    - Boundary element
        > Midpoint position (x,y)
        > Boundary condtion (p, dpdn, F, dFdn, T, dTdn)
        > Panel length (L)
        > Panel normal unit vector (nx, ny)
        > The index is arrange in sequence CCW direction
    - Internal node
        > Position (x,y)
        > Airy Stress (p)
        > Stress Properties (s_xx, s_yy, s_zz, t_xy)
        > Strain Properties (e_xx, e_yy, e_zz, e_xy)
        > Displacement Properties (u, v)
*/