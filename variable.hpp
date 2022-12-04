#ifndef INCLUDE_VARIABLE
#define INCLUDE_VARIABLE
#include "setting.hpp"

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

    // Boundary condition Dirichet and Neumann
    std::vector<double> F;      // Airy stress laplacian : Del^2(phi)
    std::vector<double> p;      // Airy stress function  : phi
    std::vector<double> dFdn;   // Airy stress laplacian normal plane derivative
    std::vector<double> dpdn;   // Airy stress function normal plane derivative
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
    std::vector<double> R;      // Distance to nearest panel

    // Boundary condition Dirichet and Neumann
    std::vector<double> phi;      // Airy stress function (phi)

    // Stress properties
    std::vector<double> s_xx;   // Normal stress in x-direction
    std::vector<double> s_yy;   // Normal stress in y-direction
    std::vector<double> s_zz;   // Normal stress in z-direction
    std::vector<double> t_xy;   // Shear stress in xy-plane

    // Strain properties
    std::vector<double> e_xx;   // Normal strain in x-direction
    std::vector<double> e_yy;   // Normal strain in y-direction
    std::vector<double> e_zz;   // Normal strain in z-direction
    std::vector<double> e_xy;   // Shear strain in xy-plane

    // Displacement properties
    std::vector<double> u;      // Displacement in x direction
    std::vector<double> v;      // Displacement in y direction
};

#endif

/* DATA STORAGE
    - Geometry
        > Boundary node position (x,y)
        > Boundary condition <?> (u,v,Tx,Ty)
        > The index is arrange in sequence CCW direction
    - Boundary element
        > Midpoint position (x,y)
        > Boundary condtion (p, dpdn, F, dFdn)
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