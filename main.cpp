#ifndef INCLUDE_VARIABLE
#include "variable.hpp"
#endif

#ifndef SETTINGs
#include "setting.hpp"
#endif

#ifndef INITIALIZATION
#include "src/Initialization/initialization.hpp"
#endif

#ifndef SAVE_DATA
#include "src/Saving/save_data.hpp"
#endif

int main(){
    // Initialize the storage data for boundary and internal element
    element PanelElement;
    std::vector<element> InnerElement;
    intElement InternalNode;
    
    // Initialize the method
    initialization init;
    dataSaving save;

    // Initialization generate the element
    init.generate_boundary_element(PanelElement, InnerElement);
    init.generate_internal_node(InternalNode, PanelElement, InnerElement);

    // // Calculate the initial condition
    // init.calculate_initial_condition(PanelElement);

    // Write the element data
    save.simulation_log();
    save.write_internal_data(InternalNode);
    save.write_BEM_data(PanelElement, InnerElement);

    // INITIALIZATION
    // (1) Element Generation (Internal and Boundary)
    // (2) Boundary Value Calculation (traction & displacement)
    // (3) Matrix Built Up

    // SOLVER
    // (1) Solving Laplace of F
    // (2) Solving Poison of phi

    // POST PROCESSING
    // (1) Calculate the phi inside domain
    // (2) Calculate the properties (sigma_x, sigma_y, tau_xy, u, v, ...) inside domain
    // (3) Write Data

    return 0;
}

/* >>>> READ ME!! <<<<
   Architecture:
   * Program code flow
     - Input Geometry
     - Create initial element distribution (the boundary element is indexed sequentialy)
     - Calculate the initial condition
     - Solve the Laplace of F
     - Solve the Poisson of p (phi - Airy Stress Function)
     - Calculate the p value at internal domain
     - Calculate the properties
     - Write data
   * Data Storage
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
   * Method
     - Geometry + Boundary condition generation (G)
     - Boundary and internal element generation (I)
     - Element boundary condition parameter calculation (I)
     - Laplacian Calculation (B)
     - Poisson Calculation (B)
     - Internal Airy Stress Calculation (B)
     - Properties from Airy Stress Calculation (P)
     - Write data (S)
     - Write parameter (S)
   * Class
     - Geometry (G)
     - Initialization (I)
     - Boundary Element Method (B)
     - Internal Domain Property (P)
     - LSMPS (L)
     - Saving (S)
*/

/*  >>>>> WRITING <<<<<
    SINGLE WORD TYPE
    * type 0: Capitalize all        e.g. GROUP   : Const variable
    * type 1: Normal case           e.g. group   : general variable
    * type 2: initial underscore    e.g. _group  : temp variable
    * type 3: numbering             e.g. group1  : sequence variable
    * type 4: shorten               e.g. grp     : temp variable, long variable name
    WORD SEPARATOR
    * type 0: Underscore                  e.g. group_element
    * type 1: Camel case                  e.g. groupElement
    * type 2: Capitalize each word (CEW)  e.g. GroupElement
    The function, class, and method naming system:
    * Class name        : normal(1) + camel(1)
    * Class method      : normal(1) + underscore(0)
    * Class variable    : normal(1) + underscore(0)
*/