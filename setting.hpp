#ifndef SETTINGs
#define SETTINGs

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>

namespace Par{
// #==================================================#
// +-------------- [SIMULATION OPTION] ---------------+
// #==================================================#
    
    // Simulation type
    // 1:= Biharmonic solver,
    // 2:= Temperature solver
    const extern int opt_sim_type;
    
    // The type of biharmonic solution
    // 1:= Plane strain,
    // 2:= Plane stress
    const extern int opt_biharmonic_type;
                    
    // The initialization option for internal node
    // 1:= Regular distribution,
    // 2:= Finer near boundary
    const extern int opt_int_init; 
    
    // The option for BEM calculation
    // 1:= Type 1 calculation -> Calculate A, B, C, and D,
    // 2:= Type 2 calculation -> Calculate G, dGdn, W, and dWdn (NOT WORK, may be deleted)
    const extern int opt_BEM;
    
    // // The flag of additional cylindrical coordinate evaluation
    // // 0:= Basic cartesian coordinate,
    // // 1:= Add the cylindrical coordinate
    // const extern int opt_cylinder;
                                    
// #==================================================#
// +--------------- [PROGRAM PARAMETER] --------------+
// #==================================================#
    
    // Property Calculation Parameter
    const extern bool flag_cylinder;        // Flag to save simulation log

    // Saving Parameter
    const extern bool flag_save_log;        // Flag to save simulation log
    const extern bool flag_save_BEM;        // Flag to save boundary element data
    const extern bool flag_save_Int_Node;   // Flag to save internal node data
    const extern bool flag_save_Neigh;      // Flag to save neighbor ID
    
// #==================================================#
// +-------------- [GEOMETRY PARAMETER] --------------+
// #==================================================#
    /* The base geometry visualization
      RECTANGULAR :
       _______________________________
      |                               |   ^
      |                               |   |
      |           Origin at           |   |
      |              (+)              |   | dom_Ly
      |         Domain Center         |   |
      |                               |   |
      |_______________________________|   v
      <------------dom_Lx------------->
      
      CIRCULAR :
                 _ _ _ _    
              *           *          ^
           *                 *       |
         *                     *     |
        *                       *    |
       *        Origin at        *   |
       *           (+)           *   | dom_Ly
       *      Domain Center      *   |
        *                       *    |
         *                     *     |
           *                 *       |
              *  _ _ _ _  *          v                                
       <---------dom_Lx---------->
    */

    // ---------------------------------------
    // Parameter of the BASE Geometry Boundary
    // ---------------------------------------
    const extern int G_type;        // Type of geometry:
                    // 1 := Rectangular
                    // 2 := Circular/Oval
    const extern double dom_Lx;     // Base geometry x length
    const extern double dom_Ly;     // Base geometry y length
    
    // Boundary Value Parameter for Rectangular geometry\
       -> traction is constant along the surface (in Pascal)\
       -> temperature of constant dirichelt or neumann (in Kelvin)
    // Bottom surface
    const extern double trac_b_x;   // Bottom traction in x direction
    const extern double trac_b_y;   // Bottom traction in y direction
    const extern double temp_b;     // Bottom Temperature value
    const extern bool temp_type_b;  // BC type
    // Right surface
    const extern double trac_r_x;   // Right traction in x direction
    const extern double trac_r_y;   // Right traction in y direction
    const extern double temp_r;     // Right temperature value
    const extern bool temp_type_r;  // BC type
    // Top surface
    const extern double trac_t_x;   // Top traction in x direction
    const extern double trac_t_y;   // Top traction in y direction
    const extern double temp_t;     // Top temperature value
    const extern bool temp_type_t;  // BC type
    // Left surface
    const extern double trac_l_x;   // Left traction in x direction
    const extern double trac_l_y;   // Left traction in y direction
    const extern double temp_l;     // Left temperature value
    const extern bool temp_type_l;  // BC type

    // Boundary Value Parameter for Circular geometry\
       -> traction is only a pressure\
       -> temperature is only a constant
    const extern double trac_press;
    const extern double Temp;
    const extern bool Temp_type;

    // ----------------------------------------
    // Parameter of the INNER Geometry Boundary
    // ----------------------------------------
    const extern int N_Gin;                         // Number of geometry inside the domain (multiply connected)

    // Parameter List of Geometry Inside the Domain 
    const extern std::vector<int> Gin_type;         // Type of geometry
                            // 1 := Rectangular
                            // 2 := Circular/Oval
    const extern std::vector<double> Gin_Xlen;      // Geometry length in x direction
    const extern std::vector<double> Gin_Ylen;      // Geometry length in y direction
    const extern std::vector<double> Gin_Xcen_pos;  // Geometry center x position
    const extern std::vector<double> Gin_Ycen_pos;  // Geometry center y position
    const extern std::vector<double> Gin_Rot;       // Geometry rotation in CCW direction by degree
    const extern std::vector<double> In_pressure;   // The value of internal pressure\
                                                       -> Traction for internal boundary still limited to internal pressure
    const extern std::vector<double> In_temp;       // The value of boundary temperature value 
    const extern std::vector<bool> In_temp_type;    // The type of boundary temperature value 

// #==================================================#
// +------------- [SIMULATION PARAMETER] -------------+
// #==================================================#
    // Panel Element Parameter
    const extern double len;      // Panel length

    // Internal Node Parameter
    const extern double spc;      // Internal node spacing
    const extern double dist_fac; // The spacing factor of finer region

    // Neighbor Parameter
    const extern double R_s;      // Support domain factor size

// #==================================================#
// +------------- [PHYSICAL PROPERTIES] --------------+
// #==================================================#
    // Solid properties
    const extern double E;          // Material elasticity in GPa
    const extern double nu;         // Material poisson ratio
    const extern double mu;         // Lame constant (shear modulus)
    const extern double lambda;     // Lame constant
}

#endif