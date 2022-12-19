#include "setting.hpp"

namespace Par{
// #==================================================#
// +-------------- [SIMULATION OPTION] ---------------+
// #==================================================#
    // Simulation type
    const int opt_sim_type = 1;
                    // 1:= Biharmonic solver,
                    // 2:= Temperature solver
    
    // The type of biharmonic solution
    const int opt_biharmonic_type = 1;
                    // 1:= Plane strain
                    // 2:= Plane stress

    // The initialization option for internal node
    const int opt_int_init = 2;
                    // 1:= Regular distribution,
                    // 2:= Finer near boundary
    
    // The option for BEM calculation
    const int opt_BEM = 1;
                    // 1:= Type 1 calculation -> Calculate A, B, C, and D,
                    // 2:= Type 2 calculation -> Calculate G, dGdn, W, and dWdn (NOT WORK, may be deleted)

// #==================================================#
// +--------------- [PROGRAM PARAMETER] --------------+
// #==================================================#
    // Property Calculation Parameter
    const bool flag_cylinder = opt_biharmonic_type == 1? true : false;

    // Saving flag parameter
    const bool flag_save_log = true;       // Flag to save simulation log
    const bool flag_save_BEM = true;        // Flag to save boundary element data
    const bool flag_save_Int_Node = true;   // Flag to save internal node data
    const bool flag_save_Neigh = false;     // Flag to save neighbor ID

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
    const int G_type = 2;           // Type of geometry:
                    // 1 := Rectangular
                    // 2 := Circular/Oval
    const double dom_Lx = 3.0e0;    // Base geometry x length
    const double dom_Ly = 3.0e0;    // Base geometry y length
    
    // Boundary Value Parameter for Rectangular geometry\
       -> traction is constant along the surface (in pascal)\
       -> temperature of constant dirichelt or neumann (in Kelvin)
    // Bottom surface
    const double trac_b_x = 0.0e3;   // Bottom traction in x direction
    const double trac_b_y = -10.0e3;   // Bottom traction in y direction
    const double temp_b = 300.0e0;   // Bottom Temperature value
    const bool temp_type_b = false;  // BC type
    // Right surface
    const double trac_r_x = 1.0e3;   // Right traction in x direction
    const double trac_r_y = 0.0e3;   // Right traction in y direction
    const double temp_r = 300.0e0;   // Right temperature value
    const bool temp_type_r = false;  // BC type
    // Top surface
    const double trac_t_x = 0.0e3;   // Top traction in x direction
    const double trac_t_y = 10.0e3;   // Top traction in y direction
    const double temp_t = 300.0e0;   // Top temperature value
    const bool temp_type_t = false;  // BC type
    // Left surface
    const double trac_l_x = -1.0e3;  // Left traction in x direction
    const double trac_l_y = 0.0e3;   // Left traction in y direction
    const double temp_l = 500.0e0;   // Left temperature value
    const bool temp_type_l = false;  // BC type
    
    // Boundary Value Parameter for Circular geometry\
       -> traction is only a pressure\
       -> temperature is only a constant
    const double trac_press = 0.0e3;
    const double Temp = 300.0;
    const bool Temp_type = false;

    // ----------------------------------------
    // Parameter of the INNER Geometry Boundary
    // ----------------------------------------
    const int N_Gin = 1;                                // Number of geometry inside the domain (multiply connected)

    // Parameter List of Geometry Inside the Domain 
    const std::vector<int> Gin_type = {2, 2, 2, 2, 2};  // Type of geometry
                            // 1 := Rectangular
                            // 2 := Circular/Oval
    const std::vector<double> Gin_Xlen = {1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0};           // Geometry length in x direction
    const std::vector<double> Gin_Ylen = {1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0};           // Geometry length in y direction
    const std::vector<double> Gin_Xcen_pos = {0.0e0, 1.0e0, 1.0e0, -1.0e0, -1.0e0};     // Geometry center x position
    const std::vector<double> Gin_Ycen_pos = {0.0e0, 1.0e0, -1.0e0, -1.0e0, 1.0e0};     // Geometry center y position
    const std::vector<double> Gin_Rot = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};           // Geometry rotation in CCW direction (in degree)
    
    // Boundary value
    const std::vector<double> In_pressure = {5.0e3, 0.0e0, 0.0e0, 0.0e0, 0.0e0};        // The value of internal pressure\
                                                                                        -> Traction for internal boundary still limited to internal pressure
    const std::vector<double> In_temp = {500.0e0, 200.0e0, 200.0e0, 200.0e0, 200.0e0};  // The value of boundary temperature value
    const extern std::vector<bool> In_temp_type = {false, false, false, false, false};  // The type of boundary temperature value

// #==================================================#
// +------------- [SIMULATION PARAMETER] -------------+
// #==================================================#
    // Panel Element Parameter
    const double len = 0.04e0;      // Panel length

    // Internal Node Parameter
    const double spc = 0.04e0;      // Internal node spacing
    const double dist_fac = 4.0e0;  // The spacing factor of finer region

    // Neighbor Parameter
    const double R_s = 3.5e0;       // Support domain factor size

// #==================================================#
// +------------- [PHYSICAL PROPERTIES] --------------+
// #==================================================#
    // Solid properties
    const double E = 69.0e9;                        // Material elasticity in GPa
    const double nu = 0.3;                          // Material poisson ratio
    const double mu = (0.5 * E) / (1 + nu);         // Lame constant (shear modulus)
    const double lambda = (2*mu*nu) / (1 - 2*nu);   // Lame constant
}
