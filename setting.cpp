#include "setting.hpp"

namespace Par{
// #==================================================#
// +-------------------- [OPTION] --------------------+
// #==================================================#
    // Basic Setting
    const int opt_sim_type = 1;  // The simulation type;\
                                    1:= Plane strain,\
                                    2:= Plane stress

    const int opt_int_init = 1;  // The init. opt. for internal node;\
                                    1:= Regular,\
                                    2:= Finer near boundary

    const int opt_saving = 1;    // The saving setting;\
                                    1:= Type 1,\
                                    2:= Continued...
    
    const int opt_size = 1;

    const bool flag_save_log = true;
    const bool flag_a = true;
    const bool flag_b = true;

// #==================================================#
// +--------------- [PROGRAM PARAMETER] --------------+
// #==================================================#
    // Saving Parameter
    // ...

    // Boundary Element Calculation Parameter
    const double tes1 = 1.0e0;
    const double tes2 = 1.2e0;

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

    // Base Geometry Parameter
    const int G_type = 1;    // Type of geometry: \
                                    * 1 := Rectangular\
                                    * 2 := Circular/Oval

    const double dom_Lx = 5.0e0;    // Base geometry x length
    const double dom_Ly = 4.0e0;    // Base geometry y length
    
    // Traction Parameter for Rectangular geometry\
       -> traction is constant along the surface
    // Bottom surface
    const double trac_b_x = 0.0e0;
    const double trac_b_y = -1.0e0;
    // Right surface
    const double trac_r_x = 1.0e0;
    const double trac_r_y = 0.0e0;
    // Top surface
    const double trac_t_x = 0.0e0;
    const double trac_t_y = 1.0e0;
    // Left surface
    const double trac_l_x = -1.0e0;
    const double trac_l_y = 0.0e0;

    // Traction Parameter for Circular geometry\
       -> traction is only a pressure
    const double trac_press = 0.0e0;

    // Parameter of Geometry inside Domain
    // ***********************************
    const int N_Gin = 2;        // Number of geometry inside the domain (multiply connected)

    // Parameter List of Geometry Inside the Domain 
    const std::vector<int> Gin_type = {2, 1};   // Type of geometry: \
                                                    * 1 := Rectangular\
                                                    * 2 := Circular/Oval
    const std::vector<double> Gin_Xlen = {2.0e0, 0.5e0};        // Geometry length in x direction
    const std::vector<double> Gin_Ylen = {1.0e0, 0.5e0};        // Geometry length in y direction
    const std::vector<double> Gin_Xcen_pos = {0.0e0, 2.0e0};    // Geometry center x position
    const std::vector<double> Gin_Ycen_pos = {0.0e0, 1.0e0};    // Geometry center y position
    const std::vector<double> Gin_Rot = {30.0e0, 45.0e0};         // Geometry rotation in CCW direction (in degree)
    const std::vector<double> In_pressure = {0.0e0, 0.0e0};     // The value of internal pressure\
                                                                    -> Traction for internal boundary still limited to internal pressure

// #==================================================#
// +------------- [SIMULATION PARAMETER] -------------+
// #==================================================#
    // Panel Element Parameter
    const double len = 0.08e0;     // Panel length

    // Internal Node Parameter
    const double spc = 0.1e0;     // Internal node spacing
}
