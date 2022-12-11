#ifndef SETTINGs
#define SETTINGs

#include <iostream>
#include <vector>
#include <cmath>

namespace Par{
// #==================================================#
// +-------------------- [OPTION] --------------------+
// #==================================================#
    // Basic Setting
    const extern int opt_sim_type;  // The simulation type;\
                                        1:= Plane strain,\
                                        2:= Plane stress

    const extern int opt_int_init;  // The init. opt. for internal node;\
                                        1:= Regular,\
                                        2:= Finer near boundary

    const extern int opt_saving;    // The saving setting;\
                                        1:= Type 1,\
                                        2:= Continued...
    
    const extern int opt_size;

    const extern bool flag_save_log;
    const extern bool flag_a;
    const extern bool flag_b;

// #==================================================#
// +--------------- [PROGRAM PARAMETER] --------------+
// #==================================================#
    // Saving Parameter
    // ...

    // Boundary Element Calculation Parameter
    const extern double tes1;
    const extern double tes2;
    
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
    const extern int G_type;       // Type of geometry: \
                                    * 1 := Rectangular\
                                    * 2 := Circular/Oval
    const extern double dom_Lx;    // Base geometry x length
    const extern double dom_Ly;    // Base geometry y length
    
    // Traction Parameter for Rectangular geometry\
       -> traction is constant along the surface
    // Bottom surface
    const extern double trac_b_x;
    const extern double trac_b_y;
    // Right surface
    const extern double trac_r_x;
    const extern double trac_r_y;
    // Top surface
    const extern double trac_t_x;
    const extern double trac_t_y;
    // Left surface
    const extern double trac_l_x;
    const extern double trac_l_y;

    // Traction Parameter for Circular geometry\
       -> traction is only a pressure
    const extern double trac_press;

    // -----------------------------------
    // Parameter of Geometry inside Domain
    // -----------------------------------
    const extern int N_Gin;                         // Number of geometry inside the domain (multiply connected)

    // Parameter List of Geometry Inside the Domain 
    const extern std::vector<int> Gin_type;         // Type of geometry: \
                                                       * 1 := Rectangular\
                                                       * 2 := Circular/Oval
    const extern std::vector<double> Gin_Xlen;      // Geometry length in x direction
    const extern std::vector<double> Gin_Ylen;      // Geometry length in y direction
    const extern std::vector<double> Gin_Xcen_pos;  // Geometry center x position
    const extern std::vector<double> Gin_Ycen_pos;  // Geometry center y position
    const extern std::vector<double> Gin_Rot;       // Geometry rotation in CCW direction by degree
    const extern std::vector<double> In_pressure;   // The value of internal pressure\
                                                       -> Traction for internal boundary still limited to internal pressure

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
}

#endif