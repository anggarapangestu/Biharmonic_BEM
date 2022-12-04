#include "setting.hpp"

namespace Par{
// #==================================================#
// +-------------------- [OPTION] --------------------+
// #==================================================#
    // Basic Setting
    const int opt_sim_type = 1;  // The simulation type;\
                                        1:= Plane strain,\
                                        2:= Plane stress

    const int opt_geom = 1;      // The geometry type;\
                                        1:= Rectangle,\
                                        2:= Continued...

    const int opt_int_init = 1;  // The init. opt. for internal node;\
                                        1:= Regular,\
                                        2:= Finer near boundary

    const int opt_saving = 1;    // The saving setting;\
                                        1:= Type 1,\
                                        2:= Continued...
    
    const int opt_size = 1;
    
    const bool flag_save_log = true;  // [FLAG] Save simulation log
    const bool flag_a;
    const bool flag_b;

// #==================================================#
// +--------------- [PROGRAM PARAMETER] --------------+
// #==================================================#
    // Saving Parameter
    // ...

    // Boundary Element Calculation Parameter
    const double tes1 = 1.0e0;
    const double tes2 = 1.2e0;
    
// #==================================================#
// +------------- [SIMULATION PARAMETER] -------------+
// #==================================================#
    // Geometry Parameter (length in [m])
    const double dom_Lx = 5.0e0;
    const double dom_Ly = 4.0e0;

    /* The geometry visualization
       _______________________________
      |                               |   ^
      |                               |   |
      |           Origin at           |   |
      |              (+)              |   | dom_Ly
      |         Domain Center         |   |
      |                               |   |
      |_______________________________|  _|_
      <------------dom_Lx------------->
    */

    // Initialization Parameter
    const double traction_1;
    const double traction_2;
    const double traction_3;
    const double traction_4;

    // Element Parameter
    const double length = 0.05e0;     // Basic panel length

    // Element Parameter
    const double ls;         // Basic panel length
}