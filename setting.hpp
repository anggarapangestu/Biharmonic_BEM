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

    const extern int opt_geom;      // The geometry type;\
                                        1:= Rectangle,\
                                        2:= Continued...

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
// +------------- [SIMULATION PARAMETER] -------------+
// #==================================================#
    // Geometry Parameter
    const extern double dom_Lx;
    const extern double dom_Ly;

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
    const extern double traction_1;
    const extern double traction_2;
    const extern double traction_3;
    const extern double traction_4;

    // Panel Element Parameter
    const extern double len;     // Basic panel length

    // Internal Node Parameter
    const extern double spc;    // Basic panel length
}

#endif