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

#ifndef BOUNDARY_ELEMENT_CALCULATION
#include "src/BEM/BEM.hpp"
#endif

#ifndef PHYSICAL_PROPERTIES_CALCULATION
#include "src/PropertyCalc/physical_prop.hpp"
#endif

int main(){
/* The Flow in this Biharmonic Program Code
    INITIALIZATION
    (1) Element Generation :
        * Boundary element for base geometry boundary (BASE)
        * Boundary element for internal boundary (INNER)
        * Nodes at internal domain (INTERNAL)
    (2) Boundary Value Calculation 
        I.  Biharmonic (sim type 1)
            * Traction
            * Displacement (still not developed)
        II. Temperature (sim type 2)
            * Temperature
    (3) Boundary element initialization for BEM calculation
        -> Combine each boundary element (BASE and INNER)

    SOLVER
    I.  Biharmonic (sim type 1)
        (1) Solving Laplace of F
        (2) Solving Poison of phi
    II. Temperature (sim type 2)
        (1) Solving Laplace of T

    POST PROCESSING
    (1) Calculate the domain properties
        I.  Biharmonic (sim type 1)
            -> Calculate the phi inside domain
            -> Calculate the properties (sigma_x, sigma_y, tau_xy, u, v, ...) inside domain
        II. Temperature (sim type 2)
            -> Calculate the T inside domain
    (2) Write Data
*/

    // The computational time counter
    clock_t _time = clock();
    
    // ======================================
    // ====== INTERNAL VARIABLE REGION ======
    // ======================================

    // Initialize the storage data for boundary element and internal node
    element PanelElement;               // BASE element
    std::vector<element> InnerElement;  // INNER element
    intElement InternalNode;            // INTERNAL node
    
    // Initialize the function method
    initialization init;
    calcBEM BEMstep;
    propertyCalc prop_step;
    dataSaving save;

    // Print simulation log
    save.simulation_log();


    // ============================
    // ====== INITIALIZATION ======
    // ============================
    // Initialization HEADER
    std::cout << std::endl;
    std::cout << "#=================================================#\n";
    std::cout << "+---------------- INITIALIZATION -----------------+\n";
    std::cout << "#=================================================#\n";

    // Initialization generate the panel element and intenral node
    init.generate_boundary_element(PanelElement, InnerElement);
    init.generate_internal_node(InternalNode, PanelElement, InnerElement);

    // Calculate the boundary value
    init.calculate_boundary_condition(PanelElement, InnerElement);
    
    
    // ===============================
    // ======= BEM CALCULATION =======
    // ===============================
    // BEM calculation HEADER
    std::cout << std::endl;
    std::cout << "#=================================================#\n";
    std::cout << "+---------------- BEM SOLVER LOG -----------------+\n";
    std::cout << "#=================================================#\n";

    // Initialize the BEM paramter
    BEMstep.Define_BEM(PanelElement, InnerElement);
    
    // // Testing
    // BEMstep.CALC_theta();
    // BEMstep.TEST_BEM(PanelElement, InnerElement);

    // Biharmonic simulation solver
    if (Par::opt_sim_type == 1){
        // Calculate the other boundary element value by BEM
        BEMstep.solve_F(PanelElement, InnerElement);
        BEMstep.solve_phi(PanelElement, InnerElement);

        // Calculate the phi value at the internal domain node
        BEMstep.calculate_internal_phi(InternalNode, PanelElement, InnerElement);
        
        // Additional of phi analytical solution
        prop_step.phi_analytic_biaxial(InternalNode);
    }
    // Heat transfer simulation solver
    else if (Par::opt_sim_type == 2){
        // Calculate the other boundary element value by BEM
        BEMstep.solve_T(PanelElement, InnerElement);
        
        // Calculate the T value at the internal domain node
        BEMstep.calculate_internal_T(InternalNode, PanelElement, InnerElement);
    }

    // =======================================
    // ======= CALCULATING PROPERTIES ========
    // =======================================
    if (Par::opt_sim_type == 1){
        // BEM calculation HEADER
        std::cout << std::endl;
        std::cout << "#=================================================#\n";
        std::cout << "+------------ PROPERTIES CALCULATION -------------+\n";
        std::cout << "#=================================================#\n";
        // Calculate all properties for biharmonic at the domain
        prop_step.calculate_property(InternalNode);
    }
    

    // ============================
    // ======= SAVING DATA ========
    // ============================
    // Write file HEADER
    std::cout << std::endl;
    std::cout << "+---------------- SAVING DATA LOG ----------------+\n";
    // Write the element data
    
    // Saving biharmonic data
    if (Par::opt_sim_type == 1){
        save.write_internal_data(InternalNode);
        save.write_BEM_data(PanelElement, InnerElement);
    }
    
    // Saving temperature data
    else if (Par::opt_sim_type == 2){
        save.write_internal_data_temp(InternalNode);
        save.write_BEM_data_temp(PanelElement, InnerElement);
    }

    // Displaying the computational time
    _time = clock() - _time;
	printf("\nTOTAL COMPUTATIONAL TIME               [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
    printf("\n#============== SIMULATION FINISHED ==============#\n\n");

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

/*  >>>>> WRITING SYSTEM <<<<<
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