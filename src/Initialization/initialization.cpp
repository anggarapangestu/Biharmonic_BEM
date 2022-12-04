#include "initialization.hpp"

// Generation of internal node
void initialization::generate_internal_node(intElement& intElm, const element& elm){
    // initialization generate internal node starting log
    printf("\nGenerate the internal node ...\n");

    // Procedure: \
       1. Generate the regular node distribution \
       2. Eleminate the node outside the domain\
       3. Save the position data into intElm\
       4. Neighbor searching

    // Generate the internal particle data
    if (Par::opt_int_init == 1){
        printf("<+> Regular internal node\n");
        this->internal_regular(intElm, elm);
    }else if (Par::opt_int_init == 2){
        printf("<+> Finer near panel internal node\n");
        this->internal_finer_near_panel(intElm, elm);
    }
    
    // Calculate the neighbor
    // The code for neighbor evaluation
}

// ===================================================
// Generation of boundary panel element
void initialization::generate_boundary_element(element& elm){
    // initialization generate boundary panel starting log
    printf("\nGenerate the boundary panel ...\n");

    // Procedure: \
       1. Determine the number of panel at each rectangle\
       2. Create the panel midpoint, normal, length sequencially in CCW direction\

    // Internal variable
    int nx = std::ceil(Par::dom_Lx/Par::len);
    int ny = std::ceil(Par::dom_Ly/Par::len);
    int num = 2 * (nx + ny);
    double lx = Par::dom_Lx/(double)nx;
    double ly = Par::dom_Ly/(double)ny;
    double x_piv = -Par::dom_Lx/2.0;
    double y_piv = -Par::dom_Ly/2.0;
    
    // Resize each boundary element variable
    elm.num = num;
    elm.xm.resize(num,0.0e0);
    elm.ym.resize(num,0.0e0);
    elm.L.resize(num,0.0e0);
    elm.xn.resize(num,0.0e0);
    elm.yn.resize(num,0.0e0);

    // Allocating the panel data (start from left-bottom corner)
    int _id = 0;
    // Bottom side -> moving to right
    for (int i = 0; i < nx; i++){
        // Midpoint coordinate position
        elm.xm[_id] = x_piv + (i + 0.5) * lx;
        elm.ym[_id] = y_piv;
        elm.L[_id] = lx;
        
        // Normal to bottom direction
        elm.xn[_id] = 0;
        elm.yn[_id] = -1;
        
        _id++;
    }

    // Right side -> moving upward
    for (int i = 0; i < ny; i++){
        // Midpoint coordinate position
        elm.xm[_id] = x_piv + Par::dom_Lx;
        elm.ym[_id] = y_piv + (i + 0.5) * ly;
        elm.L[_id] = ly;
        
        // Normal to right direction
        elm.xn[_id] = 1;
        elm.yn[_id] = 0;
        
        _id++;
    }

    // Top side -> moving to left
    for (int i = 0; i < nx; i++){
        // Midpoint coordinate position
        elm.xm[_id] = x_piv + Par::dom_Lx - (i + 0.5) * lx;
        elm.ym[_id] = y_piv + Par::dom_Ly;
        elm.L[_id] = lx;
        
        // Normal to right direction
        elm.xn[_id] = 0;
        elm.yn[_id] = 1;
        
        _id++;
    }
    
    // Left side -> moving downward
    for (int i = 0; i < ny; i++){
        // Midpoint coordinate position
        elm.xm[_id] = x_piv;
        elm.ym[_id] = y_piv + Par::dom_Ly - (i + 0.5) * ly;
        elm.L[_id] = ly;
        
        // Normal to right direction
        elm.xn[_id] = -1;
        elm.yn[_id] = 0;
        
        _id++;
    }

    // Resize the other variable and set as 0.0
    elm.F.resize(num,0.0e0);
    elm.p.resize(num,0.0e0);
    elm.dFdn.resize(num,0.0e0);
    elm.dpdn.resize(num,0.0e0);
}

// ===================================================
// Calculate the initial condition at each panel
void initialization::calculate_initial_condition(element& elm){
    // The code lies here ...
}