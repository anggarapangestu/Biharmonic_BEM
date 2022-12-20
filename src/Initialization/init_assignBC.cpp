#include "initialization.hpp"

// ================================================================================
// ================================================================================
// The calculation of F and/or dFdn boundary value
void initialization::F_val_calc(element& elm, int ID){
    // Define the boundary value of F or dFdn assignment type
    int F_calc_type = 1;
    
    // Constant neumann value dFdn = 0
    if (F_calc_type == 1){
        double _guess_value = 1.0e-12;            // Still can be changed
        for (int i = 0; i < elm.num; i++){
            elm.dFdn[i] = _guess_value;
        }
    }
    
    // Method by Dewangga [2021]
    else if (F_calc_type == 2){
        for (int i = 0; i < elm.num; i++){
            elm.dFdn[i] = elm.Tx[i] * elm.xn[i] + elm.Ty[i] * elm.yn[i];
        }
    }
    
    // The case of dirichlet F = Tx + Ty
    else if(F_calc_type == 3){
        double _guess_value = Par::trac_t_y + Par::trac_r_x;
        for (int i = 0; i < elm.num; i++){
            elm.F[i] = _guess_value;
            elm.F_type[i] = false;
        }
    }
    
    // The case of full given F and dFdn value
    else if(F_calc_type == 4){
        for (int i = 0; i < elm.num; i++){
            elm.F[i] = 0.0e0;
            elm.dFdn[i] = 0.0e0;
        }
    }
    
}

// ================================================================================
// ================================================================================
// The calculation of phi and/or dpdn boundary value
void initialization::p_val_calc(element& elm, int ID){
    // ========= Type 1 Calculation ========= 
    // Calculating the dpdn

    // Initial value of phi gradient
    double dpdx_base = 0.0e0;
    double dpdy_base = 0.0e0;
    
    // Marking whether a base or inner geometry
    int dir;    // The element direction
    bool base;  // The marker for element level (BASE/INNER)
    if (ID >= 0){
        base = false;
        dir = -1;
    }else{
        base = true;
        dir = 1;
    }

    // Initialize the internal variable
    std::vector<double> dpdx_node(elm.num + 1,0.0e0);
    std::vector<double> dpdy_node(elm.num + 1,0.0e0);
    dpdx_node[0] = {dpdx_base};
    dpdy_node[0] = {dpdy_base};
    double _dpdx;
    double _dpdy;
    
    // [CALC 1] Calculate the integration of gradient phi at each node from given traction
    for (int i = 0; i < elm.num; i++){
        dpdx_node[i+1] = dpdx_node[i] - (dir * elm.L[i] * elm.Ty[i]);
        dpdy_node[i+1] = dpdy_node[i] + (dir * elm.L[i] * elm.Tx[i]);
    }

    // Update for symmetric extremes (work as well, but still need to change the value of F and phi)
    double max_x = 0, max_y = 0, min_x = 0, min_y = 0;
    for (int i = 0; i < dpdx_node.size(); i++){
        max_x = dpdx_node[i] > max_x ? dpdx_node[i] : max_x;
        max_y = dpdy_node[i] > max_y ? dpdy_node[i] : max_y;
        min_x = dpdx_node[i] < min_x ? dpdx_node[i] : min_x;
        min_y = dpdy_node[i] < min_y ? dpdy_node[i] : min_y;
    }
    
    // [CALC 2] Create a symmetric phi gradient along the panel
    for (int i = 0; i < dpdx_node.size(); i++){
        dpdx_node[i] += (max_x + min_x)/(-2.0);
        dpdy_node[i] += (max_y + min_y)/(-2.0);
    }

    // Assign the value of dpdn
    for (int i = 0; i < elm.num; i++){
        _dpdx = 0.5 * (dpdx_node[i+1] + dpdx_node[i]);
        _dpdy = 0.5 * (dpdy_node[i+1] + dpdy_node[i]);
        elm.dpdn[i] = _dpdx * elm.xn[i] + _dpdy * elm.yn[i];
        
        // In case the traction is zero then phi = 0 -> change to dirichlet (Need a validation)
        if (elm.dpdn[i] == 0){
            // Only done for internal boundary
            // if(base == false)
            {
                elm.p[i] = 0;
                elm.p_type[i] = false;
            }
        }
    }

    
    // ========= Type 2 Calculation ========= 
    // Calculate the phi (If the algorithm is known)
    bool given_P = false;
    if (given_P == true){
        // Internal variables definition
        double A, B, C;
        A = Par::trac_t_y / 2.0;
        B = 0.0;
        C = Par::trac_r_x / 2.0;
        
        for (int i = 0; i < elm.num; i++){
            // Hard coded from a given solution
            elm.p[i] = A * elm.xm[i] * elm.xm[i] + B * elm.xm[i] * elm.ym[i] + C * elm.ym[i] * elm.ym[i];
            elm.p_type[i] = false;
        }
    }
}


// ================================================================================
// ================================================================================
// The calculation of T and/or dTdn boundary value
void initialization::T_val_calc(element& elm, int ID){
    // Marking whether a base or inner geometry
    bool base;
    if (ID >= 0){
        base = false;
    }else{
        base = true;
    }

    // Define the boundary value of T or dTdn assignment type (based on the boundary panel location)
    int T_calc_type;
    if (base == true){
        // T_calc_type = 1;
        if (Par::G_type == 1){
            T_calc_type = 3;
        }else if (Par::G_type == 2){
            T_calc_type = 2;
        }
    }else{
        T_calc_type = 2;
    }
    
    // Constant neumann value dTdn = 0
    if (T_calc_type == 1){
        double _guess_value = 0.0e0;
        for (int i = 0; i < elm.num; i++){
            elm.dTdn[i] = _guess_value;
        }
    }
    
    // Definition of boundary value prescribed from setting
    else if (T_calc_type == 2){
        if (base == true){
            for (int i = 0; i < elm.num; i++){
                elm.T[i] = Par::Temp;
                elm.T_type[i] = Par::Temp_type;
            }
        }else{
            for (int i = 0; i < elm.num; i++){
                elm.T[i] = Par::In_temp[ID];
                elm.T_type[i] = Par::In_temp_type[ID];
            }
        }
    }
    
    // Definition of boundary value for rectangular base domain
    else if (T_calc_type == 3){
        // Find the geometry extremes position
        double max_x = 0, max_y = 0, min_x = 0, min_y = 0;
        for (int i = 0; i < elm.num; i++){
            max_x = elm.xm[i] > max_x ? elm.xm[i] : max_x;
            max_y = elm.ym[i] > max_y ? elm.ym[i] : max_y;
            min_x = elm.xm[i] < min_x ? elm.xm[i] : min_x;
            min_y = elm.ym[i] < min_y ? elm.ym[i] : min_y;
        }

        // Assign the boundary value provided from setting
        for (int i = 0; i < elm.num; i++){
            if (elm.xm[i] == max_x){
                // Right boundary
                elm.T[i] = Par::temp_r;
                elm.T_type[i] = Par::temp_type_r;
            }
            else if (elm.xm[i] == min_x){
                // Left boundary
                elm.T[i] = Par::temp_l;
                elm.T_type[i] = Par::temp_type_l;
            }
            else if (elm.ym[i] == min_y){
                // Bottom boundary
                elm.T[i] = Par::temp_b;
                elm.T_type[i] = Par::temp_type_b;
            }
            else if (elm.ym[i] == max_y){
                // Top boundary
                elm.T[i] = Par::temp_t;
                elm.T_type[i] = Par::temp_type_t;
            }   
        }
    }

}

// ================================================================================
// ================================================================================
// The calculation of F and p boundary value for biharmonic: phi(x,y) = x^4 - 12(xy)^2 + y^4
void initialization::F_bihar_val_calc(element& elm, int ID){
    // Define the boundary value of F or dFdn assignment type
    // F(x,y) = - 12(x^2 + y^2)
    // Del F(x,y) = - 24x i - 24y j
    int BCtype = 1; // 1:= Neumann, 2:= Dirichlet
    if (ID >= 0){
        BCtype = 2;
    }
    double x, y;
    double nx, ny, Fx, Fy;
    
    // All Neumann value given
    if (BCtype == 1){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            nx = elm.xn[i];
            ny = elm.yn[i];
            Fx = -24.0 * x;
            Fy = -24.0 * y;
            elm.dFdn[i] = Fx*nx + Fy*ny;
        }
    }
    
    // All Dirichlet value given
    else if (BCtype == 2){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            elm.F[i] = -12.0 * (x*x + y*y);
            elm.F_type[i] = false;
        }
    }
}

void initialization::p_bihar_val_calc(element& elm, int ID){
    // Define the boundary value of p or dpdn assignment type
    // phi(x,y) = x^4 - 12(xy)^2 + y^4
    // Del phi(x,y) = (4x^3 - 24xy^2)i + (-24x^2y + 4y^3)j
    int BCtype = 2; // 1:= Neumann, 2:= Dirichlet
    if (ID >= 0){
        BCtype = 1;
    }
    double x, y;
    double nx, ny, px, py;
    
    // All Neumann value given
    if (BCtype == 1){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            nx = elm.xn[i];
            ny = elm.yn[i];
            px = (4.0*x*x*x - 24.0*x*y*y);
            py = (-24.0*x*x*y + 4.0*y*y*y);
            elm.dpdn[i] = px*nx + py*ny;
        }
    }
    
    // All Dirichlet value given
    else if (BCtype == 2){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            elm.p[i] = x*x*x*x - 12.0*x*x*y*y + y*y*y*y;
            elm.p_type[i] = false;
        }
    }
}

// ================================================================================
// ================================================================================
// The calculation of p boundary value for laplace: phi(x,y) = x^2 - 4xy + y^2
void initialization::p_lap_val_calc(element& elm, int ID){
    // Define the boundary value of p or dpdn assignment type
    // phi(x,y) = x^2 - 4xy + y^2
    // Del phi(x,y) = (2x - 4y)i + (-4x + 2y)j
    int BCtype = 1; // 1:= Neumann, 2:= Dirichlet
    if (ID >= 0){
        BCtype = 2;
    }
    double x, y;
    double nx, ny, px, py;
    
    // All Neumann value given
    if (BCtype == 1){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            nx = elm.xn[i];
            ny = elm.yn[i];
            px = (2.0*x - 4.0*y);
            py = (-4.0*x + 2.0*y);
            elm.dpdn[i] = px*nx + py*ny;
        }
    }
    
    // All Dirichlet value given
    else if (BCtype == 2){
        for (int i = 0; i < elm.num; i++){
            x = elm.xm[i];
            y = elm.ym[i];
            elm.p[i] = x*x - 4.0*x*y + y*y;
            elm.p_type[i] = false;
        }
    }
}