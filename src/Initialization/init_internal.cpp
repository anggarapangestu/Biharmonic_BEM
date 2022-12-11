#include "initialization.hpp"

// Regular internal node distribution
void initialization::internal_regular(intElement& intElm, const element& elm){
    // Under construction
    
    // Procedure:\
       1. Find the extremes boundary location\
       2. Generate a rectangular uniform distribution\
       3. Check the location whether the node is inside the domain\
       4. Delete the variable outside the domain\
       5. Assign the node position into the intElm
    
    // Internal variable
    double max[2] = {0, 0};     // Store the maximum position of x and y respectively
    double min[2] = {0, 0};     // Store the minimum position of x and y respectively
    std::vector<double> x_pos;  // Temporary x_variable
    std::vector<double> y_pos;  // Temporary y_variable
    std::vector<double> R_pos;  // Temporary nearest distance
    int nx, ny, num;
    double x_piv, y_piv;

    // Determine the extremes boundary location
    for (int i = 0; i < elm.num; i++){
        max[0] = max[0] > elm.xm[i] ? max[0] : elm.xm[i];
        max[1] = max[1] > elm.ym[i] ? max[1] : elm.ym[i];
        min[0] = min[0] < elm.xm[i] ? min[0] : elm.xm[i];
        min[1] = min[1] < elm.ym[i] ? min[1] : elm.ym[i];
    }

    // Determine the node number and spacing
    nx = 5 + std::ceil((max[0] - min[0])/Par::spc);
    ny = 5 + std::ceil((max[1] - min[1])/Par::spc);
    num = nx * ny;
    x_piv = (max[0] + min[0])/2.0 - nx * Par::spc/2.0;
    y_piv = (max[1] + min[1])/2.0 - ny * Par::spc/2.0;

    // Allocating the node position
    x_pos.resize(num,0.0e0);
    y_pos.resize(num,0.0e0);
    R_pos.resize(num,0.0e0);
    int _id = 0;
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            x_pos[_id] = x_piv + (j + 0.5) * Par::spc;
            y_pos[_id] = y_piv + (i + 0.5) * Par::spc;
            _id++;
        }
    }
    
    // Eleminate for node outside the domain
    for (int i = 0; i < num; i++){
        // Take the position of current evaluated node
        double _x = x_pos[i];
        double _y = y_pos[i];
        
        // Find the nearest panel
        // Initialize by using the first panel parameter
        int nst_id = 0;
        double _R;
        double min_R;
        min_R = std::sqrt(std::pow(_x - elm.xm[0], 2) + std::pow(_y - elm.ym[0], 2));
        for (int j = 0; j < elm.num; j++){
            // Temporary distance to the consecutive panel
            _R = std::sqrt(std::pow((elm.xm[j] - _x), 2) + std::pow((elm.ym[j] - _y ), 2));
            
            // Update the id of nearest panel and the minimum distance
            if (_R < min_R){
                nst_id = j;
                min_R = _R;
            }
        }
        
        // Assign the minimum distance into the list
        R_pos[i] = min_R;

        // Determine whether the node lies inside or outside the domain
        double R_n;  // The normal direction which is vec(r) dot hat(n)
        R_n = elm.xn[nst_id] * (_x - elm.xm[nst_id]) + elm.yn[nst_id] * (_y - elm.ym[nst_id]);
        if (R_n > 0){
            // Node lies outside the domain, then it need to be deleted
            x_pos.erase(x_pos.begin() + i);
            y_pos.erase(y_pos.begin() + i);
            R_pos.erase(R_pos.begin() + i);
            i--; num--;
        }
    }
    
    // Assign the node position data into the internal element variable
    intElm.x = x_pos;
    intElm.y = y_pos;
    intElm.R = R_pos;
    intElm.num = num;
}

// Finer near panel internal node distribution
void initialization::internal_finer_near_panel(intElement& intElm, const element& elm){
    // Under construction
}