#include "initialization.hpp"

// ==================================================================================
// ==================================================================================
// Regular internal node distribution
void initialization::internal_regular(intElement& intElm, const element& elm, const std::vector<element>& in_elm){
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
    std::vector<double> _size;  // Temporary spacing size
    std::vector<double> R_pos;  // Temporary nearest distance
    std::vector<bool> _assign;  // Temporary flag for assign the node data
    int nx, ny, num, _num;
    double x_piv, y_piv;

    // ========= PROCESS 1 =========
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
    _num = num;
    x_piv = (max[0] + min[0])/2.0 - nx * Par::spc/2.0;
    y_piv = (max[1] + min[1])/2.0 - ny * Par::spc/2.0;

    // ========= PROCESS 2 =========
    // Allocating the node position
    x_pos.resize(num,0.0e0);
    y_pos.resize(num,0.0e0);
    _size.resize(num,Par::spc);
    R_pos.resize(num,0.0e0);
    _assign.resize(num,true);
    int _id = 0;
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            x_pos[_id] = x_piv + (j + 0.5) * Par::spc;
            y_pos[_id] = y_piv + (i + 0.5) * Par::spc;
            _id++;
        }
    }
    
    // ========= PROCESS 3 & 4 =========
    // Eleminate for node outside the base geometry
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
        if (R_n > -Par::spc/2.0){
            // Node lies outside the domain, then it need to be deleted
            _assign[i] = false;
            _num--;
        }
    }

    // Eleminate for node inside the inner geometry
    for (int ID = 0; ID < Par::N_Gin; ID ++){
        for (int i = 0; i < num; i++){
            // Do not re-evaluate the erased flagged element
            if (_assign[i] == false){
                continue;
            }
            
            // Take the position of current evaluated node
            double _x = x_pos[i];
            double _y = y_pos[i];
            
            // Find the nearest panel
            // Initialize by using the first panel parameter
            int nst_id = 0;
            double _R;
            double min_R;
            min_R = std::sqrt(std::pow(_x - in_elm[ID].xm[0], 2) + std::pow(_y - in_elm[ID].ym[0], 2));
            for (int j = 0; j < in_elm[ID].num; j++){
                // Temporary distance to the consecutive panel
                _R = std::sqrt(std::pow((in_elm[ID].xm[j] - _x), 2) + std::pow((in_elm[ID].ym[j] - _y ), 2));
                
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
            R_n = in_elm[ID].xn[nst_id] * (_x - in_elm[ID].xm[nst_id]) + in_elm[ID].yn[nst_id] * (_y - in_elm[ID].ym[nst_id]);
            if (R_n > -Par::spc/2.0){
                // Node lies outside the domain, then it need to be deleted
                _assign[i] = false;
                _num--;
            }
        }
    }

    // ========= PROCESS 5 =========
    // Initialize each internal Element variable
    intElm.x.resize(_num, 0.0e0);
    intElm.y.resize(_num, 0.0e0);
    intElm.s.resize(_num, 0.0e0);
    intElm.R.resize(_num, 0.0e0);
    intElm.num = _num;

    // Assign the node position data into the internal element variable
    _id = 0;
    for (int i = 0; i < num; i++){
        if(_assign[i] == true){
            intElm.x[_id] = x_pos[i];
            intElm.y[_id] = y_pos[i];
            intElm.s[_id] = _size[i];
            intElm.R[_id] = R_pos[i];
            _id++;
        }
    }
}

// ==================================================================================
// ==================================================================================
// Finer near panel internal node distribution
void initialization::internal_finer_near_panel(intElement& intElm, const element& elm, const std::vector<element>& in_elm){
    // Under construction
    // Procedure:\
       1. Find the extremes boundary location\
       2. Generate a rectangular uniform distribution\
       3. Refine the element near the boundary\
       4. Check the location whether the node is inside the domain\
       5. Delete the variable outside the domain\
       6. Assign the node position into the intElm
    
    // Internal variable
    double max[2] = {0, 0};     // Store the maximum position of x and y respectively
    double min[2] = {0, 0};     // Store the minimum position of x and y respectively
    std::vector<double> x_pos;  // Temporary x_variable
    std::vector<double> y_pos;  // Temporary y_variable
    std::vector<double> _size;  // Temporary internal spacing
    std::vector<double> R_pos;  // Temporary nearest distance
    std::vector<bool> _assign;  // Temporary flag for assign the data into list
    std::vector<bool> _devide;  // Temporary flag for devided element
    int nx, ny, num, _num;
    double x_piv, y_piv;

    // ========= PROCESS 1 =========
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
    _num = num;
    x_piv = (max[0] + min[0])/2.0 - nx * Par::spc/2.0;
    y_piv = (max[1] + min[1])/2.0 - ny * Par::spc/2.0;

    // ========= PROCESS 2 =========
    // Allocating the node position
    x_pos.resize(num,0.0e0);
    y_pos.resize(num,0.0e0);
    _size.resize(num,Par::spc);
    R_pos.resize(num,0.0e0);
    _assign.resize(num,true);
    _devide.resize(num,false);
    int _id = 0;
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            x_pos[_id] = x_piv + (j + 0.5) * Par::spc;
            y_pos[_id] = y_piv + (i + 0.5) * Par::spc;
            _id++;
        }
    }
    
    // ========= PROCESS 3 - 5 =========
    // Eleminate for node outside the base geometry
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

        // Performing divide and conqueror
        if (R_pos[i] < Par::spc * Par::dist_fac){
            // Delete current node then update by new divided node
            _assign[i] = false;
            _num--;

            int _DnCx[4] = {-1,1,1,-1};
            int _DnCy[4] = {-1,-1,1,1};

            for (int j = 0; j < 4; j++){
                // Evaluate the distance from boundary again
                double __x = _x + Par::spc/4.0 * _DnCx[j];
                double __y = _y + Par::spc/4.0 * _DnCy[j];
                
                // Find the nearest panel
                // Initialize by using the first panel parameter
                nst_id = 0;
                min_R = std::sqrt(std::pow(__x - elm.xm[0], 2) + std::pow(__y - elm.ym[0], 2));
                for (int j = 0; j < elm.num; j++){
                    // Temporary distance to the consecutive panel
                    _R = std::sqrt(std::pow((elm.xm[j] - __x), 2) + std::pow((elm.ym[j] - __y ), 2));
                    
                    // Update the id of nearest panel and the minimum distance
                    if (_R < min_R){
                        nst_id = j;
                        min_R = _R;
                    }
                }

                // Determine whether the node lies inside or outside the domain
                double R_n;  // The normal direction which is vec(r) dot hat(n)
                R_n = elm.xn[nst_id] * (__x - elm.xm[nst_id]) + elm.yn[nst_id] * (__y - elm.ym[nst_id]);
                if (R_n < 0){
                    // If lies inside the domain assign into the new
                    x_pos.push_back(__x);
                    y_pos.push_back(__y);
                    _size.push_back(Par::spc / 2.0);
                    R_pos.push_back(min_R);
                    _assign.push_back(true);
                    _devide.push_back(true);
                    _num++;
                }
            }
            
        }else{
            // Determine whether the node lies inside or outside the domain
            double R_n;  // The normal direction which is vec(r) dot hat(n)
            R_n = elm.xn[nst_id] * (_x - elm.xm[nst_id]) + elm.yn[nst_id] * (_y - elm.ym[nst_id]);
            if (R_n > -Par::spc/2.0){
                // Node lies outside the domain, then it need to be deleted
                _assign[i] = false;
                _num--;
            }
        }
    }
    // Update the total number of element
    num = _assign.size();

    // Eleminate for node inside the inner geometry
    for (int ID = 0; ID < Par::N_Gin; ID ++){
        for (int i = 0; i < num; i++){
            // Do not re-evaluate the erased flagged element
            if (_assign[i] == false){
                continue;
            }

            // Take the position of current evaluated node
            double _x = x_pos[i];
            double _y = y_pos[i];
            
            // Find the nearest panel
            // Initialize by using the first panel parameter
            int nst_id = 0;
            double _R;
            double min_R;
            min_R = std::sqrt(std::pow(_x - in_elm[ID].xm[0], 2) + std::pow(_y - in_elm[ID].ym[0], 2));
            for (int j = 0; j < in_elm[ID].num; j++){
                // Temporary distance to the consecutive panel
                _R = std::sqrt(std::pow((in_elm[ID].xm[j] - _x), 2) + std::pow((in_elm[ID].ym[j] - _y ), 2));
                
                // Update the id of nearest panel and the minimum distance
                if (_R < min_R){
                    nst_id = j;
                    min_R = _R;
                }
            }
            
            // Assign the minimum distance into the list
            R_pos[i] = min_R;

            // Performing divide and conqueror
            if ((R_pos[i] < Par::spc * Par::dist_fac) && (_devide[i] == false)){
                // Delete current node then update by new divided node
                _assign[i] = false;
                _num--;

                int _DnCx[4] = {-1,1,1,-1};
                int _DnCy[4] = {-1,-1,1,1};

                for (int j = 0; j < 4; j++){
                    // Evaluate the distance from boundary again
                    double __x = _x + Par::spc/4.0 * _DnCx[j];
                    double __y = _y + Par::spc/4.0 * _DnCy[j];
                    
                    // Find the nearest panel
                    // Initialize by using the first panel parameter
                    nst_id = 0;
                    min_R = std::sqrt(std::pow(__x - in_elm[ID].xm[0], 2) + std::pow(__y - in_elm[ID].ym[0], 2));
                    for (int j = 0; j < in_elm[ID].num; j++){
                        // Temporary distance to the consecutive panel
                        _R = std::sqrt(std::pow((in_elm[ID].xm[j] - __x), 2) + std::pow((in_elm[ID].ym[j] - __y ), 2));
                        
                        // Update the id of nearest panel and the minimum distance
                        if (_R < min_R){
                            nst_id = j;
                            min_R = _R;
                        }
                    }

                    // Determine whether the node lies inside or outside the domain
                    double R_n;  // The normal direction which is vec(r) dot hat(n)
                    R_n = in_elm[ID].xn[nst_id] * (__x - in_elm[ID].xm[nst_id]) + in_elm[ID].yn[nst_id] * (__y - in_elm[ID].ym[nst_id]);
                    if (R_n < 0){
                        // If lies inside the domain assign into the new
                        x_pos.push_back(__x);
                        y_pos.push_back(__y);
                        _size.push_back(Par::spc / 2.0);
                        R_pos.push_back(min_R);
                        _assign.push_back(true);
                        _devide.push_back(true);
                        _num++;
                    }
                }
                
            }else{
                // Determine whether the node lies inside or outside the domain
                double R_n;  // The normal direction which is vec(r) dot hat(n)
                R_n = in_elm[ID].xn[nst_id] * (_x - in_elm[ID].xm[nst_id]) + in_elm[ID].yn[nst_id] * (_y - in_elm[ID].ym[nst_id]);
                if (R_n > -Par::spc/2.0){
                    // Node lies outside the domain, then it need to be deleted
                    _assign[i] = false;
                    _num--;
                }
            }
        }
        
        // Update the total number of element
        num = _assign.size();
    }
    
    // ========= PROCESS 5 =========
    // Initialize each internal Element variable
    intElm.x.resize(_num, 0.0e0);
    intElm.y.resize(_num, 0.0e0);
    intElm.s.resize(_num, 0.0e0);
    intElm.R.resize(_num, 0.0e0);
    intElm.num = _num;

    // Assign the node position data into the internal element variable
    _id = 0;
    for (int i = 0; i < num; i++){
        if(_assign[i] == true){
            intElm.x[_id] = x_pos[i];
            intElm.y[_id] = y_pos[i];
            intElm.s[_id] = _size[i];
            intElm.R[_id] = R_pos[i];
            _id++;
        }
    }
}
