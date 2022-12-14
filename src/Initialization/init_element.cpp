#include "initialization.hpp"

// ==================================================================================
// ==================================================================================
// Boundary element initialization for rectangular geometry
void initialization::element_rectangular(element& elm, int ID){
    // Procedure: \
       1. Determine the number of panel at each rectangle\
       2. Create the panel midpoint, normal, length sequencially in CCW direction

    // Internal variables
    bool base;      // The flag of base geometry

    if (ID >= 0){    // Internal geometry
        base = false;
    }else{          // Base geometry
        base = true;
    }

    // Internal variable
    int nx, ny, num;        // Number of panel
    double xLen, yLen;      // Geometry length
    double lx, ly;          // Panel length
    double x_piv, y_piv;    // Pivot position
    double x_cen, y_cen;    // Center position

    // Update the geometry length and center position
    if (base == true){  // Base geometry
        // Geometry length
        xLen = Par::dom_Lx;
        yLen = Par::dom_Ly;
        
        // Center position
        x_cen = 0.0e0;
        y_cen = 0.0e0;
    }else{              // Inner geometry
        // Geometry length
        xLen = Par::Gin_Xlen[ID];
        yLen = Par::Gin_Ylen[ID];
        
        // Center position
        x_cen = Par::Gin_Xcen_pos[ID];
        y_cen = Par::Gin_Ycen_pos[ID];
    }

    // Calculate the internal variable
    nx = std::ceil(xLen/Par::len);  // Number of panel
    ny = std::ceil(yLen/Par::len);
    num = 2 * (nx + ny);
    lx = xLen/(double)nx;        // Length of panel
    ly = yLen/(double)ny;
    x_piv = -xLen/2.0;           // Pivot coordinate
    y_piv = -yLen/2.0;

    // ==========================================
    // +----- ELEMENT PROPERTY CALCULATION -----+
    // ==========================================

    // Partiate into two parts
    if (base == true){      // The base geometry
        // Resize each boundary element variable
        elm.num = num;
        elm.xm.resize(num,0.0e0);
        elm.ym.resize(num,0.0e0);
        elm.L.resize(num,0.0e0);
        elm.xn.resize(num,0.0e0);
        elm.yn.resize(num,0.0e0);
        elm.Tx.resize(num,0.0e0);
        elm.Ty.resize(num,0.0e0);

        // Allocating the panel data (start from left-bottom corner) (CCW Direction)
        int _id = 0;
        // Bottom side -> moving to right
        for (int i = 0; i < nx; i++){
            // Midpoint coordinate position
            elm.xm[_id] = x_cen + x_piv + (i + 0.5) * lx;
            elm.ym[_id] = y_cen + y_piv;
            elm.L[_id] = lx;
            
            // Normal to bottom direction
            elm.xn[_id] = 0;
            elm.yn[_id] = -1;
            
            // Assign the traction stress
            elm.Tx[_id] = Par::trac_b_x;
            elm.Ty[_id] = Par::trac_b_y;
            
            _id++;
        }

        // Right side -> moving upward
        for (int i = 0; i < ny; i++){
            // Midpoint coordinate position
            elm.xm[_id] = x_cen + x_piv + xLen;
            elm.ym[_id] = y_cen + y_piv + (i + 0.5) * ly;
            elm.L[_id] = ly;
            
            // Normal to right direction
            elm.xn[_id] = 1;
            elm.yn[_id] = 0;

            // Assign the traction stress
            elm.Tx[_id] = Par::trac_r_x;
            elm.Ty[_id] = Par::trac_r_y;
            
            _id++;
        }

        // Top side -> moving to left
        for (int i = 0; i < nx; i++){
            // Midpoint coordinate position
            elm.xm[_id] = x_cen + x_piv + xLen - (i + 0.5) * lx;
            elm.ym[_id] = y_cen + y_piv + yLen;
            elm.L[_id] = lx;
            
            // Normal to top direction
            elm.xn[_id] = 0;
            elm.yn[_id] = 1;

            // Assign the traction stress
            elm.Tx[_id] = Par::trac_t_x;
            elm.Ty[_id] = Par::trac_t_y;
            
            _id++;
        }
        
        // Left side -> moving downward
        for (int i = 0; i < ny; i++){
            // Midpoint coordinate position
            elm.xm[_id] = x_cen + x_piv;
            elm.ym[_id] = y_cen + y_piv + yLen - (i + 0.5) * ly;
            elm.L[_id] = ly;
            
            // Normal to left direction
            elm.xn[_id] = -1;
            elm.yn[_id] = 0;

            // Assign the traction stress
            elm.Tx[_id] = Par::trac_l_x;
            elm.Ty[_id] = Par::trac_l_y;
            
            _id++;
        }
    }
    else{                   // The internal geometry
        // Create a temporary element
        element _elm;

        // Resize each boundary element variable
        _elm.num = num;
        _elm.xm.resize(num,0.0e0);
        _elm.ym.resize(num,0.0e0);
        _elm.L.resize(num,0.0e0);
        _elm.xn.resize(num,0.0e0);
        _elm.yn.resize(num,0.0e0);
        _elm.Tx.resize(num,0.0e0);
        _elm.Ty.resize(num,0.0e0);

        // Allocating the panel data (start from left-bottom corner) (CW Direction)
        int _id = 0;
        // Left side -> moving upward
        for (int i = 0; i < ny; i++){
            // Midpoint coordinate position
            _elm.xm[_id] = x_piv;
            _elm.ym[_id] = y_piv + (i + 0.5) * ly;
            _elm.L[_id] = ly;
            
            // Normal to left direction
            _elm.xn[_id] = 1;
            _elm.yn[_id] = 0;

            // Assign the traction stress
            _elm.Tx[_id] = -Par::In_pressure[ID];
            _elm.Ty[_id] = 0.0e0;
            
            _id++;
        }

        // Top side -> moving to right
        for (int i = 0; i < nx; i++){
            // Midpoint coordinate position
            _elm.xm[_id] = x_piv + (i + 0.5) * lx;
            _elm.ym[_id] = y_piv + yLen;
            _elm.L[_id] = lx;
            
            // Normal to top direction
            _elm.xn[_id] = 0;
            _elm.yn[_id] = -1;
            
            // Assign the traction stress
            _elm.Tx[_id] = 0.0e0;
            _elm.Ty[_id] = Par::In_pressure[ID];
            
            _id++;
        }

        // Right side -> moving downward
        for (int i = 0; i < ny; i++){
            // Midpoint coordinate position
            _elm.xm[_id] = x_piv + xLen;
            _elm.ym[_id] = y_piv + yLen - (i + 0.5) * ly;
            _elm.L[_id] = ly;
            
            // Normal to right direction
            _elm.xn[_id] = -1;
            _elm.yn[_id] = 0;

            // Assign the traction stress
            _elm.Tx[_id] = Par::In_pressure[ID];
            _elm.Ty[_id] = 0.0e0;
            
            _id++;
        }

        // Bottom side -> moving to left
        for (int i = 0; i < nx; i++){
            // Midpoint coordinate position
            _elm.xm[_id] = x_piv + xLen - (i + 0.5) * lx;
            _elm.ym[_id] = y_piv;
            _elm.L[_id] = lx;
            
            // Normal to bottom direction
            _elm.xn[_id] = 0;
            _elm.yn[_id] = 1;

            // Assign the traction stress
            _elm.Tx[_id] = 0.0e0;
            _elm.Ty[_id] = -Par::In_pressure[ID];
            
            _id++;
        }

        // ==============================
        // +----- ELEMENT ROTATION -----+
        // ==============================

        // Perform rotation here
        double _rot_angle = Par::Gin_Rot[ID] * M_PI / 180;
        if (_rot_angle != 0.0){
            // Temporary element
            element _rot_elm;

            // Resize each boundary element variable
            _rot_elm.num = _elm.num;
            _rot_elm.L = _elm.L;
            _rot_elm.xm.resize(num,0.0e0);
            _rot_elm.ym.resize(num,0.0e0);
            _rot_elm.xn.resize(num,0.0e0);
            _rot_elm.yn.resize(num,0.0e0);
            _rot_elm.Tx.resize(num,0.0e0);
            _rot_elm.Ty.resize(num,0.0e0);

            for (int i = 0; i < _elm.num; i++){
                _rot_elm.xm[i] = _elm.xm[i] * std::cos(_rot_angle) - _elm.ym[i] * std::sin(_rot_angle);
                _rot_elm.ym[i] = _elm.ym[i] * std::cos(_rot_angle) + _elm.xm[i] * std::sin(_rot_angle);
                _rot_elm.xn[i] = _elm.xn[i] * std::cos(_rot_angle) - _elm.yn[i] * std::sin(_rot_angle);
                _rot_elm.yn[i] = _elm.yn[i] * std::cos(_rot_angle) + _elm.xn[i] * std::sin(_rot_angle);
                _rot_elm.Tx[i] = _elm.Tx[i] * std::cos(_rot_angle) - _elm.Ty[i] * std::sin(_rot_angle);
                _rot_elm.Ty[i] = _elm.Ty[i] * std::cos(_rot_angle) + _elm.Tx[i] * std::sin(_rot_angle);
            }
            
            // Update the element
            _elm = _rot_elm;
        }

        // =================================
        // +----- ELEMENT TRANSLATION -----+
        // =================================
        // Perform translation from the center
        for (int i = 0; i < _elm.num; i++){
            _elm.xm[i] += x_cen;
            _elm.ym[i] += y_cen;
        }
        
        // Update the current element
        // **************************
        elm = _elm;
    }

}

// ==================================================================================
// ==================================================================================
// Boundary element initialization for circular geometry
void initialization::element_circular(element& elm, int ID){
    // Procedure: \
       1. Determine the number of panel at each rectangle\
       2. Create the panel midpoint, normal, length sequencially in CCW direction (base) or CW direction (inner)

    // Base vs Internal Geometry variable
    bool base;      // The flag of base geometry
    int _dir;       // The direction of base particle
    double _press;  // The pressure traction
    
    // Update the geometry length and center position
    if (ID >= 0){   // Internal geometry
        base = false;
    }else{          // Base geometry
        base = true;
    }

    // Internal variable
    int num;            // Number of panel
    double _a, _b;      // Geometry parameter
    double x_cen, y_cen;// Center position
    bool isCircular;    // Eliptical or Circular
    // Circular parameter
    double _dt;         // Base theta variation    
    // Eliptical parameter
    std::vector<double> _alpha;   // List of alpha variation

    // ============================================
    // +----- INTERNAL PARAMETER CALCULATION -----+
    // ============================================

    if (base == true){  // Base geometry
        // Geometry length
        _a = Par::dom_Lx/2;
        _b = Par::dom_Ly/2;

        // Calculate the internal variable
        _dir = 1;
        _press = Par::trac_press;
        
        // Center position
        x_cen = 0.0e0;
        y_cen = 0.0e0;
    }else{              // Inner geometry
        // Geometry length
        _a = Par::Gin_Xlen[ID]/2;
        _b = Par::Gin_Ylen[ID]/2;

        // Calculate the internal variable
        _dir = -1;
        _press = Par::In_pressure[ID];
        
        // Center position
        x_cen = Par::Gin_Xcen_pos[ID];
        y_cen = Par::Gin_Ycen_pos[ID];
    }

    // Determine wheter circular or eliptical
    if (_a == _b){
        isCircular = true;
        
        // Calculate circular parameter
        num = std::ceil((M_PI * (_a + _b))/Par::len);
        _dt = _dir * 2 * M_PI / (double)num;            // in radian
    }else{
        isCircular = false;
        
        // Calculate eliptical parameter
        double _Rmin = _a < _b ? _a : _b;
        double _Rmax = _a > _b ? _a : _b;
        bool _add90deg = _a > _b ? true : false;

        // Base alpha variation
        double _alp = Par::len/_Rmin;
        _alpha = {0.0};
        
        // Generate the 1/4 rotation
        int _num = 0;
        while(true){
            double Tfac = std::pow((1 + (_Rmin - _Rmax)/_Rmax * std::pow(std::sin(_alpha[_num]),2)),1.5);
            double _alp_temp = _alpha[_num] + _alp * Tfac;
            if (_alp_temp > M_PI/2){
                if (_alpha[_num] + _alp * Tfac / 2.0 > M_PI/2){
                    _alpha[_num] = M_PI/2;
                    break;
                }
                _alpha.push_back(M_PI/2);
                _num++;
                break;
            }
            _alpha.push_back(_alp_temp);
            _num++;
        }

        // Generate the 1/2 rotation
        for (int i = 0; i < _num; i++){
            _alpha.push_back(M_PI - _alpha[_num - (i + 1)]);
        }
        
        // Generate the 1 rotation
        _num = 2 * _num;
        for (int i = 0; i < _num; i++){
            _alpha.push_back(M_PI + _alpha[i + 1]);
        }

        // Perform rotation 90 degree only for _a > _b
        _num = 2 * _num + 1;
        if (_add90deg == true){
            std::vector<double> __alpha;   // List of alpha variation

            __alpha.push_back(0.0);
            for (int i = 0; i < _num; i++){
                // __alpha.push_back(_alpha[i] - M_PI/2);
                if (_alpha[i] - M_PI/2 > 0){
                    __alpha.push_back(_alpha[i] - M_PI/2);
                }
            }
            
            for (int i = 1; ; i++){
                __alpha.push_back(_alpha[i] + 3 * M_PI / 2);
                if (_alpha[i] - M_PI/2 >= 0){break;}
            }
            _alpha = __alpha;
        }

        // Update the number of panel
        num = _num - 1;
    }

    // ==========================================
    // +----- ELEMENT PROPERTY CALCULATION -----+
    // ==========================================
    
    // Temporary element
    element _elm;

    // Resize each boundary element variable
    _elm.num = num;
    _elm.xm.resize(num,0.0e0);
    _elm.ym.resize(num,0.0e0);
    _elm.L.resize(num,0.0e0);
    _elm.xn.resize(num,0.0e0);
    _elm.yn.resize(num,0.0e0);
    _elm.Tx.resize(num,0.0e0);
    _elm.Ty.resize(num,0.0e0);

    // Partiate into two parts
    // Allocating the panel data (start from left-bottom corner) (CCW Direction)
    int _id = 0;
    double _t;
    double _l;          // Panel length
    double r_i = _a;
    double x_i = _a;
    double y_i = 0.0e0;
    double r_f, x_f, y_f;
    
    // Bottom side -> moving to right
    for (int i = 0; i < num; i++){
        // Update the current node
        if (isCircular == true){
            _t = (i + 1) * _dt;
        }else{
            _t = _dir * _alpha[i + 1];
        }

        r_f = std::sqrt(1/(std::pow((std::cos(_t)/_a),2) + std::pow((std::sin(_t)/_b),2)));
        x_f = r_f * std::cos(_t);
        y_f = r_f * std::sin(_t);

        // Midpoint coordinate position
        _elm.xm[_id] = (x_f + x_i)/2;
        _elm.ym[_id] = (y_f + y_i)/2;
        _l = std::sqrt(std::pow((x_f - x_i),2) + std::pow((y_f - y_i),2));
        _elm.L[_id] = _l;
        
        // Normal to bottom direction
        _elm.xn[_id] = (y_f - y_i)/_l;
        _elm.yn[_id] = -(x_f - x_i)/_l;
        
        // Assign the traction stress
        _elm.Tx[_id] = -(_press) * _elm.xn[_id];
        _elm.Ty[_id] = -(_press) * _elm.yn[_id];
        
        _id++;

        // Update the node for next calculation
        r_i = r_f;
        x_i = x_f;
        y_i = y_f;
    }

    // Additional for inner element
    if (base == false){
        // ==============================
        // +----- ELEMENT ROTATION -----+
        // ==============================

        // Perform rotation here
        double _rot_angle = Par::Gin_Rot[ID] * M_PI / 180;
        if (_rot_angle != 0.0){
            // Temporary element
            element _rot_elm;

            // Resize each boundary element variable
            _rot_elm.num = _elm.num;
            _rot_elm.L = _elm.L;
            _rot_elm.xm.resize(num,0.0e0);
            _rot_elm.ym.resize(num,0.0e0);
            _rot_elm.xn.resize(num,0.0e0);
            _rot_elm.yn.resize(num,0.0e0);
            _rot_elm.Tx.resize(num,0.0e0);
            _rot_elm.Ty.resize(num,0.0e0);

            for (int i = 0; i < _elm.num; i++){
                _rot_elm.xm[i] = _elm.xm[i] * std::cos(_rot_angle) - _elm.ym[i] * std::sin(_rot_angle);
                _rot_elm.ym[i] = _elm.ym[i] * std::cos(_rot_angle) + _elm.xm[i] * std::sin(_rot_angle);
                _rot_elm.xn[i] = _elm.xn[i] * std::cos(_rot_angle) - _elm.yn[i] * std::sin(_rot_angle);
                _rot_elm.yn[i] = _elm.yn[i] * std::cos(_rot_angle) + _elm.xn[i] * std::sin(_rot_angle);
                _rot_elm.Tx[i] = _elm.Tx[i] * std::cos(_rot_angle) - _elm.Ty[i] * std::sin(_rot_angle);
                _rot_elm.Ty[i] = _elm.Ty[i] * std::cos(_rot_angle) + _elm.Tx[i] * std::sin(_rot_angle);
            }
            
            // Update the element
            _elm = _rot_elm;
        }

        // =================================
        // +----- ELEMENT TRANSLATION -----+
        // =================================
        // Perform translation from the center
        for (int i = 0; i < _elm.num; i++){
            _elm.xm[i] += x_cen;
            _elm.ym[i] += y_cen;
        }    
    }
    
    // Update the current element
    // **************************
    elm = _elm;
}