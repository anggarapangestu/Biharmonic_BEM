#include "BEM.hpp"
// ================================================================================
// ================================================================================
void calcBEM::Define_BEM(const element& elm, const std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM method initialization ...\n");

    // Initialize the size
    this->N = 0;
    N += elm.num;
    for (int ID = 0; ID < Par::N_Gin; ID++){
        N += in_elm[ID].num;
    }
    
    // Initialize the ID list
    this->elmID.resize(N,0);    // The element ID
    this->elmGIN.resize(N,0);   // The element geometry position (-1:= base, 0:= inner_1, ...)
    for (int i = 0; i < N; i++){
        int ID = i;
        int GIN = -1;
        if (ID - elm.num < 0){
            // The base geometry
            elmID[i] = ID;
            elmGIN[i] = GIN;
        }else{
            // The inner geometry
            GIN++;
            ID -= elm.num;

            for (int _id = 0; _id < Par::N_Gin; _id++){
                if (ID - in_elm[_id].num >= 0){
                    ID -= in_elm[_id].num;
                    GIN ++;
                }else{
                    break;
                }
            }

            elmID[i] = ID;
            elmGIN[i] = GIN;
        }
    }
    
    printf("<+> Initialization done\n");
}

// ================================================================================
// ================================================================================
// Calculate the other F boundary value
void calcBEM::solve_F(element& elm, std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM calculating F ...\n");
    clock_t _time = clock();

    // Initialize the matrix
    Eigen::MatrixXd A_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::MatrixXd B_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::VectorXd A_Vec = Eigen::VectorXd::Zero(this->N);       // Left : Value is calculated
    Eigen::VectorXd B_Vec = Eigen::VectorXd::Zero(this->N);       // Right: Value is given

    // ================= Fill the BEM matrix =================
    // *******************************************************
    double _a, _k, L, _Aij, _Bij;
    for (int i = 0; i < this->N; i++){
        double xi, yi;  // Current element evaluated
        double xj, yj;  // Iteration element
        double xn, yn;  // Normal vector

        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            xi = elm.xm[i];
            yi = elm.ym[i];
        }else{
            // Inner geometry
            xi = in_elm[elmGIN[i]].xm[elmID[i]];
            yi = in_elm[elmGIN[i]].ym[elmID[i]];
        }

        // Calculate all element
        for (int j = 0; j < this->N; j++){
            // Update the iterated element
            if (elmGIN[j] < 0){
                // Base geometry
                xj = elm.xm[j];
                yj = elm.ym[j];
                xn = elm.xn[j];
                yn = elm.yn[j];
                L  = elm.L[j];
            }else{
                // Inner geometry
                xj = in_elm[elmGIN[j]].xm[elmID[j]];
                yj = in_elm[elmGIN[j]].ym[elmID[j]];
                xj = in_elm[elmGIN[j]].xn[elmID[j]];
                yj = in_elm[elmGIN[j]].yn[elmID[j]];
                L  = in_elm[elmGIN[j]].L[elmID[j]];
            }

            // Calculating the matrix element
            if (Par::opt_BEM == 1){
                _a = this->calc_a(xi,yi,xj,yj,xn,yn);
                _k = this->calc_k(xi,yi,xj,yj,xn,yn);
                _Aij = this->calc_Aij(_a,_k,L);
                _Bij = this->calc_Bij(_a,_k,L);
            }else if (Par::opt_BEM == 2){
                _Aij = this->calc_dGdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Bij = this->calc_G_dL(xi,yi,xj,yj,xn,yn,L);
            }    
            
            // Matrix assignment
            if (i == j){
                A_Mat(i, j) = _Aij - 0.5;
            }else{
                A_Mat(i, j) = _Aij;
            }
            B_Mat(i, j) = _Bij;
        }
    }

    save.write_Matrix(A_Mat, "AF");
    save.write_Matrix(B_Mat, "BF");

    // ================= Fill the BEM vector =================
    // *******************************************************
    for (int i = 0; i < this->N; i++){
        double _dFdn, _F;   // The boundary value
        bool F_type;        // The neumann type flag
        
        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            _dFdn  = elm.dFdn[i];
            _F     = elm.F[i];
            F_type = elm.F_type[i];
        }else{
            // Inner geometry
            _dFdn  = in_elm[elmGIN[i]].dFdn[elmID[i]];
            _F     = in_elm[elmGIN[i]].F[elmID[i]];
            F_type = in_elm[elmGIN[i]].F_type[elmID[i]];
        }
        
        // Assign the boundary value
        if (F_type == true){
            // Neumann boundary condition
            B_Vec(i) = _dFdn;
        }
        if (F_type == false){
            // Dirichlet boundary condition
            calcBEM::swap_col(A_Mat, B_Mat, i);
            B_Vec(i) = _F;
        }
    }

    // Solve the matrix by using matrix invertion
    // ******************************************
    Eigen::VectorXd bi = B_Mat * B_Vec;
    A_Vec = A_Mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi);

    // ================= Update the boundary value on the element ================= 
    // ****************************************************************************
    for (int i = 0; i < this->N; i++){
        bool F_type;        // The neumann type flag
        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            F_type = elm.F_type[i];
            
            // Assign the boundary value
            if (F_type == true){
                // Neumann boundary condition
                elm.F[i] = A_Vec(i);
            }
            if (F_type == false){
                // Dirichlet boundary condition
                elm.dFdn[i] = A_Vec(i);
            }
        }else{
            // Inner geometry
            F_type = in_elm[elmGIN[i]].F_type[elmID[i]];

            // Assign the boundary value
            if (F_type == true){
                // Neumann boundary condition
                in_elm[elmGIN[i]].F[elmID[i]] = A_Vec(i);
            }
            if (F_type == false){
                // Dirichlet boundary condition
                in_elm[elmGIN[i]].dFdn[elmID[i]] = A_Vec(i);
            }
        }
    }
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Calculating F comp. time           [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// ================================================================================
// ================================================================================
// Calculate the other phi boundary value
void calcBEM::solve_phi(element& elm, std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM calculating phi ...\n");
    clock_t _time = clock();

    // Initialize the matrix
    Eigen::MatrixXd A_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::MatrixXd B_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::MatrixXd C_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::MatrixXd D_Mat = Eigen::MatrixXd::Zero(this->N, this->N);
    Eigen::VectorXd A_Vec = Eigen::VectorXd::Zero(this->N);       // Left : Value is calculated (typically p)
    Eigen::VectorXd B_Vec = Eigen::VectorXd::Zero(this->N);       // Right: Value is given      (typically dpdn)
    Eigen::VectorXd C_Vec = Eigen::VectorXd::Zero(this->N);       // Right: Value is given (F)
    Eigen::VectorXd D_Vec = Eigen::VectorXd::Zero(this->N);       // Right: Value is given (dFdn)

    // ================= Fill the BEM matrix =================
    // *******************************************************
    for (int i = 0; i < this->N; i++){
        double _Aij, _Bij, _Cij, _Dij;  // BEM matrix element
        double _a, _k, L;   // Matrix parameter
        double xi, yi;      // Current element evaluated
        double xj, yj;      // Iteration element
        double xn, yn;      // Normal vector

        // Update the current evaluated element + update C and D vector element
        if (elmGIN[i] < 0){
            // Base geometry
            xi = elm.xm[i];
            yi = elm.ym[i];
            C_Vec(i) = elm.F[i];
            D_Vec(i) = elm.dFdn[i];
        }else{
            // Inner geometry
            xi = in_elm[elmGIN[i]].xm[elmID[i]];
            yi = in_elm[elmGIN[i]].ym[elmID[i]];
            C_Vec(i) = in_elm[elmGIN[i]].F[elmID[i]];
            D_Vec(i) = in_elm[elmGIN[i]].dFdn[elmID[i]];
        }

        // Calculate all element
        for (int j = 0; j < this->N; j++){
            // Update the iterated element
            if (elmGIN[j] < 0){
                // Base geometry
                xj = elm.xm[j];
                yj = elm.ym[j];
                xn = elm.xn[j];
                yn = elm.yn[j];
                L  = elm.L[j];
            }else{
                // Inner geometry
                xj = in_elm[elmGIN[j]].xm[elmID[j]];
                yj = in_elm[elmGIN[j]].ym[elmID[j]];
                xj = in_elm[elmGIN[j]].xn[elmID[j]];
                yj = in_elm[elmGIN[j]].yn[elmID[j]];
                L  = in_elm[elmGIN[j]].L[elmID[j]];
            }

            // Calculating the matrix element
            if (Par::opt_BEM == 1){
                _a = this->calc_a(xi,yi,xj,yj,xn,yn);
                _k = this->calc_k(xi,yi,xj,yj,xn,yn);
                _Aij = this->calc_Aij(_a,_k,L);
                _Bij = this->calc_Bij(_a,_k,L);
                _Cij = this->calc_Cij(_a,_k,L);
                _Dij = this->calc_Dij(_a,_k,L);
            }else if (Par::opt_BEM == 2){
                _Aij = this->calc_dGdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Bij = this->calc_G_dL(xi,yi,xj,yj,xn,yn,L);
                _Cij = this->calc_dWdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Dij = this->calc_W_dL(xi,yi,xj,yj,xn,yn,L);
            }
            
            // Matrix
            if (i == j){
                A_Mat(i, j) = _Aij - 0.5;
            }else{
                A_Mat(i, j) = _Aij;
            }
            B_Mat(i, j) = _Bij;
            C_Mat(i, j) = _Cij;
            D_Mat(i, j) = _Dij;
        }
    }
    
    save.write_Matrix(A_Mat, "Aphi");
    save.write_Matrix(B_Mat, "Bphi");
    save.write_Matrix(C_Mat, "Cphi");
    save.write_Matrix(D_Mat, "Dphi");

    // ================= Fill the BEM vector =================
    // *******************************************************
    for (int i = 0; i < this->N; i++){
        double _dpdn, _p;   // The boundary value
        bool p_type;        // The neumann type flag
        
        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            _dpdn  = elm.dpdn[i];
            _p     = elm.p[i];
            p_type = elm.p_type[i];
        }else{
            // Inner geometry
            _dpdn  = in_elm[elmGIN[i]].dpdn[elmID[i]];
            _p     = in_elm[elmGIN[i]].p[elmID[i]];
            p_type = in_elm[elmGIN[i]].p_type[elmID[i]];
        }
        
        // Assign the boundary value
        if (p_type == true){
            // Neumann boundary condition
            B_Vec(i) = _dpdn;
        }
        if (p_type == false){
            // Dirichlet boundary condition
            calcBEM::swap_col(A_Mat, B_Mat, i);
            B_Vec(i) = _p;
        }
    }

    // Solve the matrix by using matrix invertion
    // ******************************************
    Eigen::VectorXd bi = (B_Mat * B_Vec) - (C_Mat * C_Vec) + (D_Mat * D_Vec);
    A_Vec = A_Mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi);

    // ================= Update the boundary value on the element ================= 
    // ****************************************************************************
    for (int i = 0; i < this->N; i++){
        bool p_type;        // The neumann type flag
        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            p_type = elm.p_type[i];
            
            // Assign the boundary value
            if (p_type == true){
                // Neumann boundary condition
                elm.p[i] = A_Vec(i);
            }
            if (p_type == false){
                // Dirichlet boundary condition
                elm.dpdn[i] = A_Vec(i);
            }
        }else{
            // Inner geometry
            p_type = in_elm[elmGIN[i]].p_type[elmID[i]];

            // Assign the boundary value
            if (p_type == true){
                // Neumann boundary condition
                in_elm[elmGIN[i]].p[elmID[i]] = A_Vec(i);
            }
            if (p_type == false){
                // Dirichlet boundary condition
                in_elm[elmGIN[i]].dpdn[elmID[i]] = A_Vec(i);
            }
        }
    }

    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Calculating phi comp. time         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// ================================================================================
// ================================================================================
// Calculate phi inside the domain region
void calcBEM::calculate_internal_phi(intElement& intElm, const element& elm, const std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM calculating internal node ...\n");
    clock_t _time = clock();

    // Initialize the group value
    double A_group, B_group, C_group, D_group;

    // Resize the internal node element phi variable
    intElm.phi.resize(intElm.num);

    // ================= Fill the BEM matrix =================
    // *******************************************************
    for (int i = 0; i < intElm.num; i++){
        // Internal variable
        double _Aij, _Bij, _Cij, _Dij;  // BEM element parameter
        double _p, _dpdn, _F, _dFdn;  // BEM element parameter
        double _a, _k, L;   // Matrix parameter
        double xi, yi;      // Current internal node evaluated
        double xj, yj;      // Iteration element
        double xn, yn;      // Normal vector

        // Update the current evaluated element + update C and D vector element
        xi = intElm.x[i];
        yi = intElm.y[i];
        A_group = 0;
        B_group = 0;
        C_group = 0;
        D_group = 0;

        // Calculate all boundary element
        for (int j = 0; j < this->N; j++){
            // Update the iterated element
            if (elmGIN[j] < 0){
                // Base geometry
                xj = elm.xm[j];
                yj = elm.ym[j];
                xn = elm.xn[j];
                yn = elm.yn[j];
                L  = elm.L[j];
                _p = elm.p[j];
                _F = elm.F[j];
                _dpdn = elm.dpdn[j];
                _dFdn = elm.dFdn[j];
            }else{
                // Inner geometry
                xj = in_elm[elmGIN[j]].xm[elmID[j]];
                yj = in_elm[elmGIN[j]].ym[elmID[j]];
                xj = in_elm[elmGIN[j]].xn[elmID[j]];
                yj = in_elm[elmGIN[j]].yn[elmID[j]];
                L  = in_elm[elmGIN[j]].L[elmID[j]];
                _p = in_elm[elmGIN[j]].p[elmID[j]];
                _F = in_elm[elmGIN[j]].F[elmID[j]];
                _dpdn = in_elm[elmGIN[j]].dpdn[elmID[j]];
                _dFdn = in_elm[elmGIN[j]].dFdn[elmID[j]];
            }

            // Calculating the matrix element
            if (Par::opt_BEM == 1){
                _a = this->calc_a(xi,yi,xj,yj,xn,yn);
                _k = this->calc_k(xi,yi,xj,yj,xn,yn);
                _Aij = this->calc_Aij(_a,_k,L);
                _Bij = this->calc_Bij(_a,_k,L);
                _Cij = this->calc_Cij(_a,_k,L);
                _Dij = this->calc_Dij(_a,_k,L);
            }else if (Par::opt_BEM == 2){
                _Aij = this->calc_dGdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Bij = this->calc_G_dL(xi,yi,xj,yj,xn,yn,L);
                _Cij = this->calc_dWdn_dL(xi,yi,xj,yj,xn,yn,L);
                _Dij = this->calc_W_dL(xi,yi,xj,yj,xn,yn,L);
            }
            
            // Matrix
            A_group += _Aij * _p;
            B_group += _Bij * _dpdn;
            C_group += _Cij * _F;
            D_group += _Dij * _dFdn;
        }

        // Update the phi value of the internal node
        intElm.phi[i] = A_group - B_group + C_group - D_group;
    }
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Internal node phi calculation\n");
    printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}
