#include "BEM.hpp"

// ================================================================================
// ============================== BEM INITIALIZATION ==============================
// ================================================================================
// Initialize the BEM internal variable
void calcBEM::Define_BEM(const element& elm, const std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM method initialization ...\n");
    printf("<+> Define the BEM Internal Data\n");

    // Initialize the matrix size
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

    // Assign the element geometrical data
    this->define_element(elm, in_elm);
    
    printf("<+> BEM Initialization done\n");
    printf("<+> Total number of boundary element  : %8d\n", N);
}

// Define the element for BEM calculation
void calcBEM::define_element(const element& elm, const std::vector<element>& in_elm){
    // Resize the BEM element data
    printf("<+> Define the BEM element\n");

    // Geometrical Data
    this->BEM_elm.xm.resize(this->N,0.0e0);
    this->BEM_elm.ym.resize(this->N,0.0e0);
    this->BEM_elm.xn.resize(this->N,0.0e0);
    this->BEM_elm.yn.resize(this->N,0.0e0);
    this->BEM_elm.L.resize(this->N,0.0e0);
    
    // Dirichlet and Neumann Data
    this->BEM_elm.F.resize(this->N,0.0e0);
    this->BEM_elm.p.resize(this->N,0.0e0);
    this->BEM_elm.dFdn.resize(this->N,0.0e0);
    this->BEM_elm.dpdn.resize(this->N,0.0e0);
    this->BEM_elm.F_type.resize(this->N,true);
    this->BEM_elm.p_type.resize(this->N,true);

    // Assign the value of the element
    for (int i = 0; i < this->N; i++){
        // Update the current evaluated element
        if (elmGIN[i] < 0){
            // Base geometry
            this->BEM_elm.xm[i] = elm.xm[i];
            this->BEM_elm.ym[i] = elm.ym[i];
            this->BEM_elm.xn[i] = elm.xn[i];
            this->BEM_elm.yn[i] = elm.yn[i];
            this->BEM_elm.L[i]  = elm.L[i];
            
            // Solid and Biharmonic simulation
            if (Par::opt_sim_type == 1 || Par::opt_sim_type == 3){
                this->BEM_elm.F[i]  = elm.F[i];
                this->BEM_elm.p[i]  = elm.p[i];
                this->BEM_elm.dFdn[i]    = elm.dFdn[i];
                this->BEM_elm.dpdn[i]    = elm.dpdn[i];
                this->BEM_elm.F_type[i]  = elm.F_type[i];
                this->BEM_elm.p_type[i]  = elm.p_type[i];
            }
            // Temperature simulation
            else if (Par::opt_sim_type == 2){
                this->BEM_elm.F[i]      = elm.T[i];
                this->BEM_elm.dFdn[i]   = elm.dTdn[i];
                this->BEM_elm.F_type[i] = elm.T_type[i];
            }
            // Laplace simulation
            else if (Par::opt_sim_type == 4){
                this->BEM_elm.F[i]      = elm.p[i];
                this->BEM_elm.dFdn[i]   = elm.dpdn[i];
                this->BEM_elm.F_type[i] = elm.p_type[i];
            }

        }else{
            // Inner geometry
            this->BEM_elm.xm[i] = in_elm[elmGIN[i]].xm[elmID[i]];
            this->BEM_elm.ym[i] = in_elm[elmGIN[i]].ym[elmID[i]];
            this->BEM_elm.xn[i] = in_elm[elmGIN[i]].xn[elmID[i]];
            this->BEM_elm.yn[i] = in_elm[elmGIN[i]].yn[elmID[i]];
            this->BEM_elm.L[i]  = in_elm[elmGIN[i]].L[elmID[i]];
            
            // Solid and Biharmonic simulation
            if (Par::opt_sim_type == 1 || Par::opt_sim_type == 3){
                this->BEM_elm.F[i]  = in_elm[elmGIN[i]].F[elmID[i]];
                this->BEM_elm.p[i]  = in_elm[elmGIN[i]].p[elmID[i]];
                this->BEM_elm.dFdn[i]    = in_elm[elmGIN[i]].dFdn[elmID[i]];
                this->BEM_elm.dpdn[i]    = in_elm[elmGIN[i]].dpdn[elmID[i]];
                this->BEM_elm.F_type[i]  = in_elm[elmGIN[i]].F_type[elmID[i]];
                this->BEM_elm.p_type[i]  = in_elm[elmGIN[i]].p_type[elmID[i]];
            }
            // Temperature simulation
            else if (Par::opt_sim_type == 2){
                this->BEM_elm.F[i]  = in_elm[elmGIN[i]].T[elmID[i]];
                this->BEM_elm.dFdn[i]    = in_elm[elmGIN[i]].dTdn[elmID[i]];
                this->BEM_elm.F_type[i]  = in_elm[elmGIN[i]].T_type[elmID[i]];
            }
            // Laplace simulation
            else if (Par::opt_sim_type == 4){
                this->BEM_elm.F[i]      = in_elm[elmGIN[i]].p[elmID[i]];
                this->BEM_elm.dFdn[i]   = in_elm[elmGIN[i]].dpdn[elmID[i]];
                this->BEM_elm.F_type[i] = in_elm[elmGIN[i]].p_type[elmID[i]];
            }
        }
    }
}

// Take the final data of each element
void calcBEM::retrieve_element(element& elm, std::vector<element>& in_elm){
    // Assign the value of the element
    for (int i = 0; i < this->N; i++){
        // Update the current evaluated element
        
        // Base geometry
        if (elmGIN[i] < 0){
            // Solid and Biharmonic simulation
            if (Par::opt_sim_type == 1 || Par::opt_sim_type == 3){
                elm.F[i] = this->BEM_elm.F[i];
                elm.p[i] = this->BEM_elm.p[i];
                elm.dFdn[i] = this->BEM_elm.dFdn[i];
                elm.dpdn[i] = this->BEM_elm.dpdn[i];
                elm.F_type[i] = this->BEM_elm.F_type[i];
                elm.p_type[i] = this->BEM_elm.p_type[i];
            }
            // Temperature simulation
            else if (Par::opt_sim_type == 2){
                elm.T[i] = this->BEM_elm.F[i];
                elm.dTdn[i] = this->BEM_elm.dFdn[i];
                elm.T_type[i] = this->BEM_elm.F_type[i];
            }
            // Laplace simulation
            else if (Par::opt_sim_type == 4){
                elm.p[i] = this->BEM_elm.F[i];
                elm.dpdn[i] = this->BEM_elm.dFdn[i];
                elm.p_type[i] = this->BEM_elm.F_type[i];
            }
        }
        
        // Inner geometry
        else{
            // Solid and Biharmonic simulation
            if (Par::opt_sim_type == 1 || Par::opt_sim_type == 3){
                in_elm[elmGIN[i]].F[elmID[i]] = this->BEM_elm.F[i];
                in_elm[elmGIN[i]].p[elmID[i]] = this->BEM_elm.p[i];
                in_elm[elmGIN[i]].dFdn[elmID[i]] = this->BEM_elm.dFdn[i];
                in_elm[elmGIN[i]].dpdn[elmID[i]] = this->BEM_elm.dpdn[i];
                in_elm[elmGIN[i]].F_type[elmID[i]] = this->BEM_elm.F_type[i];
                in_elm[elmGIN[i]].p_type[elmID[i]] = this->BEM_elm.p_type[i];
            }
            // Temperature simulation
            else if (Par::opt_sim_type == 2){
                in_elm[elmGIN[i]].T[elmID[i]] = this->BEM_elm.F[i];
                in_elm[elmGIN[i]].dTdn[elmID[i]] = this->BEM_elm.dFdn[i];
                in_elm[elmGIN[i]].T_type[elmID[i]] = this->BEM_elm.F_type[i];
            }
            // Laplace simulation
            else if (Par::opt_sim_type == 4){
                in_elm[elmGIN[i]].p[elmID[i]] = this->BEM_elm.F[i];
                in_elm[elmGIN[i]].dpdn[elmID[i]] = this->BEM_elm.dFdn[i];
                in_elm[elmGIN[i]].p_type[elmID[i]] = this->BEM_elm.F_type[i];
            }
        }
    }
}

// ================================================================================
// ============================ BEM CALCULATION MANAGER ===========================
// ================================================================================

// The main BEM calculation
void calcBEM::calc_BEM(intElement& intElm, element& elm, std::vector<element>& in_elm){
    // initialization the BEM calculation log
    printf("\nBEM calculation initialization ...\n");

    // Biharmonic simulation solver
    if (Par::opt_sim_type == 1){
        // BEM type log
        printf("<+> Biharmonic solid simulation\n");

        // // Calculate the other boundary element value by BEM
        // this->solve_F(elm, in_elm);
        // this->solve_phi(elm, in_elm);

        // // Calculate the phi value at the internal domain node
        // this->calculate_internal_phi(intElm, elm, in_elm);

        // Calculate the other boundary element value by BEM
        this->solve_laplace();
        this->solve_biharmonic();
        
        // Get the other boundary condition value
        this->retrieve_element(elm, in_elm);
        
        // Calculate phi at the internal domain node
        this->bhm_internal_calc(intElm);
    }
    // Heat transfer simulation solver
    else if (Par::opt_sim_type == 2){
        // BEM type log
        printf("<+> Steady heat transfer simulation\n");

        // // Calculate the other boundary element value by BEM
        // this->solve_T(elm, in_elm);
        
        // // Calculate the T value at the internal domain node
        // this->calculate_internal_T(intElm, elm, in_elm);

        // Calculate the other boundary element value by BEM
        this->solve_laplace();
        
        // Get the other boundary condition value
        this->retrieve_element(elm, in_elm);
        
        // Calculate phi at the internal domain node
        this->lap_internal_calc(intElm);
        intElm.T = intElm.phi;
    }
    // Biharmonic solver
    else if (Par::opt_sim_type == 3){
        // BEM type log
        printf("<+> Biharmonic function solver\n");
        
        // Calculate the other boundary element value by BEM
        this->solve_laplace();
        this->solve_biharmonic();
        
        // Get the other boundary condition value
        this->retrieve_element(elm, in_elm);
        
        // Calculate phi at the internal domain node
        this->bhm_internal_calc(intElm);
    }
    // Laplace solver
    else if (Par::opt_sim_type == 4){
        // BEM type log
        printf("<+> Laplace function solver\n");

        // Calculate the other boundary element value by BEM
        this->solve_laplace();
        
        // Get the other boundary condition value
        this->retrieve_element(elm, in_elm);
        
        // Calculate phi at the internal domain node
        this->lap_internal_calc(intElm);
    }
}

// ================================================================================
// =============================== BEM TESTING LIST ===============================
// ================================================================================

// The code below for TESTING !!!

// Calculate the other F boundary value
void calcBEM::TEST_BEM(element& elm, std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nBEM testing on calc. F ...\n");
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
                xn = in_elm[elmGIN[j]].xn[elmID[j]];
                yn = in_elm[elmGIN[j]].yn[elmID[j]];
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
                std::cout << "> A (" << i << ", " << j << ") : " << _Aij << std::endl;
            }else{
                A_Mat(i, j) = _Aij;
            }
            B_Mat(i, j) = _Bij;
        }
    }

    save.write_Matrix(A_Mat, "A_TEST");
    save.write_Matrix(B_Mat, "B_TEST");

    // Displaying the computational time
    _time = clock() - _time;
    printf("<-> Finish testing F comp. time        [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// Evaluate the theta value at each position near the evaluated panel
void calcBEM::CALC_theta(){
    // Create a dummy internal node data with regular distribution (from 0 < x,y < 1)
    intElement _test;
    double N = 1000;
    for (int i = 0; i < (int)N; i++){
        for (int j = 0; j < (int)N; j++){
            _test.x.push_back(j/N);
            _test.y.push_back(i/N);
        }
    }

    // Evaluate the theta for each internal node related to point M below
    double xm = 0.5;
    double ym = 0.5;
    for (int i = 0; i < _test.x.size(); i++){
        double a = this->calc_a(_test.x[i],_test.y[i],xm,ym,0,1);
        double k = this->calc_k(_test.x[i],_test.y[i],xm,ym,0,1);
        _test.phi.push_back(this->calc_Aij(a,k,0.2));
    }

    // Write the data
    save.save_Test(_test);
}

