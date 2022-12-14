#include "LSMPSa.hpp"

#pragma region public methods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSa::set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                       const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = x.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(x, y, s, f, neighborlist);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddx()
{
    return this->_ddx;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddy()
{
    return this->_ddy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2x()
{
    return this->_d2d2x;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2dxdy()
{
    return this->_d2dxdy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2y()
{
    return this->_d2d2y;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
#pragma endregion

#pragma region private methods
// ===========================================================================
// ===========================================================================
void LSMPSa::calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                             const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    int nparticle = x.size();
    
    // Evaluate the LSMPS A of all particle
    for (size_t i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::VectorXd bi = Eigen::VectorXd::Zero(MAT_SIZE);
        
        // Initialize H_rs Matrix
        Hi(0, 0) = std::pow(s[i], -1);
        Hi(1, 1) = std::pow(s[i], -1);
        Hi(2, 2) = std::pow(s[i], -2) * 2;
        Hi(3, 3) = std::pow(s[i], -2);
        Hi(4, 4) = std::pow(s[i], -2) * 2;

        // Evaluate at each neighbor particle
        std::vector<int> _neighborIndex = neighborlist[i];
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            // Note that neighbor particle and evaluated particle inside the same variable
            int idxi = i;                   // Particle ID
            int idxj = _neighborIndex[j];   // Neighbor Particle ID
            
            // Particle position
            double _xi = x[idxi];
            double _yi = y[idxi];
            
            // Neighbor particle position
            double _xj = x[idxj];
            double _yj = y[idxj];

            // Coordinate distance between Particle and Neighbor Particle
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            
            // LSMPS A interpolated variable value different
            double _fij = f[idxj] - f[idxi];    // Neighbor particle - evaluated particle

            // Resultant distance(_rij) and effective radius(_Rij)
            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2)); //distance between particles. 
            double _Ri = this->R_fac * s[idxi];   // Effective radius of Current particle

            std::vector<double> _p1 = this->get_p(_xij, _yij, s[idxi]);
            std::vector<double> _p2 = this->get_p(_xij, _yij, s[idxi]);
            double _wij = this->weight_function(_rij, _Ri) * std::pow(s[idxj]/s[idxi],2);

            // Calculatin moment matrix M and b
            for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
            {
                for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                {
                    // generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                }
                // generate moment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _fij);
            }
        }

        // Solve Least Square
        Eigen::VectorXd MiInv_Bi = Mi.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi);
        Eigen::VectorXd Dx = Hi * MiInv_Bi;

        // Assign to private variables
        this->_ddx[i] = Dx[0];
        this->_ddy[i] = Dx[1];
        this->_d2d2x[i] = Dx[2];
        this->_d2dxdy[i] = Dx[3];
        this->_d2d2y[i] = Dx[4];
    }
}

// ===========================================================================
// ===========================================================================

// Calculation on p LSMPS
std::vector<double> LSMPSa::get_p(const double &xij, const double &yij, const double &si)
{
    std::vector<double> _p(MAT_SIZE);

    double _xij = xij / si;
    double _yij = yij / si;

    _p[0] = _xij;
    _p[1] = _yij;
    _p[2] = _xij * _xij;
    _p[3] = _xij * _yij;
    _p[4] = _yij * _yij;

    return _p;
}

// Weight function definition
double LSMPSa::weight_function(const double &rij, const double &Rij)
{
    double _wij;
    if (rij <= Rij)
    {
        _wij = std::pow(1 - (rij / Rij), 2);
    }
    else
    {
        _wij = 0.0e0;
    }

    return _wij;
}

#pragma endregion
