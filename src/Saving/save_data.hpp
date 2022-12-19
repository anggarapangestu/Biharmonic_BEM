#ifndef SAVE_DATA
#define SAVE_DATA
#include "../../variable.hpp"
#include <iomanip>
#include <fstream>
#include <string>

#include "../Eigen/Dense"

class dataSaving
{
private:
    // The class private method belong here
    std::ofstream save;
public:
    // The class public method belong here
    void simulation_log();                                 // Displaying and saving the simulation log
    void write_internal_data(const intElement& intElm);    // Saving the internal node data
    void write_internal_data_temp(const intElement& intElm);    // Saving the internal node data temperature
    void write_BEM_data(const element& elm, const std::vector<element>& in_elm);        // Saving the boundary panel data
    void write_BEM_data_temp(const element& elm, const std::vector<element>& in_elm);   // Saving the boundary panel data
    void write_Matrix(const Eigen::MatrixXd& MAT, std::string name);  // Saving the boundary panel data
    void write_Matrix(const Eigen::VectorXd& VEC, std::string name);  // Saving the boundary panel data
    void save_Test(const intElement & elm);
};

#endif