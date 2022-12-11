#ifndef SAVE_DATA
#define SAVE_DATA
#include "../../variable.hpp"
#include <iomanip>
#include <fstream>
#include <string>

class dataSaving
{
private:
    // The class private method belong here
    std::ofstream save;
public:
    // The class public method belong here
    void simulation_log();                                 // Displaying and saving the simulation log
    void write_internal_data(const intElement& intElm);    // Saving the internal node data
    void write_BEM_data(const element& elm, const std::vector<element>& in_elm);  // Saving the boundary panel data
};

#endif