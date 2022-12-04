#ifndef INITIALIZATION
#define INITIALIZATION
#include "../../variable.hpp"

class initialization
{
private:
    // The class private method belong here
    void internal_regular(intElement& intElm, const element& elm);
    void internal_finer_near_panel(intElement& intElm, const element& elm);
public:
    // The class public method belong here
    void generate_internal_node(intElement& intElm, const element& elm);
    void generate_boundary_element(element& elm);
};

#endif

/*

class className
{
private:
    // The class private method belong here
public:
    // The class public method belong here
};

*/

