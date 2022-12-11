#ifndef INCLUDED_NEIGHBOR
#define INCLUDED_NEIGHBOR

#ifndef INCLUDE_VARIABLE
#include "../../variable.hpp"
#endif

#include <algorithm>

class neighbor
{
private:
	// The instances inside neighbor
	void create_grid(double xfixed, double hG, double xGmin, double xGmax, 
      	int nG, int &nG_new, std::vector<double> &xGi);
	void find_cell(double xPi, double yPi, int nGx, int nGy, const std::vector<double> &xGi, 
		const std::vector<double> &yGi, double hG, int &i_cell, int &j_cell);
	void find_neighborcells(int i_cell, int j_cell, double neighbor_scale, int nGx, int nGy,
		int &lx_nbc, int &rx_nbc, int &ly_nbc, int &ry_nbc);

public:
	// TODO: generate neighborhood using link-list algorithm
	std::vector<std::vector<int>> link_list(int np, const std::vector<double> &xp, const std::vector<double> &yp, double ngh_scale);
};

#endif
