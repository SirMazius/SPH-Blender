#pragma once

#include "Vec3.h"
#include "Kernels.h"
#include <vector>
#include "FluidParams.h"

using namespace std;

class ExternalForces {
public:
	static void ComputeGravity(vector<Vec3> & l_externalForce,
			const vector<float> & l_density);
    static void ComputeInwardNormal(const vector<float> &l_density, const vector<vector<int>> & l_neighbors, const vector<Vec3> &l_positions, vector<Vec3> &l_normals);
    static void ComputeColorField(const vector<float> &l_density, const vector<vector<int>> & l_neighbors, const vector<Vec3> &l_positions, vector<float> &l_color);
    //static void ComputeSurfaceTension(const vector<float> &l_color, const vector<Vec3> &l_normals, vector<Vec3> & l_externalForce);
    static void ComputeSurfaceTension(vector <float> & l_density,vector <float> & l_densityBorder, const vector<float> &l_color, const vector<Vec3> &l_normals,
    		vector<Vec3> & l_externalForce);
};
