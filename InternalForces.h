#pragma once

#include <vector>
#include "Vec3.h"
#include "Kernels.h"
#include "FluidParams.h"

using namespace std;

class InternalForces {
public:
	static void ComputeMassDensity(vector<float> & l_density,
			const vector<Vec3> & l_positions,
			const vector<vector<int>> & l_neighbors);
	static void ComputePressures(const vector<float> & l_density,
			vector<float> & l_pressures, float restDensity);

	static void ComputePressureForce(const vector<float> & l_density,
			const vector<Vec3> & l_positions, const vector<float> & l_pressures,
			vector<Vec3> & l_pressureForce,
			const vector<vector<int>> l_neighbors);
	static void ComputeViscosityForce(const vector<float> & l_density,
			const vector<Vec3> & l_velocity, vector<Vec3> & l_internalForce, const vector<vector<int>> & l_neighbors);

	static void Salut(vector<Vec3> & v) {
		v.push_back(Vec3(0, 0, 0));
	}
};
