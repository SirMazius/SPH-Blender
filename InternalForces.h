#pragma once

#include <vector>
#include "Vec3.h"
#include "Kernels.h"
#include "FluidParams.h"
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class InternalForces {
public:
	static void ComputeMassDensity(vector<float> & l_density, const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors);

	static void ComputeDensityDelta(vector<float> & l_auxDensity, const vector<Vec3> & l_auxPositions, const vector<vector<int>> & l_neighbors);

	static void ComputePressureCorrection(vector<float> & l_auxDensity, vector<float> & l_pressures, float & pError);

	static void ComputePressures(const vector<float> & l_density, vector<float> & l_pressures, float restDensity);

	static void ComputePressureForce(const vector<float> & l_density, const vector<Vec3> & l_positions, const vector<float> & l_pressures,
			vector<Vec3> & l_pressureForce, const vector<vector<int>> & l_neighbors);

	static void ComputePCI_SPH_PressureForce(const vector<float> & l_density, const vector<Vec3> & l_positions, const vector<float> & l_pressures,
			vector<Vec3> & l_pressureForce, const vector<vector<int>> & l_neighbors);

	static void ComputeViscosityForce(const vector<float> & l_density, vector<Vec3> & l_velocity, vector<Vec3> & l_internalForce,
			const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions);

	static void GetMatrix(MatrixXd & input, MatrixXd & G, Vector3d & center, double & det, Vector3d & currentPoint);

	static void ComputeAnisotropy(const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions, vector<Vector3d> & l_centers, vector<MatrixXd> & l_Gs,
			vector<double> & l_det);
	static void ComputeAnisotropyMassDensity(vector<float> & l_density, const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors,
			vector<Vector3d> l_centers, vector<MatrixXd> l_Gs, vector<double> l_det);


	static void Salut(vector<Vec3> & v) {
		v.push_back(Vec3(0, 0, 0));
	}
};
