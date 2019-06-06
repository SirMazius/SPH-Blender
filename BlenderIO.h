#pragma once

#include <vector>
#include "Vec3.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "FluidParams.h"
#include "Kernels.h"
#include <numeric>

using namespace std;
using namespace Eigen;

class BlenderIO {
public:
	static bool ReadParams(string fileName, Vec3 l_bounds[2]);
	static bool ReadPOSVEL(string fileName, vector<Vec3> & l_pos, vector<Vec3> & l_velocity);
//	static void WritePOSVEL(string name, float t, int iteration, vector<Vec3> & l_pos, vector<Vec3> & l_velocity);
	static void WritePOSVEL(string name, float t, int iteration, vector<Vec3> l_pos, vector<Vec3> l_velocity);
	static void WriteExcelData(string name, float deepDensity, float surfaceDensity, float executionTime);
	static void WriteHeightDensityData(vector<Vec3> l_positions, vector<float> l_density, int iteration);
	static void WriteDensityPerParticle(vector<Vec3> l_positions, vector<float> l_density, int iteration);
	static void WriteEigenVectorsNEigenValues(const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions, int iteration);
	static void WriteForces(vector<Vec3> & l_internalForce, vector<Vec3> & l_pressureForce,
			vector<Vec3> & l_externalForce, int iteration);

	static unordered_map<string, string> parametersMap;
	static vector<string> l_parameters;
};
