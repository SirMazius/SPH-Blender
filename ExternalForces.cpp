#include "ExternalForces.h"

void ExternalForces::ComputeGravity(vector<Vec3> &l_externalForce, const vector<float> &l_density) { // @suppress("Member declaration not found")

	int count = FluidParams::nParticles;
	Vec3 gravity(0.0, 0.0, -9.8);

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_externalForce.at(i).SetZero();
		l_externalForce.at(i) += l_density[i] * gravity;
	}

	// cout << "ComputeGravity" << endl;
}

void ExternalForces::ComputeInwardNormal(const vector<float> &l_density, const vector<vector<int>> & l_neighbors,
		const vector<Vec3> &l_positions, vector<Vec3> &l_normals) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_normals.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			Vec3 vAux;
			Vec3::vDirector(l_positions[j], l_positions[i], vAux);
			l_normals[i] += (mass / l_density[j]) * Kernels::ValueGradient(vAux);
		}
	}
	//cout << "ComputeInwardNormal" << endl;
}

void ExternalForces::ComputeColorField(const vector<float> &l_density, const vector<vector<int>> & l_neighbors,
		const vector<Vec3> &l_positions, vector<float> &l_color) {
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_color.at(i) = 0;
		for (int j : l_neighbors[i]) {
			Vec3 vAux;
			Vec3::vDirector(l_positions[j], l_positions[i], vAux);
			l_color[i] += (mass / l_density[j]) * Kernels::ValueLaplacian(vAux);
		}
	}
	//cout << "ComputeColorField" << endl;
}

void ExternalForces::ComputeSurfaceTension(vector <float> & l_density, vector <float> & l_densityBorder, const vector<float> &l_color, const vector<Vec3> &l_normals,
		vector<Vec3> & l_externalForce) {
	float threshold = FluidParams::threshold;
	float sigma = FluidParams::surfaceTension;
	int count = FluidParams::nParticles;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 AuxNormal = l_normals[i];
		float AuxNormalValue = AuxNormal.mag();
		if (AuxNormalValue > threshold && abs(AuxNormalValue > 0)) {
			l_externalForce[i] += -sigma * l_color[i] * (AuxNormal / AuxNormalValue);
			l_densityBorder[i] = l_density[i];
		}
	}
	//cout << "ComputeSurfaceTension" << endl;
}
