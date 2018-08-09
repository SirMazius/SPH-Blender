#include "InternalForces.h"
#include "Vec3.h"
#include "Kernels.h"
void InternalForces::ComputeMassDensity(vector<float> & l_density, const vector<Vec3> & l_positions,
		const vector<vector<int>> & l_neighbors) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
	// unsigned int auxCounter = 0;
	for (int i = 0; i < count; i++) {
		l_density.at(i) = 0;
		for (int j : l_neighbors.at(i)) {
			Vec3 vAux;
			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vAux);
			l_density.at(i) += mass * Kernels::Value(vAux);
			// auxCounter++;
		}
	}
	cout << "ComputeMassDensity" << endl;
}

void InternalForces::ComputePressures(const vector<float> & l_density, vector<float> & l_pressures, float restDensity) {
	int count = FluidParams::nParticles;
	float stiffness = FluidParams::stiffness;
	for (int i = 0; i < count; i++) {
		l_pressures.at(i) = stiffness * (l_density.at(i) - restDensity);
//		l_pressures.at(i) = stiffness * l_density.at(i);
	}
	cout << "ComputePressures" << endl;
}

void InternalForces::ComputePressureForce(const vector<float> & l_density, const vector<Vec3> & l_positions,
		const vector<float> & l_pressures, vector<Vec3> & l_internalForce, const vector<vector<int>> l_neighbors) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;

	for (int i = 0; i < count; i++) {
		l_internalForce.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			if (i != j && l_density.at(j) != 0) {
				Vec3 vAux;
				Vec3::vDirector(l_positions[j], l_positions[i], vAux);
				l_internalForce[i] -= ( ((l_pressures[i] + l_pressures[j]) / (2 * l_density[j])) * mass
						* Kernels::SpikyGradient(vAux));
			}
		}

	}
	cout << "ComputePressureForce" << endl;
}

void InternalForces::ComputeViscosityForce(const vector<float> & l_density, const vector<Vec3> & l_velocity,
		vector<Vec3> & l_internalForce, const vector<vector<int>> & l_neighbors) {
	float mu = FluidParams::viscosity;
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
	for (int i = 0; i < count; i++) {
		Vec3 viscosityForce;
		for (int j : l_neighbors[i]) {
			if (i != j && l_density.at(j) != 0) {
				Vec3 vAux, vDir;
				// cout << "SIZEEZ "<< l_velocity.size() << endl;
				Vec3::vDirector(l_velocity.at(i), l_velocity.at(j), vAux);
				Vec3::vDirector(l_velocity.at(j), l_velocity.at(i), vDir);
				viscosityForce += (vAux * (mass / l_density.at(j)) * Kernels::ViscosityLaplacian(vDir));
			}
		}
		l_internalForce.at(i) += viscosityForce * mu;
	}
	cout << "ComputeViscosityForce" << endl;
}
