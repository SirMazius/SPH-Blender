#include "InternalForces.h"
#include "Vec3.h"
#include "Kernels.h"
void InternalForces::ComputeMassDensity(vector<float> & l_density, const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_density.at(i) = 0;
		for (int j : l_neighbors.at(i)) {
			Vec3 vAux;
			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vAux);
//			l_density.at(i) += mass * Kernels::Value(vAux);
			l_density.at(i) += Kernels::Value(vAux);
		}
		l_density.at(i) *= mass;
	}

	//cout << "ComputeMassDensity" << endl;
}

void InternalForces::ComputeDensityDelta(vector<float> & l_auxDensity, const vector<Vec3> & l_auxPositions, const vector<vector<int>> & l_neighbors) {
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
	//float dt2 = FluidParams::dt2;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_auxDensity.at(i) = 0;
		for (int j : l_neighbors.at(i)) {
			Vec3 vDir;
			Vec3::vDirector(l_auxPositions.at(j), l_auxPositions.at(i), vDir);
			l_auxDensity.at(i) += mass * Kernels::Value(vDir);
		}
	}
}

void InternalForces::ComputePressureCorrection(vector<float> & l_auxDensity, vector<float> & l_pressures, float & pError) {

	Vec3 sum1, sum2;
	const Vec3 center;
	float mass2 = FluidParams::mass * FluidParams::mass;
	float sum3 = 0, h = FluidParams::kernelRadius, pH = FluidParams::particleRadius, beta = FluidParams::beta, restDensity = FluidParams::restDensity;
	int count = FluidParams::nParticles;
	int count2 = 0;
	float multi = 2.0;

	for (float x = -h - pH; x <= h + pH; x += multi * pH) {
		for (float y = -h - pH; y <= h + pH; y += multi * pH) {
			for (float z = -h - pH; z <= h + pH; z += multi * pH) {

				Vec3 r(x, y, z);
				Vec3 vDir;
				Vec3::vDirector(r, center, vDir);
				Vec3 grad = Kernels::SpikyGradient(r);
				sum1 -= grad;
				sum2 += grad;
				sum3 += Vec3::Dot(grad, grad);
				count2++;
			}
		}
	}

	vector<float> l_error(count);
	vector<float> l_errorAux(count);
	float factor = -1 / (beta * (Vec3::Dot(sum1, sum2) - sum3)); // (restDensity * restDensity) / (2 * FluidParams::dt2 * mass2 * (sum3));
//	cout << "FACTOR "<< factor << endl;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		float error = max(0.0f, l_auxDensity.at(i) - restDensity);
//		float errorpercent = (error * 100) / restDensity;
//		if (errorpercent < 50)
		{
			l_error.at(i) = error;
			l_pressures.at(i) += factor * (error);
		}
		l_errorAux.at(i) = error;
	}

	float bigger = 0;

	for (auto e : l_error) {
		float absE = abs(e);
//		cout << absE << endl;
		if (absE > bigger)
			bigger = absE;
	}

	pError = (bigger * 100) / restDensity;
//	cout << "MAX ERROR ->> " << (*max_element(l_error.begin(), l_error.end()) * 100) / restDensity<< endl;
//	cout << pError << endl;
}

void InternalForces::ComputePressures(const vector<float> & l_density, vector<float> & l_pressures, float restDensity) {
	int count = FluidParams::nParticles;
	float stiffness = FluidParams::stiffness;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_pressures.at(i) = stiffness * (l_density.at(i) - restDensity);
		//		l_pressures.at(i) = stiffness * l_density.at(i);
	}

//cout << "ComputePressures" << endl;
}

void InternalForces::ComputePressureForce(const vector<float> & l_density, const vector<Vec3> & l_positions, const vector<float> & l_pressures,
		vector<Vec3> & l_pressureForce, const vector<vector<int>> & l_neighbors) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 pressureAux;
		// l_pressureForce.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			if (i != j && l_density.at(j) != 0) {
				Vec3 vAux;
				Vec3::vDirector(l_positions[j], l_positions[i], vAux);
//				l_internalForce[i] -= ( ((l_pressures[i] + l_pressures[j]) / (2 * l_density[j])) * mass
//						* Kernels::SpikyGradient(vAux));

//				l_internalForce.at(i) = l_internalForce.at(i)
//						- (l_pressures.at(i) + l_pressures.at(j) / 2) * (mass / l_density.at(j))
//								* Kernels::SpikyGradient(vAux);

				pressureAux = pressureAux
						+ (l_pressures.at(i) / (l_density.at(i) * l_density.at(i)) + l_pressures.at(j) / (l_density.at(j) * l_density.at(j))) * mass
								* Kernels::SpikyGradient(vAux);

			}
		}
		l_pressureForce.at(i) = pressureAux * -l_density.at(i);
	}
//cout << "ComputePressureForce" << endl;
}

void InternalForces::ComputePCI_SPH_PressureForce(const vector<float> & l_density, const vector<Vec3> & l_positions, const vector<float> & l_pressures,
		vector<Vec3> & l_pressureForce, const vector<vector<int>> & l_neighbors) {

	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 pressureAux;
		// l_pressureForce.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			if (i != j && l_density.at(j) != 0) {
				Vec3 vAux;
				Vec3::vDirector(l_positions[j], l_positions[i], vAux);
//				l_internalForce[i] -= ( ((l_pressures[i] + l_pressures[j]) / (2 * l_density[j])) * mass
//						* Kernels::SpikyGradient(vAux));

//				l_internalForce.at(i) = l_internalForce.at(i)
//						- (l_pressures.at(i) + l_pressures.at(j) / 2) * (mass / l_density.at(j))
//								* Kernels::SpikyGradient(vAux);

				pressureAux = pressureAux
						+ (l_pressures.at(i) / (l_density.at(i) * l_density.at(i)) + l_pressures.at(j) / (l_density.at(j) * l_density.at(j))) * mass
								* Kernels::SpikyGradient(vAux);

			}
		}
		l_pressureForce.at(i) = pressureAux * -l_density.at(i);
	}
//cout << "ComputePressureForce" << endl;
}

void InternalForces::ComputeViscosityForce(const vector<float> & l_density, vector<Vec3> & l_velocity, vector<Vec3> & l_internalForce,
		const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions) {
	float mu = FluidParams::viscosity;
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 viscosityForce;
		// l_internalForce.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			if (i != j && l_density.at(j) != 0) {
				Vec3 vAux, vDir;
				Vec3::vDirector(l_velocity.at(i), l_velocity.at(j), vAux);
				Vec3::vDirector(l_positions.at(j), l_positions.at(i), vDir);
				viscosityForce += (vAux * (mass / l_density.at(j)) * Kernels::ViscosityLaplacian(vDir));
			}
		}
		l_internalForce.at(i) = viscosityForce * mu;
	}
//cout << "ComputeViscosityForce" << endl;
}
