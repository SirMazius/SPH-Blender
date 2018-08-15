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
			l_density.at(i) += mass * Kernels::Value(vAux);
		}
	}

	//cout << "ComputeMassDensity" << endl;
}

void InternalForces::ComputeDensityDelta(vector<float> & l_density, vector<float> & l_auxDensity, vector<Vec3> & l_pressureForce,
		const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors) {
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;
	float dt2 = FluidParams::dt2;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_auxDensity.at(i) = 0;
		Vec3 sum2, deltaXi;
		float sum1 = 0.0;

		for (int j : l_neighbors.at(i)) {
			Vec3 vDir, grad, deltaXj;

			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vDir);
			grad = Kernels::ValueGradient(vDir);

			deltaXj = dt2 * l_pressureForce.at(j) / mass;
			sum1 += Vec3::Dot(grad, deltaXj);
			sum2 += grad;
		}

		deltaXi = dt2 * l_pressureForce.at(i) / mass;
		l_auxDensity.at(i) = mass * (Vec3::Dot(sum2, deltaXi) - sum1);
//		cout << "DENSITY T+1 ->>>>>>>>>> " << l_auxDensity.at(i);
	}
}

void InternalForces::ComputePressureCorrection(vector<float> & l_density, vector<float> & l_auxDensity, vector<float> & l_pressures,
		const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors, float & pError) {

	vector<float> l_pError(FluidParams::nParticles);
	vector<float> l_newDensity(FluidParams::nParticles);
	int count = FluidParams::nParticles;
	float beta = FluidParams::beta;
	float restDensity = FluidParams::restDensity;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {

		Vec3 sum1, sum2;
		float sum3 = 0;

		for (int j : l_neighbors.at(i)) {
			Vec3 vDir, grad;
			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vDir);
			grad = Kernels::ValueGradient(vDir);
			sum1 -= grad; // Esto nos lo podemos ahorrar
			sum2 += grad;
			sum3 += Vec3::Dot(grad, grad);
		}
		float delta = (-1/(beta *Vec3::Dot(sum1, sum2) - sum3));
		l_pError.at(i) = std::max(0.0f,(l_density.at(i) + l_auxDensity.at(i)) - restDensity);
		l_pressures.at(i) += delta * l_pError.at(i);
	}

	pError = 0;

//	for (int j = 0; j < count; j++) {
//		pError += l_pError.at(j);
//	}
//	pError /= count;
	pError = std::accumulate(l_pError.begin(), l_pError.end(), 0) / count;

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
