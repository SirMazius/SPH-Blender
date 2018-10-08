#include "InternalForces.h"
#include  "FluidParams.h"
#include "Vec3.h"
#include "Kernels.h"
#include <cstddef>
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

void InternalForces::ComputeAnisotropyMassDensity(vector<float> & l_density, const vector<Vec3> & l_positions, const vector<vector<int>> & l_neighbors,
		vector<Vector3d> l_centers, vector<MatrixXd> l_Gs, vector<double> l_det) {
	int count = FluidParams::nParticles;
	float mass = FluidParams::mass;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {

		l_density.at(i) = 0;
		double det = l_det.at(i);
		MatrixXd G = l_Gs.at(i);
		Vector3d c = l_centers.at(i);

		for (int j : l_neighbors.at(i)) {
			Vec3 vAux;
			Vector3d vEigenAux(c[0] - l_positions.at(j).x, c[1] - l_positions.at(j).y, c[2] - l_positions.at(j).z);
			Vector3d vDirectorAux = G * vEigenAux;
			vAux.x = vDirectorAux[0];
			vAux.y = vDirectorAux[1];
			vAux.z = vDirectorAux[2];
			l_density.at(i) += mass * det * Kernels::Value(vAux);
		}
	}

	cout << "ComputeAnisotropyMassDensity" << endl;
}

void InternalForces::GetMatrix(MatrixXd & input, MatrixXd & G, Vector3d & center, double & det, Vector3d & currentPoint) {

	float particleOffset = FluidParams::particleOffset;
//	MatrixXd b(2, 3);
//	b << 2.5, 2.4, 7.0, 0.5, 0.7, 3.0;
//	input = b;
	if (input.rows() < 10) {
		G = Vector3d(1, 1, 1).asDiagonal();
//		cout << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGg" << endl <<G << endl;
		center = currentPoint;
		det = 1.0;
		return;
	}

	MatrixXd a = input.transpose();
	VectorXd mean = input.colwise().mean();
	MatrixXd var = (a).colwise() - mean;
	MatrixXd t = var * var.transpose();
	MatrixXd dst = t * (1 / (a.cols() - 1)); //S
	center = mean;
//
	EigenSolver<MatrixXd> es(dst);
	VectorXd eigenValues = es.eigenvalues().real();
	MatrixXd eigenVectors = es.eigenvectors().real();
// (2*np.sqrt(abs(L)) + 2.0*self.particle_offset) /self.h
	VectorXd ab(3);
	ab[0] = (2.0 * sqrt(abs(eigenValues[0])) + 2.0 * particleOffset) / FluidParams::kernelRadius;
	ab[1] = (2.0 * sqrt(abs(eigenValues[1])) + 2.0 * particleOffset) / FluidParams::kernelRadius;
	ab[2] = (2.0 * sqrt(abs(eigenValues[2])) + 2.0 * particleOffset) / FluidParams::kernelRadius;

	for (int i = 0; i < 3; i++)
		if (ab[i] <= 0.0)
			ab[i] = 1;

//	Vector3d h;
//	h << FluidParams::kernelRadius, FluidParams::kernelRadius, FluidParams::kernelRadius;
//	h << 0.045, 0.045, 0.045;
	// Vector3d diagG = h.array() / ab.array();
	Vector3d diagG(1.0 / ab[0], 1.0 / ab[1], 1.0 / ab[2]);
	MatrixXd diagGmatrix = diagG.asDiagonal();

	G = eigenVectors * (diagGmatrix * eigenVectors.transpose());
	det = diagG[0] * diagG[1] * diagG[2];

//	cout << "Cov -> " << endl << dst << endl << endl;
//	cout << "EigenValues -> " << endl << eigenValues<< endl << endl;
//	cout << "EigenVectors -> " << endl << eigenVectors << endl << endl;
//	cout << "Diag -> " << endl << diagG << endl << endl;
//	cout << "G -> " << endl << G << endl;
//	float jja;
//	cin >> jja;

//		cout << "abR ->> " << endl << abReal << endl;
//		cout << "DiagG ->> " << endl << diagG.asDiagonal() << endl;
//		cout << "DIMENSIONES " << ab.rows() << "   " << ab.cols() << endl;
//		std::cout << "AB -> " << endl << ab << std::endl << endl;
//		cout << "Cov -> " << endl << dst << endl;
//		cout << c << << endl;
//		cout << "PseudoEigenValues -> " << endl << pseudoValues << endl;
//		cout << "PseudoEigenVectors -> " << endl <<  pseudoVectors << endl;
//		cout << "EigenValues -> " << endl << eigenValues[0].real() << endl;
//		cout << "EigenVectors -> " << endl << eigenVectors << endl;
//		cout << "REAL ->> "<<eigenValues.real() << endl;
//		cout << "IMAG ->> "<<eigenValues.imag() << endl;
}

void InternalForces::ComputeAnisotropy(const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions, vector<Vector3d> & l_centers,
		vector<MatrixXd> & l_Gs, vector<double> & l_det) {
	int count = FluidParams::nParticles;
	float h = FluidParams::kernelRadius;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {

		vector<Vector3d> l_points;

		for (int j : l_neighbors.at(i)) {
			Vec3 vAux;
			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vAux);
			if (vAux.mag() <= h) {
				Vector3d point(l_positions.at(j).x, l_positions.at(j).y, l_positions.at(j).z);
				l_points.push_back(point);
			}
		}

		int auxCounter = l_points.size();
		MatrixXd points(auxCounter, 3);

		for (int j = 0; j < auxCounter; j++)
			points.row(j) = l_points.at(j);

		Vector3d center;
		Vector3d currentPoint(l_positions.at(i).x, l_positions.at(i).y, l_positions.at(i).z);
		double det;
		MatrixXd G;
		GetMatrix(points, G, center, det, currentPoint);
		l_centers.at(i) = center;
		l_det.at(i) = det;
		l_Gs.at(i) = G;
	}
	cout << "ComputeAnisotropy" << endl;
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
		float errorpercent = (error * 100) / restDensity;
		if (errorpercent < 30/*30*/) {
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
		float density_i = l_density.at(i);

		// l_pressureForce.at(i).SetZero();
		for (int j : l_neighbors[i]) {
			float density_j = l_density.at(j);
			if (i != j) {
				Vec3 vAux;
				Vec3::vDirector(l_positions[j], l_positions[i], vAux);

				pressureAux = pressureAux
						+ (l_pressures.at(i) / (density_i * density_i) + l_pressures.at(j) / (density_j * density_j)) * mass
								* Kernels::SpikyGradient(vAux);

			}

		}
		l_pressureForce.at(i) = pressureAux * -density_i;
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
