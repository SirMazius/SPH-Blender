//============================================================================
// Name        : HelloEclipseFinal.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <omp.h>
#include "Kernels.h"
#include "BlenderIO.h"
#include "HashTable.h"
#include "InternalForces.h"
#include "ExternalForces.h"
#include "FluidParams.h"
#include "Integrator.h"
#include "Collision.h"
#include <chrono>
#include <algorithm>
//#include "string"

using namespace std;
using namespace std::chrono;

//void Compute_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_positions, vector<Vec3> & l_pressureForce, vector<Vec3> & l_internalForce,
//		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
//		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]);

void Compute_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_positions, vector<Vec3> & l_pressureForce, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2], vector<Vector3d> & l_centers, vector<MatrixXd> & l_Gs,
		vector<double> & l_det);

void Compute_PCI_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_pressureForce, vector<Vec3> & l_positions, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]);
int mode, test;
int main(int argc, char *argv[]) {
	cout << "ATOOOOI ->>> " << atoi(argv[1]) << endl;
	//BlenderIO::WriteExcelData("ExampleFile", 10, 10, 10);
//	float jajaja;
//	cin >> jajaja;
	if (argc == 3) {
		mode = atoi(argv[1]);
		test = atoi(argv[2]);
		cout << "MODE " << mode << endl;
		switch (mode) {
		case 0:
			cout << "EXECUTING SPH" << endl;
			break;
		case 1:
			cout << "EXECUTING ASPH" << endl;
			break;
		case 2:
			cout << "EXECUTING PCI-SPH" << endl;
			break;
		default:
			cout << "WRONG MODE DETECTED" << endl;
			return 0;
			break;
		}
	} else {
		cout << "Error" << endl;
		return 0;
	}

	struct timespec start, stop;

//	if (mode == 4) {
//		cout << "MODO DE SIMULACION" << endl;
//		cout << "0) SPH" << endl << "1) A-SPH" << endl << "2) PCI-SPH" << endl;
//		cin >> mode;
//	}

	Vec3 l_bounds[2];
	vector<Vec3> l_positions, l_velocity;
	cout << "INTRODUZCA NOMBRE DE FICHERO" << endl << endl;
	string name = "Matias1";
	//cin >> name;

	if (!BlenderIO::ReadParams(name, l_bounds)) {
		cout << "ERRROROROROOROROROROROR PARAMS" << endl;
		return 0;
	}

	FluidParams::Initialize();
	Kernels::Initialize(FluidParams::kernelRadius);

	if (!BlenderIO::ReadPOSVEL(name, l_positions, l_velocity)) {
		cout << "ERRROROROROOROROROROROR POS VEL" << endl;
		return 0;
	}

	HashTable::Initialize();

	vector<Vec3> l_pressureForce(FluidParams::nParticles), l_internalForce(FluidParams::nParticles), l_externalForce(FluidParams::nParticles), l_normals(
			FluidParams::nParticles), l_acceleration(FluidParams::nParticles), l_prevV(FluidParams::nParticles);
	vector<float> l_density(FluidParams::nParticles), l_pressures(FluidParams::nParticles), l_color(FluidParams::nParticles);
	vector<vector<int>> l_neighbors(FluidParams::nParticles);

	//bool isFirst = true;
	vector<Vector3d> l_centers(FluidParams::nParticles);
	vector<double> l_det(FluidParams::nParticles);
	vector<MatrixXd> l_Gs(FluidParams::nParticles);
//
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start);
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	cout << endl;

	cout << "RADIUS ->> " << FluidParams::kernelRadius << endl;
	cout << "NRadius ->> " << FluidParams::kernelParticles << endl;
	cout << "NParticles ->> " << FluidParams::nParticles << endl;
	cout << "STEPS ->>" << FluidParams::simulationSteps << endl;
	cout << "SIZE POS->>> " << l_positions.size() << endl;
	cout << "SIZE V->>> " << l_positions.size() << endl;
	cout << "VISCOSITY ->> " << FluidParams::viscosity << endl;
	cout << "STIFF ->> " << FluidParams::stiffness << endl;
	cout << "DT ->> " << FluidParams::dt << endl;
	cout << "THRESHOLD ->> " << FluidParams::threshold << endl;
	cout << "TENSION ->> " << FluidParams::surfaceTension << endl;
	cout << "PARTICLE RADIUS ->> " << FluidParams::particleRadius << endl;
//
	if (mode == 0 || mode == 1)
		Compute_SPH(l_neighbors, l_positions, l_pressureForce, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density, l_pressures,
				l_color, name, l_bounds, l_centers, l_Gs, l_det);
//////
	else if (mode == 2)
		Compute_PCI_SPH(l_neighbors, l_pressureForce, l_positions, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density,
				l_pressures, l_color, name, l_bounds);
	else
		cout << "Not correct mode selected" << endl;

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	double accum = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec * 1e-9 - start.tv_nsec * 1e-9);
	auto duration = duration_cast<microseconds>(t2 - t1).count();

	cout << endl << endl << endl << "EXITO" << endl;
	printf("%lf\n", accum);
	cout << "NEW DURATION-->>> " << duration << " " << duration / 1000000 << endl;
	ofstream timeFile;
	timeFile.open("Time");
		timeFile << duration / 1000000;
	timeFile.close();
	return 0;
}

void Compute_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_positions, vector<Vec3> & l_pressureForce, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2], vector<Vector3d> & l_centers, vector<MatrixXd> & l_Gs,
		vector<double> & l_det) {

	vector<Vec3> l_prevPos(FluidParams::nParticles);
	l_prevPos = l_positions;
	vector<float> l_densityBorder(FluidParams::nParticles);
	for (int i = 1; i < FluidParams::simulationSteps; i++) {

		for (int j = 0; j < 10; j++) {

			HashTable::InsertParticles(l_positions);
			HashTable::RetrieveNeighbors(l_neighbors, l_positions);
////
			if (mode == 1) {
				InternalForces::ComputeAnisotropy(l_neighbors, l_positions, l_centers, l_Gs, l_det);
				InternalForces::ComputeAnisotropyMassDensity(l_density, l_positions, l_neighbors, l_centers, l_Gs, l_det);
			} else {
				InternalForces::ComputeMassDensity(l_density, l_positions, l_neighbors);
			}

//

//			float sDensity = 0;
//			for (auto d : l_density) {
//				if (d > 250)
//					sDensity += d;
//			}
//			sDensity /= FluidParams::nParticles;
//
//			cout << sDensity << endl;

			InternalForces::ComputePressures(l_density, l_pressures, FluidParams::restDensity);
			InternalForces::ComputePressureForce(l_density, l_positions, l_pressures, l_pressureForce, l_neighbors);
			InternalForces::ComputeViscosityForce(l_density, l_velocity, l_internalForce, l_neighbors, l_positions);

			ExternalForces::ComputeGravity(l_externalForce, l_density);
			ExternalForces::ComputeInwardNormal(l_density, l_neighbors, l_positions, l_normals);
			ExternalForces::ComputeColorField(l_density, l_neighbors, l_positions, l_color);
			ExternalForces::ComputeSurfaceTension(l_density, l_densityBorder, l_color, l_normals, l_externalForce);

			Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
			Integrator::LeapFrog(l_positions, l_velocity, l_acceleration, l_prevPos);
//			Integrator::EulerSemi(l_positions, l_velocity, l_acceleration);

			Collision::Collide(l_positions, l_velocity, l_bounds, test);

		}

		float sDensity = 0;
		int densityCounter = 0;
		for (auto d : l_density) {
			if (d > 600) {
				sDensity += d;
				densityCounter++;
			}

		}
		sDensity /= densityCounter;
		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);

		float densityBorder = 0;
		int borderCounter = 0;
		for (auto dB : l_densityBorder) {
			if (dB) {
				densityBorder += dB;
				borderCounter++;
			}
		}
		densityBorder /= borderCounter;
		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << " " << densityBorder << endl;

		BlenderIO::WriteExcelData("ExampleFile", sDensity, densityBorder, 0.0);

	}
}

void Compute_PCI_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_pressureForce, vector<Vec3> & l_positions, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]) {

	//vector<float> l_auxDensity(FluidParams::nParticles);
	vector<Vec3> l_prevPos(FluidParams::nParticles);
	l_prevPos = l_positions;
	vector<float> l_densityBorder(FluidParams::nParticles);
	for (int i = 1; i < FluidParams::simulationSteps; i++) {
		float MaxpError;
		for (int j = 0; j < 10; j++) {
			int iterations = 0;
			MaxpError = 0;

			HashTable::InsertParticles(l_positions);
			HashTable::RetrieveNeighbors(l_neighbors, l_positions);

			InternalForces::ComputeMassDensity(l_density, l_positions, l_neighbors);

			//InternalForces::ComputePressures(l_density, l_pressures, FluidParams::restDensity);
			//InternalForces::ComputePressureForce(l_density, l_positions, l_pressures, l_internalForce, l_neighbors);
			InternalForces::ComputeViscosityForce(l_density, l_velocity, l_internalForce, l_neighbors, l_positions);
			ExternalForces::ComputeGravity(l_externalForce, l_density);
			ExternalForces::ComputeInwardNormal(l_density, l_neighbors, l_positions, l_normals);
			ExternalForces::ComputeColorField(l_density, l_neighbors, l_positions, l_color);
			ExternalForces::ComputeSurfaceTension(l_density, l_densityBorder, l_color, l_normals, l_externalForce);

			std::fill(l_pressures.begin(), l_pressures.end(), 0.0f);
			std::fill(l_pressureForce.begin(), l_pressureForce.end(), Vec3());

			vector<Vec3> l_auxVelocity = l_velocity;
			vector<Vec3> l_auxPos = l_positions;
			vector<Vec3> l_auxPrevPos = l_positions;
			vector<float> l_auxDensity = l_density;

//			float sDensity = 0;
//			for (auto d : l_density) {
//				if (d > 250)
//					sDensity += d;
//			}
//			sDensity /= FluidParams::nParticles;
////
////			 cout << "DENSIDAD INICIAL ->> " << sDensity << endl;
			while (iterations < 3 || MaxpError > 10) {

				Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_auxDensity);
				Integrator::EulerSemi(l_auxPos, l_auxVelocity, l_acceleration);
//				Integrator::LeapFrog(l_auxPos, l_auxVelocity, l_acceleration, l_auxPrevPos);

//				Collision::Collide(l_auxPos, l_auxVelocity, l_bounds);
				HashTable::InsertParticles(l_auxPos);
				HashTable::RetrieveNeighbors(l_neighbors, l_auxPos);

				InternalForces::ComputeMassDensity(l_auxDensity, l_auxPos, l_neighbors);
				InternalForces::ComputePressureCorrection(l_auxDensity, l_pressures, MaxpError);

				InternalForces::ComputePressureForce(l_auxDensity, l_positions, l_pressures, l_pressureForce, l_neighbors);

				l_auxPos = l_positions;
				l_auxVelocity = l_velocity;

				if (iterations > 200) {
					cout << "SALIMOS" << endl;
					break;
				}

//				cout << iterations << endl;
				iterations++;
			}

			cout << iterations << endl;

			Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
			Integrator::EulerSemi(l_positions, l_velocity, l_acceleration);
//			Integrator::LeapFrog(l_positions, l_velocity, l_acceleration, l_prevPos);
			Collision::Collide(l_positions, l_velocity, l_bounds, test);

		}

		float densityBorder = 0;
		int borderCounter = 0;
		for (auto dB : l_densityBorder) {
			if (dB) {
				densityBorder += dB;
				borderCounter++;
			}
		}
		densityBorder /= borderCounter;

		float sDensity = 0;
		int densityCounter = 0;
		for (auto d : l_density) {
			if (d > 600) {
				sDensity += d;
				densityCounter++;
			}
		}
		sDensity /= densityCounter;

		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);
//		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << "  " << MaxpError << endl;
		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << " " << densityBorder << " " << MaxpError << endl;
		BlenderIO::WriteExcelData("ExampleFile", sDensity, densityBorder, 0.0);
	}
}

