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
#include <future>
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

void add(int x, int y) {
//    return x+y;
}

int main(int argc, char *argv[]) {

	cout << "LA SIMULACION USARA ESTOS HILOS ->>> " << 4 << endl;
//	omp_set_dynamic(0);
//	omp_set_num_threads(1);

	cout << "ATOOOOI ->>> " << atoi(argv[1]) << endl;
	//BlenderIO::WriteExcelData("ExampleFile", 10, 10, 10);

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

	Vec3 l_bounds[2];
	vector<Vec3> l_positions, l_velocity;
	cout << "INTRODUZCA NOMBRE DE FICHERO" << endl << endl;
	string name = "Matias1";

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

	if (mode == 0 || mode == 1) {
		Compute_SPH(l_neighbors, l_positions, l_pressureForce, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density, l_pressures,
				l_color, name, l_bounds, l_centers, l_Gs, l_det);
	} else if (mode == 2) {
		Compute_PCI_SPH(l_neighbors, l_pressureForce, l_positions, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density,
				l_pressures, l_color, name, l_bounds);
	} else {
		cout << "Not correct mode selected" << endl;
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(t2 - t1).count();

	cout << endl << endl << endl << "EXITO" << endl;
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

	high_resolution_clock::time_point taux1, taux2;

	int heightCounter = 50, auxHeightCounter = 0;
	double accumTime = 0;

	vector<Vec3> l_prevPos(FluidParams::nParticles);
	l_prevPos = l_positions;
	vector<float> l_densityBorder(FluidParams::nParticles);

	for (int i = 1; i < FluidParams::simulationSteps; i++) {
		high_resolution_clock::time_point t3 = high_resolution_clock::now();
		for (int j = 0; j < 10; j++) {

			taux1 = high_resolution_clock::now();
			HashTable::InsertParticles(l_positions);
			taux2 = high_resolution_clock::now();
			auto auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//			cout << "INSERT DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

			taux1 = high_resolution_clock::now();
			HashTable::RetrieveNeighbors(l_neighbors, l_positions);
			taux2 = high_resolution_clock::now();
			auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//			cout << "RETRIEVE DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

			if (mode == 1) {
				InternalForces::ComputeAnisotropy(l_neighbors, l_positions, l_centers, l_Gs, l_det);
				InternalForces::ComputeAnisotropyMassDensity(l_density, l_positions, l_neighbors, l_centers, l_Gs, l_det);
			} else {
				taux1 = high_resolution_clock::now();
				InternalForces::ComputeMassDensity(l_density, l_positions, l_neighbors);
				taux2 = high_resolution_clock::now();
				auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//				cout << "DENSITY DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;
			}

			taux1 = high_resolution_clock::now();
			InternalForces::ComputePressures(l_density, l_pressures, FluidParams::restDensity);
			InternalForces::ComputePressureForce(l_density, l_positions, l_pressures, l_pressureForce, l_neighbors);
			InternalForces::ComputeViscosityForce(l_density, l_velocity, l_internalForce, l_neighbors, l_positions);
			taux2 = high_resolution_clock::now();
			auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//			cout << "INTERNAL DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

			taux1 = high_resolution_clock::now();
			ExternalForces::ComputeGravity(l_externalForce, l_density);
			ExternalForces::ComputeInwardNormal(l_density, l_neighbors, l_positions, l_normals);
			ExternalForces::ComputeColorField(l_density, l_neighbors, l_positions, l_color);
			ExternalForces::ComputeSurfaceTension(l_density, l_densityBorder, l_color, l_normals, l_externalForce);
			taux2 = high_resolution_clock::now();
			auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//			cout << "EXTERNAL DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

			taux1 = high_resolution_clock::now();
			Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
			Integrator::LeapFrog(l_positions, l_velocity, l_acceleration, l_prevPos);
//			Integrator::EulerSemi(l_positions, l_velocity, l_acceleration);
			Collision::Collide(l_positions, l_velocity, l_bounds, test);
			taux2 = high_resolution_clock::now();
			auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
//			cout << "INTEGRATION DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

		}

		float sDensity = 0;
		int densityCounter = 0;
		for (auto d : l_density) {
			if (d > 600) {
				sDensity += d;
				densityCounter++;
			}

		}

		if (heightCounter == 50) {
//			BlenderIO::WriteHeightDensityData(l_positions, l_density, auxHeightCounter);
			std::async(std::launch::async, BlenderIO::WriteHeightDensityData, l_positions, l_density, auxHeightCounter);
			auxHeightCounter++;
			heightCounter = 0;
		}
		heightCounter++;
		taux1 = high_resolution_clock::now();

//		BlenderIO::WritePOSVEL( name, i * FluidParams::dt, i, l_positions, l_velocity);
		std::async(std::launch::async, BlenderIO::WritePOSVEL, name, i * FluidParams::dt, i, l_positions, l_velocity);
		taux2 = high_resolution_clock::now();
		auto auxDuration = duration_cast<microseconds>(taux2 - taux1).count();
		cout << "WRITING DURATION-->>> " << auxDuration << " " << auxDuration / 1000000 << endl;

		high_resolution_clock::time_point t4 = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(t4 - t3).count();
		cout << "NEW DURATION-->>> " << duration << " " << duration / 1000000 << endl;
		accumTime += duration;
		cout << "ACCUM DURATION -->>> " << accumTime / 1000000 << endl;
//		float exsaif;
//		cin >> exsaif;

		sDensity /= densityCounter;

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

//		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);
		BlenderIO::WriteExcelData("ExampleFile", sDensity, densityBorder, 0.0);
	}
}

void Compute_PCI_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_pressureForce, vector<Vec3> & l_positions, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]) {

	//vector<float> l_auxDensity(FluidParams::nParticles);
	int heightCounter = 50, auxHeightCounter = 0;
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
			while (iterations < 6 || MaxpError > 5) {

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

		if (heightCounter == 50) {
			//			BlenderIO::WriteHeightDensityData(l_positions, l_density, auxHeightCounter);
			std::async(std::launch::async, BlenderIO::WriteHeightDensityData, l_positions, l_density, auxHeightCounter);
			auxHeightCounter++;
			heightCounter = 0;
		}
		heightCounter++;

		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);
//		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << "  " << MaxpError << endl;
		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << " " << densityBorder << " " << MaxpError << endl;
		BlenderIO::WriteExcelData("ExampleFile", sDensity, densityBorder, 0.0);
	}
}
