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

void Compute_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_positions, vector<Vec3> & l_pressureForce, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]);

void Compute_PCI_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_pressureForce, vector<Vec3> & l_positions, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]);

int main() {
	struct timespec start, stop;

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

//	Compute_SPH(l_neighbors, l_positions, l_pressureForce, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density, l_pressures, l_color, name,
//			l_bounds);

	Compute_PCI_SPH(l_neighbors, l_pressureForce, l_positions, l_internalForce, l_externalForce, l_velocity, l_acceleration, l_normals, l_density, l_pressures,
			l_color, name, l_bounds);

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	double accum = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec * 1e-9 - start.tv_nsec * 1e-9);
	auto duration = duration_cast<microseconds>(t2 - t1).count();

	cout << endl << endl << endl << "EXITO" << endl;
	printf("%lf\n", accum);
	cout << "NEW DURATION-->>> " << duration << " " << duration / 1000000 << endl;
	return 0;
}

void Compute_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_positions, vector<Vec3> & l_pressureForce, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]) {
	for (int i = 1; i < FluidParams::simulationSteps; i++) {
		for (int j = 0; j < 10; j++) {

			HashTable::InsertParticles(l_positions);
			HashTable::RetrieveNeighbors(l_neighbors, l_positions);

			InternalForces::ComputeMassDensity(l_density, l_positions, l_neighbors);
			InternalForces::ComputePressures(l_density, l_pressures, FluidParams::restDensity);
			InternalForces::ComputePressureForce(l_density, l_positions, l_pressures, l_pressureForce, l_neighbors);
			InternalForces::ComputeViscosityForce(l_density, l_velocity, l_internalForce, l_neighbors, l_positions);

			ExternalForces::ComputeGravity(l_externalForce, l_density);
			ExternalForces::ComputeInwardNormal(l_density, l_neighbors, l_positions, l_normals);
			ExternalForces::ComputeColorField(l_density, l_neighbors, l_positions, l_color);
			ExternalForces::ComputeSurfaceTension(l_color, l_normals, l_externalForce);

			Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
//			Integrator::LeapFrog(l_positions, l_velocity, l_acceleration, l_prevV, isFirst);
			Integrator::EulerSemi(l_positions, l_velocity, l_acceleration);

			Collision::Collide(l_positions, l_velocity, l_bounds);

		}

		float sDensity = 0;
		for (auto d : l_density) {
			if (d > 250)
				sDensity += d;
		}
		sDensity /= FluidParams::nParticles;
		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);
		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << endl;

	}
}

void Compute_PCI_SPH(vector<vector<int>> & l_neighbors, vector<Vec3> & l_pressureForce, vector<Vec3> & l_positions, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<Vec3> & l_velocity, vector<Vec3> l_acceleration, vector<Vec3> l_normals, vector<float> l_density,
		vector<float> l_pressures, vector<float> l_color, string name, Vec3 l_bounds[2]) {

	vector<float> l_auxDensity(FluidParams::nParticles);

	for (int i = 1; i < FluidParams::simulationSteps; i++) {
//		for (int j = 0; j < 10; j++)
		{
			int iterations = 0;
			float pError = 0;

			vector<Vec3> l_auxVelocity = l_velocity;
			vector<Vec3> l_auxPos = l_positions;

			HashTable::InsertParticles(l_positions);
			HashTable::RetrieveNeighbors(l_neighbors, l_positions);

			InternalForces::ComputeMassDensity(l_density, l_positions, l_neighbors);

			//InternalForces::ComputePressures(l_density, l_pressures, FluidParams::restDensity);
			//InternalForces::ComputePressureForce(l_density, l_positions, l_pressures, l_internalForce, l_neighbors);
			InternalForces::ComputeViscosityForce(l_density, l_velocity, l_internalForce, l_neighbors, l_positions);
			ExternalForces::ComputeGravity(l_externalForce, l_density);
			ExternalForces::ComputeInwardNormal(l_density, l_neighbors, l_positions, l_normals);
			ExternalForces::ComputeColorField(l_density, l_neighbors, l_positions, l_color);
			ExternalForces::ComputeSurfaceTension(l_color, l_normals, l_externalForce);

			std::fill(l_pressures.begin(), l_pressures.end(), 0.0f);
			std::fill(l_pressureForce.begin(), l_pressureForce.end(), Vec3());

//			float totalDensity = std::accumulate(l_density.begin(), l_density.end(), 0) / FluidParams::nParticles;
//			totalDensity /= FluidParams::nParticles;
//			pError = totalDensity - FluidParams::restDensity;

			while (iterations < 10 || pError > 100) {

//				if (iterations > 20)
//					break;

				Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
				Integrator::EulerSemi(l_auxPos, l_auxVelocity, l_acceleration);

//				InternalForces::ComputeDensityDelta(l_density, l_auxDensity, l_pressureForce, l_auxPos, l_neighbors);
				InternalForces::ComputePressureCorrection(l_density, l_auxDensity, l_pressures, l_auxPos, l_neighbors, pError);

//				cout << pError << endl;

				InternalForces::ComputePressureForce(l_density, l_auxPos, l_pressures, l_pressureForce, l_neighbors);
				iterations++;
			}

//			Integrator::LeapFrog(l_positions, l_velocity, l_acceleration, l_prevV, isFirst);
			Integrator::ComputeAccelerations(l_acceleration, l_internalForce, l_pressureForce, l_externalForce, l_density);
			Integrator::EulerSemi(l_positions, l_velocity, l_acceleration);

			Collision::Collide(l_positions, l_velocity, l_bounds);

		}

		float sDensity = 0;
		for (auto d : l_density) {
			if (d > 250)
				sDensity += d;
		}
		sDensity /= FluidParams::nParticles;
		BlenderIO::WritePOSVEL(name, i * FluidParams::dt, i, l_positions, l_velocity);
		cout << "///////////////////////>>>>>>>>>" << i << "                     " << sDensity << endl;

	}
}

