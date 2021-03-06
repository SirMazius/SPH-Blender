#pragma once

#include "BlenderIO.h"
#include <string>
#include "Kernels.h"
using namespace std;
class FluidParams{
public:
	static int nParticles;
	static int simulationSteps;
	static float restDensity;
	static float mass;
	static float viscosity;
	static float surfaceTension;
	static float threshold;
	static float stiffness;
	static float restitutionCoef;
	static float kernelParticles;
	static float kernelRadius;
	static float dt;
	static float dt2;
	static float fluidVolume;
	static float beta;
	static float particleRadius;
	static float particleVolume;
	static float particleOffset;
	static float pciKernelFactor;

	static void Initialize(/*int _nParticles, float _restDensity, float _mass,
			float _viscosity, float _surfaceTension, float _threshold,
			float _stiffness, float _restitutionCoef, float _kernelParticles,
			float _kernelRadius*/);
};
