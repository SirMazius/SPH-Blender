#include "FluidParams.h"
#include "BlenderIO.h"
#include <string>
int FluidParams::nParticles = 0;
int FluidParams::simulationSteps = 0;
float FluidParams::dt = 0;
float FluidParams::dt2 = 0;
float FluidParams::restDensity = 0;
float FluidParams::mass = 0;
float FluidParams::viscosity = 0;
float FluidParams::surfaceTension = 0;
float FluidParams::threshold = 0;
float FluidParams::stiffness = 0;
float FluidParams::restitutionCoef = 0;
float FluidParams::kernelParticles = 0;
float FluidParams::kernelRadius = 0;
float FluidParams::fluidVolume = 0;
float FluidParams::beta = 0;

void FluidParams::Initialize(/*int _nParticles, float _restDensity, float _mass,
		float _viscosity, float _surfaceTension, float _threshold,
		float _stiffness, float _restitutionCoef, float _kernelParticles,
		float _kernelRadius*/) {

//	std::unordered_map<std::string,double>::const_iterator got =

	nParticles = stoi(BlenderIO::parametersMap.find("nparts")->second);
	dt = stof(BlenderIO::parametersMap.find("tstep")->second);
	dt2 = dt * dt;
	restDensity = stof(BlenderIO::parametersMap.find("density")->second);
	mass = 0.02;
	viscosity = stof(BlenderIO::parametersMap.find("visco")->second);
	surfaceTension = stof(BlenderIO::parametersMap.find("surften")->second);

	stiffness = stof(BlenderIO::parametersMap.find("stiff")->second);
	restitutionCoef = stof(BlenderIO::parametersMap.find("collision_restitution")->second);
	kernelParticles = stof(BlenderIO::parametersMap.find("kernel_parts")->second);
	fluidVolume = ( nParticles * mass ) / restDensity;
	kernelRadius = cbrt((3*fluidVolume*kernelParticles)/(4*nParticles*3.14159));//0.0457;
	simulationSteps = stof(BlenderIO::parametersMap.find("tfin")->second) / (dt * 1/*10*/);
	threshold = sqrt(restDensity/kernelParticles);// stof(BlenderIO::parametersMap.find("surften_threshold")->second);

	beta = dt2 * mass * mass  * 2 / (restDensity * restDensity);
}
