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
float FluidParams::particleRadius = 0;
float FluidParams::particleVolume = 0;
float FluidParams::particleOffset = 0;
float FluidParams::pciKernelFactor = 0;

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
	fluidVolume = (nParticles * mass) / restDensity;
	kernelRadius = cbrt((3 * fluidVolume * kernelParticles) / (4 * nParticles * 3.14159)); //0.0457;
	simulationSteps = stof(BlenderIO::parametersMap.find("tfin")->second) / (dt * 10);
	threshold = sqrt(restDensity / kernelParticles); // stof(BlenderIO::parametersMap.find("surften_threshold")->second);

	beta = dt2 * mass * mass * 2 / (restDensity * restDensity);
	particleRadius = cbrt(pow(kernelRadius, 3) / kernelParticles);
	particleVolume = 4.0 / 3.0 * 3.14159 * pow(particleRadius,3);
	//self.particle_offset = (3.0/(4.0*np.pi)*self.particle_volume)**(1.0/3.0)
	particleOffset = pow(3.0 / (4.0*3.14159) * particleVolume,1.0/3.0);
//	Vec3 sum1, sum2;
//
//	float sum3 = 0, h = kernelRadius, pH = particleRadius;
//
//	float multi = 2.0;
//
//	for (float x = -h - pH; x <= h + pH; x += multi * pH) {
//		for (float y = -h - pH; y <= h + pH; y += multi * pH) {
//			for (float z = -h - pH; z <= h + pH; z += multi * pH) {
//
//				Vec3 r(x, y, z);
//				Vec3 vDir;
//				// Vec3::vDirector(r, center, vDir);
//				Vec3 grad = Kernels::SpikyGradient(r);
//				sum1 -= grad;
//				sum2 += grad;
//				sum3 += Vec3::Dot(grad, grad);
//			}
//		}
//	}
//	pciKernelFactor = -1 / (beta * (Vec3::Dot(sum1, sum2) - sum3));
//	cout << "FACTOR" << pciKernelFactor << endl;
}
