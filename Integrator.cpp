#include "Integrator.h"

void Integrator::LeapFrog(vector<Vec3> & l_pos, vector<Vec3> & l_velocity,
		vector<Vec3> & l_acceleration) {

	Vec3 previousV, laterV, newPos, newV, auxAcc;
	float dt = FluidParams::dt; //TEMPORAL
	float dt2 = dt * dt;
	int count = FluidParams::nParticles;

	for (int i = 0; i < count; i++) {
		laterV = l_velocity[i] + l_acceleration[i] * dt / 2;
		newPos = l_pos[i] + laterV * dt;

		float distance = Vec3::dist(newPos, l_pos[i]);

		Vec3 auxV = (2 * l_velocity[i]) / dt;
		float auxDistance = (2 * distance) / dt2;
		auxAcc.x = auxDistance - auxV.x;
		auxAcc.y = auxDistance - auxV.y;
		auxAcc.z = auxDistance - auxV.z;

		newV = laterV + auxAcc * dt / 2;

		l_pos[i] = newPos;
		l_velocity[i] = newV;
	}
	cout << "LeapFrog" << endl;
}

void Integrator::EulerSemi(vector<Vec3> & l_position, vector<Vec3> & l_velocity,
		vector<Vec3> & l_acceleration) {

	float dt = FluidParams::dt;
	int count = FluidParams::nParticles;

	for (int i = 0; i < count; i++) {
		l_velocity.at(i) = l_velocity.at(i) + l_acceleration.at(i) * dt;
		l_position.at(i) = l_position.at(i) + l_velocity.at(i) * dt;
	}
	cout << "EulerSemi" << endl;

}

void Integrator::ComputeAccelerations(vector<Vec3> & l_acceleration,
		 vector<Vec3> & l_internalForce,
		 vector<Vec3> & l_externalForce,  vector<float> & l_density)
{
	int count = FluidParams::nParticles;
	for (int i = 0; i < count; i++) {
		if (l_density.at(i) !=0) {
			l_acceleration[i] = (l_internalForce[i] + l_externalForce[i]) / l_density[i];
		}

//		l_acceleration[i] = (l_internalForce[i] + l_externalForce[i]) / mass;
	}
	cout << "ComputeAccelerations" << endl;
}
