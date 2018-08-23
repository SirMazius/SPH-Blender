#include "Integrator.h"

void Integrator::LeapFrog(vector<Vec3> & l_position, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration, vector<Vec3> & l_prevPos) {

	int count = FluidParams::nParticles;
	float dt = FluidParams::dt;
	float dt2 = FluidParams::dt2;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
//		Vec3 last_position = l_position.at(i);
//
//		Vec3 pastVelocity = l_prevV.at(i);
//		Vec3 futureVelocity = pastVelocity + l_acceleration.at(i) * dt;
//		l_prevV.at(i) = futureVelocity;
//		l_position.at(i) = l_position.at(i) + futureVelocity * dt;
//		l_velocity.at(i) = (l_position.at(i) - last_position) / dt;

		Vec3 position = l_position.at(i);
		l_position.at(i) = 2 * position - l_prevPos.at(i) + l_acceleration.at(i) * dt2;
		l_prevPos.at(i) = position;
		l_velocity.at(i) = (l_position.at(i) - position) / dt;
	}
}

void Integrator::EulerSemi(vector<Vec3> & l_position, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration) {

	float dt = FluidParams::dt;
	int count = FluidParams::nParticles;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		l_velocity.at(i) = l_velocity.at(i) + l_acceleration.at(i) * dt;
		l_position.at(i) = l_position.at(i) + l_velocity.at(i) * dt;
	}

// VERLET
//#pragma omp parallel for
//	for (int i = 0; i < count; i++) {
//		l_position.at(i) = l_position.at(i)+ l_velocity.at(i) * dt + 0.5f * l_acceleration.at(i)* dt * dt;
//		l_velocity.at(i) = (l_position.at(i) - l_position.at(i)) / dt;
//	}

//#pragma omp parallel for
//	for (int i = 0; i < count; i++) {
//		Vec3 oldVel = l_velocity.at(i);
//		l_velocity.at(i) = l_velocity.at(i) + l_acceleration.at(i) * dt;
//		l_position.at(i) = l_position.at(i) + (oldVel + l_velocity.at(i)) * 0.5 * dt;
//	}

// cout << "EulerSemi" << endl;

}

void Integrator::ComputeAccelerations(vector<Vec3> & l_acceleration, vector<Vec3> & l_internalForce, vector<Vec3> & l_pressureForce,
		vector<Vec3> & l_externalForce, vector<float> & l_density) {
	int count = FluidParams::nParticles;
	//float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		//if (l_density.at(i) > 0.0001)
		{
			l_acceleration[i] = ((l_internalForce[i] + l_externalForce[i]) + l_pressureForce[i]) / l_density[i];
			//l_acceleration[i] = (l_internalForce[i] + l_externalForce[i]) / mass;
		}

	}
	// cout << "ComputeAccelerations" << endl;
}
