#include "Integrator.h"

void Integrator::LeapFrog(vector<Vec3> & l_position, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration,
		vector<Vec3> & l_prevV, bool & isFirst) {

	int count = FluidParams::nParticles;
	float dt = FluidParams::dt;
//	if (isFirst) {
//#pragma omp parallel for
//		for (int i = 0; i < count; i++) {
//			Vec3 previousV, laterV;
//			previousV = l_velocity.at(i) - 0.5 * dt * l_acceleration.at(i);
//			laterV = previousV + l_acceleration.at(i) * dt;
//			l_prevV.at(i) = laterV;
//
//			l_position.at(i) = l_position.at(i) + laterV * dt;
//			l_velocity.at(i) = (previousV + laterV) * 0.5 + l_acceleration.at(i) * dt;
//		}
//		cout << "LeapFrog->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
//		isFirst = false;
//		return;
//	}

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 previousV, laterV, AuxPos;
		previousV = l_velocity.at(i) - 0.5 * dt * l_acceleration.at(i);
		laterV = previousV + l_acceleration.at(i) * dt;
//		l_prevV.at(i) = laterV;
		AuxPos = l_position.at(i) + laterV * dt;

		//l_velocity.at(i) = (previousV + laterV) * 0.5 + l_acceleration.at(i) * dt;
		l_velocity.at(i) =  (Vec3::sub(AuxPos, l_position.at(i))) / dt;
		l_position.at(i) = AuxPos;
	}

//#pragma omp parallel for
//	for (int i = 0; i < count; i++) {
//		Vec3 auxV, auxAcc, auxPos;
//		auxV = l_velocity.at(i) + l_acceleration.at(i) * dt * 0.5;
//		auxPos = l_position.at(i) + auxV * dt;
//		auxAcc = (auxPos - l_position.at(i)) / dt;
//		l_velocity.at(i) = auxV + auxAcc * dt * 0.5;
//		l_position.at(i) = auxPos;
//	}
	//cout << "LeapFrog" << endl;
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

void Integrator::ComputeAccelerations(vector<Vec3> & l_acceleration, vector<Vec3> & l_internalForce,
		vector<Vec3> & l_externalForce, vector<float> & l_density) {
	int count = FluidParams::nParticles;
	//float mass = FluidParams::mass;
#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		//if (l_density.at(i) > 0.0001)
		{
			l_acceleration[i] = ((l_internalForce[i] + l_externalForce[i]))/ l_density[i];
			//l_acceleration[i] = (l_internalForce[i] + l_externalForce[i]) / mass;
		}

	}
	// cout << "ComputeAccelerations" << endl;
}
