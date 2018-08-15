#pragma once
#include <vector>
#include <cmath>
#include "Vec3.h"
#include "FluidParams.h"
using namespace std;

class Integrator {
public:
	static void ComputeAccelerations(vector<Vec3> & l_acceleration, vector<Vec3> & l_internalForce, vector<Vec3> & l_pressureForce,
			vector<Vec3> & l_externalForce, vector<float> & l_density);

	static void LeapFrog(vector<Vec3> & l_pos, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration,
			vector<Vec3> & l_prevV, bool & isFirst);

	static void EulerSemi(vector<Vec3> & l_pos, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration);

	static void RK4(vector<Vec3> & l_pos, vector<Vec3> & l_velocity, vector<Vec3> & l_acceleration);

	static Vec3 f(const Vec3 & v1, const Vec3 & v2);
private:
};
