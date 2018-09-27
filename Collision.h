#pragma once
#include <iostream>
#include <vector>
#include "Vec3.h"
#include "FluidParams.h"

using namespace std;

class Collision {
public:
	static void Collide(vector<Vec3> & l_positions, vector<Vec3> & l_velocity, Vec3 box[8], int testNumber);
	static void CollideAgainstSpheres(vector<Vec3> & l_positions, vector<Vec3> & l_velocity);
	static void CollideAgainstCapsules(vector<Vec3> & l_positions, vector<Vec3> & l_velocity);
};
