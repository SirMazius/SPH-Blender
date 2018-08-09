#pragma once

#include <vector>
#include "Vec3.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include "FluidParams.h"
using namespace std;

class BlenderIO {
public:
	static bool ReadParams(string fileName, Vec3 l_bounds[2]);
	static bool ReadPOSVEL(string fileName, vector<Vec3> & l_pos, vector<Vec3> & l_velocity);
	static void WritePOSVEL(string name, float t, int iteration, vector<Vec3> & l_pos, vector<Vec3> & l_velocity);

	static unordered_map<string, string> parametersMap;
	static vector<string> l_parameters;
};
