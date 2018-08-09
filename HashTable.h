#pragma once

#include <vector>
#include "Vec3.h"
#include "math.h"
#include "FluidParams.h"

using namespace std;

class HashTable {
public:
	static void Initialize();
	static int Hash(const Vec3 & pos);
	static void Discretize(const Vec3 & pos, Vec3 & r);
	static Vec3 Discretize(const Vec3 & pos);
	static void InsertParticles(const vector<Vec3> & l_pos);
	static void ClearTable();
	static void RetrieveNeighbors(vector<vector<int>> & l, const vector<Vec3> & l_position);
private:
	static int NextPrime(int);
	static Vec3 bbMin, bbMax;
	static vector<vector<int>> buckets;
	static int size;
	static float l;

	const static int p1 = 73856093;
	const static int p2 = 19349663;
	const static int p3 = 83492791;

};
