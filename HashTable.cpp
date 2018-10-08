#include "HashTable.h"

int HashTable::size = 0;
float HashTable::l = 0;
vector<vector<int>> HashTable::buckets;
Vec3 HashTable::bbMin;
Vec3 HashTable::bbMax;

void HashTable::Initialize() {
	size = NextPrime(FluidParams::nParticles * 2); //TEMPORAL
	l = FluidParams::kernelRadius;
	buckets = vector<vector<int>>(size);
//	buckets.at(10).push_back(2);
}

int HashTable::Hash(const Vec3 & pos) {
//	int rX, rY, rZ;
//
//	rX = (int) pos.x * p1;
//	rY = (int) pos.y * p2;
//	rZ = (int) pos.z * p3;

	unsigned int h = (((int) pos.x) * 73856093) ^ (((int) pos.y) * 19349663) ^ (((int) pos.z) * 83492791);
	h %= size;
	return h;
}

void HashTable::Discretize(const Vec3 & pos, Vec3 & r) {
	r.x = floor(pos.x / l);
	r.y = floor(pos.y / l);
	r.z = floor(pos.z / l);
}

Vec3 HashTable::Discretize(const Vec3 & pos) {
	Vec3 aux;
	aux.x = floor(pos.x / l);
	aux.y = floor(pos.y / l);
	aux.z = floor(pos.z / l);

	return aux;
}

void HashTable::InsertParticles(const vector<Vec3> & l_pos) {

	Vec3 auxV;
	int count = 0;
	ClearTable();
	for (const Vec3 & p : l_pos) {

		Discretize(p, auxV);
		int index = Hash(auxV);
		buckets.at(index).push_back(count);
		count++;
	}

//	for (auto & bucket : buckets)
//		for (auto index : bucket)
//			cout << index << endl;
}

void HashTable::ClearTable() {
	for (vector<int> & v : buckets)
		v.clear();
}

void HashTable::ClearTableX() { // MAS LENTA
#pragma omp parallel for
	for (int i = 0; i < size; i++) {
		buckets.at(i).clear();
	}
}

int HashTable::NextPrime(int n) { // HAY QUE PASARLE EL SIZE QUE TOQUE
	bool found = false;
	int prime = n;

	while (!found) {
		bool searching = true;
		prime++;

		if (prime == 2 || prime == 3) {
			found = true;
			return prime;
		}

		if (prime % 2 == 0 || prime % 3 == 0) {
			found = false;
			searching = false;
		}

		int divisor = 6;

		while (searching && divisor * divisor - 2 * divisor + 1 <= prime) {
			if (prime % (divisor - 1) == 0 || prime % (divisor + 1) == 0) {
				found = false;
				searching = false;
			}

			divisor += 6;
		}

		if (searching) {
			found = true;
		}

	}

	return prime;
}

void HashTable::RetrieveNeighbors(vector<vector<int>> & l_neighbors, const vector<Vec3> & l_position) {

	int count = FluidParams::nParticles;
	double hLenght = FluidParams::kernelRadius;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 bbMin, bbMax, vAux;
		const Vec3 *pos = &l_position[i];

		bbMin.x = pos->x - hLenght;
		bbMin.y = pos->y - hLenght;
		bbMin.z = pos->z - hLenght;
		bbMin = Discretize(bbMin);

		bbMax.x = pos->x + hLenght;
		bbMax.y = pos->y + hLenght;
		bbMax.z = pos->z + hLenght;
		bbMax = Discretize(bbMax);

		//ClearTable();
		l_neighbors.at(i).clear();

		for (vAux.x = bbMin.x; vAux.x <= bbMax.x; vAux.x++) {
			for (vAux.y = bbMin.y; vAux.y <= bbMax.y; vAux.y++) {
				for (vAux.z = bbMin.z; vAux.z <= bbMax.z; vAux.z++) {

					int key = Hash(vAux);
					for (int j : buckets[key])
						l_neighbors.at(i).push_back(j);

				}
			}
		}
	}
}

