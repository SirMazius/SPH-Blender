#pragma once
#include <math.h>
#include <iostream>
using namespace std;
class Vec3 {
public:
	Vec3();
	Vec3(float, float, float);
	~Vec3();

	inline Vec3 normalize() {
		return Vec3(this->x, this->y,this->z) / mag();
	}
	inline static float dist(Vec3 & v1, Vec3 & v2) {
		return sqrt(pow(v2.x - v1.x, 2) + pow(v2.y - v1.y, 2) + pow(v2.z - v1.z, 2));
	}

	inline void SetZero() {
		x = y = z = 0;
	}

	inline static void vDirector(const Vec3 & v1, const Vec3 & v2, Vec3 & v3) {
		v3.x = v2.x - v1.x;
		v3.y = v2.y - v1.y;
		v3.z = v2.z - v1.z;
	}

	inline float norm_squared() {
		return (x * x) + (y * y) + (z * z);
	}

	inline float mag() {
		return sqrt((x * x) + (y * y) + (z * z));
	}

	Vec3 operator*(float value) {
		return Vec3(x * value, y * value, z * value);
	}

	friend ostream& operator<<(ostream& os, const Vec3 & v) {
		os << v.x << " " << v.y << " " << v.z << "\n";
		return os;
	}

	friend Vec3 operator*(float value, Vec3 vector) {
		return Vec3(vector.x * value, vector.y * value, vector.z * value);
	}

	Vec3 operator/(float value) {
		return Vec3(x / value, y / value, z / value);
	}

	friend Vec3 operator/(float value, Vec3 vector) {
		return Vec3(vector.x / value, vector.y / value, vector.z / value);
	}

	Vec3 operator+(const Vec3& v) {
		return Vec3(x + v.x, y + v.y, z + v.z);
	}

	Vec3 operator-(const Vec3& v) {
		return Vec3(x - v.x, y - v.y, z - v.z);
	}

	Vec3& operator+=(const Vec3& vAux) {

		this->x += vAux.x;
		this->y += vAux.y;
		this->z += vAux.z;
		return *this;
	}

	Vec3& operator-=(const Vec3& vAux) {

		this->x -= vAux.x;
		this->y -= vAux.y;
		this->z -= vAux.z;
		return *this;
	}

	inline static Vec3 Sum(const Vec3 & v1, const Vec3 & v2) {
		Vec3 vAux;

		vAux.x = v1.x + v2.x;
		vAux.y = v1.y + v2.y;
		vAux.z = v1.z + v2.z;

		return vAux;
	}

	inline static Vec3 sub(const Vec3 & v1, const Vec3 & v2) {
		Vec3 vAux;

		vAux.x = v1.x - v2.x;
		vAux.y = v1.y - v2.y;
		vAux.z = v1.z - v2.z;

		return vAux;
	}

	inline static Vec3 Div(const Vec3 & v1, float i) {
		Vec3 vAux;

		vAux.x = v1.x / i;
		vAux.y = v1.y / i;
		vAux.z = v1.z / i;

		return vAux;
	}

	inline static Vec3 Mult(const Vec3 & v1, float i) {
		Vec3 vAux;

		vAux.x = v1.x * i;
		vAux.y = v1.y * i;
		vAux.z = v1.z * i;

		return vAux;
	}

	inline static float Dot(const Vec3 & v1, const Vec3 & v2) {
		return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
	}

	double x, y, z;
};
