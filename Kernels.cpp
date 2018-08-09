#include "Kernels.h"

float Kernels::h = 0;
float Kernels::h2 = 0;
float Kernels::h6 = 0;
float Kernels::h9 = 0;

float Kernels::valueFactor = 0;
float Kernels::valueGradientFactor = 0;
float Kernels::valueLaplacianFactor = 0;
float Kernels::spikyGradientFactor = 0;
float Kernels::viscosityLaplacianFactor = 0;

void Kernels::SayYEy() {
	cout << "HOLA" << endl;
}

void Kernels::Smen() {
	cout << Kernels::h9 << endl;
}

void Kernels::Initialize(float _h) {
	h = _h;
	h2 = h * h;
	h6 = pow(h, 6);
	h9 = pow(h, 9);
	valueFactor = 315 / (64 * PI * h9);
	valueGradientFactor = -945 / (32 * PI * h9);
	valueLaplacianFactor = -945 / (32 * PI * h9);
	spikyGradientFactor = -45 / (PI * h6);
	viscosityLaplacianFactor = -45 / (PI * h6);
	viscosityLaplacianFactor = 45 / (PI * h6);
}

float Kernels::Value(Vec3 & r) {
	if (r.mag() <= h) {
		return valueFactor * pow(h2 - r.norm_squared(), 3);
	}
	return 0;
}

Vec3 Kernels::ValueGradient(Vec3 & r) {
	if (r.mag() <= h) {
		return valueGradientFactor * r * pow(h2 - r.norm_squared(), 2);
	}
	return Vec3();
}

float Kernels::ValueLaplacian(Vec3 & r) {
	if (r.mag() <= h) {
		return valueLaplacianFactor * (h2 - r.norm_squared())
				* (3 * h2 - 7 * r.norm_squared());
	}
	return 0;
}

Vec3 Kernels::SpikyGradient(Vec3 & r) {
	float mag = r.mag();
	if (mag <= h && abs(mag) > 0.001) {
		return spikyGradientFactor * (r / mag) * pow(h - mag, 2);
	}
	return Vec3();
}

float Kernels::ViscosityLaplacian(Vec3 & r) {
	if (r.mag() <= h) {
		return viscosityLaplacianFactor * (h - r.mag());
	}
	return 0;
}

