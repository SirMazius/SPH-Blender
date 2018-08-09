#pragma once

#include <iostream>
#include <math.h>
#include "Vec3.h"

using namespace std;

class Kernels {

public:
	static void Initialize(float);
	/*
		VALUE
	*/
	static float Value(Vec3&);
	static Vec3 ValueGradient(Vec3&);
	static float ValueLaplacian(Vec3&);
	/*
		Spiky
	*/
	static Vec3 SpikyGradient(Vec3&);
//
//	/*
//		Viscosity
//	*/
	static float ViscosityLaplacian(Vec3&);
	static void SayYEy();
	static void Smen();

private:
	static float h, h2;
	static double h6, h9;
	static double valueFactor, valueGradientFactor, valueLaplacianFactor;
	static double spikyGradientFactor;
	static double viscosityLaplacianFactor;
	static constexpr float PI = 3.1415926;
};
