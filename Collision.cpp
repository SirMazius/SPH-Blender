#include "Collision.h"
#include "Vec3.h"
void Collision::Collide(vector<Vec3> & l_positions, vector<Vec3> & l_velocity, Vec3 box[2], int testNumber) {
	int count = FluidParams::nParticles;
	float restitution = FluidParams::restitutionCoef;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {

		Vec3 *pos = &l_positions[i];
		Vec3 *vel = &l_velocity[i];
		const float b = 0.9999;

		if (pos->x < box[0].x) {
			pos->x = box[0].x;
			if (vel->x < 0.0) {
				vel->x *= -restitution;
			}
			vel->y *= b;
			vel->z *= b;
		} else if (pos->x > box[1].x) {
			pos->x = box[1].x;
			if (vel->x > 0.0) {
				vel->x *= -restitution;
			}
			vel->y *= b;
			vel->z *= b;
		}

		if (pos->y < box[0].y) {
			pos->y = box[0].y;
			if (vel->y < 0.0) {
				vel->y *= -restitution;
			}
			vel->x *= b;
			vel->z *= b;
		} else if (pos->y > box[1].y) {
			pos->y = box[1].y;
			if (vel->y > 0.0) {
				vel->y *= -restitution;
			}
			vel->x *= b;
			vel->z *= b;
		}

		if (pos->z < box[0].z) {
			pos->z = box[0].z;
			if (vel->z < 0.0) {
				vel->z *= -restitution;
			}
			vel->x *= b;
			vel->y *= b;
		} else if (pos->z > box[1].z) {
			pos->z = box[1].z;
			if (vel->z > 0.0) {
				vel->z *= -restitution;
			}
			vel->x *= b;
			vel->y *= b;
		}
	}

	if (testNumber == 4)
		CollideAgainstSpheres(l_positions, l_velocity);
	else if (testNumber == 5)
		CollideAgainstCapsules(l_positions, l_velocity);
}

void Collision::CollideAgainstSpheres(vector<Vec3> & l_positions, vector<Vec3> & l_velocity) {
	int count = FluidParams::nParticles;
	Vec3 sphereCenter(-2.011, -2.97, 2.13);
	float radius = 0.520f;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {
		Vec3 vDir = l_positions.at(i) - sphereCenter;
		float distance = vDir.mag();
		if (distance < radius) {
			Vec3 normal = vDir.normalize();
			l_positions.at(i) = normal * radius + sphereCenter;
			l_velocity.at(i) = l_velocity.at(i) - (1 + 0.95f) * Vec3::Dot(l_velocity.at(i), normal) * normal;
		}
	}
}

inline void CrossProduct(const Vec3 & u, const Vec3 & v, Vec3 & result) {
	result.x = u.y * v.z - u.z * v.y;
	result.y = -1 * (u.x * v.z - u.z * v.x);
	result.z = u.x * v.y - u.y - v.x;
}

void Collision::CollideAgainstCapsules(vector<Vec3> & l_positions, vector<Vec3> & l_velocity) {
	Vec3 line1(-2, -2.9, 1.8), line2(-2, -2.9, 4);
	float radius = 0.120f;
	int count = FluidParams::nParticles;

#pragma omp parallel for
	for (int i = 0; i < count; i++) {

		Vec3 pos = l_positions[i];

		if (pos.z > line1.z || pos.z < line2.z) {

			Vec3 vDir, vDir2, vCross;
			Vec3::vDirector(line1, pos, vDir);
			Vec3::vDirector(line1, line2, vDir2);
			CrossProduct(vDir, vDir2, vCross);
			float distance = vCross.mag() / vDir2.mag();

			if (distance < radius) {
				Vec3 AuxNormal;
				Vec3 sphereCenter = vDir2.normalize() * sqrt(pow(vDir.mag(), 2) + pow(distance, 2)) + line1;
				Vec3::vDirector(sphereCenter, pos, AuxNormal);
				Vec3 normal = AuxNormal.normalize();
				l_positions.at(i) = normal * radius + sphereCenter;
				l_velocity.at(i) = l_velocity.at(i) - (1 + 0.95f) * Vec3::Dot(l_velocity.at(i), normal) * normal;
			}
		}

	}

}

