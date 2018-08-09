#include "Collision.h"

void Collision::Collide(vector<Vec3> & l_positions, vector<Vec3> & l_velocity, Vec3 box[2]) {
	int count = FluidParams::nParticles;
	float restitution = FluidParams::restitutionCoef;
	for (int i = 0; i < count; i++) {

		Vec3 *pos = &l_positions[i];
		Vec3 *vel = &l_velocity[i];
		const float b = 0.95;

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

}
