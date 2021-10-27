#pragma once
#include "vec3.h"
#include <vector>

class WorleyNoise {
public:

	WorleyNoise(){}

	virtual ~WorleyNoise() {};

	WorleyNoise(int points, Vec3f min, Vec3f max) :nb_feat_point(points) {
		feat_points.resize(points);
		for (int i = 0; i < points; i++) {
			float dx11 = max[0] - min[0];
			float dy11 = max[1] - min[1];
			float dz11 = max[2] - min[2];
			float jitterX = (float)(rand() % 2000) / 2000 * dx11;
			float jitterY = (float)(rand() % 2000) / 2000 * dy11;
			float jitterZ = (float)(rand() % 2000) / 2000 * dz11;
			feat_points[i] = Vec3f{min[0] + jitterX, min[1] + jitterY, min[2] + jitterZ};
		}
	}

	float Noise(Vec3f position) {
		float dist1 = 999999.f;
		float dist2 = 999999.f;
		for (int i = 0; i < nb_feat_point; i++) {
			float d = (position - feat_points[i]).squaredLength();
			if (d <= dist1) {
				dist2 = dist1; dist1 = d;
			}
			else if (d <= dist2)
				dist2 = d;
		}
		return sqrt(dist2);
	}

private:
	int nb_feat_point;
	std::vector<Vec3f> feat_points;

};