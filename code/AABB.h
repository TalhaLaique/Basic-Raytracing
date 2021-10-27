#pragma once
#include "Ray.h"
#include "Vec3.h"


class AABB {
public:
	AABB() {
		min_corner = Vec3f(0.f,0.f,0.f);
		max_corner = Vec3f(0.f, 0.f, 0.f);
	}
	AABB(Vec3f min, Vec3f max) {
		min_corner = min;
		max_corner = max;
	}
	virtual ~AABB() {};

	Vec3f get_min_corner() {
		return min_corner;
	}

	Vec3f get_max_corner() {
		return max_corner;
	}

	void swap(float& min, float& max) {
		float temp = min;
		min = max;
		max = temp;
	}

	bool ray_intersect(Ray& ray) {
		Vec3f org = ray.origin();
		Vec3f dir = ray.direction();

		//Intersection along x-axis
		float tmin = (min_corner[0] - org[0]) / dir[0];
		float tmax = (max_corner[0] - org[0]) / dir[0];
		if (tmin > tmax) swap(tmin, tmax);
		Vec3f interSXMIN = org + tmin * dir;
		Vec3f interSXMAX = org + tmax * dir;
		if (interSXMIN[1] > min_corner[1] && interSXMIN[1] < max_corner[1] && interSXMIN[2]>min_corner[2] && interSXMIN[2] < max_corner[2]) return true;		
		if (interSXMAX[1] > min_corner[1] && interSXMAX[1] < max_corner[1] && interSXMAX[2]>min_corner[2] && interSXMAX[2] < max_corner[2]) return true;
		

		//Intersection along y-axis
		float tymin = (min_corner[1] - org[1]) / dir[1];
		float tymax = (max_corner[1] - org[1]) / dir[1];
		if (tymin > tymax) swap(tymin, tymax);
		Vec3f interSYMIN = org + tymin * dir;
		Vec3f interSYMAX = org + tymax * dir;
		if (interSYMIN[0] > min_corner[0] && interSYMIN[0] < max_corner[0] && interSYMIN[2]>min_corner[2] && interSYMIN[2] < max_corner[2]) return true;			
		else if (interSYMAX[0] > min_corner[0] && interSYMAX[0] < max_corner[0] && interSYMAX[2]>min_corner[2] && interSYMAX[2] < max_corner[2]) return true;			

		//Intersection along z-axis
		float tzmin = (min_corner[2] - org[2]) / dir[2];
		float tzmax = (max_corner[2] - org[2]) / dir[2];
		if (tzmin > tzmax) swap(tzmin, tzmax);
		Vec3f interSZMIN = org + tzmin * dir;
		Vec3f interSZMAX = org + tzmax * dir;
		if (interSZMIN[0] > min_corner[0] && interSZMIN[0] < max_corner[0] && interSZMIN[1]>min_corner[1] && interSZMIN[1] < max_corner[1]) return true;		
		else if (interSZMAX[0] > min_corner[0] && interSZMAX[0] < max_corner[0] && interSZMAX[1]>min_corner[1] && interSZMAX[1] < max_corner[1]) return true;
				
		return false;
	}
private:
	Vec3f min_corner;
	Vec3f max_corner;
};