#pragma once
#include "Vec3.h"
#include <random>
#include <ctime>
#include<cstdlib>
#include<map>

using namespace std;

class LightSource {
public:
	inline LightSource(const Vec3f& position, const Vec3f color,float intensity) :
		m_position(position), m_color(color), m_intensity(intensity) {
	}

	inline Vec3f Position() {
		//Uncomment for below for light area source configuration
		/*Vec3f corner = Vec3f(-1.f, -1.f, -1.f);
		Vec3f u = normalize(cross(corner, Vec3f(0.f, 1.f, 0.f)));
		Vec3f v = normalize(cross(u, corner));
		float x = (float)((rand() % 1000) / 500.f) - 1;
		float y = (float)((rand() % 1000) / 500.f) - 1;
		Vec3f uvec = 0.2f * x * u;
		Vec3f vvec = 0.2f * y * v;
		return m_position = m_position + uvec + vvec;*/
		return m_position;// +uvec + vvec;
	};

	inline float Intensity() const {
		return m_intensity;
	};
private:
	Vec3f m_position;
	Vec3f m_color;
	float m_intensity;
};