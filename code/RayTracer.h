#pragma once

#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include<omp.h>

#include "Vec3.h"
#include "Image.h"
#include "Camera.h"
#include "Scene.h"
#include "LightSource.h"
#include "AABB.h"

using namespace std;

class RayTracer {
public:
	RayTracer () {}
	virtual ~RayTracer() {}
	   
	inline bool rayTrace (const Ray & ray, 
						  const Scene & scene, 
						  size_t & meshIndex, 
						  size_t & triangleIndex, 
						  float & u, 
						  float & v, 
						  float & d) {
		const auto& meshes = scene.meshes();
		float closest = std::numeric_limits<float>::max();
		bool intersectionFound = false;
		for (size_t mIndex = 0; mIndex < meshes.size(); mIndex++) {
			const auto& P = meshes[mIndex].vertexPositions();
			const auto& T = meshes[mIndex].indexedTriangles();
			for (size_t tIndex = 0; tIndex < T.size(); tIndex++) {
				const Vec3i& triangle = T[tIndex];
				float ut, vt, dt;
				if (ray.triangleIntersect(P[triangle[0]], P[triangle[1]], P[triangle[2]], ut, vt, dt) == true) {
					if (dt > 0.f && dt < closest) {
						intersectionFound = true;
						closest = dt;
						meshIndex = mIndex;
						triangleIndex = tIndex;
						u = ut;
						v = vt;
						d = dt;
					}
				}
			}
		}
		return intersectionFound;
	}

	bool isvisible(Ray& ray, const Scene& scene) {
		const auto& meshes = scene.meshes();
		float closest = std::numeric_limits<float>::max();
		for (size_t mIndex = 0; mIndex < meshes.size(); mIndex++) {
			const auto& P = meshes[mIndex].vertexPositions();
			const auto& T = meshes[mIndex].indexedTriangles();
			for (size_t tIndex = 0; tIndex < T.size(); tIndex++) {
				const Vec3i& triangle = T[tIndex];
				float ut, vt, dt;
				if (ray.triangleIntersect(P[triangle[0]], P[triangle[1]], P[triangle[2]], ut, vt, dt) == true) {
					closest = dt;
					return false;
				}
			}
		}
		return true;
	}

	inline Vec3f shade (Scene & scene, size_t meshIndex, size_t triangleIndex, float u, float v, LightSource lightsource) {
		const auto& mesh = scene.meshes()[meshIndex];
		const auto& P = mesh.vertexPositions();
		const auto& N = mesh.vertexNormals();
		const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);

		//Shadding
		Vec3f intersection = (1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];
		Vec3f wo = normalize(scene.camera().position() - intersection);
		Material material = scene.material();
		float intensity = lightsource.Intensity();
		Vec3f light_position = lightsource.Position();
		float distance = length(light_position - intersection);
		Vec3f wi = normalize(light_position - intersection);

		//Radiance Responses of different BRDFs
		Vec3f response_MFM = material.evaluateColorResponse_MFM(hitNormal, wi, wo, intersection);
		//Vec3f response_MFM = material.evaluateColorResponse_Diff(hitNormal, wi);
		//Vec3f response_DIFF = material.evaluateColorResponse_Diff(hitNormal, wi);
		//Vec3f response_PHNG = material.evaluateColorResponse_Phong(hitNormal, wi, wo);
		//Vec3f response_BPHNG = material.evaluateColorResponse_Blinn_Phong(hitNormal, wi, wo);

		//Vec3f wAvg = (0.5f * response_MFM) + (0.5f * response_BPHNG);
		//Vec3f response = (response_MFM + response_BPHNG) / wAvg;

		return (intensity * response_MFM) / pow(distance, 2); //Intensity of light for every single pixel
	}

	inline void render (Scene& scene, Image& image) {
		size_t w = image.width();
		size_t h = image.height();
		Vec3f color, color1; Vec3f accum_color, accum_color1;
		const Camera& camera = scene.camera();	
		for (int y = 0; y < h; y++) {
			static int progress = 0;
			progress++;
			if (progress % 10 == 0)
				std::cout << ".";
#pragma omp parallel for
			for (int x = 0; x < w; x++) {
				color = Vec3f(0.f, 0.f, 0.f);		
				accum_color = Vec3f(0.f, 0.f, 0.f);
				std::vector<bool> bb_interesect;
				bb_interesect.resize(scene.meshes().size());
				//for (int n = 0; n < 3; n++) {
					for (int l = 0; l < scene.lightSource().size(); l++) {
						LightSource lightS = scene.lightSource()[l];
						//random sampling of rays on light sources
						//float xp = (float)(rand() % 250) / 250 - 0.10;
						//float yp = (float)(rand() % 250) / 250 - 0.10;
						//Ray ray = camera.rayAt(float(x + xp) / w, 1.f - float(y + yp) / h);
						Ray ray = camera.rayAt(float(x) / w, 1.f - float(y) / h);
						for (int m = 0; m < scene.meshes().size(); m++) {
							Mesh c_mesh = scene.meshes()[m];
							AABB c_AABB = c_mesh.get_AABB();
							bb_interesect[m] = c_AABB.ray_intersect(ray);
						}
						for (int m_meshes = 0; m_meshes < scene.meshes().size(); m_meshes++) {
							size_t meshIndex, triangleIndex;
							float u, v, d;
							if (bb_interesect[m_meshes]) {
								bool intersectionFound = rayTrace(ray, scene, meshIndex, triangleIndex, u, v, d);
								if (intersectionFound && d > 0.f) {
									Vec3f intersectP = d * ray.direction();
									Vec3f lightPos = lightS.Position();
									Vec3f direc = lightPos - intersectP;
									Ray raytolight = Ray(intersectP + 0.05f * normalize(direc), normalize(direc));
									if (isvisible(raytolight, scene)) {
										image.set_pixel(x, y, shade(scene, meshIndex, triangleIndex, u, v, lightS));
										//color = shade(scene, meshIndex, triangleIndex, u, v, lightS);
									}
									else {
										image.set_pixel(x, y, Vec3f(0.f, 0.f, 0.f));
										//color = Vec3f(0.f, 0.f, 0.f);
									}
								}
							}
						}
					}
					//accum_color += color;
				//}
				//accum_color /= 3;
				//image.set_pixel(y, x, accum_color);		
			}
		}
	}

	void createCoordinateSystem(const Vec3f& N, Vec3f& Nt, Vec3f& Nb)
	{
		if (std::fabs(N[0]) > std::fabs(N[1]))
			Nt = Vec3f(N[2], 0, -N[0]) / sqrtf(N[0] * N[0] + N[2] * N[2]);
		else
			Nt = Vec3f(0, -N[2], N[1]) / sqrtf(N[1] * N[1] + N[2] * N[2]);
		Nb = cross(N, Nt);
	}

	Vec3f uniformSampleHemisphere(const float& r1, const float& r2)
	{
		float sinTheta = sqrtf(1 - r1 * r1);
		float phi = 2 * M_PI * r2;
		float x = sinTheta * cosf(phi);
		float z = sinTheta * sinf(phi);
		return Vec3f(x, 1, z);
	}

	//Vec3f pathtracing(int x, int y, int depth, Vec3f position, Vec3f N, Scene& scene, const Camera& camera, Image& image) {
	//	std::vector<bool> bb_interesect;
	//	Vec3f directLight = { 0.0f,0.0f,0.0f };
	//	Vec3f indirectLight = { 0.0f,0.0f,0.0f };
	//	Vec3f Normal = { 0.0f,0.0f,0.0f };
	//	Ray random_ray = Ray(position, N);
	//	bb_interesect.resize(scene.meshes().size());
	//	for (int l = 0; l < scene.lightSource().size(); l++) {
	//		LightSource lightS = scene.lightSource()[l];
	//		float xp = (float)(rand() % 250) / 250 - 0.10;
	//		float yp = (float)(rand() % 250) / 250 - 0.10;
	//		Ray ray = camera.rayAt(float(x + xp) / image.width(), 1.f - float(y + yp) / image.height());

	//		for (int m = 0; m < scene.meshes().size(); m++) {
	//			Mesh c_mesh = scene.meshes()[m];
	//			AABB c_AABB = c_mesh.get_AABB();
	//			bb_interesect[m] = c_AABB.ray_intersect(ray);
	//		}
	//		for (int m_meshes = 0; m_meshes < scene.meshes().size(); m_meshes++) {
	//			size_t meshIndex, triangleIndex;
	//			float u, v, d;
	//			if (bb_interesect[m_meshes]) {
	//				bool intersectionFound = rayTrace(ray, scene, meshIndex, triangleIndex, u, v, d);
	//				if (intersectionFound && d > 0.f) {
	//					Vec3f intersectP = d * ray.direction();
	//					Vec3f lightPos = lightS.Position();
	//					Vec3f direc = lightPos - intersectP;
	//					Ray raytolight = Ray(intersectP + 0.05f * normalize(direc), normalize(direc));
	//					if (isvisible(raytolight, scene)) {
	//						directLight = shade(scene, meshIndex, triangleIndex, u, v, lightS);
	//						const auto& mesh = scene.meshes()[meshIndex];
	//						const auto& P = mesh.vertexPositions();
	//						const auto& N = mesh.vertexNormals();
	//						const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
	//						Normal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
	//					}
	//					Vec3f Nt, Nb;
	//					createCoordinateSystem(Normal, Nt, Nb);
	//					float pdf = 1 / (2 * M_PI);
	//					for (uint32_t n = 0; n < 16; ++n) {
	//						float r1 = (float)(rand()%200)/200;
	//						float r2 = (float)(rand()%200)/200;
	//						Vec3f sample = uniformSampleHemisphere(r1, r2);
	//						Vec3f sampleWorld(
	//							sample[0] * Nb[0] + sample[1] * Normal[0] + sample[2] * Nt[0],
	//							sample[0] * Nb[1] + sample[1] * Normal[1] + sample[2] * Nt[1],
	//							sample[0] * Nb[2] + sample[1] * Normal[2] + sample[2] * Nt[2]);	
	//						depth += 1;
	//						scene.set_depth(depth);
	//						indirectLight += r1 * pathtracing(x, y, depth , intersectP + sampleWorld,	sampleWorld, scene, camera, image) / pdf;
	//					}	
	//					// divide by N
	//					indirectLight /= (float)16;
	//				}
	//			}
	//		}
	//	}
	//	if (scene.get_depth() > 20) return directLight + indirectLight;
	//	return directLight + indirectLight;
	//}

	//With global illumination - not working
	/*inline void render2(Scene& scene, Image& image) {
		size_t w = image.width();
		size_t h = image.height();
		Vec3f color;
		const Camera& camera = scene.camera();
		for (int y = 0; y < h; y++) {
			static int progress = 0;
			progress++;
			if (progress % 10 == 0)
				std::cout << ".";
#pragma omp parallel for
			for (int x = 0; x < w; x++) {
				color = pathtracing(x, y, scene.get_depth(), Vec3f{0.0f,0.0f,0.0f}, Vec3f{ 0.0f,0.0f,0.0f },scene, camera, image);
				image.set_pixel(x, y, color);
			}
		}
	}*/

};