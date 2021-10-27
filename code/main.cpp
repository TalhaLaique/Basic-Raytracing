#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <chrono>

#include "Image.h"
#include "Ray.h"
#include "Camera.h"
#include "Mesh.h"
#include "Scene.h"
#include "RayTracer.h"

using namespace std;

typedef Vec3i Triangle;

int main() {
	int width = 300;
	int height = 300;
	Image img(width, height);
	Image img1(width, height);
	for (size_t y = 0; y < height; y++) {
		for (size_t x = 0; x < width; x++) {
			img.set_pixel(y, x, Vec3f(0.f, 0.f, 1.f));
			//img1.set_pixel(y, x, Vec3f(0.f, 0.f, 1.f));
		}
	}

	Scene scene;

	Camera camera(Vec3f(0.f, 0.f, 1.2f),
		Vec3f(),
		Vec3f(0.f, 1.f, 0.f),
		60.f,
		float(width / height));
	scene.camera() = camera;

	LightSource lightSource1(Vec3f(2.0, 2.0, 2.0), Vec3f(1.0, 1.0, 1.0), 15.f);
	/*LightSource lightSource2(Vec3f(2.0, 1.5, 2.0), Vec3f(1.0, 1.0, 1.0), 15.0);
	LightSource lightSource3(Vec3f(3.0, 2.0, 2.0), Vec3f(1.0, 1.0, 1.0), 20.0);*/
	scene.lightSource().push_back(lightSource1);
	/*scene.lightSource().push_back(lightSource2);
	scene.lightSource().push_back(lightSource3);*/
	Mesh mesh, ground;
	mesh.loadOFF("lowresol.off");
	//ground.loadOFF("ground.off");
	auto& T = ground.indexedTriangles();
	auto& V = ground.vertexPositions();
	V.push_back(Vec3f(-5.0, -4.0, -5.0)); V.push_back(Vec3f(-5.0, -4.0, 5.0));
	V.push_back(Vec3f(5.0, -4.0, -5.0)); V.push_back(Vec3f(5.0, -4.0, 5.0));
	T.push_back(Triangle(0, 1, 3));
	T.push_back(Triangle(0, 3, 2));
	ground.recomputeNormals();
	scene.meshes().push_back(mesh);
	scene.meshes().push_back(ground);
	scene.set_AABB();


	//RayTracing
	RayTracer rayTracer;
	img.fillBackground();
	std::cout << "Ray tracing: starts";
	auto start = std::chrono::high_resolution_clock::now();
	rayTracer.render(scene, img);
	//rayTracer.render2(scene, img);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cout << "TimeLapsed:" << diff.count() << std::endl;
	std::cout << "ends." << std::endl;
	img.savePPM("testing_AABB.ppm");
	return 0;
};