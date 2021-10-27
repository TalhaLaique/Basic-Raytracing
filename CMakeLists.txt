cmake_minimum_required(VERSION 3.9)

project(MyRayTracer LANGUAGES CXX)

find_package(OpenMP REQUIRED)

add_executable (
	MyRayTracer
	AABB.h
	Bvh.h
 	Camera.h
	Image.h
	LightSource.h
	main.cpp
	Material.h
	Mesh.h
	Ray.cpp
	Ray.h
	RayTracer.h
	Scene.h
	Vec3.h
	WorleyNoise.h
)

set_target_properties(MyRayTracer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_link_libraries(MyRayTracer PRIVATE OpenMP::OpenMP_CXX)

