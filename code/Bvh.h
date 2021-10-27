#pragma once
#include"AABB.h"
class Bvh {
public:
	Bvh(){}
	virtual ~Bvh() {};
	/*
	
	TODO
	
	*/
private:
	Bvh* left_node;
	Bvh* right_node;
	AABB aabb;
};