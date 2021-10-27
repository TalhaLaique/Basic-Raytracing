#pragma once

#include <vector>

#include "Camera.h"
#include "Mesh.h"
#include "LightSource.h"
#include "Material.h"

class Scene {
public:
	inline Scene() { m_depth = 0; }

	inline const Camera& camera() const { return m_camera; }

	inline Camera& camera() { return m_camera; }
	
	inline const std::vector<Mesh> & meshes () const { return m_meshes;  }
	
	inline std::vector<Mesh> & meshes () { return m_meshes; }

	//inline const LightSource& lightSource() const { return m_lightsource; }

	inline std::vector <LightSource> & lightSource() { return m_lightsource; }

	inline const Material& material() const { return m_material; }

	inline Material& material() { return m_material; }

	inline int get_depth() { return m_depth; }

	inline void set_depth(int d) { m_depth = d; }

	void set_AABB() {
		for (int i = 0; i < m_meshes.size(); i++) {
			m_meshes[i].AABB_Mesh();
		}
	}

private:
	Camera m_camera;
	std::vector<Mesh> m_meshes;
	std::vector <LightSource> m_lightsource;
	Material m_material;
	int m_depth;
};