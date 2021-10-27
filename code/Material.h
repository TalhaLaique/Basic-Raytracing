#pragma once
#include "Vec3.h"
#include "WorleyNoise.h"
#include <algorithm>

#define M_PI 3.1416
class Material {
public:
	Material(const float diffuse = 0.05f, const Vec3f& color = Vec3f(0.9f, 0.4f, 0.5f),
		const float alpha = 0.05f, const Vec3f F0 = Vec3f(0.31f, 0.31f, 0.31f), const int shine=1) :
		m_diffuse(diffuse), m_color(color), m_alpha(alpha), m_F0(F0), m_shine(shine) {};

	virtual ~Material() {};


	float Cook_Torrance_GeoMetric(float val, float shadowing, float masking) {
		if (val < shadowing && val < masking) return val;
		else if (shadowing < val && shadowing < masking) return shadowing;
		else return masking;
	}

	void CreateNoise(int pts, Vec3f min, Vec3f max) {
		m_wn_rough = WorleyNoise(pts, min, max);
		m_wm_albedo = WorleyNoise(pts, min, max);;
		m_wn_metal = WorleyNoise(pts, min, max);;
	}

	Vec3f evaluateColorResponse_MFM(const Vec3f& normal, const Vec3f& wi, const Vec3f& wo, Vec3f& pos) {
		/*
						Microfacet BRDF Model

		f(wi,wh) = GGX(wi,wo)F(wi,wh)G(wi,wo,wh)/4.(n.wi)(n.wo)
		wi: incoming light direction
		wo: Outgoing observer direction
		n: normal
		wh: average of wi and wo

		GGX: Normal Distribution: D(wi,wo)= alphasqr/M_PI(1+(alphasqr-1).(n.wh)^2)^2
		F: Fresenel Term (Schlick Appx): F(wi,wh)=F0 + (1-F0)(1-max(0,wi.wh))^5
		G: Geometric Shadowing: G(wi,wo,wh)=min(1, 2(n.wh)(n.wi)/(wo.wh), 2(n.wh)(n.wo)/(wo.wh))
		*/

		//Diffuse + Specular components
		Vec3f wh = normalize((wi + wo) / 2.f);
		float dot_wi_wh = dot(wi,wh);
		float dot_wo_wh = dot(wo, wh);
		float dot_n_wh = dot(normal,wh);
		float dot_n_wi = dot(normal, wi);
		if (dot_n_wi < 0) return Vec3f(0.f, 0.f, 0.f);
		float dot_n_wo = dot(normal,wo);

		//Deriving roughness, metallicness, and albedo from worley Noise
		CreateNoise(500, Vec3f{-2.0f,-2.0f,-2.0f}, Vec3f{2.0f,2.0f,2.0f});
		m_alpha = m_wn_rough.Noise(pos);
		float F0 = m_wn_metal.Noise(pos);
		float albedo = m_wm_albedo.Noise(pos);


		float GGX_Distrit = pow(m_alpha,2) / (M_PI * pow(1 + (pow(m_alpha, 2) - 1) * pow(dot_n_wh, 2), 2));
		//Vec3f Frensel = Vec3f(m_F0[0] + (1 - m_F0[0])*pow(1 - (dot_wi_wh > 0.f) ? dot_wi_wh : 0.f, 5), m_F0[1] + (1 - m_F0[1]) * pow(1 - (dot_wi_wh > 0.f) ? dot_wi_wh : 0.f, 5), m_F0[2] + (1 - m_F0[2]) * pow(1 - (dot_wi_wh > 0.f) ? dot_wi_wh : 0.f, 5));
		float Frensel = F0 + (1 - F0) * pow(1 - (dot_wi_wh > 0.0f) ? dot_wi_wh : 0.0f, 5);
		//float denom_term = 1 / (4 * dot_n_wi * dot_n_wo);

		float shadowing = (2 * (dot_n_wh) * (dot_n_wi)) / dot_wo_wh;
		float masking = (2 * (dot_n_wh) * (dot_n_wo)) / dot_wo_wh;
		float GeoMetric = Cook_Torrance_GeoMetric(1.f, shadowing, masking);

		//Final Response = Diffuse Term + Specular Term
		//Vec3f response = m_diffuse + GGX_Distrit * Frensel * GeoMetric * denom_term * dot_n_wi;
		//Vec3f response = m_diffuse + GGX_Distrit * Frensel * GeoMetric * dot_n_wi;
		float response = (m_diffuse + GGX_Distrit * Frensel * GeoMetric) * dot_n_wi;
		m_color[0] -= albedo;
		m_color[1] += albedo;
		m_color[2] -= albedo;

		return response * m_color;
	}

	Vec3f evaluateColorResponse_Diff(const Vec3f& normal, const Vec3f& wi) {
		//Diffuse Component
		float dot_n_wi = dot(normal, wi);
		//float response = (m_diffuse / M_PI) * dot_n_wi;
		float response = m_diffuse * (dot_n_wi>0)? dot_n_wi: 0.0f;
		return response * m_color;
	}

	//Vec3f evaluateColorResponse_Phong(const Vec3f& normal, const Vec3f& wi, const Vec3f& wo) {
	//	//Glossy Component or Specular reflection
	//	float dot_n_wi = dot(normal, wi);
	//	Vec3f r = 2.f * normal * dot(wi, normal) - wi;
	//	float response = m_F0 * pow(dot(r, wo), m_shine) * dot_n_wi;
	//	return response * m_color;
	//}

	//Vec3f evaluateColorResponse_Blinn_Phong(const Vec3f& normal, const Vec3f& wi, const Vec3f& wo) {
	//	//Glossy Component or Specular reflection
	//	float dot_n_wi = dot(normal, wi);
	//	Vec3f wh = normalize((wi + wo) / 2);
	//	float response = m_F0 * pow(dot(normal, wh), m_shine) * dot_n_wi;
	//	return response * m_color;

	//}

private:
	Vec3f m_color;
	float m_diffuse; //For diffusing the albedo
	float m_alpha; //For roughness
	Vec3f m_F0; //For Reflectance (Diaelectric, Metal, Insulator)
	float m_shine; //shininess
	WorleyNoise m_wn_rough;
	WorleyNoise m_wm_albedo;
	WorleyNoise m_wn_metal;
};