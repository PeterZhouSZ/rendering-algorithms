#include <dirt/parser.h>
#include <dirt/microfacet.h>
#include <dirt/surface.h>
#include <dirt/vec.h>
#include <dirt/onb.h>

Vec3f GGX::sample(const float alpha, const Vec3f &N) const
{
    // float x = randf();
    // float y = randf();
    // float phi = xi.y()*M_PI;
    // float tanThetaSq = alpha * alpha * xi.x()/(1.0f - xi.x());
    // float cosTheta = 1.0f / std::sqrt(1.0f + tanThetaSq);
        
    // float sinTheta = std::sqrt(max(1.0f - cosTheta*cosTheta, 0.0f));
    // return Vec3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    // onb uvw;
    // uvw.build_from_w(hit.sn);
    // Vec3f w_h = uvw.local(randomCosineHemisphere());
    // return w_h;

    float rng1 = randf();
    float rng2 = randf();

    // float theta_m = atanf(alpha*sqrtf(rng1) / sqrtf(1 - rng1));
    // float phi_m = 2 * M_PI * rng2;

    // Vec3f direction = Vec3f(sin(theta_m) * cos(phi_m), sin(theta_m) * sin(phi_m), cos(theta_m));

    float phi = rng1 * M_PI * 2;
    float cosTheta = 0.0f;
    float tanThetaSq = alpha * alpha * rng2/(1.0f - rng2);
    cosTheta = 1.0f/std::sqrt(1.0f + tanThetaSq);

    float r = std::sqrt(max(1.0f - cosTheta*cosTheta, 0.0f));
    // cout << alpha << endl;
    Vec3f direction = normalize(Vec3f(std::cos(phi)*r, std::sin(phi)*r, cosTheta));

	onb uvw;
    uvw.build_from_w(N);
    Vec3f t_normal = uvw.local(direction);
    return normalize(t_normal);
}

float GGX::D(const Vec3f& m, const HitInfo & hit, const float alpha) const
{
    
    float cos_theta = dot(normalize(m), hit.sn);

    if (dot(normalize(m), hit.sn) <= 0)
		return 0.0f;
    
    
    // float cos_theta4 = powf(cos_theta, 4.0f);
    float tan_theta = tan(acosf(cos_theta));
    
    float term = alpha * alpha + tan_theta * tan_theta;
    // float term2 = term * term;

    return alpha * alpha / (M_PI * powf(cos_theta, 4.0f) * term * term);
}

float GGX::pdf(const Vec3f &m, const HitInfo & hit, const float alpha) const
{
    // cout << "a:" << dot(hit.sn, m) << endl;
    return D(m, hit, alpha) * dot(hit.sn, m);
}

float GGX::G(const Vec3f& i, const Vec3f& o, const Vec3f& m, const HitInfo & hit, const float alpha) const
{
    return G1(i, m, hit, alpha) * G1(o, m, hit, alpha);
}

float GGX::G1(const Vec3f& v, const Vec3f& m, const HitInfo & hit, const float alpha) const
{
    float cos_theta = dot(normalize(v), hit.sn);
    float tan_theta = tan(acosf(cos_theta)); // angle between v, hit.sn

	// Perpendicular incidence 
	if (tan_theta == 0.0f)
		return 1.0f;

	// X+ term
	if (dot(v, m) * dot(v, hit.sn) <= 0)
		return 0.0f;

	return 2.0f / (1.0f + sqrtf(1.0f + alpha * alpha * tan_theta * tan_theta));
}