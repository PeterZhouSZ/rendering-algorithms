/*
    This file is part of Dirt, the Dartmouth introductory ray tracer.

    Copyright (c) 2017-2019 by Wojciech Jarosz

    Dirt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Dirt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <dirt/material.h>
#include <dirt/parser.h>
#include <dirt/scene.h>
#include <dirt/surface.h>
#include <dirt/onb.h>
#include <dirt/microfacet.h>


bool refract(const Vec3f &in, const Vec3f &normal, Vec3f &out, const float ratio)
{	
	// air into material
	
	Vec3f tmp_normal = normalize(normal);
	Vec3f tmp_in = normalize(in);
	float angle_in = dot(tmp_in, tmp_normal);
	float k = 1.0f - pow(ratio, 2) * (1 - angle_in * angle_in);
	if(k > 0) {
		out = ratio * (tmp_in - angle_in*tmp_normal) - tmp_normal * sqrt(k);
		return true;
	}
	else{
		return false;
	}
	
}

// Vec3f bump_N(const HitInfo &hit, shared_ptr<const Texture> BumpMap, float offset)
// {
// 	Color3f shift = BumpMap->value(hit);

// 	onb uvw;
//     uvw.build_from_w(hit.sn);
//     Vec3f t_normal = uvw.local(shift * offset);
    
//     return normalize(t_normal + hit.sn);

// 	// return normalize(shift * offset + hit.sn);
// }


// Vec3f BumpMap::bump_N(const HitInfo &hit, BumpMap img)
// {
// 	if(!img._bump) 
// 		return hit.sn;
	
// 	int i = hit.uv.x * img.content.width();
//     int j = (1 - hit.uv.y) * img.content.height() - 0.001;

// 	int i1 = i - 1;
// 	int i2 = i + 1;
// 	int j1 = j - 1;
// 	int j2 = j + 1;

// 	if (i1 < 0) i1 = 0;
// 	if (i1 > img.content.width()-1) i1 = img.content.width() - 1;
// 	if (i2 < 0) i2 = 0;
// 	if (i2 > img.content.width()-1) i2 = img.content.width() - 1;

// 	if (j1 < 0) j1 = 0;
// 	if (j1 > img.content.height()-1) j1 = img.content.height() - 1;
// 	if (j2 < 0) j2 = 0;
// 	if (j2 > img.content.height()-1) j2 = img.content.height() - 1;

//     // clamp
//     if (i < 0) i = 0;
//     if (j < 0) j = 0;
//     if (i > img.content.width() - 1) i = img.content.width() - 1;
//     if (j > img.content.height() - 1) j = img.content.height() - 1;

// 	float u_gradient = luminance(img.content(i1, j)) - luminance(img.content(i2, j));
// 	float v_gradient = luminance(img.content(i, j1)) - luminance(img.content(i, j2));

// 	onb uvw;
// 	uvw.build_from_w(hit.sn);
// 	Vec3f tangent = normalize(uvw.local(Vec3f(0, 0, 1)) - Vec3f(u_gradient, v_gradient, 0));
// 	return hit.sn + 20 * u_gradient * tangent + 20 * v_gradient * (-hit.sn);
// }


float schlick(float cos, float ior)
{
	float r0 = (1 - ior) / (1 + ior);
	r0 = r0*r0;
	return r0 + (1 - r0) * pow((1 - cos) , 5);
}
Vec3f reflect(Vec3f dir, Vec3f normal)
{
	Vec3f tmp_d = normalize(dir);
	Vec3f tmp_sn = normalize(normal);
	
	Vec3f reflected = tmp_d - 2*dot(tmp_d, tmp_sn)*tmp_sn;
	return reflected;
}

float fresnel(float eta, float cosI, float &cosT)
{
	// eta = eta_i / eta_t
	if(eta == 1.0f) return 0;

	if (cosI < 0.0f) {
        eta = 1.0f / eta;
        cosI *= -1.0f;
    }
    float sinSq = eta * eta * (1.0f - cosI * cosI);
    if (sinSq > 1.0f) {
        cosT = 0.0f;
        return 1.0f;
    }
    cosT = sqrt(max(1.0f - sinSq, 0.0f));

    float Rs = (eta * cosI - cosT) / (eta * cosI + cosT);
    float Rp = (eta * cosT - cosI) / (eta * cosT + cosI);

    return (Rs * Rs + Rp * Rp) / 2.0;
	
}

float sign(const Vec3f a, const Vec3f b)
{
	return dot(a, b) < 0.0f ? -1.0f : 1.0f;
}


namespace
{

auto g_defaultMaterial = make_shared<Lambertian>(json{ {"albedo", 0.8} });

} // namespace


shared_ptr<const Material> Material::defaultMaterial()
{
	return g_defaultMaterial;
}




Lambertian::Lambertian(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	stored_photons = true;
	isDelta = false;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
}

bool Lambertian::scatter(const Ray3f &ray, const HitInfo &hit, Color3f &attenuation, Ray3f &scattered) const
{
	// TODO: Implement Lambertian reflection
	//       You should assign the albedo to ``attenuation'', and
	//       you should assign the scattered ray to ``s cattered''
	//       The origin of the scattered ray should be at the hit point,
	//       and the scattered direction is the shading normal plus a random
	//       point on a sphere (please look at the text book for this)

	//       You can get the hit point using hit.p, and the shading normal using hit.sn

	//       Hint: You can use the function randomInUnitSphere() to get a random
	//       point in a sphere. IMPORTANT: You want to add a random point *on*
	//       a sphere, not *in* the sphere (the text book gets this wrong)
	//       If you normalize the point, you can force it to be on the sphere always, so
	//       add normalize(randomInUnitSphere()) to your shading normal
	
	
	// construct unit sphere, center = hit.p + hit.sn
	
	// Vec3f scatteredPoint =  hit.sn + hit.p + normalize(randomInUnitSphere());
	
	scattered = Ray3f(hit.p, hit.sn + normalize(randomInUnitSphere()));
	attenuation = albedo->value(hit);
	return true;
	
}

bool Lambertian::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	srec.isSpecular = false;
	// srec.scattered = randomCosineHemisphere();
	// srec.attenuation = albedo->value(hit);
	// return true;

	Vec3f tmp_N = hit.sn;
	if(albedo->_bump)
		tmp_N = albedo->bump_N(hit);

	onb uvw;
    uvw.build_from_w(tmp_N);
    Vec3f direction = uvw.local(randomCosineHemisphere());
    srec.scattered = normalize(direction);
    srec.attenuation = albedo->value(hit);

	if (dot(srec.scattered, tmp_N) < 0.0f) return false;

    
    return true;
	// bool result = scatter(Ray3f(Vec3f(randf(), randf(), randf()), dirIn), hit, srec.attenuation, srec.scattered);
}

Color3f Lambertian::get_albedo(const HitInfo & hit) const
{
	return albedo->value(hit);
}

Color3f Lambertian::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	// cout << max(0.0f, dot(scattered, hit.sn)) << endl;
	Vec3f tmp_N = hit.sn;
	if(albedo->_bump)
		tmp_N = albedo->bump_N(hit);
	return albedo->value(hit) * max(0.0f, dot(scattered, tmp_N))/ M_PI;
}

float Lambertian::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	Vec3f tmp_N = hit.sn;
	if(albedo->_bump)
		tmp_N = albedo->bump_N(hit);
	return max(0.0f, dot(scattered, tmp_N))/M_PI;
}


Metal::Metal(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	roughness = parseTexture(j.at("roughness"));
	stored_photons = false;
	isDelta = true;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
	// roughness = clamp(j.value("roughness", roughness), 0.f, 1.f);
	
}

bool Metal::scatter(const Ray3f &ray, const HitInfo &hit, Color3f &attenuation, Ray3f &scattered) const
{
	// TODO: Implement metal reflection
	//       This function proceeds similar to the lambertian material, except that the
	//       scattered direction is different.
	//       Instead of adding a point on a sphere to the normal as before, you should add the point
	//       to the *reflected ray direction*.
	//       You can reflect a vector by the normal using reflect(vector, hit.sn); make sure the vector is normalized.
	//       Different to before you can't just use randomInUnitSphere directly; the sphere should be scaled by roughness.
	//       (see text book). In other words, if roughness is 0, the scattered direction should just be the reflected direction.
	//       
	//       This procedure could produce directions below the surface. Handle this by returning false if the scattered direction and the shading normal
	//       point in different directions (i.e. their dot product is negative)
	
	
	Vec3f reflected = reflect(ray.d, hit.sn);
	
	
	scattered = Ray3f(hit.p, normalize(normalize(reflected) + roughness->value(hit) * (randomInUnitSphere())));
	attenuation = albedo->value(hit);
	return (dot(scattered.d, hit.sn) > 0);
}

bool Metal::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	Vec3f reflected = reflect(dirIn, hit.sn);
	srec.isSpecular = true;
	srec.attenuation = albedo->value(hit);
	srec.scattered = normalize(normalize(reflected) + roughness->value(hit) * (randomInUnitSphere()));

	return (dot(srec.scattered, hit.sn) > 0);
}

Color3f Metal::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	return albedo->value(hit);
}

float Metal::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	return -1.0f;
}


Dielectric::Dielectric(const json & j)
{
	ior = j.value("ior", ior); // index of refraction
	stored_photons = false;
	isDelta = true;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
}

bool Dielectric::scatter(const Ray3f &ray, const HitInfo &hit, Color3f &attenuation, Ray3f &scattered) const
{
	// TODO: Implement dielectric scattering
	
	Vec3f tmp_d = normalize(ray.d);
	Vec3f reflected = reflect(tmp_d, hit.sn);
	Vec3f refracted;
	attenuation = Vec3f(1.0, 1.0, 1.0);
	Vec3f tmp_normal;
	float reflect_prob;
	float cos; // angle for refraction
	float ratio;
	
	if (dot(tmp_d, hit.sn) > 0){
		tmp_normal = -hit.sn;
		cos = ior * dot(tmp_d, normalize(hit.sn)) / length(tmp_d);
		ratio = ior;
	}
	else {
		cos = dot(-tmp_d, normalize(hit.sn)) / length(tmp_d);
		ratio = 1 / ior;
		tmp_normal = hit.sn;
	}
	if (refract(tmp_d, tmp_normal, refracted, ratio)) {
		reflect_prob = schlick(cos, ior);
	}
	else reflect_prob = 1.0f;
	if (drand48() < reflect_prob)
		scattered = Ray3f(hit.p, reflected);
	else scattered = Ray3f(hit.p, refracted);
	
	
	return true;
	
	
}

bool Dielectric::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	srec.isSpecular = true;
	srec.attenuation = Color3f(1.0f);

	Vec3f tmp_d = normalize(dirIn);
	Vec3f reflected = reflect(tmp_d, hit.sn);
	Vec3f refracted;
	
	Vec3f tmp_normal;
	float reflect_prob;
	float cos; // angle for refraction
	float ratio;
	
	
	if (dot(tmp_d, hit.sn) > 0){
		tmp_normal = -hit.sn;
		cos = ior * dot(tmp_d, normalize(hit.sn)) / length(tmp_d);
		ratio = ior;
	}
	else {
		cos = dot(-tmp_d, normalize(hit.sn)) / length(tmp_d);
		ratio = 1 / ior;
		tmp_normal = hit.sn;
	}
	if (refract(tmp_d, tmp_normal, refracted, ratio)) {
		reflect_prob = schlick(cos, ior);
	}
	else reflect_prob = 1.0f;
	if (drand48() < reflect_prob){
		// cout << "wtf1" << endl;
		srec.scattered = reflected;
	}
		
	else {
		// cout << "wtf2" << endl;
		srec.scattered = refracted;
	}

	

	return true;
}

Color3f Dielectric::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	cout << "fuck" << endl;
}

float Dielectric::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	return -2.0f;
}

DiffuseLight::DiffuseLight(const json & j)
{
	emit = j.value("emit", emit);
}

Color3f DiffuseLight::emitted(const Ray3f &ray, const HitInfo &hit) const
{
	// only emit from the normal-facing side
	if (dot(ray.d, hit.sn) > 0)
		return Color3f(0,0,0);
	else
		return emit;
}

BlendMaterial::BlendMaterial(const json & j )
{
	mat_A = parseMaterial(j["a"]);
	mat_B = parseMaterial(j["b"]);
	blend_ratio = parseTexture(j.at("amount"));
	stored_photons = true;
	
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
}

bool BlendMaterial::scatter(const Ray3f &ray, const HitInfo &hit, Color3f &attenuation, Ray3f &scattered) const
{
	double rand_num = drand48();
	if (rand_num <= luminance(blend_ratio->value(hit))){
		return mat_B->scatter(ray, hit,attenuation, scattered);
	}
	else{
		return mat_A->scatter(ray, hit,attenuation, scattered);
	}
}

bool BlendMaterial::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	srec.isSpecular = false;
}

Color3f BlendMaterial::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{

}

float BlendMaterial::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	return -1.0f;
}

Phong::Phong(const json & j )
{
	albedo = parseTexture(j.at("albedo"));
	
	exponent = (j.value("exponent", exponent));
	// if(exponent >= 500)
		stored_photons = false;
	// else stored_photons = true;
	isDelta = false;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
}



Color3f Phong::get_albedo(const HitInfo & hit) const
{
	return albedo->value(hit);
}

bool Phong::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	srec.isSpecular = true;
	
	Vec3f direction = randomCosinePowerHemisphere(exponent);

	onb uvw;
    uvw.build_from_w(reflect(dirIn, hit.sn));
    direction = uvw.local(direction);
    srec.scattered = normalize(direction);
    srec.attenuation = albedo->value(hit);

	if (dot(srec.scattered, hit.sn) < 0.0f) return false;

	return true;
}

Color3f Phong::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	// if (exponent > 2000)
	// 	cout << pdf(dirIn, scattered, hit) << endl;
	return albedo->value(hit) * pdf(dirIn, scattered, hit);
}

float Phong::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	Vec3f mirrorDir = normalize(reflect(dirIn, hit.sn));
	float cosine = max(dot(normalize(scattered), mirrorDir), 0.0f);
	float constant = (exponent + 1)/(2*M_PI);
	return constant * pow(cosine, exponent);
}

BlinnPhong::BlinnPhong(const json & j )
{
	albedo = parseTexture(j.at("albedo"));
	
	exponent = (j.value("exponent", exponent));
	if(exponent >= 500)
		stored_photons = false;
	else stored_photons = true;
	isDelta = false;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
}

Color3f BlinnPhong::get_albedo(const HitInfo & hit) const
{
	return albedo->value(hit);
}

bool BlinnPhong::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	srec.isSpecular = false;
	
	Vec3f direction = randomCosinePowerHemisphere(exponent);

	onb uvw;
    uvw.build_from_w(hit.sn);
    Vec3f t_normal = uvw.local(direction);
    srec.scattered = normalize(reflect(dirIn, t_normal));
    srec.attenuation = albedo->value(hit);

	if (dot(srec.scattered, hit.sn) < 0.0f) return false;

	return true;
}

Color3f BlinnPhong::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	return albedo->value(hit) * pdf(dirIn, scattered, hit);
}

float BlinnPhong::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	if (dot(dirIn, hit.sn) > 0.0f || dot(scattered, hit.sn) < 0.0f){
		
		return 0.0f;
	}
		
	Vec3f randomNormal = normalize(-normalize(dirIn) + normalize(scattered));
	float cosine = max(dot(randomNormal, hit.sn), 0.0f);
	float normalPdf = (exponent + 1) / (2*M_PI) * powf(cosine, exponent);
	float finalPDF = normalPdf/(4 * dot(-dirIn, randomNormal));
	// if (finalPDF == 0)
	// 	cout << "lul:" << normalPdf << endl;

	return finalPDF;
}



RoughDielectric::RoughDielectric(const json & j)
{
	ior = j.value("ior", ior); // index of refraction
	alpha = j.value("alpha", 0.5); //default for glass

	stored_photons = false;
	isDelta = false;
	// _bump = j.value("bump", false);
	// if(_bump) BumpMap = parseTexture(j.at("bumpmap"));
	
	micro_type = j.value("microfacet", "GGX");
	if (micro_type == "GGX")
		distribution = new GGX();
}

bool RoughDielectric::sample(const Vec3f & dirIn, const HitInfo &hit, ScatterRecord &srec) const
{
	
	srec.isSpecular = true;
	srec.attenuation = Color3f(1.0f);

	Vec3f w_h = distribution->sample(alpha, hit.sn);

	Vec3f tmp_d = normalize(dirIn);
	Vec3f reflected = reflect(tmp_d, w_h);
	Vec3f refracted;
	
	Vec3f tmp_normal;
	float reflect_prob;
	float cos; // angle for refraction
	float ratio;
	
	
	if (dot(tmp_d, w_h) > 0){
		tmp_normal = -w_h;
		cos = ior * dot(tmp_d, normalize(w_h)) / length(tmp_d);
		ratio = ior;
	}
	else {
		cos = dot(-tmp_d, normalize(w_h)) / length(tmp_d);
		ratio = 1 / ior;
		tmp_normal = w_h;
	}
	if (refract(tmp_d, tmp_normal, refracted, ratio)) {
		reflect_prob = schlick(cos, ior);
	}
	else reflect_prob = 1.0f;
	if (drand48() < reflect_prob){
		
		srec.scattered = reflected;
		_refr = false;
	}
		
	else {
		
		srec.scattered = refracted;
		_refr = true;
	}

	return true;
}



Color3f RoughDielectric::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	

    bool reflect = dot(-dirIn, hit.sn) * dot(hit.sn, scattered) >= 0.0f;

    float eta = dot(-dirIn, hit.sn) < 0.0f ? ior : 1.0f / ior;

	if(reflect){
		Vec3f m = normalize(sign(-dirIn, hit.sn)*(-dirIn + scattered));
		float cos;

		float F_r = fresnel(1.0f / ior, dot(-dirIn, m), cos);
		float G_r = distribution->G(-dirIn, scattered, m, hit, alpha);
		float D_r = distribution->D(m, hit, alpha);
		float fr = (F_r * G_r * D_r * 0.25f)/abs(dot(-dirIn, hit.sn));
        return Color3f(fr);
	}
	else{
		Vec3f m = -normalize(-dirIn * eta + scattered);
		float cos;

		float F_t = fresnel(1.0f/ior, dot(-dirIn, m), cos);
		float G_t = distribution->G(-dirIn, scattered, m, hit, alpha);
		float D_t = distribution->D(m, hit, alpha);
		float denom = (eta * dot(-dirIn, m) + dot(scattered, m)) * (eta * dot(-dirIn, m) + dot(scattered, m)) * abs(dot(-dirIn, hit.sn));

		float fs = abs(dot(-dirIn, m) * dot(scattered, m)) * (1.0f - F_t) * G_t * D_t / denom;
        return Color3f(fs);
	}
	
}

float RoughDielectric::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	if(alpha == 0.0f) return -1;

	bool reflect = dot(-dirIn, hit.sn) * dot(hit.sn, scattered) >= 0.0f;
    

	float eta = dot(-dirIn, hit.sn) < 0.0f ? ior : 1.0f/ior;
    


	if(reflect){
		Vec3f m = normalize(sign(-dirIn, hit.sn)*(-dirIn + scattered));
		float cos;

		float F = fresnel(1.0f / ior, dot(-dirIn, m), cos);
		
		float pm = distribution->pdf(m, hit, alpha);
		float pdf = F * pm * 0.25f / abs(dot(-dirIn, m));
		return pdf;
	}
	else{
		Vec3f m = -normalize(-dirIn * eta + scattered);
		float cos;

		float F = fresnel(1.0f / ior, dot(-dirIn, m), cos);
		
		float pm = distribution->pdf(m, hit, alpha);
		float denom = eta * dot(-dirIn, m) + dot(scattered, m);
		float pdf = (1.0 - F) * pm * fabsf(abs(dot(scattered, m))) / (denom * denom);
		return pdf;
	}
	
}