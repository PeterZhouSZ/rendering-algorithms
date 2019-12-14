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
		out = ratio * (tmp_in - angle_in*tmp_normal) - tmp_normal*sqrt(k);
		return true;
	}
	else{
		return false;
	}
	
}

float sqr(float a)
{
	return a*a;
}

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

float fresnel(float eta, float cosThetaI, float &cosThetaT)
{
	// if(ior == 1.0f) return 0;

	
    // float cosi = dot(I, N); 
    // float etai = 1, etat = ior; 
    // if (cosi > 0) 
	// 	std::swap(etai, etat);  
    
    // float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); 

    // // Total internal reflection
    // if (sint >= 1) { 
    //     return 1.0f; 
    // } 
    // else { 
    //     float cost = sqrtf(std::max(0.f, 1 - sint * sint)); 
    //     cosi = fabsf(cosi); 
    //     float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
    //     float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
    //     return (Rs * Rs + Rp * Rp) / 2; 
    // }

	if(eta == 1.0f) return 0;

	if (cosThetaI < 0.0f) {
        eta = 1.0f/eta;
        cosThetaI = -cosThetaI;
    }
    float sinThetaTSq = eta*eta*(1.0f - cosThetaI*cosThetaI);
    if (sinThetaTSq > 1.0f) {
        cosThetaT = 0.0f;
        return 1.0f;
    }
    cosThetaT = sqrt(max(1.0f - sinThetaTSq, 0.0f));

    float Rs = (eta*cosThetaI - cosThetaT)/(eta*cosThetaI + cosThetaT);
    float Rp = (eta*cosThetaT - cosThetaI)/(eta*cosThetaT + cosThetaI);

    return (Rs*Rs + Rp*Rp)*0.5f;
	
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

// Vec3f Material::bump_N(const HitInfo &hit) const
// {
// 	if(!_bump) 
// 		return hit.sn;
	
// 	Color3f shift = BumpMap->value(hit);
// 	return normalize(shift * 3 + hit.sn);

// }


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
	// if(_bump)
	// 	tmp_N = bump_N(hit);

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
	// if(_bump)
	// 	tmp_N = bump_N(hit);
	return albedo->value(hit) * max(0.0f, dot(scattered, tmp_N))/ M_PI;
}

float Lambertian::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	Vec3f tmp_N = hit.sn;
	// if(_bump)
	// 	tmp_N = bump_N(hit);
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
	// Vec3f w_h = distribution->sample(alpha, hit.sn); //microface normal

	// cout << "w_h1:" << w_h << endl;
	// // w_h = hit.sn;

	// cout << "angle:" << dot(w_h, hit.sn) << endl;
	// // cout << "?" << endl;
	srec.isSpecular = true;
	srec.attenuation = Color3f(1.0f);

	// Vec3f tmp_normal;

	// if (dot(dirIn, hit.sn) > 0){
	// 	tmp_normal = -hit.sn;
	// 	cos = ior * dot(dirIn, normalize(hit.sn)) / length(dirIn);
	// 	ratio = ior;
	// }
	// else {
	// 	cos = dot(-dirIn, normalize(hit.sn)) / length(dirIn);
	// 	ratio = 1 / ior;
	// 	tmp_normal = hit.sn;
	// }
	// if (refract(dirIn, tmp_normal, refracted, ratio)) {
	// 	reflect_prob = schlick(cos, ior);
	// }
	// else reflect_prob = 1.0f;
	
    // // Vec3f w_h = hit.sn;
	// float F = fresnel(dirIn, w_h, ior);
	// // cout << F << endl;

	// if(randf() < F){
	// 	_refr = false;
	// 	srec.scattered = reflect(normalize(dirIn), w_h);
	// 	// cout << "wtf" << endl;
	// 	return true;
	// }
	// else {
		
	// 	_refr = true;
	// 	float cosi = dot(dirIn, hit.sn); 
	// 	float etai = 1, etat = ior;

	// 	if (cosi > 0.0f)
	// 		std::swap(etai, etat);

	// 	float eta = etai / etat;
	// 	float c = dot(-dirIn, w_h);

		
	// 	float discriminant = 1.0f + eta * (c * c - 1.0f);
	// 	if (discriminant < 0.0f)
	// 		return false;

	// 	srec.scattered = (eta * c - sign(-dirIn, hit.sn) * sqrtf(discriminant)) * w_h - eta * -dirIn;
	// 	// cout << "wtf2" << endl;

	// 	return true;
	// }

	// float Fr = fresnel(dirIn, w_h, ior);
	// // cout << Fr << endl;
	// if(randf() < Fr){
		
	// 	_refr = false;
	// 	srec.scattered = reflect(normalize(dirIn), w_h);

	// 	// cout << dot(normalize(-dirIn + srec.scattered), w_h) << endl;
	// 	return true;
	// }
	// else{
	// 	// cout << "???" << endl;
	// 	_refr = true;
	// 	float cos = dot(-dirIn, hit.sn);

	// 	// bool entering = Frame::cosTheta(bRec.wi) > 0.0f;dot(-dirIn, hit.sn)
	// 	float eta_i = 1;
	// 	float eta_t = ior;

	// 	if(cos < 0){
	// 		std::swap(eta_i, eta_t);
	// 	}


	// 	float eta = eta_i / eta_t;
	// 	float c = dot(w_h, -dirIn);

	// 	// float sign = bRec.wi.z() > 0.0f ? 1.0f : -1.0f;
	// 	float sign = cos > 0.0f ? 1.0f : -1.0f;
	// 	float discriminant = 1.0f + eta * (c * c - 1.0f);
	// 	if (discriminant < 0.0f)
	// 	{
	// 		return false;
	// 	}

	// 	srec.scattered = (eta * c - sign * sqrtf(1.0f + eta * (c * c - 1.0f)))*w_h - eta * -dirIn;
	// 	return true;
	// }


	// new

	// float wiDotN = event.wi.z();

    // float eta = wiDotNdot < 0.0f ? ior : 1.0f/ior;

	// float eta = dot(-dirIn, hit.sn) < 0.0f ? ior : 1.0f/ior;

    // // float sampleRoughness = (1.2f - 0.2f*std::sqrt(std::abs(wiDotN)))*roughness;
    // // float alpha = Microfacet::roughnessToAlpha(distribution, roughness);
    // // float sampleAlpha = Microfacet::roughnessToAlpha(distribution, sampleRoughness);

    // // Vec3f m = Microfacet::sample(distribution, sampleAlpha, event.sampler->next2D());
	// Vec3f m = distribution->sample(alpha, hit.sn);
	// // m = hit.sn;
    // // float pm = Microfacet::pdf(distribution, sampleAlpha, m);
	// float pm = distribution->pdf(m, hit, alpha);


    // if (pm < 1e-10f)
    //     return false;

    // // float wiDotM = event.wi.dot(m);
    // float cosThetaT = 0.0f;
    // float F = fresnel(1.0f/ior, dot(-dirIn, m), cosThetaT);
    // // float etaM = wiDotM < 0.0f ? ior : 1.0f/ior;

	// float etaM = dot(-dirIn, m) < 0.0f ? ior : 1.0f/ior;

    // bool _reflect = randf() < F ? true : false;

	// if (_reflect)
    //     srec.scattered = reflect(dirIn, m);
    // else srec.scattered = (etaM * dot(-dirIn, m) - sign(-dirIn, m) * cosThetaT) * m - etaM * -dirIn;
        
	// // float woDotN = event.wo.z();

    // bool reflected = dot(-dirIn, hit.sn) * dot(hit.sn, srec.scattered) > 0.0f;
    // if (reflected != _reflect)
    //     return false;

	// return true;

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
		// cout << "wtf1" << endl;
		srec.scattered = reflected;
	}
		
	else {
		// cout << "wtf2" << endl;
		srec.scattered = refracted;
	}

	

	return true;
}



Color3f RoughDielectric::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	// Color3f fr = Color3f(0.0f);
	// Color3f ft = Color3f(0.0f);

	// //reflection part
	// if (!_refr){
	// 	Vec3f w_h = sign(-dirIn, hit.sn) * normalize(-dirIn + scattered);

	// 	float D_reflect = distribution->D(w_h, hit, alpha);
	// 	float F = fresnel(dirIn, hit.sn, ior);
	// 	float G_reflect = distribution->G(-dirIn, scattered, w_h, hit, alpha);

	// 	fr = Color3f(F * D_reflect * G_reflect / (4 * abs( dot(hit.sn, scattered))*dot(hit.sn, -dirIn) ));
	// 	if ((4 * abs( dot(hit.sn, scattered))*dot(hit.sn, -dirIn) ) == 0.0f) 
	// 		fr = Color3f(0.0f);	
	// 	return fr;
	// }
	// else {
	// 	// transmission part
	// 	float cosi = dot(dirIn, hit.sn);
	// 	float eta_i = 1, eta_t = ior; 
	// 	if (cosi > 0.0f)
	// 		std::swap(eta_i, eta_t);

	// 	Vec3f w_ht = -normalize(eta_i * -dirIn + eta_t * scattered);

	// 	float term = fabsf(dot(-dirIn, w_ht)) * fabsf(dot(scattered, w_ht)) / (fabsf(dot(-dirIn, hit.sn)) * fabsf(dot(scattered, hit.sn)));
	// 	float fr_t = 1.0f - fresnel(dirIn, w_ht, ior);
	// 	float denom = eta_i * dot(-dirIn, w_ht) + eta_t * dot(w_ht, scattered);

	// 	float G_t = distribution->G(-dirIn, scattered, w_ht, hit, alpha);
	// 	float D_t = distribution->D(w_ht, hit, alpha);

	// 	ft = Color3f(term * eta_t * eta_t * fr_t * G_t * D_t / (denom * denom));
	// 	return ft;
	// }
	
	// Color3f fr = Color3f(0.0f);
	// Color3f ft = Color3f(0.0f);

	// //reflection part
	// Vec3f w_h = sign(-dirIn, hit.sn) * normalize(-dirIn + scattered);

	// float D_reflect = distribution->D(w_h, hit, alpha);
	// float F = fresnel(dirIn, hit.sn, ior);
	// float G_reflect = distribution->G(-dirIn, scattered, w_h, hit, alpha);

	// fr = Color3f(F * D_reflect * G_reflect / (4 * abs( dot(hit.sn, scattered)*dot(hit.sn, -dirIn)) ));
	// if ((4 * abs( dot(hit.sn, scattered))*dot(hit.sn, -dirIn) ) == 0.0f) 
	// 	fr = Color3f(0.0f);	

	// // transmission part
	// float cosi = dot(dirIn, hit.sn);
	// float eta_i = 1, eta_t = ior; 
	// if (cosi > 0.0f)
	// 	std::swap(eta_i, eta_t);

	// Vec3f w_ht = -normalize(eta_i * -dirIn + eta_t * scattered);

	// float term = fabsf(dot(-dirIn, w_ht)) * fabsf(dot(scattered, w_ht)) / (fabsf(dot(-dirIn, hit.sn)) * fabsf(dot(scattered, hit.sn)));
	// float fr_t = 1.0f - fresnel(dirIn, w_ht, ior);
	// float denom = eta_i * dot(-dirIn, w_ht) + eta_t * dot(w_ht, scattered);

	// float G_t = distribution->G(-dirIn, scattered, w_ht, hit, alpha);
	// float D_t = distribution->D(w_ht, hit, alpha);

	// ft = Color3f(term * eta_t * eta_t * fr_t * G_t * D_t / (denom * denom));
	// // cout << fr+ft << endl;
	// return fr + ft;


	
	// Color3f fr, ft;
				
	// float sign = dot(-dirIn, hit.sn) < 0.0f ? -1.0f : 1.0f;
	// Vec3f w_h = sign * normalize(-dirIn + scattered);
	// Vec3f half_vector = -dirIn + scattered;
	// if (length(half_vector) < 1e-3f)
	// {
	// 	fr = Color3f(0.0f);
	// }

	// // if (w_h.isZero()) fr = 0.0f;
	// float D_r = distribution->D(w_h, hit, alpha);
	// float F_r = fresnel(dirIn, w_h, ior);;
	// float G_r = distribution->G(-dirIn, scattered, w_h, hit, alpha);
	
	// fr = Color3f(F_r * D_r * G_r / (4 * fabsf(dot(-dirIn, hit.sn) * dot(scattered, hit.sn))));
	// if (std::isnan(fr.x) || std::isnan(fr.y) || std::isnan(fr.z)) fr = Color3f(0.0f);		// might arise due 0/0


	


	// // Transmission
	// float eta_i = 1, eta_t = ior;
	// if (dot(-dirIn, hit.sn) < 0.0f)
	// 	std::swap(eta_i, eta_t);

	// Vec3f w_ht = -normalize(eta_i * -dirIn + eta_t * scattered);

	// // float term = fabsf(bRec.wi.dot(w_ht)) * fabsf(bRec.wo.dot(w_ht)) / (fabsf(bRec.wi.z()) * fabsf(bRec.wo.z()));
	// float term = fabsf(dot(-dirIn, w_ht)) * fabsf(dot(scattered, w_ht)) / fabsf(dot(-dirIn, hit.sn)) / fabsf(dot(scattered, hit.sn));
	// // float fr_t = 1.0f - fresnel(w_ht.dot(bRec.wi), m_extIOR, m_intIOR);
	// float fr_t = 1.0f - fresnel(dirIn, w_ht, ior);
	// // float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
	
	// float denom = eta_i * dot(-dirIn, w_ht) + eta_t * dot(w_ht, scattered);
	// // cout << "t_denom" << fr_t << endl;
	// // float G_t = m_distribution.G(bRec.wi, bRec.wo, w_ht);
	// float G_t = distribution->G(-dirIn, scattered, w_ht, hit, alpha);
	// // float D_t = m_distribution.D(w_ht);
	// float D_t = distribution->D(w_ht, hit, alpha);
	// ft = Color3f(term * eta_t * eta_t * fr_t * G_t * D_t / (denom * denom));
	// if (std::isnan(ft.x) || std::isnan(ft.y) || std::isnan(ft.z)) ft = Color3f(0.0f);
	
	// cout << "fr: " << fr << endl;
	// cout << "ft: "<< ft << endl;

	// if(luminance(fr) == 0.0f){
	// 	cout << "F_r:" << F_r << endl;
	// 	cout << "D_r:" << D_r << endl;
	// 	cout << "G_r:" << G_r << endl;
	// 	cout << "denom:" << (4 * fabsf(dot(-dirIn, hit.sn) * dot(scattered, hit.sn))) << endl;
	// }
	// if(luminance(ft) == 0.0f){
	// 	cout << "F_rt:" << fr_t << endl;
	// 	cout << "D_r:" << D_t << endl;
	// 	cout << "G_r:" << G_t << endl;
	// 	cout << "term:" << term << endl;
	// 	cout << "denom:" << denom << endl;
	// }


	// return fr + ft;


	//new

	// float wiDotN = event.wi.z();
    // float woDotN = event.wo.z();

    bool reflect = dot(-dirIn, hit.sn) * dot(hit.sn, scattered) >= 0.0f;
    // if ((reflect && !sampleR) || (!reflect && !sampleT))
    //     return Vec3f(0.0f);

    

    float eta = dot(-dirIn, hit.sn) < 0.0f ? ior : 1.0f/ior;
    Vec3f m;
    if (reflect)
        // m = sgnE(wiDotN)*(-dirIn + scattered).normalized();
		m = normalize(sign(-dirIn, hit.sn)*(-dirIn + scattered));
    else
        m = -normalize(-dirIn * eta + scattered);
    // float wiDotM = event.wi.dot(m);
    // float woDotM = event.wo.dot(m);
	float cos;

    float F = fresnel(1.0f/ior, dot(-dirIn, m), cos);
    float G = distribution->G(-dirIn, scattered, m, hit, alpha);
    float D = distribution->D(m, hit, alpha);

    if (reflect) {
		
        float fr = (F*G*D*0.25f)/abs(dot(-dirIn, hit.sn));
        return Color3f(fr);
    } 
	else {
		// cout << "refr" << endl;
        float fs = abs(dot(-dirIn, m)*dot(scattered, m)) * (1.0f - F)*G*D / (sqr(eta*dot(-dirIn, m) + dot(scattered, m))*abs(dot(-dirIn, hit.sn)));
        return Color3f(fs);
    }



	
}

float RoughDielectric::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	
	// if (!_refr){
	// 	//reflect
	// 	Vec3f w_h = sign(-dirIn, hit.sn) * normalize(-dirIn + scattered);
	// 	float pdf = distribution->pdf(w_h, hit, alpha);
	// 	Vec3f half_vector = -dirIn + scattered;
	// 	float jacobian_reflect = 1 / 4.0f / fabsf(dot(scattered, w_h));
	// 	// cout << "refl_pdf:" << pdf * jacobian_reflect << endl;

	// 	return pdf * jacobian_reflect;
	// }
	// else {
	// 	//refract
	// 	float cosi = dot(dirIn, hit.sn);
	// 	float etai = 1, etat = ior; 

	// 	if (cosi > 0.0f)
	// 		std::swap(etai, etat);


	// 	Vec3f w_ht = -normalize(etai * -dirIn + etat * scattered);
	// 	float denom = etai * dot(-dirIn, w_ht) + etat * dot(scattered, w_ht);
	// 	float jacobian_t = etat * etat * fabsf(dot(scattered, w_ht)) / (denom * denom);
	// 	float pdf_t = distribution->D(w_ht, hit, alpha) * jacobian_t;
	// 	if (pdf_t == 0.0f) return 0.0f;

	// 	// cout << "refr_pdf:" << pdf_t << endl;
	// 	return pdf_t;
	// }

	// Vec3f w_h = sign(-dirIn, hit.sn) * normalize(-dirIn + scattered);
	// float pdf_r = distribution->pdf(w_h, hit, alpha);

	// // cout << pdf_r << endl;

	// if(!_refr){
	// 	float jacobian = 0.25f / fabsf(dot(w_h, scattered));
	// 	pdf_r *= jacobian;
	// 	return pdf_r;
	// }
	// else{
	// 	float eta_i = 1, eta_t = ior;
	// 	if (dot(-dirIn, hit.sn) < 0.0f)
	// 		std::swap(eta_i, eta_t);

	// 	Vec3f w_ht = -normalize(eta_i * -dirIn + eta_t * scattered);
	// 	float denom = eta_i * dot(w_ht, -dirIn) + eta_t * dot(w_ht, scattered);
	// 	float jacobian_t = eta_t * eta_t * fabsf(dot(w_ht, scattered)) / (denom * denom);
	// 	float pdf_t = distribution->D(w_ht, hit, alpha) * jacobian_t;
	// 	if (pdf_t == 0.0f) return 0.0f;
		

	// 	return pdf_t;
	// }


	//new

	// float wiDotN = event.wi.z();
    // float woDotN = event.wo.z();

    // bool reflect = wiDotN*woDotN >= 0.0f;
	bool reflect = dot(-dirIn, hit.sn) * dot(hit.sn, scattered) >= 0.0f;
    // if ((reflect && !sampleR) || (!reflect && !sampleT))
    //     return 0.0f;

    // float sampleRoughness = (1.2f - 0.2f*std::sqrt(std::abs(wiDotN)))*roughness;
    // float sampleAlpha = Microfacet::roughnessToAlpha(distribution, sampleRoughness);

    // float eta = wiDotN < 0.0f ? ior : 1.0f/ior;
    // Vec3f m;
    // if (reflect)
    //     m = sgnE(wiDotN)*(event.wi + event.wo).normalized();
    // else
    //     m = -(event.wi*eta + event.wo).normalized();


	float eta = dot(-dirIn, hit.sn) < 0.0f ? ior : 1.0f/ior;
    Vec3f m;
    if (reflect)
        // m = sgnE(wiDotN)*(-dirIn + scattered).normalized();
		m = normalize(sign(-dirIn, hit.sn)*(-dirIn + scattered));
    else
        m = -normalize(-dirIn * eta + scattered);
    // float wiDotM = event.wi.dot(m);
    // float woDotM = event.wo.dot(m);
    // float F = Fresnel::dielectricReflectance(1.0f/ior, wiDotM);
	float cos;

    float F = fresnel(1.0f/ior, dot(-dirIn, m), cos);
    // float pm = Microfacet::pdf(distribution, sampleAlpha, m);
	float pm = distribution->pdf(m, hit, alpha);

    float pdf;
    if (reflect){
		pdf = pm * 0.25f/abs(dot(-dirIn, m));
		// cout << "refl:" << pm << endl;
	}  
    else{
		pdf = pm * fabsf(abs(dot(scattered, m)))/sqr(eta*dot(-dirIn, m) + dot(scattered, m));
		// cout << "refr:" << pdf << endl;
	}
        

	if(alpha == 0.0f) return -1;
    
	if (reflect)
		pdf *= F;
	else
		pdf *= 1.0f - F;
    return pdf;
}