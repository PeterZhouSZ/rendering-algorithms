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
	Vec3f scatteredPoint =  hit.sn + hit.p + normalize(randomInUnitSphere());
	
	scattered = Ray3f(hit.p, hit.sn + normalize(randomInUnitSphere()));
	attenuation = albedo->value(hit);
	return true;
}


Metal::Metal(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	roughness = parseTexture(j.at("roughness"));
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


Dielectric::Dielectric(const json & j)
{
	ior = j.value("ior", ior); // index of refraction
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

