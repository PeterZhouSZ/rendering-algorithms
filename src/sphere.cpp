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

#include <dirt/sphere.h>
#include <dirt/scene.h>
#include <dirt/onb.h>


Vec2f get_sphere_uv( const Vec3f p){
    
    // float v = acos(clamp(p.z, -1.0f, 1.0f))/ M_PI; 
    // v = (1 - v);
    
    float phi = atan2(p.y, p.x);
    float theta = acos(p.z);
    float u = 0.5 + phi / (2*M_PI);
    float v = 1 - theta / (M_PI);
    
    return Vec2f(u, v);
}

Sphere::Sphere(float radius,
               shared_ptr<const Material> material,
               const Transform & xform)
    : Surface(xform), m_radius(radius), m_material(material)
{

}

Sphere::Sphere(const Scene & scene, const json & j)
    : Surface(scene, j)
{
	m_radius = j.value("radius", m_radius);
    m_material = scene.findOrCreateMaterial(j);
}

Box3f Sphere::localBBox() const
{
    return Box3f(Vec3f(-m_radius), Vec3f(m_radius));
}



bool Sphere::intersect(const Ray3f &ray, HitInfo &hit) const
{
    INCREMENT_INTERSECTION_TESTS;
    // TODO: Assignment 1: Implement ray-sphere intersection
    
    auto tray = m_xform.inverse().ray(ray);

    putYourCodeHere("Assignment 1: Insert your ray-sphere intersection code here");
    Vec3f ray2center = tray.o - center;
    float a = dot(tray.d, tray.d);
    float b = 2.0f * dot(ray2center, tray.d);
    float c = dot(ray2center, ray2center) - m_radius * m_radius;
    
    //cout << "b2: " << b*b  << endl;
    //cout << "c:" << c << endl;
    // TODO: If the ray misses the sphere, you should return false
    // TODO: If you successfully hit something, you should compute the hit point, 
    //       hit distance, and normal and fill in these values
    float hitT = 0.0f;
    Vec3f hitPoint;
    Vec3f geometricNormal;

    // For this assignment you can leave these values as is
    Vec3f shadingNormal;
    Vec2f uvCoordinates;

    // You should only assign hit and return true if you successfully hit something
    if ((b*b - 4*a*c) < 0)
        return false;
    else {
        float t1 = (-b - sqrt(b*b - 4*a*c)) / (2*a);
        float t2 = (-b + sqrt(b*b - 4*a*c)) / (2*a);


        if (t1 < tray.maxt && t1 > tray.mint) {
            hitT = t1;
            hitPoint = tray.o + tray.d * hitT;
            hitPoint *= m_radius / length(hitPoint - center);
            
            geometricNormal = hitPoint - center;
            geometricNormal = normalize(geometricNormal);
            shadingNormal = geometricNormal;
            uvCoordinates = get_sphere_uv(hitPoint);
            
            hit = HitInfo(hitT,
            m_xform.point(hitPoint),
            normalize(m_xform.normal(geometricNormal)),
            normalize(m_xform.normal(shadingNormal)),
            uvCoordinates,
            m_material.get(),
            this);
            //cout << "point:" << hit.p <<endl;

            return true;
        }
        
        if (t2 < tray.maxt && t2 > tray.mint) {
            hitT = t2;
            hitPoint = tray.o + tray.d * hitT;
            hitPoint *= m_radius / length(hitPoint - center);

            
            geometricNormal = (hitPoint - center);
            geometricNormal = normalize(geometricNormal);
            shadingNormal = geometricNormal;
            uvCoordinates = get_sphere_uv(hitPoint);
            
            hit = HitInfo(hitT,
            m_xform.point(hitPoint),
            normalize(m_xform.normal(geometricNormal)),
            normalize(m_xform.normal(shadingNormal)),
            uvCoordinates,
            m_material.get(),
            this);
            //cout << "point:" << hit.p <<endl;


            return true;
        }
    }
    return false;
    
}

Vec3f Sphere::sample(const Vec3f& o) const
{
    
    // Vec3f p = sampleRect(m_xform.point(Vec3f(0.f, 0.f, 0.f)), v0, v1);
    // float d = length(o - center);
    float d = length(m_xform.point(center) - o);
    float cosThetaMax = sqrt(d * d - m_radius * m_radius) / d;

    if (d * d - m_radius * m_radius <= 0){
        Vec3f p = sampleSphere(m_xform.point(center), m_radius);
        return normalize(p - o); 
    }
    else{
        onb k;
        k.build_from_w(m_xform.point(center) - o);
        return k.local(randomSphericalCap(cosThetaMax));
    }
    

    
}

float Sphere::pdf(const Vec3f& o, const Vec3f& v) const
{
    HitInfo hit;
    Ray3f ray = Ray3f(o, v);

    
    if (intersect(ray, hit)){
        float d = length(o - m_xform.point(center));
        float cosThetaMax = sqrt(d * d - m_radius * m_radius) / d;
        if(d * d - m_radius * m_radius <= 0){
            
            return 1 / (4*M_PI*m_radius*m_radius) ; 
        }
        return 1.0f / (2 * M_PI * (1 - cosThetaMax));
    }
    else return 0.0f;
}
