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

#include <dirt/quad.h>
#include <dirt/scene.h>

Vec2f get_quad_uv(const Vec3f p, const Vec3f botleft, const Vec3f topright){
    float u = (p.x-botleft.x) / (topright.x - botleft.x);
    float v = (p.y-botleft.y) / (topright.y - botleft.y);
    
    //cout << "x: " << p.x << "y: " << p.y << endl;
    //cout << "u = " << u << " v = " << v << endl;
    // if(u <0 || u>1.0) cout << "wrong!!" << endl;
    // if(v <0 || v>1.0) cout << "wrong!!!!!!" << endl;
    return Vec2f(u, v);
}

Quad::Quad(const Vec2f & size,
           shared_ptr<const Material> material,
           const Transform & xform)
	: Surface(xform), m_size(size*0.5f), m_material(material)
{
    v0 = m_xform.vector({m_size.x, 0, 0});
    v1 = m_xform.vector({0, m_size.y, 0});
}

Quad::Quad(const Scene & scene, const json & j)
    : Surface(scene, j)
{
    m_size = j.value("size", m_size);
	m_size /= 2.f;
    
    m_material = scene.findOrCreateMaterial(j);
    v0 = m_xform.vector({m_size.x, 0, 0});
    v1 = m_xform.vector({0, m_size.y, 0});
}

bool Quad::intersect(const Ray3f &ray, HitInfo &hit) const
{
    INCREMENT_INTERSECTION_TESTS;

    // compute ray intersection (and ray parameter), continue if not hit
    auto tray = m_xform.inverse().ray(ray);
    if (tray.d.z == 0)
        return false;
    auto t = -tray.o.z / tray.d.z;
    auto p = tray(t);

    if (m_size.x < p.x || -m_size.x > p.x || m_size.y < p.y || -m_size.y > p.y)
        return false;

    // check if computed param is within ray.mint and ray.maxt
    if (t < tray.mint || t > tray.maxt)
        return false;

	// project hitpoint onto plane to reduce floating-point error
	p.z = 0;

    Vec3f gn = normalize(m_xform.normal({0,0,1}));
    Vec2f uv = get_quad_uv((p), (Vec3f(-m_size.x, -m_size.y, 0)), (Vec3f(m_size.x, m_size.y, 0)));

    // if hit, set intersection record values
    hit = HitInfo(t, m_xform.point(p), gn, gn,
                         uv, /* TODO: Compute proper UV coordinates */
                         m_material.get(), this);
    return true;
}


Box3f Quad::localBBox() const
{
    return Box3f(-Vec3f(m_size.x,m_size.y,0) - Vec3f(1e-4f), Vec3f(m_size.x,m_size.y,0) + Vec3f(1e-4f));
}

Vec3f Quad::sample(const Vec3f& o) const
{
    
    // Vec3f p = sampleRect(m_xform.point(Vec3f(0.f, 0.f, 0.f)), v0, v1);
    Vec3f p = sampleRect(Vec3f(0.f, 0.f, 0.f), {m_size.x, 0, 0}, {0, m_size.y, 0});

    return normalize(m_xform.point(p) - o);
}

float Quad::pdf(const Vec3f& o, const Vec3f& v) const
{
    HitInfo hit;
    Ray3f ray = Ray3f(o, v);

    
    if (intersect(ray, hit)){
        float area = length(cross(v0, v1)) * 4;
        float cos = abs(dot(normalize(-v), normalize(hit.sn)));
        float geometry_factor = length2(hit.p - o) / cos; 
        
        return 1.0f / area * geometry_factor;
    }
    else return 0.0f;
}


void Quad::SampleFromEmit(Vec3f &pos, Vec3f &dir) const
{
    Vec3f o = sampleRect(Vec3f(0.f, 0.f, 0.f), {m_size.x, 0, 0}, {0, m_size.y, 0}) ;
    Vec3f p = randomOnUnitHemisphere() + o;
    
    
    pos = m_xform.point(o);
    dir = normalize(m_xform.vector(p-o));
}