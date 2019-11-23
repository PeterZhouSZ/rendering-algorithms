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

#include <dirt/mesh.h>
#include <dirt/scene.h>

// Ray-Triangle intersection
// p0, p1, p2 - Triangle vertices
// n0, n1, n2 - optional per vertex normal data
// t0, t1, t2 - optional per vertex texture coordinates
bool singleTriangleIntersect(const Ray3f& ray,
	                         const Vec3f& p0, const Vec3f& p1, const Vec3f& p2,
	                         const Vec3f* n0, const Vec3f* n1, const Vec3f* n2,
                             const Vec2f* uv0, const Vec2f* uv1, const Vec2f* uv2,
	                         HitInfo& hit,
	                         const Material * material,
	                         const SurfaceBase * surface)
{
    // TODO: Implement ray-triangle intersection
    // TODO: If the ray misses the triangle, you should return false
    //       You can pick any ray triangle intersection routine you like.
    //       I recommend you follow the Moller-Trumbore algorithm

    putYourCodeHere("Insert your ray-triangle intersection code here");
    

    // First, check for intersection and fill in the hitT
    float hitT = 0.0f;
    // You should also compute the u/v (i.e. the alpha/beta barycentric coordinates) of the hit point
    // (Moller-Trumbore gives you this for free)
    float u, v;

    const float epsilon = 0.0000001;
    Vec3f edge1 = p1 - p0;
    Vec3f edge2 = p2 - p0;
    Vec3f pvec = cross(ray.d, edge2);

    float det = dot(edge1, pvec);
    if(abs(det) < epsilon) 
        return false;
    
    float invDet = 1 / det;
    Vec3f tvec = ray.o - p0;
    u = dot(tvec, pvec) * invDet;

    if(u < 0 || u > 1)
        return false;

    Vec3f qvec = cross(tvec, edge1); 
    v = dot(ray.d, qvec) * invDet; 
    if (v < 0 || (u + v) > 1) 
        return false; 

    hitT = invDet * dot(edge2, qvec);

    if (hitT < ray.mint || hitT > ray.maxt)
        return false;
    // TODO: If you successfully hit the triangle, you should check if the hitT lies
    //       within the ray's tmin/tfar, and return false if it does not

    // TODO: Fill in the geometric normal with the geometric normal of the triangle (i.e. normalized cross product of the sides)
    Vec3f gn = cross(edge1, edge2);
    gn = normalize(gn);

    // Compute the shading normal
    Vec3f sn;
    if (n0 != nullptr && n1 != nullptr && n2 != nullptr) { // Do we have per-vertex normals available?
        // We do -> dereference the pointers
        Vec3f normal0 = *n0;
        Vec3f normal1 = *n1;
        Vec3f normal2 = *n2;

        // TODO: You should compute the shading normal by
        //       doing barycentric interpolation of the per-vertex normals (normal0/1/2)
        //       Make sure to normalize the result
        sn = normal0 * (1 - u - v) + normal1 * u + normal2 * v;
        sn = normalize(sn);
    } else {
        // We don't have per-vertex normals - just use the geometric normal
        sn = gn;
    }

    Vec2f interp_uv = Vec2f(u, v);
    if(uv0 != nullptr && uv1 != nullptr && uv2 != nullptr){
        interp_uv = *uv0 * (1 - u - v) + *uv1 * u + *uv2 * v;
    }
    

    // Because we've hit the triangle, fill in the intersection data
    hit = HitInfo(hitT, ray(hitT), gn, sn, interp_uv, material, surface);
    return true;
}

bool TestsingleTriangleIntersect(const Ray3f& ray,
	                         const Vec3f& p0, const Vec3f& p1, const Vec3f& p2,
	                         const Vec3f* n0, const Vec3f* n1, const Vec3f* n2,
	                         HitInfo& hit,
	                         const Material * material,
	                         const SurfaceBase * surface)
{
    // TODO: Implement ray-triangle intersection
    // TODO: If the ray misses the triangle, you should return false
    //       You can pick any ray triangle intersection routine you like.
    //       I recommend you follow the Moller-Trumbore algorithm

    putYourCodeHere("Insert your ray-triangle intersection code here");
    

    // First, check for intersection and fill in the hitT
    float hitT = 0.0f;
    // You should also compute the u/v (i.e. the alpha/beta barycentric coordinates) of the hit point
    // (Moller-Trumbore gives you this for free)
    float u, v;

    const float epsilon = 0.0000001;
    Vec3f edge1 = p1 - p0;
    Vec3f edge2 = p2 - p0;
    Vec3f pvec = cross(ray.d, edge2);

    float det = dot(edge1, pvec);
    if(abs(det) < epsilon) 
        return false;
    
    float invDet = 1 / det;
    Vec3f tvec = ray.o - p0;
    u = dot(tvec, pvec) * invDet;

    if(u < 0 || u > 1)
        return false;

    Vec3f qvec = cross(tvec, edge1); 
    v = dot(ray.d, qvec) * invDet; 
    if (v < 0 || (u + v) > 1) 
        return false; 

    hitT = invDet * dot(edge2, qvec);

    if (hitT < ray.mint || hitT > ray.maxt)
        return false;
    // TODO: If you successfully hit the triangle, you should check if the hitT lies
    //       within the ray's tmin/tfar, and return false if it does not

    // TODO: Fill in the geometric normal with the geometric normal of the triangle (i.e. normalized cross product of the sides)
    Vec3f gn = cross(edge1, edge2);
    gn = normalize(gn);

    // Compute the shading normal
    Vec3f sn;
    if (n0 != nullptr && n1 != nullptr && n2 != nullptr) { // Do we have per-vertex normals available?
        // We do -> dereference the pointers
        Vec3f normal0 = *n0;
        Vec3f normal1 = *n1;
        Vec3f normal2 = *n2;

        // TODO: You should compute the shading normal by
        //       doing barycentric interpolation of the per-vertex normals (normal0/1/2)
        //       Make sure to normalize the result
        sn = normal0 * (1 - u - v) + normal1 * u + normal2 * v;
        sn = normalize(sn);
    } else {
        // We don't have per-vertex normals - just use the geometric normal
        sn = gn;
    }
    
    Vec2f interp_uv = Vec2f(u, v);
    

    // Because we've hit the triangle, fill in the intersection data
    hit = HitInfo(hitT, ray(hitT), gn, sn, interp_uv, material, surface);
    return true;
}



Triangle::Triangle(const Scene & scene, const json & j, shared_ptr<const Mesh> mesh, uint32_t triNumber)
    : m_mesh(mesh), m_faceIdx(triNumber)
{
    
}

bool Triangle::intersect(const Ray3f &ray, HitInfo &hit) const
{
    INCREMENT_INTERSECTION_TESTS;

    auto i0 = m_mesh->F[m_faceIdx].x,
         i1 = m_mesh->F[m_faceIdx].y,
         i2 = m_mesh->F[m_faceIdx].z;
    auto p0 = m_mesh->V[i0],
         p1 = m_mesh->V[i1],
         p2 = m_mesh->V[i2];
    const Vec3f * n0 = nullptr, *n1 = nullptr, *n2 = nullptr;
    if (!m_mesh->N.empty())
    {
        n0 = &m_mesh->N[i0];
        n1 = &m_mesh->N[i1];
        n2 = &m_mesh->N[i2];
    }
    const Vec2f * t0 = nullptr, *t1 = nullptr, *t2 = nullptr;
    if (!m_mesh->UV.empty())
    {
        t0 = &m_mesh->UV[i0];
        t1 = &m_mesh->UV[i1];
        t2 = &m_mesh->UV[i2];
    }

    return singleTriangleIntersect(ray,
                                   p0, p1, p2,
                                   n0, n1, n2,
                                   t0, t1, t2,
                                   hit, m_mesh->material.get(), this);
}

Box3f Triangle::localBBox() const
{
	// all mesh vertices have already been transformed to world space,
	// so we need to transform back to get the local space bounds
    Box3f result;
    result.enclose(m_mesh->m_xform.inverse().point(vertex(0)));
    result.enclose(m_mesh->m_xform.inverse().point(vertex(1)));
    result.enclose(m_mesh->m_xform.inverse().point(vertex(2)));
    
    // if the triangle lies in an axis-aligned plane, expand the box a bit
    auto diag = result.diagonal();
    for (int i = 0; i < 3; ++i)
    {
        if (diag[i] < 1e-4f)
        {
            result.pMin[i] -= 5e-5f;
            result.pMax[i] += 5e-5f;
        }
    }
    return result;
}

Box3f Triangle::worldBBox() const
{
    // all mesh vertices have already been transformed to world space,
    // so just bound the triangle vertices
    Box3f result;
    result.enclose(vertex(0));
    result.enclose(vertex(1));
    result.enclose(vertex(2));

    // if the triangle lies in an axis-aligned plane, expand the box a bit
    auto diag = result.diagonal();
    for (int i = 0; i < 3; ++i)
    {
        if (diag[i] < 1e-4f)
        {
            result.pMin[i] -= 5e-5f;
            result.pMax[i] += 5e-5f;
        }
    }
    return result;
}

Vec3f Triangle::sample(const Vec3f& o) const
{
    
    // Vec3f p = sampleRect(m_xform.point(Vec3f(0.f, 0.f, 0.f)), v0, v1);
    auto i0 = m_mesh->F[m_faceIdx].x,
         i1 = m_mesh->F[m_faceIdx].y,
         i2 = m_mesh->F[m_faceIdx].z;

    auto p0 = m_mesh->V[i0],
         p1 = m_mesh->V[i1],
         p2 = m_mesh->V[i2];

    
    Vec3f p = sampleTriangle(p0, p1, p2);

    return normalize(p - o);
}

float Triangle::pdf(const Vec3f& o, const Vec3f& v) const
{
    HitInfo hit;
    Ray3f ray = Ray3f(o, v);

    Vec3f edge1 = vertex(1) - vertex(0);
    Vec3f edge2 = vertex(2) - vertex(0);

    if (intersect(ray, hit)){
        float area = length(cross(edge2, edge1)) / 2; //why
        
        float cos = abs(dot(normalize(-v), normalize(hit.sn)));
        float geometry_factor = length2(hit.p - o) / cos; 
        
        return 1.0f / area * geometry_factor;
    }
    else return 0.0f;
}