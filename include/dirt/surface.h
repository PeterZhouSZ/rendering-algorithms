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

#pragma once

#include <dirt/fwd.h>
#include <dirt/material.h>
#include <dirt/transform.h>
#include <dirt/parser.h>


/**
    Contains information about a ray-surface intersection hit point.

    Used by surface intersection routines to return more than just a single
    value. Includes the position, traveled ray distance, uv coordinates, the
    geometric and interpolated shading normals, and a pointer to the intersected
    surface and underlying material.
 */
struct HitInfo
{
	float t;                                ///< Ray parameter for the hit
	Vec3f p;                                ///< Hit position
	Vec3f gn;                               ///< Geometric normal
	Vec3f sn;                               ///< Interpolated shading normal
	Vec2f uv;                               ///< UV texture coordinates
	const Material * mat = nullptr;         ///< Material at the hit point
	const SurfaceBase * surface = nullptr;  ///< Surface at the hit point

	/// Default constructor that leaves all members uninitialized
	HitInfo() = default;

	/// Parameter constructor that initializes all data members
	HitInfo(float t,
            const Vec3f &p,
            const Vec3f &gn,
            const Vec3f &sn,
            const Vec2f &uv,
            const Material * m = nullptr,
            const SurfaceBase * s = nullptr) :
		t(t), p(p), gn(gn), sn(sn), uv(uv), mat(m), surface(s)
	{

	}
};

/// This is the abstract superclass for all surfaces.
class SurfaceBase
{
public:
    virtual ~SurfaceBase() {}

    /// Return the surface's local-space AABB.
    virtual Box3f localBBox() const = 0;

    /**
        Return the surface's world-space AABB.

        By default just returns the local bounding box.
    */
    virtual Box3f worldBBox() const {return localBBox();}

    /**
        Add a child surface (if this is an aggregate).

        This function will become useful once we create groups of objects. The
        base class implementation just throws an error.

        This function should only be used before \ref build() is called.
    */
    virtual void addChild(shared_ptr<SurfaceBase> surface)
    {
        throw DirtException("This surface does not support children.");
    }

    /**
        Perform any necessary precomputation before ray tracing.
       
        If a surface requires some precomputation (e.g. building an 
        acceleration structure) overload this function. This will be called 
        after the last child has been added and before any call to intersect().

        The base class implementation just does nothing.
    */
    virtual void build() {};

    /**
        Ray-Surface intersection test.

        Intersect a ray against this surface and return detailed intersection
        information.
        
        \param ray
             A 3-dimensional ray data structure with minimum/maximum
             extent information
        \param hit
             A detailed intersection record, which will be filled by the
             intersection query
        \return  \c true if an intersection was found
     */
    virtual bool intersect(const Ray3f &ray, HitInfo &hit) const = 0;


    /// Return whether or not this Surface's Material is emissive.
    virtual bool isEmissive() const {return false;}
};


/**
    Surfaces represent the geometry of the scene. A Surface could be an
    individual primitive like a \ref Sphere, or it could be composed of
    many smaller primitives, e.g. the triangles composing a \ref Mesh.
    Each surface currently stores a transformation matrix which
    positions/orients the surface in the scene and pointer to a single \ref
    Material which specifies the light reflectance properties.
 */
class Surface : public SurfaceBase
{
public:
	Surface(const Transform & xform = Transform()) : m_xform(xform) {}

    Surface(const Scene & scene, const json & j = json::object());
    virtual ~Surface() {}

    /**
        Return the surface's world-space AABB.

        By default just returns the local bounding box transformed to world 
        space by the Surface's transformation.
    */
    Box3f worldBBox() const override;

protected:
    Transform m_xform = Transform();        ///< Local-to-world Transformation
};