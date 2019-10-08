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

#include <dirt/transform.h>
#include <dirt/vec.h>
#include <dirt/parser.h>

/**
    This class represents a virtual pinhole camera.
   
    The camera is responsible for generating primary rays. It is positioned
    using a Transform and points along the -z axis of the local coordinate
    system. It has an image plane positioned a z = -dist with size
    (width, height).
   
    We currently only support pinhole perspective cameras. This class could
    be made into a virtual base class to support other types of cameras
    (e.g. an orthographic camera, or omni-directional camera).
 */




class Camera
{
public:
    /// Construct a camera from json parameters.
    Camera(const json & j = json())
    {
		m_xform = j.value("transform", m_xform);
	    m_resolution = j.value("resolution", m_resolution);
	    m_focalDistance = j.value("fdist", m_focalDistance);
	    m_apertureRadius = j.value("aperture", m_apertureRadius);
        
		
		// Default vfov value. Override this with the value from json
        // TODO: Assignment 1: read the vertical field-of-view from j ("vfov"),
        // and compute the width and height of the image plane. Remember that
        // the "vfov" parameter is specified in degrees, but C++ math functions
        // expect it in radians. You can use deg2rad() from common.h to convert
        // from one to the other
        putYourCodeHere("Assignment 1: Compute the image plane size.");
        fdist = j.value("fdist", 1.0f);
        vfov = j.value("vfov", vfov); // vertical fov
        m_size = Vec2f(2*(fdist * tan(deg2rad(vfov / 2.0f))) * (m_resolution[0] / m_resolution[1]), 2*(fdist * tan(deg2rad(vfov / 2.0f))));
        //cout << "size: " << m_size.x << "," << m_size.y << endl;
        /*
        if (!key_exists(j, "transform"))
            cout << "wtf" << endl;
        */
        /*
        if (!key_exists(j["transform"], "up"))
            cout << "default up vector = (0,1,0)" << endl;
        else {
            up = j["transform"]["up"];
            cout << "LUL" << endl;
        }
        */
       /*
        if (key_exists(j, "transform")){
            if (key_exists(j["transform"], "up"))
                up = j["transform"]["up"];
        }
        */
            
        //up = j["transform"]["up"];

        unit_horizontal = normalize(cross(look, up));
        
        unit_vertical = normalize(up);
        
        
    }

	/// Return the camera's image resolution
	Vec2d resolution() const {return m_resolution;}

    /**
        Generate a ray going through image-plane location (u,v).

        (\c u,\c v) range from 0 to m_resolution.x() and m_resolution.y() along
       	the x- and y-axis of the rendered image, respectively

        \param u 	The horizontal position within the image
        \param v  	The vertical position within the image
        \return 	The \ref Ray3f data structure filled with the
       				appropriate position and direction
     */
    Ray3f generateRay(float u, float v) const
    {
        // TODO: Assignment 1: Implement camera ray generation
        putYourCodeHere("Assignment 1: Insert your camera ray generation code here");
        Vec3f midPoint = origin + normalize(look) * fdist;

        Vec3f topleft = midPoint + unit_vertical * m_size[1] / 2 - unit_horizontal * m_size[0] / 2;
        

        return Ray3f(m_xform.point(origin), m_xform.vector(topleft + u * m_size.x * unit_horizontal / m_resolution[0] - v * unit_vertical * m_size.y / m_resolution[1] - origin)); 
        
    }
    bool key_exists(const json& j, const string& key)
    {
        return j.find(key) != j.end();
    }

private:
	//
	// The camera setup looks something like this, where the
	// up vector points out of the screen:
	//
	//         top view                         side view
	//            ^                    up
	//            |                     ^
	//            |                     |             _,-'
	//          width                   |         _,-'   |
	//       +----|----+     +          |     _,-'       | h
	//        \   |   /    d |        e | _,-'           | e
	//         \  |  /     i |        y +'---------------+-i----->
	//          \ | /      s |        e  '-,_   dist     | g
	//           \|/       t |               '-,_        | h
	//            +          +                   '-,_    | t
	//           eye                                 '-,_|
	//


	Transform m_xform = Transform();      ///< Local coordinate system
	Vec2f m_size = Vec2f(1,1);            ///< Physical size of the image plane
	float m_focalDistance = 1.f;          ///< Distance to image plane along local z axis
	Vec2d m_resolution = Vec2d(512,512);  ///< Image resolution
	float m_apertureRadius = 0.f;         ///< The size of the aperture for depth of field
    float vfov = 90.0f;
    Vec3f origin = Vec3f(0.f, 0.f, 0.f);
    Vec3f look  = Vec3f(0.f, 0.f, -1.0f);
    Vec3f up = Vec3f(0.f, 1.f, 0.f);
    float fdist = 1.0f;
    Vec3f unit_horizontal;
    Vec3f unit_vertical;
};
