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

#include <dirt/scene.h>
#include <dirt/progress.h>
#include <fstream>

/// Construct a new scene from a json object
Scene::Scene(const json & j)
{
    parseFromJSON(j);
}

/// Read a scene from a json file
Scene::Scene(const string & filename)
{
    // open file
    std::ifstream stream(filename, std::ifstream::in);
    if (!stream.good())
    	throw DirtException("Cannot open file: %s.", filename);

    json j;
    stream >> j;
    parseFromJSON(j);
}

Scene::~Scene()
{
    m_materials.clear();
}


shared_ptr<const Material> Scene::findOrCreateMaterial(const json & jp, const string& key) const
{
    auto it = jp.find(key);
    if (it == jp.end())
        return Material::defaultMaterial();
    
    auto j = it.value();
    if (j.is_string())
    {
        string name = j.get<string>();
        // find a pre-declared material
        auto i = m_materials.find(name);
        if (i != m_materials.end())
            return i->second;
        else
            throw DirtException("Can't find a material with name '%s' here:\n%s", name, jp.dump(4));
    }
    else if (j.is_object())
    {
	    // create a new material
        return parseMaterial(j);
    }
    else
        throw DirtException("Type mismatch: Expecting either a material or material name here:\n%s", jp.dump(4));
}

// compute the color corresponding to a ray by raytracing
Color3f Scene::recursiveColor(const Ray3f &ray, int depth) const
{
    putYourCodeHere("Assignment 1: Insert your recursiveColor() code here");
    

    // TODO: Recursively raytrace the scene, similar to the code you wrote in 01_dirt_tutorial
    //       Different to before, you should also take into account surfaces that are self-emitting
    // Pseudo-code:
    //
	// if scene.intersect:
    //      get emitted color (hint: you can use hit.mat->emitted)
	// 		if depth < MaxDepth and hit_material.scatter(....) is successful:
	//			recursive_color = call this function recursively with the scattered ray and increased depth
	//          return emitted color + attenuation * recursive_color
	//		else
	//			return emitted color;
	// else:
	// 		return background color (hint: look at m_background)
    //cout << "dir= " << ray.d << endl;
    const int MaxDepth = 64;
	HitInfo hit;

    Ray3f scattered =  Ray3f(Vec3f(0.0f), Vec3f(0.0f, 0.0f, 0.0f));
	Color3f attenuation;
	//cout << "depth:" << depth << endl;
	if (intersect(ray, hit)){
		Vec3f emit = hit.mat->emitted(ray, hit);
		if(depth < MaxDepth && hit.mat->scatter(ray, hit, attenuation, scattered)){
			return emit + attenuation * recursiveColor(scattered, ++depth);
		}
		else return emit;
	}
	else return m_background;
}

// raytrace an image
Image3f Scene::raytrace() const
{
    // allocate an image of the proper size
    auto image = Image3f(m_camera->resolution().x, m_camera->resolution().y);

    putYourCodeHere("Assignment 1: insert your raytrace() code here");
    const int m_imageSamples = 140;

    // TODO: Render the image, similar to the tutorial
    // Pseudo-code:
    //
        // foreach image row (go over image height)
            // foreach pixel in the row (go over image width)
                // init accumulated color to zero
                // repeat m_imageSamples times:
                    // compute a random point within the pixel (you can just add a random number between 0 and 1
                    //                                          to the pixel coordinate. You can use randf() for this)
                    // compute camera ray
                    // accumulate color raytraced with the ray (by calling recursiveColor)
                // divide color by the number of pixel samples

    // Hint: you can create a Progress object (progress.h) to provide a 
    // progress bar during rendering.
    {
		Progress progress("Rendering", image.size());
		// Generate a ray for each pixel in the ray image
		for (auto y : range(image.height()))
		{
			for (auto x : range(image.width()))
			{
				INCREMENT_TRACED_RAYS;
				
				Color3f total = Color3f(0.f);
                
				for (int i = 0; i < m_imageSamples; i++){
                    auto ray = m_camera->generateRay(x + randf(), y + randf());
					Color3f init = recursiveColor(ray, 0);
					total += init;
				}
				total /= m_imageSamples;
				
				Color3f color = total;

				// TODO: Call recursiveColorFunction ``NumSamples'' times and average the
				// results. Assign the average color to ``color''

				image(x, y) = color;
				++progress;
			}
		}
	}	// progress reporter goes out of scope here

	// return the ray-traced image
    return image;
}

