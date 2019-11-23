
#include <dirt/parser.h>
#include <dirt/scene.h>
#include <dirt/surface.h>
#include <dirt/integrator.h>



Color3f NormalIntegrator::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    HitInfo hit;
    if (scene.intersect(ray, hit))  
        return Color3f(abs(hit.sn));
    else return Color3f(0.0f);
    
}



Color3f AmbientOcclusionIntegrator::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    HitInfo hit;
    ScatterRecord srec;
    if (!scene.intersect(ray, hit))  
        return Color3f(0.0f);
    bool scattered_result = hit.mat->sample(ray.d, hit, srec);
    Ray3f shadowRay = Ray3f(hit.p, srec.scattered);
    if (scene.intersect(shadowRay, hit))
        return Color3f(0.0f);
    else return Color3f(1.0f);
    
}



Color3f PathTracerMaterials::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    HitInfo hit;
    ScatterRecord srec;
    
    const int MaxDepth = 64;

	if (scene.intersect(ray, hit)){
		Vec3f emit = hit.mat->emitted(ray, hit);
        // cout << "a" << endl;
		if(depth < MaxDepth && hit.mat->sample(ray.d, hit, srec)){
            // cout << "b" << endl;
            if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0){
                
                return emit + hit.mat->eval(ray.d, srec.scattered, hit) / hit.mat->pdf(ray.d, srec.scattered, hit) * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);
            }
			    
            else return emit + srec.attenuation * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);

		}
		else return emit;
	}
	else return scene.get_background();
}




Color3f DirectMats::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    HitInfo hit;
    ScatterRecord srec;
    // depth = 0;
    Color3f L(0.0f);
    Color3f throughput(1.0f);

    Ray3f traced_ray = ray;
    
    while (depth < 2){
        if(!scene.intersect(traced_ray, hit)){
            L += throughput * scene.get_background();
            break;
        }
        else{
            Vec3f emit = hit.mat->emitted(traced_ray, hit);
            L += throughput * emit;
        
            if(hit.mat->sample(traced_ray.d, hit, srec)){
                if (!srec.isSpecular && hit.mat->pdf(traced_ray.d, srec.scattered, hit) != 0){
                    
                    throughput *= hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit) ;
                }
                else throughput *= srec.attenuation ; //abs(dot(normalize(hit.sn), srec.scattered))??

            }
            traced_ray = Ray3f(hit.p, srec.scattered);
            depth++;
        }
    }
    return L;
    /*
    recursuve way
    */

	// if (scene.intersect(ray, hit)){
	// 	Vec3f emit = hit.mat->emitted(ray, hit);
        
	// 	if(depth < max_depth && hit.mat->sample(ray.d, hit, srec)){
    //         if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0){
                
    //             return emit + hit.mat->eval(ray.d, srec.scattered, hit) / hit.mat->pdf(ray.d, srec.scattered, hit) * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);
    //         }
			    
    //         else return emit + srec.attenuation * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);

	// 	}
	// 	else return emit;
	// }
	// else return scene.get_background();
    
}

Color3f DirectNEE::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    HitInfo hit;
    ScatterRecord srec;
    
    Color3f L(0.0f);
    Color3f throughput(1.0f);
    depth = 0;

    Ray3f traced_ray = ray;
    
    while (depth < 2){
        if(!scene.intersect(traced_ray, hit)){
            L += throughput * scene.get_background();
            break;
        }
        else{
            Vec3f light_sample = scene.emitters().sample(hit.p);
            Vec3f emit = hit.mat->emitted(traced_ray, hit);
            bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);
            srec.scattered = light_sample;
            L += throughput * emit;
        
            
            if (scene.emitters().pdf(hit.p, light_sample) != 0 && !srec.isSpecular){
                
                throughput *= hit.mat->eval(traced_ray.d, light_sample, hit) / scene.emitters().pdf(hit.p, light_sample) ;
            }
            else throughput *= srec.attenuation ;

            
            traced_ray = Ray3f(hit.p, srec.scattered);
            depth++;
        }
    }
    return L;

    /*
    recursuve way
    */
    

	// if (scene.intersect(ray, hit)){
	// 	Vec3f emit = hit.mat->emitted(ray, hit);
        
	// 	if(depth < max_depth){
    //         Vec3f light_sample = scene.emitters().sample(hit.p);
            
    //         srec.scattered = light_sample;
    //         if ( scene.emitters().pdf(hit.p, light_sample) != 0){
                
    //             return emit + hit.mat->eval(ray.d, light_sample, hit) / scene.emitters().pdf(hit.p, light_sample) * Li(scene, Ray3f(hit.p, light_sample), ++depth);
    //         }
	// 		else{
                
    //             return emit + srec.attenuation * Li(scene, Ray3f(hit.p, light_sample), ++depth);
    //         }

	// 	}
	// 	else return emit;
	// }
	// else return scene.get_background();
    
}

Color3f DirectMIS::Li(const Scene &scene, const Ray3f &ray, int depth) const
{  
    HitInfo hit;
    ScatterRecord srec;
    // depth = 0;
    // Color3f L(0.0f);
    // Color3f throughput(1.0f);

    Ray3f traced_ray = ray;
    
    // while (depth < 2){
    //     if(!scene.intersect(traced_ray, hit)){
    //         L += throughput * scene.get_background();
    //         break;
    //     }
    //     else{
    //         float prob = randf();
        
    //         Vec3f light_sample = scene.emitters().sample(hit.p);
    //         Vec3f emit = hit.mat->emitted(traced_ray, hit);
    //         bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);
    //         L += throughput * emit;
            
    //         // if (hit.mat->pdf(traced_ray.d, srec.scattered, hit) == 0) average_pdf = scene.emitters().pdf(hit.p, light_sample);
    //         // if (scene.emitters().pdf(hit.p, light_sample) == 0) average_pdf = hit.mat->pdf(traced_ray.d, srec.scattered, hit);
        
    //         if (prob > 0.5){
    //             if (scene.emitters().pdf(hit.p, light_sample) != 0){
    //                 float average_pdf = (hit.mat->pdf(traced_ray.d, light_sample, hit) + scene.emitters().pdf(hit.p, light_sample)) / 2;
    //                 throughput *= hit.mat->eval(traced_ray.d, light_sample, hit) / average_pdf ;
    //             }
    //             else throughput *= srec.attenuation ;
    //             traced_ray = Ray3f(hit.p, light_sample);
    //         }
    //         else{
    //             if(sample_result){
    //                 if (!srec.isSpecular && hit.mat->pdf(traced_ray.d, srec.scattered, hit) != 0){
    //                     float average_pdf = (hit.mat->pdf(traced_ray.d, srec.scattered, hit) + scene.emitters().pdf(hit.p, srec.scattered)) / 2;

    //                     throughput *= hit.mat->eval(traced_ray.d, srec.scattered, hit) / average_pdf ;
    //                 }
    //                 else throughput *= srec.attenuation ; 
            
    //             }
    //             traced_ray = Ray3f(hit.p, srec.scattered);
    //         }
            
            
    //         depth++;
            
    //     }
        
    // }
    // return L;
    
    float prob = randf();
	if (scene.intersect(ray, hit)){
		Vec3f emit = hit.mat->emitted(ray, hit);
        Vec3f light_sample = scene.emitters().sample(hit.p);
        bool m_sample_result = hit.mat->sample(ray.d, hit, srec); // material sample result

        if (prob > 0.5){
            if(depth < max_depth && m_sample_result){
                
                if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0){
                    float average_pdf = (hit.mat->pdf(ray.d, srec.scattered, hit) + scene.emitters().pdf(hit.p, srec.scattered)) / 2;
                    return emit + hit.mat->eval(ray.d, srec.scattered, hit) / average_pdf * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);
                }
                else {
                    return emit + srec.attenuation * Li(scene, Ray3f(hit.p, srec.scattered), ++depth);
                }
            }
            else return emit;
        }
        else {
            if(depth < max_depth){
                
                // srec.scattered = light_sample;
                if ( scene.emitters().pdf(hit.p, light_sample) != 0){
                    float average_pdf = (hit.mat->pdf(ray.d, light_sample, hit) + scene.emitters().pdf(hit.p, light_sample)) / 2;
                    // cout << scene.emitters().pdf(hit.p, light_sample) << endl;

                    return emit + hit.mat->eval(ray.d, light_sample, hit) / average_pdf * Li(scene, Ray3f(hit.p, light_sample), ++depth);
                }
                else{
                    return emit + srec.attenuation * Li(scene, Ray3f(hit.p, light_sample), ++depth);
                }

		    }
		    else return emit;
        }
		
	}
	else return scene.get_background();
    
}

Color3f PathTracerNEE::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    return Li(scene, ray, depth, true);
}

Color3f PathTracerNEE::Li(const Scene &scene, const Ray3f &ray, int depth, bool include) const
{
    HitInfo hit;
    ScatterRecord srec;
    
    if (scene.intersect(ray, hit)){
		Vec3f emit = hit.mat->emitted(ray, hit) * include;
        if(depth == max_depth) return emit;

        Vec3f light_sample = scene.emitters().sample(hit.p);
        bool m_sample_result = hit.mat->sample(ray.d, hit, srec); // material sample result
        
        Color3f L_dir = Color3f(0.0f);
        
    
        float pdf_e = scene.emitters().pdf(hit.p, light_sample);
        // cout << "e:" << pdf_e << endl;
        // cout << "m:" << pdf_m << endl;
        // cout << mis << endl;
        
        if (depth < max_depth && pdf_e > 0){
            if(!srec.isSpecular)
                L_dir = hit.mat->eval(ray.d, light_sample, hit) / pdf_e * Li(scene, Ray3f(hit.p, light_sample), max_depth, true);
        }
        

        if(depth < max_depth && m_sample_result){
                
            if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0){
                
                return emit + hit.mat->eval(ray.d, srec.scattered, hit) / hit.mat->pdf(ray.d, srec.scattered, hit) * Li(scene, Ray3f(hit.p, srec.scattered), ++depth, false)  + L_dir;
            }
            else {
                return emit + srec.attenuation * Li(scene, Ray3f(hit.p, srec.scattered), ++depth, true) + L_dir;
            }
        }
        else return emit ;
        
    }
    else return scene.get_background();
    
}


Color3f PathTracerMis::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    // HitInfo hit;
    // ScatterRecord srec;
    
    // if (scene.intersect(ray, hit)){
	// 	Vec3f emit = hit.mat->emitted(ray, hit);
    //     if(depth == max_depth) return emit;

    //     Vec3f light_sample = scene.emitters().sample(hit.p);
    //     bool m_sample_result = hit.mat->sample(ray.d, hit, srec); // material sample result
        
    //     Color3f L_dir = Color3f(0.0f);
        
    
    //     float pdf_m = 0; //bug with metal,dielectric
    //     float pdf_e = scene.emitters().pdf(hit.p, light_sample);

    //     if(!srec.isSpecular) pdf_m = hit.mat->pdf(ray.d, srec.scattered, hit);
        
    //     float mis = pdf_e / (pdf_m + pdf_e);
    //     // cout << "e:" << pdf_e << endl;
    //     // cout << "m:" << pdf_m << endl;
    //     // cout << mis << endl;
        
    //     if (depth < max_depth && pdf_e > 0){
    //         if(!srec.isSpecular)
    //             L_dir = hit.mat->eval(ray.d, light_sample, hit) / pdf_e * Li(scene, Ray3f(hit.p, light_sample), max_depth) * mis ;
    //     }
        

    //     if(depth < max_depth && m_sample_result){
                
    //         if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0){
                
    //             return emit * (1-mis)+ hit.mat->eval(ray.d, srec.scattered, hit) / hit.mat->pdf(ray.d, srec.scattered, hit) * Li(scene, Ray3f(hit.p, srec.scattered), ++depth)  + L_dir;
    //         }
    //         else {
                
    //             return emit * (1-mis)+ srec.attenuation * Li(scene, Ray3f(hit.p, srec.scattered), ++depth)   + L_dir;
    //         }

    //     }
    //     else return emit ;
        
    // }
    // else return scene.get_background();
    
    

    HitInfo hit;
    float rr_start;

    if (depth == max_depth){
        if (scene.intersect(ray, hit)){
            Vec3f emit = hit.mat->emitted(ray, hit);
            return emit;
        }
    }

    ScatterRecord srec;
    depth = 0;
    Color3f L(0.0f);
    Color3f throughput(1.0f);

    Ray3f traced_ray = ray;
    
    while (depth < max_depth){
        if(!scene.intersect(traced_ray, hit)){
            L += throughput * scene.get_background();
            break;
        }
        else{

            Color3f L_dir = Color3f(0.0f);
        
            Vec3f light_sample = scene.emitters().sample(hit.p);
            Vec3f emit = hit.mat->emitted(traced_ray, hit);
            bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);
            L += throughput * emit;

            float pdf_m = 0; 
            float pdf_e = scene.emitters().pdf(hit.p, light_sample);

            pdf_m = hit.mat->pdf(traced_ray.d, light_sample, hit);
            float mis_1;
            if (pdf_m < 0) mis_1 = 0; // delta pdf_m
            else mis_1 = pdf_e / (pdf_m + pdf_e);

            
            if(!hit.mat->isDelta && pdf_e > 0) {
                L_dir = hit.mat->eval(traced_ray.d, light_sample, hit) / pdf_e * Li(scene, Ray3f(hit.p, light_sample), max_depth) * mis_1;
                
            }
            
            L += L_dir * throughput;

            // if(depth == 0 && hit.mat->get_name() == "RoughDielectric") cout << sample_result << endl;

            if(sample_result){
                // if(hit.mat->get_name() == "RoughDielectric") cout << hit.mat->pdf(traced_ray.d, srec.scattered, hit) << endl;
                if (!hit.mat->isDelta && hit.mat->pdf(traced_ray.d, srec.scattered, hit) > 0){
                    pdf_e = scene.emitters().pdf(hit.p, srec.scattered);
                    pdf_m = hit.mat->pdf(traced_ray.d, srec.scattered, hit);
                    float mis_2 = pdf_m / (pdf_m + pdf_e);
                    throughput *= hit.mat->eval(traced_ray.d, srec.scattered, hit) / pdf_m * mis_2;

                    if(depth == 0 && hit.mat->get_name() == "RoughDielectric"){
                        // cout <<"pdf:" << pdf_m << endl;
                        // cout <<"eval:" <<  hit.mat->eval(traced_ray.d, srec.scattered, hit) << endl;
                        // cout << "throughput" << throughput << endl;
                    }
                }
                else throughput *= srec.attenuation * 1.0f; // delta pdf_m
        
            }
            traced_ray = Ray3f(hit.p, srec.scattered);

            // russian roulette

            // if(depth > rr_start){
            //     float rr_prob = luminance(throughput);
            //     if (randf() > rr_prob)
            //         break;
            //     else throughput /= rr_prob;
            // }
            
            
            depth++;
            
        }
        
    }
    return L;

}




void PPM::preprocess(const Scene &scene, int count)
{
    search_radius = GetSearchRadius(count);
    // generatePhotonMap(scene);

    // m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        // m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
    tree.clearTree();
		



	
    int stored_photons = 0;
    int stored_caustics = 0;
    int emitted_photons = 0;
    int n_lights = scene.getAllEmitters().GetSize();
    

    while (stored_photons < max_photon_count && stored_caustics < max_caustics){
        // First choose a light
        // const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
        
        Vec3f pos, dir;
        int rand_idx = rand() % n_lights;

        scene.getAllEmitters().GetSurface(rand_idx)->SampleFromEmit(pos, dir);
        Ray3f photon_ray = Ray3f(pos, dir);
        // Color3f photon_power = emitter->samplePhoton(photon_ray, sampler->next2D(), sampler->next2D(), sampler->next1D());
        Color3f photon_power = Color3f(200000.0f);
        photon_power = Color3f(15.0 * M_PI * 130*130 );

        bool isCaustics = false;
        
        if (luminance(photon_power) != 0.0f)
        {
            emitted_photons++;			// keep track of how many photons we shot to divide the contrib of all stored photons finally

            // start the trace of the photon
            // Intersection isect;
            HitInfo hit;
            int depth = 0;
            ScatterRecord srec;

            while (depth < max_depth && stored_photons < max_photon_count)
            {
                // Check for intersection
                // If not intersection, don't do anything
                if (!scene.intersect(photon_ray, hit)) break;
                
                
                // const BSDF* bsdf = isect.mesh->getBSDF();
                bool sample_result = hit.mat->sample(photon_ray.d, hit, srec);

                if (hit.mat->get_name() == "Dielectric" || hit.mat->get_name() == "RoughDielectric")
                    isCaustics = true;

                
                if (hit.mat->stored_photons && sample_result)
                {
                    // Store photon
                    // m_photonMap->push_back(Photon(hit.p, photon_ray.d, photon_power));
                    tree.addPhoton(Photon(hit.p, dir, photon_power));
                    // cout << "global:" << hit.p << endl;
                    if(isCaustics){
                        caustics_map.addPhoton(Photon(hit.p, dir, photon_power));
                        stored_caustics++;
                        // cout << hit.p << endl;
                    }
                        
                    stored_photons++;
                }

                // Now sample next direction
                // BSDFQueryRecord bRec(isect.toLocal(-photon_ray.d));
                // Color3f f = bsdf->sample(bRec, sampler->next2D());
                // cout << "b" << endl;
                if (hit.mat->pdf(photon_ray.d, srec.scattered, hit) > 0 && !srec.isSpecular){
                    Color3f f = hit.mat->eval(photon_ray.d, srec.scattered, hit) / hit.mat->pdf(photon_ray.d, srec.scattered, hit);
                    // cout << "?" << endl;
                
                // cout << "a" << endl;
                // Vector3f reflected_dir = isect.toWorld(bRec.wo);

                //update the photon power
                // pdf is included in the f term
                    Color3f incoming_power = photon_power;
                    // photon_power *= f * fabsf(Frame::cosTheta(bRec.wo));
                    
                    photon_power *= f;
                    
                    // Check for zero bsdf
                    if (luminance(f) == 0.0f)
                        break;
                    
                    // Check for russian roulette
                    if (depth > rr_start)
                    {
                        float p = 1.0f - luminance(photon_power) / luminance(incoming_power);
                        if (randf() < p)
                            break;
                        else photon_power /= (1.0f - p);
                    }
                    if (depth > max_depth ) break;
                    

                    // Generate the next bounce direction
                    photon_ray = Ray3f(hit.p, srec.scattered);
                    depth++;
                }
                else {
                    Color3f f = srec.attenuation;
                    Color3f incoming_power = photon_power;
                    // photon_power *= f * fabsf(Frame::cosTheta(bRec.wo));
                    
                    photon_power *= f;
                    
                    // Check for zero bsdf
                    if (luminance(f) == 0.0f)
                        break;
                    
                    // Check for russian roulette
                    if (depth > rr_start)
                    {
                        float p = 1.0f - luminance(photon_power) / luminance(incoming_power);
                        if (randf() < p)
                            break;
                        else photon_power /= (1.0f - p);
                    }
                    if (depth > max_depth ) break;
                    

                    // Generate the next bounce direction
                    photon_ray = Ray3f(hit.p, srec.scattered);
                    depth++;
                }
            }
        }
        
    }

    // Divide all photons by the total emitted
    stored = emitted_photons;
    // m_photonMap->scale(emitted_photons);

    /* Build the photon map */
    // m_photonMap->build();
    tree.buildTree();
    caustics_map.buildTree();
    cout << "Done building photon map" << endl;
}

float PPM::GetSearchRadius(const int iter)
{
    float radius_sqr = search_radius * search_radius;
    if (iter == 0) return radius_sqr;
    float alpha = 0.667; //decrease parameter
    float factor = (alpha + (float)(iter - 1.0f)) / (1.0f + (float)(iter - 1.0f));

    return factor * radius_sqr;
}

void PPM::generatePhotonMap(const Scene &scene)
{
    cout << "total " << max_photon_count << " photons per iteration" << endl;
    // while (!isFull){
    //     // Vec3f pos = Vec3f((randf() - 0.5f)*100.0f, (randf() - 0.5f)*100.0f, (randf() - 0.5f)*100.0f);
    //     // Vec3f dir = Vec3f(0.0f, 1.0f, 0.0f);
    //     Vec3f pos, dir;
    //     Vec3f power = Vec3f(10.0f);
    //     Ray3f photon_ray = Ray3f(pos, dir);

    //     // sample light sources propotional to the area size of light 
    //     int rand_idx = rand() % scene.getAllEmitters().GetSize();
    //     float prob = 1 / scene.getAllEmitters().GetSize();
    //     scene.getAllEmitters().GetSurface(rand_idx)->SampleFromEmit(pos, dir);
    //     power = power / prob / max_photon_count; // divide by number of photons
    //     tracePhoton(scene, pos, dir, power, 0);
    // }
    // tree.buildTree();

    


}

void PPM::tracePhoton(const Scene &scene, const Vec3f pos, const Vec3f dir, const Vec3f power, int bounces)
{
    // if (stored >= max_photon_count) {
    //     cout << "Photon Map is full" << endl;
    //     isFull = true;
    //     return;
	// }


    // Ray3f ray = Ray3f(pos, dir);
    // HitInfo hit;
    // ScatterRecord srec;

    // if(scene.intersect(ray, hit)){
        
    //     if(!hit.mat->isEmissive()){
    //         bool sample_result = hit.mat->sample(ray.d, hit, srec);
    //         Vec3f photon_power = power;
    //         if (sample_result){
    //             if (!srec.isSpecular && hit.mat->pdf(ray.d, srec.scattered, hit) != 0)
    //                 photon_power *= hit.mat->eval(ray.d, srec.scattered, hit) / hit.mat->pdf(ray.d, srec.scattered, hit) ;
    //             tree.addPhoton(Photon(hit.p, dir, power ));
    //             // visualize_map.push_back(hit.p);
    //             stored++;
    //         }
            
    //         if (bounces > rr_start){
    //             float p = min(luminance(photon_power) / luminance(power), 1.0f);
    //             if (randf() > p) 
    //                 return;
    //             else photon_power /= p;
    //         }
    //         // cout << photon_power << endl;

    //         tracePhoton(scene, hit.p, srec.scattered, photon_power, ++bounces);
    //     }
    // }
}

Color3f PPM::Li(const Scene &scene, const Ray3f &ray, int depth) const
{
    // HitInfo hit;

    

    // ScatterRecord srec;
    // depth = 0;
    // Color3f L(0.0f);
    // Color3f throughput(1.0f);

    // Ray3f traced_ray = ray;

    // while (depth < 10){
    //     if(!scene.intersect(traced_ray, hit)){
    //         L += throughput * scene.get_background();
    //         break;
    //     }
    //     else{
    //         Vec3f emit = hit.mat->emitted(traced_ray, hit);
    //         if (hit.mat->isEmissive()) {
    //             L += emit;
    //             break;
    //         }
        
    //         if(hit.mat->sample(traced_ray.d, hit, srec)){
    //             std::vector<const Photon *> photons;
    //             std::vector<float> distances;
                
    //             tree.nearestNeighbours(hit.p, photons, distances, search_radius);
    //             // cout << photons.size() << endl;
    //             // for (int i = 0; i < photons.size(); ++i)
    //             //     cout << "Photon " << (i + 1) << ": Pos " << photons[i]->pos << " Distance to query point: " << distances[i] << endl;
    
    //             if(photons.size() > 0){
    //                 for (int i = 0 ; i < photons.size() ; i++){
    //                     Vec3f photon_dir = photons[i]->dir;
    //                     float area = M_PI * search_radius * search_radius;
    //                     // cout << throughput * hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit) * photons[i]->power / area << endl;
    //                     // cout << hit.mat->eval(traced_ray.d, -photon_dir, hit) / hit.mat->pdf(traced_ray.d, -photon_dir, hit) << endl;
    //                     L += throughput * photons[i]->power  / area;
    //                 }
                    
    //             }

    //         }
    //         traced_ray = Ray3f(hit.p, srec.scattered);
    //         depth++;
    //     }
    // }
    // return L;

    
    // while (depth < max_depth){
    //     if(!scene.intersect(traced_ray, hit)){
    //         L += throughput * scene.get_background();
    //         break;
    //     }
    //     else{
            
    //         Vec3f emit = hit.mat->emitted(traced_ray, hit);

    //         //no bounce off from light source
            // if (hit.mat->isEmissive()) {
            //     L += emit;
            //     break;
            // }

    //         // cout << hit.p << endl;
    //         bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);
    //         // L += throughput * emit;

            
    //         // only count for non specular surface
    //         if(sample_result){
    //             if (!srec.isSpecular && hit.mat->pdf(traced_ray.d, srec.scattered, hit) > 0){

    //                 std::vector<const Photon *> photons;
    //                 std::vector<float> distances;
                    
    //                 tree.nearestNeighbours(hit.p, photons, distances, search_radius);
                    // cout << photons.size() << endl;
                    // for (int i = 0; i < photons.size(); ++i)
                    //     cout << "Photon " << (i + 1) << ": Pos " << photons[i]->pos << " Distance to query point: " << distances[i] << endl;

    //                 if(photons.size() > 0){
    //                     for (int i = 0 ; i < photons.size() ; i++){
    //                         Vec3f photon_dir = photons[i]->dir;
    //                         float area = M_PI * search_radius * search_radius;
    //                         // cout << throughput * hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit) * photons[i]->power / area << endl;
    //                         // cout << hit.mat->eval(traced_ray.d, photon_dir, hit) / hit.mat->pdf(traced_ray.d, photon_dir, hit) << endl;
    //                         L += throughput * hit.mat->eval(traced_ray.d,  -photon_dir, hit) / hit.mat->pdf(traced_ray.d,  -photon_dir, hit) * photons[i]->power / area;
    //                     }
    //                     // break;
    //                 }

    //             }
        
    //         }

    //         if(!srec.isSpecular && hit.mat->pdf(traced_ray.d, srec.scattered, hit) > 0)
    //             throughput *= hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit);
    //         else throughput *= srec.attenuation;
            
    //         traced_ray = Ray3f(hit.p, srec.scattered);


    //         // russian roulette

    //         // if(depth > rr_start){
    //         //     float rr_prob = luminance(throughput);
    //         //     if (randf() > rr_prob)
    //         //         break;
    //         //     else throughput /= rr_prob;
    //         // }
            
    //         depth++;
            
    //     }
        
    // }
    // return L;

    HitInfo hit;
    max_depth = 4;
    

    ScatterRecord srec;
    // depth = 0;
    Color3f L(0.0f);
    Color3f L_dir(0.0f);
    Color3f throughput(1.0f);

    Ray3f traced_ray = ray;
    // cout << depth << endl;
    // cout << max_depth << endl;
    if (depth == max_depth){
        // cout << "?" << endl;
        if (scene.intersect(ray, hit)){
            // cout << "lol" << endl;
            Vec3f emit = hit.mat->emitted(ray, hit);
            return emit;
        }
    }
    
    if (!final_gather){
        while (depth < max_depth){
            if(!scene.intersect(traced_ray, hit)){
                return scene.get_background();
                
            }


            if (hit.mat->isEmissive()) {
                Vec3f emit = hit.mat->emitted(traced_ray, hit);
                L += emit * throughput;
                break;
            }
            // const BSDF* bsdf = isect.mesh->getBSDF();
            bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);

            Vec3f light_sample = scene.emitters().sample(hit.p);
                
            float pdf_e = scene.emitters().pdf(hit.p, light_sample);
            if(hit.mat->stored_photons && pdf_e > 0 && !hit.mat->isDelta) {
                // cout << depth << endl;
                L_dir = hit.mat->eval(traced_ray.d, light_sample, hit) / pdf_e * Li(scene, Ray3f(hit.p, light_sample), max_depth);
            }

            L += L_dir * throughput;

            // const BSDF* bsdf = isect.mesh->getBSDF();
            // bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);

            // Compute indirect contribution only from diffuse surfaces
            if (hit.mat->stored_photons){
                // std::vector<uint32_t> results;
                std::vector<const Photon *> photons;
                std::vector<float> distances;
                // m_photonMap->search(isect.p, m_photonRadius, results);
                tree.nearestNeighbours(hit.p, photons, distances, 30, search_radius);
                float area = M_PI * search_radius * search_radius;
                
                float max_dis = distances[0];
                for (int i = 0 ; i < distances.size() ; i++){
                    if(max_dis < distances[i])
                        max_dis = distances[i];
                }
                area = M_PI * max_dis * max_dis;

                // The uint32_t makes the size() - 1 wrap around. Subtle bug.
                if (photons.size() > 0)
                {
                    for (uint32_t i = 0; i < photons.size() - 1; i++)
                    {
                        // const Photon &photon = (*m_photonMap)[results[i]];
                        Vec3f photon_dir = photons[i]->dir;

                        // Compute the integral equation
                        // BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d), isect.toLocal(-photon.getDirection()), ESolidAngle);
                        float pdf_m = hit.mat->pdf(traced_ray.d, -photon_dir, hit);
                        if (pdf_m >= 0){
                            Color3f f = hit.mat->eval(traced_ray.d, -photon_dir, hit) / pdf_m;
                            if (pdf_m == 0.0f && luminance(hit.mat->eval(traced_ray.d, -photon_dir, hit)) == 0.0f) 
                                f = hit.mat->get_albedo(hit);
                            // cout << throughput * f * photons[i]->power / area  << endl;
                            L += throughput * f * photons[i]->power / area / stored;
                        }
                            
                    }
                }
                
                break;
            }

            // Sample a reflection ray
            // BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
            // Color3f f = bsdf->sample(bRec, sampler->next2D());
            // Vector3f reflected_dir = isect.toWorld(bRec.wo);
            Color3f f;
            // float cos_theta = fabsf(Frame::cosTheta(bRec.wo));
            if (hit.mat->pdf(traced_ray.d, srec.scattered, hit) > 0 && !srec.isSpecular){
                f = hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit);

                throughput *= f;
            }
            else throughput *= srec.attenuation;

            // Check if we've fa
            if (luminance(throughput) == 0.0f)
                break;

            // Check for russian roulette
            // if (depth > rr_start)
            // {
            //     if (randf() < 0.5f)
            //         break;
            //     else throughput *= 2.0f;
            // }
            // else if (depth > max_depth)
            // {
            //     // forcibly terminate
            //     break;
            // }

            // Propogate
            traced_ray = Ray3f(hit.p, srec.scattered);
            depth++;
        }
        return L;
    }

    if (final_gather){

        while(depth < max_depth){
            if(!scene.intersect(traced_ray, hit)){
                return scene.get_background();
            }

            if (hit.mat->isEmissive()) {
                Vec3f emit = hit.mat->emitted(traced_ray, hit);
                L += emit * throughput;
                break;
            }
            
            Vec3f light_sample = scene.emitters().sample(hit.p);
            bool sample_result = hit.mat->sample(traced_ray.d, hit, srec);

                        
            float pdf_e = scene.emitters().pdf(hit.p, light_sample);
            if(hit.mat->stored_photons && pdf_e > 0 && !hit.mat->isDelta) {
                L_dir = hit.mat->eval(traced_ray.d, light_sample, hit) / pdf_e * Li(scene, Ray3f(hit.p, light_sample), max_depth);
            }

            L += L_dir * throughput;


            if (srec.isSpecular){
                Color3f f;
                if (hit.mat->pdf(traced_ray.d, srec.scattered, hit) > 0 && !srec.isSpecular){
                    f = hit.mat->eval(traced_ray.d, srec.scattered, hit) / hit.mat->pdf(traced_ray.d, srec.scattered, hit);

                    throughput *= f;
                }
                else throughput *= srec.attenuation;

                // Check if we've fa
                if (luminance(throughput) == 0.0f)
                    break;

                // Propogate
                traced_ray = Ray3f(hit.p, srec.scattered);
                depth++;
            }
            else{
                for(int n = 1 ; n < FG_num ; n++){
                    HitInfo tmp_hit = hit;
                    auto factor = throughput;
                    
                    // Compute indirect contribution only from diffuse surfaces
                    int max_bounces = 20;
                    //second bounce
                    int count_bounce = 0;

                    std::vector<const Photon *> photons;
                    std::vector<float> distances;

                    caustics_map.nearestNeighbours(tmp_hit.p, photons, distances, 30, search_radius);
                    float max_dis = distances[0];
                    for (int i = 0 ; i < distances.size() ; i++){
                        if(max_dis < distances[i])
                            max_dis = distances[i];
                    }
                    float area = M_PI * max_dis * max_dis;
                    if(photons.size() == 0) cout << "???" << endl;
                    if (photons.size() > 0){
                        for (uint32_t i = 0; i < photons.size() - 1; i++){
                            // const Photon &photon = (*m_photonMap)[results[i]];
                            Vec3f photon_dir = photons[i]->dir;

                            // Compute the integral equation
                            
                            float pdf_m = tmp_hit.mat->pdf(traced_ray.d, -photon_dir, tmp_hit);
                            // cout << pdf_m << endl;

                            if(pdf_m < 0) cout << srec.isSpecular << endl;

                            if (pdf_m >= 0){
                                // cout << "a" << endl;
                                Color3f f = tmp_hit.mat->eval(traced_ray.d, -photon_dir, tmp_hit) / pdf_m;
                                // cout << "b" << endl;
                                if (pdf_m == 0.0f && luminance(tmp_hit.mat->eval(traced_ray.d, -photon_dir, tmp_hit)) == 0.0f) 
                                    f = tmp_hit.mat->get_albedo(hit);
                                // cout << throughput * f * photons[i]->power / area / stored / FG_num  << endl;
                                L += throughput * f * photons[i]->power / area / stored / FG_num;
                            }
                                
                        }
                    }


                    while(true){
                        sample_result = tmp_hit.mat->sample(traced_ray.d, tmp_hit, srec); // sample from current hit point
                        
                        if (sample_result){
                            HitInfo prev_hit = tmp_hit;
                            Vec3f _in = traced_ray.d;
                            Vec3f _out = srec.scattered;
                            auto attenuation = srec.attenuation;
                            // update next hit point
                            
                            if (scene.intersect(Ray3f(tmp_hit.p, srec.scattered), tmp_hit)){


                                auto pdf_m = prev_hit.mat->pdf(_in, _out, prev_hit);
                                traced_ray.d = srec.scattered; // next hit point incidance ray
                                // bool indirect_sample = tmp_hit.mat->sample(srec.scattered, tmp_hit, srec);

                                // throughput *= prev_hit.mat->eval(traced_ray.d, srec.scattered, prev_hit) / pdf_m;
                                if(prev_hit.mat->isDelta)
                                    throughput *= attenuation;
                                else{
                                    if(pdf_m == 0.0f && luminance(prev_hit.mat->eval(_in, _out, prev_hit)) < 1e-6f) {
                                        throughput *= hit.mat->get_albedo(hit);
                                    }
                                    else {
                                        throughput *= (prev_hit.mat->eval(_in, _out, prev_hit) / pdf_m);                                    
                                    }
                                }
                                if(tmp_hit.mat->get_name() == "Dielectric" || hit.mat->get_name() == "RoughDielectric") throughput *= 0;
                                
                                if(tmp_hit.mat->stored_photons){
                                    break;
                                }
                                
                            }
                        }
                        if(count_bounce > max_bounces){
                            Vec3f emit = tmp_hit.mat->emitted(ray, hit);
                            return emit;
                        }
                        count_bounce++;

                        
                        // sample_result = hit.mat->sample(traced_ray.d, tmp_hit, srec);

                        // Ray3f shadowRay = Ray3f(hit.p, srec.scattered);
                        // auto pdf_m = hit.mat->pdf(traced_ray.d, srec.scattered, hit);
                        // // auto eval = hit.mat->eval(traced_ray.d, srec.scattered, hit);
                        // // cout << pdf_m << endl;
                        // // cout << eval << endl;

                        // if (scene.intersect(shadowRay, tmp_hit)){
                        //     bool indirect_sample = tmp_hit.mat->sample(shadowRay.d, tmp_hit, srec);
                        //     if(tmp_hit.mat->stored_photons) {
                                
                        //         if(pdf_m == 0.0f && luminance(hit.mat->eval(traced_ray.d, shadowRay.d, hit)) < 1e-6f) {
                        //             throughput *= hit.mat->get_albedo(hit);
                        //             // cout << "c: " << throughput << endl;
                        //         }
                        //         else {

                        //             // if(pdf_m == 0.0f && luminance(hit.mat->eval(traced_ray.d, shadowRay.d, hit)) == 0.0f)
                        //             //     cout << "fuck" << endl;

                        //             throughput *= (hit.mat->eval(traced_ray.d, shadowRay.d, hit) / pdf_m);                                    
                        //         }
                        //         if(luminance(throughput) == 0.0f)
                        //             cout << "?? " << endl;
                        //         break;
                        //     }
                        // }

                        // if(count_bounce > max_bounces){
                        //     Vec3f emit = hit.mat->emitted(ray, hit);
                        //     return emit;
                        // }
                        // count_bounce++;
                        
                    }
                    // cout << "????:"  <<  throughput << endl;

                    Vec3f tmp_traced_ray = traced_ray.d;

                    // cout << "after:" << tmp_hit.p << endl;

                    photons.clear();
                    distances.clear();
                    
                    tree.nearestNeighbours(tmp_hit.p, photons, distances, 30, search_radius);
                    // float area = M_PI * search_radius * search_radius;
                    
                    max_dis = distances[0];
                    for (int i = 0 ; i < distances.size() ; i++){
                        if(max_dis < distances[i])
                            max_dis = distances[i];
                    }
                    area = M_PI * max_dis * max_dis;

                    // The uint32_t makes the size() - 1 wrap around. Subtle bug.
                    if (photons.size() > 0){
                        for (uint32_t i = 0; i < photons.size() - 1; i++){
                            // const Photon &photon = (*m_photonMap)[results[i]];
                            Vec3f photon_dir = photons[i]->dir;

                            // Compute the integral equation
                            
                            float pdf_m = tmp_hit.mat->pdf(tmp_traced_ray, -photon_dir, tmp_hit);
                            // cout << pdf_m << endl;

                            if(pdf_m < 0) cout << srec.isSpecular << endl;

                            if (pdf_m >= 0){
                                // cout << "a" << endl;
                                Color3f f = tmp_hit.mat->eval(tmp_traced_ray, -photon_dir, tmp_hit) / pdf_m;
                                // cout << "b" << endl;
                                if (pdf_m == 0.0f && luminance(tmp_hit.mat->eval(tmp_traced_ray, -photon_dir, tmp_hit)) == 0.0f) 
                                    f = tmp_hit.mat->get_albedo(hit);
                                // cout << throughput * f * photons[i]->power / area / stored / FG_num  << endl;
                                L += throughput * f * photons[i]->power / area / stored / FG_num;
                            }
                                
                        }
                    }

                    //caustics map calculation
                    

                    
                    
                }
                break;
            }
        }

        return L;
    }

    
}