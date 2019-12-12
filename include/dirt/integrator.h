#pragma once

#include <dirt/fwd.h>
#include <dirt/parser.h>
#include <stdlib.h>
#include <dirt/kdtree.h>


class Integrator
{
    public:
        Integrator(const json & j = json::object()) {max_depth = j.value("max_bounces", 64);};
        virtual Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const { return Color3f(1.0f, 0.0f, 1.0f);}
        virtual void preprocess(const Scene &scene, int count){};
        int iter;

    private:
        int max_depth;
        
};


class NormalIntegrator: public Integrator
{
    public:
        NormalIntegrator(const json & j = json::object()) {max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
    protected:
        int max_depth;
};

class AmbientOcclusionIntegrator: public Integrator
{
    public:
        AmbientOcclusionIntegrator(const json & j = json::object()) {max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;

    protected:
        int max_depth;
};

class PathTracerMaterials: public Integrator
{
    public:
        PathTracerMaterials(const json & j = json::object()) {max_depth = j.value("max_bounces", 64);};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;

    protected:
        int max_depth;

};

class DirectMats: public Integrator
{
    public:
        DirectMats(const json & j = json::object()) {max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
    protected:
        int max_depth;
};

class DirectNEE: public Integrator
{
    public:
        DirectNEE(const json & j = json::object()) { max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
    protected:
        int max_depth;
};

class DirectMIS: public Integrator
{
    public:
        DirectMIS(const json & j = json::object()) { max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;

    protected:
        int max_depth;
};

class PathTracerNEE : public Integrator
{
    public:
        PathTracerNEE(const json & j = json::object()) { max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth, bool include) const;
    protected:
        int max_depth;
};

class PathTracerMis: public Integrator
{
    public:
        PathTracerMis(const json & j = json::object()) { max_depth = j.value("max_bounces", 64); iter = 1;};
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
        // Color3f Li(const Scene &scene, const Ray3f &ray, int depth, float Le_weight) const ;
    protected:
        int max_depth;
};

class PPM: public Integrator
{
    public:
        PPM(const json & j = json::object()) { 
            max_depth = j.value("max_bounces", 64);
            search_radius = j.value("search_radius", 10.0f);
            iter = j.value("iter", 1);
            final_gather = j.value("FG", false);
            FG_num = j.value("FG_num", 2);
        };
        Color3f Li(const Scene &scene, const Ray3f &ray, int depth) const override;
        void tracePhoton(const Scene &scene, const Vec3f pos, const Vec3f dir, const Vec3f power, int bounces);
        void generatePhotonMap(const Scene &scene);
        float GetSearchRadius(const int iter);
        void preprocess(const Scene &scene, int count);

    protected:
        mutable int max_depth; // for photon tracing
        float search_radius;
        int max_photon_count = 200000; //max for number in photonMap
        int max_caustics = 50000;
        // vector <const Photon *> photonMap;
        KdTree tree;
        KdTree caustics_map;
        bool final_gather;
        int FG_num;
        int rr_start = 5.0;
        bool isFull = false;
        int stored = 0;
        vector<Vec3f> visualize_map;
};