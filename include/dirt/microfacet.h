#pragma once

#include <dirt/fwd.h>
#include <dirt/parser.h>
#include <dirt/vec.h>
#include <stdlib.h>
#include <dirt/surface.h>
#include <dirt/onb.h>


class Microfacet
{
    public:
        Microfacet() {};
        // sample w_h
        virtual Vec3f sample(const float alpha, const Vec3f& N) const {return Vec3f(0.0f);};
        // distribution
        virtual float D(const Vec3f& m, const HitInfo & hit, const float alpha) const {return 0;};
        virtual float pdf(const Vec3f &m, const HitInfo & hit, const float alpha) const {return 0;};
        // shadowing term
        virtual float G(const Vec3f& v, const Vec3f& o, const Vec3f& m, const HitInfo & hit, const float alpha) const {return 0;};
        virtual float G1(const Vec3f& v, const Vec3f& m, const HitInfo & hit, const float alpha) const {return 0;};

    private:
};

class GGX : public Microfacet
{
    public:
        GGX() {};
        Vec3f sample(const float alpha, const Vec3f& N) const override;
        float D(const Vec3f& m, const HitInfo & hit, const float alpha) const override;
        float pdf(const Vec3f &m, const HitInfo & hit, const float alpha) const override;
        float G(const Vec3f& v, const Vec3f& o, const Vec3f& m, const HitInfo & hit, const float alpha) const override;
        // m is microface normal
        float G1(const Vec3f& v, const Vec3f& m, const HitInfo & hit, const float alpha) const override; 
        
    protected:
        
};