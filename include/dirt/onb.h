#pragma once


#include <stdlib.h>
#include <dirt/fwd.h>
#include <dirt/vec.h>
#include <dirt/common.h>

class onb
{
    public:
        onb() {}
        inline Vec3f operator[](int i) const { return axis[i]; }
        Vec3f u() const { return axis[0]; }
        Vec3f v() const { return axis[1]; }
        Vec3f w() const { return axis[2]; }
        Vec3f local(float a, float b, float c) const { return a*u() + b*v() + c*w(); }
        Vec3f local(const Vec3f& a) const { return a.x * u() + a.y * v() + a.z * w(); }
        void build_from_w(const Vec3f&);
        
    private:
        Vec3f axis[3];
};

