#include <stdlib.h>
#include <dirt/onb.h>

void onb::build_from_w(const Vec3f& n) {
    axis[2] = normalize(n);
    Vec3f a;
    if (abs(w().x) > 0.9)
        a = Vec3f(0, 1, 0);
    else
        a = Vec3f(1, 0, 0);
    axis[1] = normalize( cross( w(), a ) );
    axis[0] = cross(w(), v());
}