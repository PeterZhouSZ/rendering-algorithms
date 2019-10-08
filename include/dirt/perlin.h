#pragma once

#include <stdlib.h>
#include <dirt/vec.h>
#include <dirt/common.h>


inline float perlin_interp(Vec3f c[2][2][2], float u, float v, float w) {
    float uu = u*u*(3-2*u);
    float vv = v*v*(3-2*v);
    float ww = w*w*(3-2*w);
    float accum = 0;
    for (int i=0; i < 2; i++){
        for (int j=0; j < 2; j++){
            for (int k=0; k < 2; k++) {
                Vec3f weight_v(u-i, v-j, w-k);
                accum += (i*uu + (1-i)*(1-uu)) * (j*vv + (1-j)*(1-vv))* (k*ww + (1-k)*(1-ww))*dot(c[i][j][k], weight_v);
            }
        }
    }
        
    return accum;
}



static Vec3f* perlin_generate(){
    Vec3f * p = new Vec3f[256];
    for (int i = 0 ; i < 256 ; ++i){
        double x_random = 2*drand48() - 1;
        double y_random = 2*drand48() - 1;
        double z_random = 2*drand48() - 1;
        p[i] = normalize(Vec3f(x_random, y_random, z_random));
    }
        
    return p;
}

static void permute(int *p, const int n){
    for(int i = n-1 ; i > 0 ; i--){
        int target = int(drand48() * (i+1));
        int tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
}

static int* perlin_generate_permute(){
    int *p = new int[256];
    for (int i = 0 ; i < 256 ; ++i)
        p[i] = i;
    permute(p, 256);
    return p;
}

class perlin {
    public:
        float noise(const Vec3f &p) const{
            float u = p.x - floor(p.x);
            float v = p.y - floor(p.y);
            float w = p.z - floor(p.z);

            int i = floor(p.x);
            int j = floor(p.y);
            int k = floor(p.z);
            Vec3f c[2][2][2];
            for (int di=0; di < 2; di++)
                for (int dj=0; dj < 2; dj++)
                    for (int dk=0; dk < 2; dk++)
                        c[di][dj][dk] = ranvec[perm_x[(i+di) & 255] ^ perm_y[(j+dj) & 255] ^ perm_z[(k+dk) & 255]];
            return perlin_interp(c, u, v, w);
        }

        float turb(const Vec3f& p, int depth=7) const {
            float accum = 0;
            Vec3f tmp_p = p;
            float weight = 1.0;
            for (int i = 0; i < depth; i++) {
                accum += weight*noise(tmp_p);
                weight /= 2;
                tmp_p *= 2;
            }
            return abs(accum);
        }

        Vec3f *ranvec = perlin_generate();;
        int *perm_x = perlin_generate_permute();
        int *perm_y = perlin_generate_permute();
        int *perm_z = perlin_generate_permute();
    private:
        
};



// Vec3f *perlin::ranvec = perlin_generate();
// int *perm_x = perlin_generate_permute();
// int *perm_y = perlin_generate_permute();
// int *perm_z = perlin_generate_permute();
