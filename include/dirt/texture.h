#pragma once

#include <dirt/fwd.h>

#include <dirt/surface.h>
#include <dirt/transform.h>
#include <dirt/perlin.h>
#include <dirt/image.h>
#include <stdlib.h>


class Texture{
    public:
        virtual Vec3f value (const HitInfo& hit) const = 0;
    private:
};

class ConstantTexture : public Texture {
    public:
        ConstantTexture() { }
        ConstantTexture(const json &j) : color(j.value("color", color)) { }
        Color3f value(const HitInfo& hit) const {
            return color;
        }
        Color3f color;
    private:
        
};

class CheckerTexture : public Texture{
    public:
        CheckerTexture() {}
        CheckerTexture(const json &j){
            json tmp;
            tmp["color"] = j["even"];
            
            even = make_shared<ConstantTexture>(tmp);

            tmp["color"] = j["odd"];
            
            odd = make_shared<ConstantTexture>(tmp);
            
            scale = j.value("scale", scale);
        }
        Color3f value(const HitInfo& hit) const;
        
        shared_ptr<ConstantTexture> even;
        shared_ptr<ConstantTexture> odd;
        protected:
            float scale = 1.0f;
};

class MarbleTexture : public Texture{
    public:
        MarbleTexture() {}
        MarbleTexture(const json &j) {
            json tmp;
            tmp["color"] = j["veins"];
            
            veins = make_shared<ConstantTexture>(tmp);

            tmp["color"] = j["base"];
            
            base = make_shared<ConstantTexture>(tmp);
            
            scale = j.value("scale", scale);
        }
        Color3f value(const HitInfo& hit) const;

        perlin noise;
        shared_ptr<ConstantTexture> veins;
        shared_ptr<ConstantTexture> base;
    protected:
        float scale = 1.0f;
};


class ImageTexture : public Texture {
    public:
        ImageTexture(){}
        ImageTexture(const json &j) ;
        Color3f value(const HitInfo& hit) const;
        Image3f img;
    protected:
        string path;
        
};