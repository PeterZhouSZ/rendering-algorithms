#pragma once

#include <dirt/fwd.h>

#include <dirt/surface.h>
#include <dirt/transform.h>
#include <dirt/perlin.h>
#include <dirt/image.h>
#include <stdlib.h>
#include <filesystem/resolver.h>


class Texture{
    public:
        virtual Vec3f value (const HitInfo& hit) const = 0;
        bool _bump;
        float b_offset;
        Image3f content; // for bumpmap

        Vec3f bump_N(const HitInfo& hit) const;
    private:
};

class ConstantTexture : public Texture {
    public:
        ConstantTexture() { }
        ConstantTexture(const json &j){
            color = j.value("color", color);
            _bump = j.value("bump", false);
            b_offset = j.value("b_offset", 1.0);

            // if(_bump){
            //     bool result = content.load(j.value("bump_file", "../scenes/04_mc_integration2/Normal.png"));
            //     cout << b_offset <<endl;
            //     cout << "result:" << result << endl;
            // } 
            if(_bump){
                string b_filename = j.value("bump_file", "./Normal.png");
                filesystem::path path(b_filename);
                string par_path = getFileResolver().resolve(b_filename).str();
                content.load(par_path);
            } 
        }
        Color3f value(const HitInfo& hit) const;
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
            _bump = j.value("bump", false);
            b_offset = j.value("b_offset", 1.0);

            

            if(_bump){
                string b_filename = j["filename"];
                filesystem::path path(b_filename);
                string par_path = getFileResolver().resolve(b_filename).str();
                content.load(j.value("bump_file", ""));
            } 
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
            _bump = j.value("bump", false);
            b_offset = j.value("b_offset", 1.0);
            if(_bump) content.load(j.value("bump_file", ""));
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