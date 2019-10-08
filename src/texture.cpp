#include <dirt/texture.h>
#include <dirt/parser.h>
#include <dirt/surface.h>
#include <filesystem/resolver.h>
#include <dirt/image.h>

Color3f CheckerTexture::value(const HitInfo& hit) const{ 
    // hit nan
    float sines = sin(scale*(hit.p.x)) * sin(scale*(hit.p.y)) * sin(scale*(hit.p.z));
    //cout << hit.sn << endl;
    if (sines < 0){
        return odd->value(hit);
    }
        
    else {
        return even->value(hit);
    }
            
}

Color3f MarbleTexture::value(const HitInfo& hit) const {

    return lerp(base->value(hit), veins->value(hit), 0.5 * (1 + sin(scale*hit.p.z + 10*noise.turb(hit.p))));
}

ImageTexture::ImageTexture(const json &j){
    string filename = j["filename"];
    filesystem::path path(filename);
    string par_path = getFileResolver().resolve(filename).str();

    
    if(img.load(par_path))
        cout << "image loaded successfully" << endl;
    else cout << "fail to load image" << endl;
}

Color3f ImageTexture::value(const HitInfo& hit) const {
    
    int i = hit.uv.x * img.width();
    int j = (1 - hit.uv.y) * img.height() - 0.001;
    // clamp
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > img.width() - 1) i = img.width() - 1;
    if (j > img.height() - 1) j = img.height() - 1;

    return img(i, j);
}
