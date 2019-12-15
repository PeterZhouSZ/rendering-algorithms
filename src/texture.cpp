#include <dirt/texture.h>
#include <dirt/parser.h>
#include <dirt/surface.h>
#include <filesystem/resolver.h>
#include <dirt/image.h>


Vec3f Texture::bump_N(const HitInfo &hit) const
{
	if(!_bump) 
		return hit.sn;

	
	int i = hit.uv.x * content.width();
    int j = (1 - hit.uv.y) * content.height() - 0.001;

	int i1 = i - 1;
	int i2 = i + 1;
	int j1 = j - 1;
	int j2 = j + 1;

	if (i1 < 0) i1 = 0;
	if (i1 > content.width()-1) i1 = content.width() - 1;
	if (i2 < 0) i2 = 0;
	if (i2 > content.width()-1) i2 = content.width() - 1;

	if (j1 < 0) j1 = 0;
	if (j1 > content.height()-1) j1 = content.height() - 1;
	if (j2 < 0) j2 = 0;
	if (j2 > content.height()-1) j2 = content.height() - 1;

    // clamp
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > content.width() - 1) i = content.width() - 1;
    if (j > content.height() - 1) j = content.height() - 1;

	float u_gradient = luminance(content(i1, j)) - luminance(content(i2, j));
	float v_gradient = luminance(content(i, j1)) - luminance(content(i, j2));

	onb uvw;
	uvw.build_from_w(hit.sn);
	Vec3f tangent = normalize(uvw.local(Vec3f(0, 0, 1)) - Vec3f(u_gradient, v_gradient, 0));
	return hit.sn + b_offset * u_gradient * tangent + b_offset * v_gradient * (-hit.sn);

    // return hit.sn;
}

Color3f ConstantTexture::value(const HitInfo& hit) const{
    return color;
}

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

    _bump = j.value("bump", false);
    b_offset = j.value("b_offset", 1.0);
    if(_bump) {
        // cout << "yeap" << endl;
        
    }
}

Color3f ImageTexture::value(const HitInfo& hit) const {
    
    // int i = hit.uv.x * img.width();
    // int j = (1 - hit.uv.y) * img.height() - 0.001;



    int i = (hit.uv.x *4 - int(hit.uv.x*4 )) * img.width();
    int j = ((1 - hit.uv.y)*4  - int((1 - hit.uv.y)*4) ) * img.height();
    // clamp
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > img.width() - 1) i = img.width() - 1;
    if (j > img.height() - 1) j = img.height() - 1;

    return img(i, j);
}
