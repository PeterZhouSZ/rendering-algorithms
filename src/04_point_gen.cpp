#include<dirt/common.h>
#include<ctime>
using namespace std;

int main(int argc, char ** argv)
{
    
    srand(time(0));
    
    cout << "x,y,z" << endl;
    // cout << "0,0,0" <<endl;
    Vec3f center = Vec3f(0, 1, 0);
    float radius = 3.0f;
    Vec3f v0 = Vec3f(0,0,0);
    Vec3f v1 = Vec3f(1,1,1);
    Vec3f v2 = Vec3f(0.5, 1, 2);

    float cosThetaMax = cos(M_PI / 4);
    cout << cosThetaMax << endl;

    for (int i = 0 ; i < 500 ; i++){
        Vec3f p = randomSphericalCap( cosThetaMax);
        cout << p.x << "," << p.y << "," << p.z << "," << endl;
        
    }

    return 0;
}