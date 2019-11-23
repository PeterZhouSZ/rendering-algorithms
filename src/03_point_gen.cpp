#include<dirt/common.h>
#include<ctime>
using namespace std;

int main(int argc, char ** argv)
{
    
    srand(time(0));
    
    cout << "x,y,z" << endl;
    cout << "0,0,0" <<endl;
    
    for (int i = 0 ; i < 500 ; i++){
        Vec3f p = randomInUnitSphere();
        cout << p.x << "," << p.y << "," << p.z << "," << endl;
        
    }

    return 0;
}