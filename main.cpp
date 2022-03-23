//
//  main.cpp
//  SistemaSolar
//
//  Created by Luis Sánchez Travesí on 21/3/22.
//

#include <iostream>
#include <cmath>
using namespace std;

#define G  6.67384*pow(10, -11)
#define c  1.496*pow(10, 11)
#define Ms  1.99*pow(10, 30)



double mod(float x, float y);

int main() {
    //Declaramos las variables
    float h = 0.5, t = 0;
    float m[10], r[10][2], v[10][2], a[10][2];
    
    //Reescalamiento
    for (int i=0; i<10; i++) {
        r[i][0]*=1./c;
        r[i][1]*=1./c;
        m[i]*=1./Ms;
    }
    
    
    for (t=0; t<=100; t+=h) {
        for (int i = 0; i<=9; i++) {
            for (int j=0; j<=9; j++) {
                if (i!=j) {
                    a[i][0] = -m[j]*(r[i][0] - r[j][0])/pow(mod(r[i][0] - r[j][0], r[i][1] - r[j][1]), 3);
                    a[i][1] = -m[j]*(r[i][1] - r[j][1])/pow(mod(r[i][0] - r[j][0], r[i][1] - r[j][1]), 3);
                }
                
                r[i][0]+= h*v[i][0] + h*h/2*a[i][0];
                r[i][1]+= h*v[i][1] + h*h/2*a[i][1];
                v[i][0]+= h/2*
                
                
            }
        }
    }
    
    
    
    return 0;
}



double mod(double x, double y){
    return sqrt(x*x + y*y);
}
