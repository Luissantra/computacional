//
//  Simulacion del sistema solar mediante el algoritmo de Verlet
//
//
//
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


//Definimos las constantes
#define c  147.6e9 //Distancia de la Tierra al Sol (10^9m)
#define Ms  1.989e30 //Masa del sol (10e30kg)
#define G 6.6738e-11 //gravitacion universal
#define cte 3.35695e-5 //sqrt(c/GMs)
#define h 0.05 //amplitud de paso

double mod(double x, double y); //devuelve el modulo de un vector v=(x, y)
void iniciarDatos(double m[], double r[10][2], double v[10][2]); //lee los datos de un fichero y los guarda en m, r, v además de reescalarlos
void verlet(int N, double m[], double r[10][2], double v[10][2], double a[10][2], double LE[10][2], bool geo); //realiza el algoritmo n pasos

int main() {
    // Variables
    double t; // Tiempo de la simulacion en dias
    int n = 0; // Numero de iteraciones
    double m[10], r[10][2], v[10][2], a[10][2], LE[10][2];
    char ans;
    bool geo = false;
    
    cout << "#################" << endl;
    cout << "# SISTEMA SOLAR #" << endl;
    cout << "#################" << endl << endl;
    cout << "Introduzca el numero de dias a simular: "<<endl;
    cin >> t;
    cout << "Modelo geocentrico (y/n): " <<endl;
    cin >> ans;
    if (ans == 'y') {
        geo = true;
    }

    n = t/h/58.1 + 3; // convertimos los dias a numero de pasos, +3 para no quedar cortos
    
    iniciarDatos(m, r, v);
    verlet(n, m, r, v, a, LE, geo);
    
    return 0;
}

void verlet(int N, double m[], double r[10][2], double v[10][2], double a[10][2], double LE[10][2], bool geo){
    
    double sumax = 0, sumay = 0, sumaL = 0, sumaE =  0;
    double U[10], w[10][2];
    bool T[10];
    
    
    
    //iniciamos a 0 para evitar posibles problemas
    for (int i=0; i<10; i++) {
        LE[i][0] = LE[i][1] = U[i] = w[i][0] = w[i][1] = 0;
    }
    for (int i=0; i<10; i++) {
        T[i] = false;
    }
    
    ofstream planets_data, dataE, dataL, dataT;
    planets_data.open("/Users/Luissantra/Desktop/SistemaSolarDef/SistemaSolarDef/planets_data.txt");
    dataE.open("/Users/Luissantra/Desktop/SistemaSolarDef/SistemaSolarDef/dataE.txt");
    dataL.open("/Users/Luissantra/Desktop/SistemaSolarDef/SistemaSolarDef/dataL.txt");
    dataT.open("/Users/Luissantra/Desktop/SistemaSolarDef/SistemaSolarDef/dataT.txt");
    
    
    //Escribimos las cond inciales en el output
    if (geo) {
        
        for (int i=0; i<10; i++) {
            planets_data<< r[i][0] - r[3][0] << ", " << r[i][1] - r[3][1] << endl;
        }
    } else {
        
        for (int i=0; i<10; i++) {
            planets_data<< r[i][0] << ", " << r[i][1] << endl;
        }
    }
    planets_data<<endl;
    
    //Paso 1
    //Calculamos una primera aceleracion
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            if (j!=i) {
                sumax += m[j]*(r[i][0]-r[j][0])/pow(mod(r[i][0]-r[j][0], r[i][1]-r[j][1]), 3);
                sumay += m[j]*(r[i][1]-r[j][1])/pow(mod(r[i][0]-r[j][0], r[i][1]-r[j][1]), 3);
            }
        }
        a[i][0] = -sumax;
        a[i][1] = -sumay;
        sumax = sumay = 0;
    }
    
    for (int i=1; i<=N; i++) {
        for (int j=0; j<10; j++) {
            //Paso 2
            //Calculamos r(t+h) y wi
            r[j][0] += h*v[j][0] + 0.5*h*h*a[j][0];
            r[j][1] += h*v[j][1] + 0.5*h*h*a[j][1];
            w[j][0] = v[j][0] + 0.5*h*a[j][0];
            w[j][1] = v[j][1] + 0.5*h*a[j][1];
            
            //Paso 3
            //Calculamos a(t+h) y la energia potencial U
            for (int k = 0; k<10; k++) {
                if (k!=j) {
                    sumax += m[k]*(r[j][0]-r[k][0])/pow(mod(r[j][0]-r[k][0], r[j][1]-r[k][1]), 3);
                    sumay += m[k]*(r[j][1]-r[k][1])/pow(mod(r[j][0]-r[k][0], r[j][1]-r[k][1]), 3);
                    U[j] -= m[j]*m[k]/mod(r[j][0] - r[k][1], r[j][1] - r[k][1]); //REVISAR
                }
            }
            a[j][0] = -sumax;
            a[j][1] = -sumay;
            sumax = sumay = 0;
            
            //Paso 4
            //Calculamos v(t+h)
            v[j][0] = w[j][0] + 0.5*h*a[j][0];
            v[j][1] = w[j][1] + 0.5*h*a[j][1];
            
            //Energia y  momento angular
            LE[j][1] = 0.5*m[j]*mod(v[j][0], v[j][1])*mod(v[j][0], v[j][1]) + U[j];
            LE[j][0] = m[j]*(r[j][0]*v[j][1]-r[j][1]*v[j][0]);
            U[j] = 0;
            //escribimos el output
            if (geo) {
                planets_data<< r[j][0] - r[3][0] << ", " << r[j][1] - r[3][1] <<endl;
            } else planets_data<< r[j][0] << ", " << r[j][1] <<endl;
            }
            
            
        planets_data<<endl;
        
        // Calculamos energÌa y momento total
        for(int l=0; l<10; l++) {
            sumaL += LE[l][0];
            sumaE += LE[l][1];
        }
        
        if (i>1) {
            dataE<<i-1<<" "<<sumaE<<" "<<'0'<<endl;
            dataL<<i-1<<" "<<sumaL<<" "<<'0'<<endl;
            
        }
        sumaL = sumaE = 0;
        
        //calculamos el periodo
        if (i>50) {
            for (int p=1; p<10; p++) {
                if (r[p][1]>-0.2&&r[p][1]<0.1&&r[p][0]>0&&T[p]==false) {
                    dataT<< p << " " << i*h*58.1 << endl;
                    T[p] = true;
                }
            }
        }
      
        
    }
    
    
    planets_data.close();
    dataL.close();
    dataE.close();
    dataT.close();
    

    return;
}



void iniciarDatos(double m[], double r[10][2], double v[10][2]){
    
    //leemos los datos del fichero
    ifstream datos;
    datos.open("/Users/Luissantra/Desktop/SistemaSolarDef/SistemaSolarDef/condiciones.txt");
    
     if (datos.is_open()) {
             int i=0;
             while (!datos.eof()) {
                 datos>>m[i];
                 datos>>r[i][0];
                 datos>>r[i][1];
                 datos>>v[i][0];
                 datos>>v[i][1];
                 i++;
             }
                 
         }
     
    datos.close();

    //Reescalamiento
    for (int i=0; i<10; i++) {
    r[i][0]/=c;
    r[i][1]/=c;
    v[i][0] *= cte;
    v[i][1] *= cte;
    m[i]/=Ms;
    }
    
    return;
}


double mod(double x, double y) {
    return sqrt(x*x + y*y);
}
