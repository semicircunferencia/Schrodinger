/*SIMULACIÓN DE UNA PARTÍCULA CON UNA BARRERA DE POTENCIAL FINITO 1D

DISCRETIZACIÓN DEL TIEMPO Y EL ESPACIO
Consideramos x_j=j*h, t_n=n*s, donde h y s son los steps espacial y temporal

ALGORITMO
1) Dar los parámetros iniciales, N (número de intervalitos en los que dividimos el espacio), n_ciclos
(tiene que ver con el momento de la partícula, o, equivalentemente, el k de la onda), y lambda
(elevación de la barrera).

2) Generar sgorro (paso s en otras unidades), k0gorro (k en otras unidades),
Vgorro (potencial igual), función de onda inicial, y calcular alfa.

3) Calcular beta usando la recurrencia del guion

4) Calcular chi usando la recurrencia del guion

5) n=n+1. Ir al paso 2
*/

#include <iostream>
#include <fstream> // Ficheros
#include <cmath>
#include <iomanip> // Formato en los ficheros
#include <complex>

#define N 100
#define h 0.01
#define itertemp 1e5
#define pi 3.1415926536

using namespace std;

void calcpotencialgorro(double potencialgorro[], double lambda, double k0gorro);
void inicializarfunconda(complex<double> funconda[], double k0gorro);
void calcalfa(complex<double> alfa[], double potencialgorro[], double sgorro);
void iteracion(complex<double> alfa[], complex<double> funconda[], double sgorro);
double norma(complex<double> funconda[]);

/**************************************************** FUNCIÓN PRINCIPAL ****************************************************/
int main(void) {
    // Recuerda que los vectores posicion[], potencial[]... tienen dimensión N+1
    // Defino variables
    int nciclos; // Tiene que ver con el momento de la onda, toma valores enteros de 1 a N/4
    double lambda; // Elevación del obstáculo
    double k0gorro;
    double sgorro;

    double potencialgorro[N+1];
    complex<double> funconda[N+1];
    complex<double> alfa[N];

    // Pido variables
    cout << "Elija un valor entero de 1 a " << N/4 << " para n_ciclos\n";
    cin >> nciclos;

    cout << "Escoja la elevación del obstáculo lambda\n";
    cin >> lambda;

    // Calculo k0gorro, sgorrro
    k0gorro=2*pi*nciclos/(1.0*N);
    sgorro=1/(4*k0gorro*k0gorro);

    // Inicializo el potencial y la función de onda
    calcpotencialgorro(potencialgorro, lambda, k0gorro);
    inicializarfunconda(funconda, k0gorro);

    // Calculo los alfas
    calcalfa(alfa, potencialgorro, sgorro);

    // Abro el archivo donde escribo tanto la norma al cuadrado como el potencial en cada instante
    ofstream fichero;
    ofstream fich_norma;
    fichero.open("schrodinger_data.dat");
    fich_norma.open("norma.dat");

    // Escribo la situación inicial
    for(int j=0; j<=N; j++) fichero << j << "," << norm(funconda[j]) << "," << potencialgorro[j] << "\n";
    fich_norma << norma(funconda);

    // Itero en el tiempo
    for(int k=1; k<itertemp; k++) {

        // Calculo la siguiente función de onda
        iteracion(alfa, funconda, sgorro);

        // Escribo los datos, solo algunas veces
        if(k%100==0) {
            fichero << "\n";
            fich_norma << "\n";

            for(int j=0; j<=N; j++) fichero << j << "," << norm(funconda[j]) << "," << potencialgorro[j] << "\n";
            fich_norma << norma(funconda);
        }

    }

    fichero.close();
    fich_norma.close();

    return 0;
}
/***************************************************************************************************************************/

/*Función calcpotencial. Inicializa el potencial, en este caso, a una barrera cuadrada de centro N/2 y anchura N/5*/
void calcpotencialgorro(double potencialgorro[], double lambda, double k0gorro) {
    int extremo1=2*N/5;
    int extremo2=3*N/5;

    // Calculo el valor del potencial en el punto j
    for(int j=0; j<extremo1; j++) potencialgorro[j]=0;
    for(int j=extremo1; j<=extremo2; j++) potencialgorro[j]=lambda*k0gorro*k0gorro;
    for(int j=extremo2+1; j<N; j++) potencialgorro[j]=0;

    return;
}

/*Función inicializarfunconda. Inicializa la función de onda a una gaussiana*/
void inicializarfunconda(complex<double> funconda[], double k0gorro) {

    complex<double> aux;
    double denominador;

    // Condiciones de contorno
    funconda[0]=0; funconda[N]=0;

    // Resto de puntos
    for(int j=1; j<N; j++) {
        aux={cos(k0gorro*j),sin(k0gorro*j)};
        funconda[j]=exp(-8*(4*j-N)*(4*j-N)/(1.0*N*N))*aux;
    }

    // Normalizo
    denominador=sqrt(norma(funconda));
    for(int j=1; j<N; j++) funconda[j]=funconda[j]/denominador;

    return;
}

/*Función calcalfa. Calcula los alfa_j a partir de los datos, usando la recurrencia:
alfa_(j-1)=-1/(-2+2*i/sgorro-Vgorro_j+alfa_j)*/
void calcalfa(complex<double> alfa[], double potencialgorro[], double sgorro) {

    // Defino el número i
    complex<double> i={0,1};

    // alfa[N-1]=0
    alfa[N-1]=0;

    // Resto de alfas, hasta alfa[0]
    for(int j=N-1; j>0; j--) alfa[j-1]=-1.0/(-2.0+2.0*i/sgorro-potencialgorro[j]+alfa[j]);

    return;
}

/*Función iteraacion. Calcula la función de onda en el siguiente step temporal. Para ello, calcula los beta_jn,
usando los datos: beta_(j-1)n=-alfa_(j-1)*(b_jn - beta_jn), siendo b_jn=4*i*funconda_jn/sgorro*/
void iteracion(complex<double> alfa[], complex<double> funconda[], double sgorro) {

    // Defino el número i
    complex<double> i={0,1};

    // Defino los vectores beta y chi
    complex<double> beta[N];
    complex<double> chi[N+1];

    // beta[N-1]=0
    beta[N-1]=0;

    // Resto de betas, hasta beta[0]
    for(int j=N-1; j>0; j--) beta[j-1]=alfa[j-1]*(beta[j]-4.0*i*funconda[j]/sgorro);

    // Calculo los chi
    chi[0]=0;
    for(int j=0; j<N; j++) chi[j+1]=alfa[j]*chi[j]+beta[j];

    // Finalmente
    for(int j=0; j<=N; j++) funconda[j]=chi[j]-funconda[j];

    return;
}

/*Función norma. Devuelve la norma de la función de onda. Asume que funconda[0]=funconda[N]=0*/
double norma(complex<double> funconda[]) {
    
    double suma=0;
    for(int j=1; j<N; j++) suma+=norm(funconda[j]);

    return h*suma;
}
