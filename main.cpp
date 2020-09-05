#include <bits/stdc++.h>
#include "Tools.hpp"
#include "Solvers.hpp"
using namespace std;
void Ejercicio1()
{
   //1) Resolver una matriz Diagonal A
    vector< vector<double> > D;
    vector<double> x, b;
    ReadMatrix(D, "MatricesT02/M_DIAG.txt");
    ReadVector(b, "MatricesT02/V_DIAG.txt");
    x.assign((int)b.size(), 0.0);
    Diagonal(D, b, x);
    //cout.precision(2) << x <<endl; //imprimir en pantalla
    WriteVector(x, "out/Punto1_X_SMALL.txt");

}
void Ejercicio2()
{
   //2) Resolver una matriz triangula superior U
    vector< vector<double> > U;
    vector<double> x, b;
    ReadMatrix(U, "MatricesT02/M_TSUP.txt");
    ReadVector(b, "MatricesT02/V_TSUP.txt");
    x.assign((int)b.size(), 0.0);
    Triangular_Superior(U, b, x);
    cout.precision(2); 
    cout << x <<endl; //imprimir en pantalla
    WriteVector(x, "out/Punto2_X_SMALL.txt");
}
void Ejercicio3()
{
   //3) Resolver una matriz triangula inferior L
    vector< vector<double> > L;
    vector<double> x, b;
    ReadMatrix(L, "MatricesT02/M_TINF.txt");
    ReadVector(b, "MatricesT02/V_TINF.txt");
    x.assign((int)b.size(), 0.0);
    Triangular_Inferior(L, b, x);
    cout << x<<endl;
    WriteVector(x, "out/Punto3_X_SMALL.txt");
}
void Ejercicio4()
{

   //4) Eliminación Guassiana
    vector< vector<double> > A;
    vector<double> x, b;
    ReadMatrix(A, "MatricesT02/M_LARGE.txt");
    ReadVector(b, "MatricesT02/V_LARGE.txt");
    x.assign((int)b.size(), 0.0);
    Eliminacion_Gaussiana(A, b, x);
    cout << x<<endl;
    WriteVector(x, "out/Punto4_X_SMALL.txt");
}
void Ejercicio5()
{
   //5) Eliminación Gaussiana con pivoteo (opcional)
    vector< vector<double> > A;
    vector<double> x, b;
    ReadMatrix(A, "MatricesT02/M_LARGE.txt");
    ReadVector(b, "MatricesT02/V_LARGE.txt");
    x.assign((int)b.size(), 0.0);
    Eliminacion_Gaussiana_Pivoteo(A, b, x);
    cout << x<<endl;
    WriteVector(x, "out/Punto5_X_SMALL.txt");
}
void Ejercicio6()
{
   //6) Descomposición LU
   vector< vector<double> > A;
   vector<vector<double> > L, U;
   ReadMatrix(A, "MatricesT02/M_SMALL.txt");
   L.assign((int)A.size(), vector<double> ((int)A[0].size(), 0.0));
   U.assign((int)A.size(), vector<double> ((int)A[0].size(), 0.0));
   Descomposicion_LU_Pivot(A, L, U); 
   cout << "Lower:"<<endl;
   cout << L <<endl;
   cout << "Upper:"<<endl;
   cout << U <<endl;
   WriteMatrix(L, "out/L_6_SMALL.txt");
   WriteMatrix(U, "out/U_6_SMALL.txt");
}
int main()
{

  //1) Resolver una matriz Diagonal A
   // Ejercicio1();
  //2) Resolver una matriz triangula superior U
    //Ejercicio2();
  //3) Resolver una matriz triangula inferior L
     //Ejercicio3();
  //4) Eliminación Guassiana
     //Ejercicio4();
  //5) Eliminación Gaussiana con pivoteo (opcional)
     //Ejercicio5();
  //6) Descomposición LU
     Ejercicio6();

  return 0;
}
