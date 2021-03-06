#include <bits/stdc++.h>
#include "Tools.hpp"
#include "Solvers.hpp"
using namespace std;
void Ejercicio1()
{
   //1) Resolver una matriz Diagonal A
    cout << endl;
    cout << "****************************************** RESOLVER MATRIZ DIAGONAL ******************************************" << endl;
    vector< vector<double> > D;
    vector<double> x, b;
    ReadMatrix(D, "MatricesT02/M_DIAG.txt");
    ReadVector(b, "MatricesT02/V_DIAG.txt");
    x.assign((int)b.size(), 0.0);
    Diagonal(D, b, x);
    cout << "#####################################################  VECTOR SOLUCION ################################################### " << endl;
    cout << x <<endl; //imprimir en pantalla
    WriteVector(x, "out/Punto1_X_SMALL.txt");

}
void Ejercicio2()
{
   //2) Resolver una matriz triangula superior U
    cout << endl;
    cout << "****************************************** RESOLVER MATRIZ TRIANGULAR SUPERIOR ******************************************" << endl;
    vector< vector<double> > U;
    vector<double> x, b;
    ReadMatrix(U, "MatricesT02/M_TSUP.txt");
    ReadVector(b, "MatricesT02/V_TSUP.txt");
    x.assign((int)b.size(), 0.0);
    Triangular_Superior(U, b, x);
    cout << "#####################################################  VECTOR SOLUCION ################################################### " << endl;
    cout << x <<endl; //imprimir en pantalla
    // prueba solucion
    ReadMatrix(U, "MatricesT02/M_TSUP.txt");
    ReadVector(b, "MatricesT02/V_TSUP.txt");
    Try_Sol(U, b, x);
    
    WriteVector(x, "out/Punto2_X_SMALL.txt");
}
void Ejercicio3()
{
   //3) Resolver una matriz triangula inferior L
    cout << endl;
    cout << "****************************************** RESOLVER MATRIZ TRIANGULAR INFERIOR ******************************************" << endl;
    vector< vector<double> > L;
    vector<double> x, b;
    ReadMatrix(L, "MatricesT02/M_TINF.txt");
    ReadVector(b, "MatricesT02/V_TINF.txt");
    x.assign((int)b.size(), 0.0);
    Triangular_Inferior(L, b, x);
    cout << "#####################################################  VECTOR SOLUCION ################################################### " << endl;
    cout << x <<endl;
    // prueba solucion
    ReadMatrix(L, "MatricesT02/M_TINF.txt");
    ReadVector(b, "MatricesT02/V_TINF.txt");
    Try_Sol(L, b, x);
    
    WriteVector(x, "out/Punto3_X_SMALL.txt");
}
void Ejercicio4()
{
   //4) Eliminación Guassiana
    cout << endl;
    cout << "****************************************** RESOLVER ELIMINACION GAUSSIANA ******************************************" << endl;
    vector< vector<double> > A;
    vector<double> x, b;
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    x.assign((int)b.size(), 0.0);
    Eliminacion_Gaussiana(A, b, x);
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++ MATRIZ RESULTANTE +++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    cout << A << endl;
    cout <<endl;
    cout << "##################################################  VECTOR SOLUCION ############################################### " << endl;
    cout << x <<endl;
    // prueba solucion
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    Try_Sol(A, b, x);
    
    WriteVector(x, "out/Punto4_X_SMALL.txt");
}
void Ejercicio5()
{
   //5) Eliminación Gaussiana con pivoteo (opcional)
    cout << endl;
    cout << "****************************************** RESOLVER ELIMINACION GAUSSIANA PIVOTEO COMPLETO ******************************************" << endl;
    vector< vector<double> > A;
    vector<double> x, b;
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    x.assign((int)b.size(), 0.0);
    Eliminacion_Gaussiana_Pivoteo(A, b, x);
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++ MATRIZ RESULTANTE +++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    cout << A << endl;
    cout <<endl;
    cout << "########################################################  VECTOR SOLUCION ##################################################### " << endl;
    cout << x<<endl;
    //prueba solucion
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    Try_Sol(A, b, x);
    
    WriteVector(x, "out/Punto5_X_SMALL.txt");
}
void Ejercicio6()
{
   //6) Descomposición LU
   cout << endl; 
   cout << "************************************************ FACTORIZACION LU **************************************************" << endl; 
   vector< vector<double> > A;
   vector<vector<double> > L, U;
   ReadMatrix(A, "MatricesT02/M_SMALL.txt");
   L.assign((int)A.size(), vector<double> ((int)A[0].size(), 0.0));
   U.assign((int)A.size(), vector<double> ((int)A[0].size(), 0.0));
   Descomposicion_LU(A, L, U); 
   cout << "------------------------------------------------ Lower: -------------------------------------------------------------"<<endl;
   cout << L <<endl;
   cout << "------------------------------------------------ Upper: -------------------------------------------------------------"<<endl;
   cout << U <<endl;
   WriteMatrix(L, "out/L_6_SMALL.txt");
   WriteMatrix(U, "out/U_6_SMALL.txt");
}

void Ejercicio7()
{
   //4) Resolver SISTEMA con LU
    cout << endl;
    cout << "****************************************** RESOLVER POR FACTORIZACION LU ******************************************" << endl;
    vector< vector<double> > A;
    vector<double> x, b;
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    x.assign((int)b.size(), 0.0);
    LU_Solve(A,b,x);
    cout << "########################################################  VECTOR SOLUCION ##################################################### " << endl;
    cout << x <<endl;
    //prueba solucion
    ReadMatrix(A, "MatricesT02/M_SMALL.txt");
    ReadVector(b, "MatricesT02/V_SMALL.txt");
    Try_Sol(A, b, x);
    
    WriteVector(x, "out/Punto7_X_SMALL.txt");
}


int main()
{

  //1) Resolver una matriz Diagonal A
    Ejercicio1();
  //2) Resolver una matriz triangula superior U
    Ejercicio2();
  //3) Resolver una matriz triangula inferior L
    Ejercicio3();
  //4) Eliminación Guassiana
    Ejercicio4();
  //5) Eliminación Gaussiana con pivoteo (opcional)
    Ejercicio5();
  //6) Descomposición LU
    Ejercicio6();
  //7) Resolver por LU
    Ejercicio7();   
   

  return 0;
}
