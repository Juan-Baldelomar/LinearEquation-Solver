#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <math.h>
using namespace std;
/*
   A: vector de vectores nxn
   b: es un vector de n
   x: es un vector de n
   Nota: 
     Todos los parámetros son paso por referencia, esto quiere decir que si se modifican dentro de la función también serán modificados fuera de la función, por lo tanto no es necesrio regresar ningun valor.
*/
void Diagonal(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Triangular_Superior(vector<vector<double> > &U, vector<double> &b, vector<double> &x);
void Triangular_Inferior(vector<vector<double> > &L, vector<double> &b, vector<double> &x);
void Eliminacion_Gaussiana(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Eliminacion_Gaussiana_Pivoteo(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Descomposicion_LU(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U);
void pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index) ;
#endif
