#include <bits/stdc++.h>
#include "Solvers.hpp"
#include "Tools.hpp"

void Diagonal(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {
    int n = A.size();
    for (int i = 0; i<n; i++){
        if (A[i][i] == 0)
            cout << "WARNING: Elemento nulo en posicion [" << i << "][" << i <<"]" <<endl;
        
        x[i] = b[i]/A[i][i];
    }
    
}

void Triangular_Superior(vector<vector<double> > &U, vector<double> &b, vector<double> &x) {
    int n = U.size();
    double X[n];

    for (int i = n - 1; i >= 0; i--) {
        double acc = 0;
        for (int j = n - 1; j > i; j--) {
            acc += U[i][j] * x[j];
        }
        x[i] = (b[i] - acc) / U[i][i];
    }
}

void Triangular_Inferior(vector<vector<double> > &L, vector<double> &b, vector<double> &x) {
    int n = L.size();

    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = 0; j < i; j++) {
            acc += L[i][j] * x[j];
        }
        x[i] = (b[i] - acc) / L[i][i];
    }
}

void Eliminacion_Gaussiana(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {

    int n = A.size();
    //indice k indica la fila con la cual estamos transformando a una Triangular Superior
    for (int k = 0; k < n - 1; k++) {
        if (A[k][k]== 0)
            cout << "WARNING: Elemento nulo en posicion [" << k << "][" << k <<"]" <<endl;
        
        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }
    }
    Triangular_Superior(A, b, x);
}

void Eliminacion_Gaussiana_Pivoteo(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {
    int n = A.size();
    vector<int> index ;
    vector<double> orderedX;
    
    //inicializar
    orderedX.assign(n, 0);
    index.assign(n, 0);
    
    // llenar de las posiciones correspondientes de cada variable
    for (int i = 0; i < n; i++) {
        index[i] = i;
    }

    //indice k indica la fila con la cual estamos transformando a una Triangular Superior
    for (int k = 0; k < n-1; k++) {
        pivot(A, k, b, index);
        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }

    }
    Triangular_Superior(A, b, x);
    
    //ordenar arreglo x para dar la solucion
    for (int i = 0; i<n; i++){
        orderedX[index[i]] =  x[i];
    }
    
    //copiar data a arreglo X
    for (int i = 0; i<n; i++){
        x[i] = orderedX[i];
    }
    
}

void Descomposicion_LU(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U) {
    int n = A.size();
    for (int i = 0; i<n; i++){
        for (int j =0; j<n; j++){
            
            double acc = 0;
            if(j<=i){
                for(int k = 0; k<j; k++){
                    acc+= A[i][k]*A[k][j];   
                }
                
                A[i][j] = A[i][j] - acc;
                
            }else{
                for(int k = 0; k<i; k++){
                    acc+= A[i][k]*A[k][j];   
                }
                A[i][j] = (A[i][j] - acc)/A[i][i];
            }
            
        }
    }
    
    // COPIA a MATRIZ L, U (No es necesario, pero para respetar la plantilla)
    for (int i = 0; i<n; i++){
        for (int j = 0; j<=i; j++){
            L[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i<n; i++){
        for (int j = i; j<n; j++){
            if (i==j)
                U[i][j] = 1;
            else
                U[i][j] = A[i][j];
        }
    }
    
    vector<vector<double>> C;
    C.assign(n, vector<double>(n,0.0));
    
    Matrix_Mult(L, U, C);
    cout << "MULT RES" <<endl;
    cout << C<< endl;
    
    cout << A << endl;
}


void Descomposicion_LU_Pivot(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U) {
    int n = A.size();
    
    for (int k = 0; k<n; k++){
        double acc;
        LU_pivot(A,k);
        for (int j = k; j<n; j++){
            acc = 0;
            //multiplicar   
            int km = 0;
            while(km < k && km < j){
                if(k==1)
                    cout << endl;
                acc += A[k][km]*A[km][j];
                km++;
            }
            
            if (j<=k)
                A[k][j] = A[k][j]- acc;
            else
                A[k][j] = (A[k][j]- acc)/A[k][k];
            
        }
        for (int i = k+1; i<n; i++){
            acc = 0;
            //multiplicar   
            int km = 0;
            while(km < k && km < i){
                acc += A[i][km]*A[km][k];
                km++;
            }
            
            if (k <=i)
                A[i][k] = (A[i][k]- acc);
            else
                A[i][k] = (A[i][k]- acc)/A[i][i]; 
        }  
    }
    // COPIA a MATRIZ L, U (No es necesario, pero para respetar la plantilla)
    for (int i = 0; i<n; i++){
        for (int j = 0; j<=i; j++){
            L[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i<n; i++){
        for (int j = i; j<n; j++){
            if (i==j)
                U[i][j] = 1;
            else
                U[i][j] = A[i][j];
        }
    }
    
    vector<vector<double>> C;
    C.assign(n, vector<double>(n,0.0));
    
    Matrix_Mult(L, U, C);
    cout << "MULT R" <<endl;
    cout << C<< endl;
    cout << " LU" << endl;
    cout << A << endl;
}





void pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index) {
    int n = A.size();
    int i_max = k, j_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k; i < n; i++) {
        for (int j = k; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                i_max = i;
                j_max = j;
                max = fabs(A[i][j]);
            }
        }
    }

    double temp;

    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }

    // swap column
    if (j_max != k) {
        for (int i = 0; i < n; i++) {
            temp = A[i][k];
            A[i][k] = A[i][j_max];
            A[i][j_max] = temp;
        }
        //swap variable order
        int at_j = index[j_max];
        int at_k = index[k];
        index[j_max] = at_k;
        index[k] = at_j;
    }
}

void LU_pivot(vector<vector<double> > &A, int k) {
    int n = A.size();
    int i_max = k, j_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k; i < n; i++) {
        for (int j = k; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                i_max = i;
                j_max = j;
                max = fabs(A[i][j]);
            }
        }
    }

    double temp;

    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
    }

    // swap column
    if (j_max != k) {
        for (int i = 0; i < n; i++) {
            temp = A[i][k];
            A[i][k] = A[i][j_max];
            A[i][j_max] = temp;
        }
    }
}

void Matrix_Mult(vector< vector<double > > &A, vector<vector<double>>&B, vector<vector<double>>&C){
    int n = A.size();
    
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            double acc = 0;
            for (int k =0; k<n; k++){
                acc+= A[i][k]*B[k][j];
            }
            C[i][j] = acc;
        }
    }
}



