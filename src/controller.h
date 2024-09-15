
#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"

using namespace BLA;

namespace LMPC {

template<int n, int m, int hx>
class ControllerWeights {
    public:
        ControllerWeights(const Matrix<n>& Q_, 
                          const Matrix<n>& Qf_, 
                          const Matrix<m>& R_);
        void print_weights();
    private:
        static constexpr int nhx = static_cast<int>(n*hx);
        static constexpr int mhx = static_cast<int>(m*hx);
        // Tracking and control weights   
        Matrix<nhx,nhx> Q;
        Matrix<mhx,mhx> R;
};

template<int n, int m, int hx>
class PredictionMatrices {
    public:
        PredictionMatrices(const Matrix<n,n>& A_, 
                           const Matrix<n,m>& B_);

        Matrix<static_cast<int>(n*hx),n> get_predmat_A();
        Matrix<static_cast<int>(n*hx),static_cast<int>(m*hx)> get_predmat_B();

        // Tracking and control weights   
    private:
        static constexpr int nhx = static_cast<int>(n*hx);
        static constexpr int mhx = static_cast<int>(m*hx);
        
        Matrix<nhx,  n> Ap;
        Matrix<nhx,mhx> Bp;
};

template<int m, int n, int k, int p>
void insert_at(Matrix<m,n>& M, Matrix<k,p> SubM, int row, int col){
    for (int i = 0; i < k; i++){
        for (int j = col; j < p; j++){
            M(row+i,col+j) = SubM(i,j); 
        }
    }
}

}

#include "controller.tpp"