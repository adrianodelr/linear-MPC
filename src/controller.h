
#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"

using namespace BLA;

namespace LMPC {

template<int m, int n, int k, int p, typename DType = float>
void insert_at(Matrix<m,n, DType>& M, Matrix<k,p, DType> SubM, int row, int col);

template <int n, typename DType = float>
struct DiagonalMatrix : public MatrixBase<DiagonalMatrix<n>, n, n, DType> {
    Matrix<n, 1, DType> diagonal;

    DiagonalMatrix();

    DiagonalMatrix(const Matrix<n, 1, DType>& diag);

    DType operator()(int row, int col) const;
};

template<int n, int m, int hx, typename DType = float>
class ControllerWeights {
    public:
        ControllerWeights(const Matrix<n, 1, DType>& Q_, 
                          const Matrix<n, 1, DType>& Qf_, 
                          const Matrix<m, 1, DType>& R_);
        
        void print_weights();

    private:
        static constexpr int nhx = static_cast<int>(n*hx);
        static constexpr int mhx = static_cast<int>(m*hx);

        // Tracking and control weights   
        Matrix<nhx, nhx, DType> Q;
        Matrix<mhx, mhx, DType> R;
};

template<int n, int m, int hx, typename DType = float>
class PredictionMatrices {
    public:
        PredictionMatrices(const Matrix<n,n, DType>& A_, 
                           const Matrix<n,m, DType>& B_);

        Matrix<static_cast<int>(n*hx),n, DType> get_predmat_A();
        Matrix<static_cast<int>(n*hx),static_cast<int>(m*hx), DType> get_predmat_B();

        // Tracking and control weights   
    private:
        static constexpr int nhx = static_cast<int>(n*hx);
        static constexpr int mhx = static_cast<int>(m*hx);
        
        Matrix<nhx,   n, DType> Ap;
        Matrix<nhx, mhx, DType> Bp;
};

}

#include "controller.tpp"