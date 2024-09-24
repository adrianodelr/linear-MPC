
#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"
#include "ALQPSolver.h"

using namespace BLA;
using namespace ALQPS;

namespace LMPC {

template<int m, int n, int k, int p, typename DType = float>
void insert_at(Matrix<m, n, DType>& M, const Matrix<k, p, DType>& SubM, int row, int col);
 
template <int n, typename DType = float>
struct DiagonalMatrix : public MatrixBase<DiagonalMatrix<n>, n, n, DType> {
    Matrix<n, 1, DType> diagonal;

    DiagonalMatrix();

    DiagonalMatrix(const Matrix<n, 1, DType>& diag);

    BLA::Matrix<n, n, DType> toMatrix() const; 

    DType operator()(int row, int col) const;  
};

template<int n, int m, int hx, typename DType = float>
class CondensedMPC {
    public: 
        static constexpr int na = static_cast<int>(n+m);
        static constexpr int nahx = static_cast<int>(na*hx);
        static constexpr int mhx = static_cast<int>(m*hx);
 
        CondensedMPC(const LTIModel<n, m, DType>& model_,
                     const Matrix<n, 1, DType>& Q_,
                     const Matrix<n, 1, DType>& Qf_,
                     const Matrix<m, 1, DType>& R_,
                     const Matrix<n, 1, DType>& xref_); 

        // destructor
        ~CondensedMPC(){
            if(alloc) {
                delete QuadProg;
            }
        }

        Matrix<mhx,mhx,DType> get_P();
        Matrix<mhx,1,DType> get_q();        

    private: 
        void update_QP_matrices(const Matrix<na,     1, DType>& x,
                                const Matrix<nahx,  na, DType>& Ap,
                                const Matrix<nahx, mhx, DType>& Bp,
                                const DiagonalMatrix<nahx, DType>& Q, 
                                const DiagonalMatrix<mhx, DType>& R);                             

        void build_prediction_matrices(Matrix<nahx,  na, DType>& Ap, 
                                       Matrix<nahx, mhx, DType>& Bp, 
                                       const Matrix<na, na, DType>& Aaug, 
                                       const Matrix<na,  m, DType>& Baug);

        void build_weight_matrices(DiagonalMatrix<nahx, DType>& Q, 
                                   DiagonalMatrix<mhx, DType>& R, 
                                   const Matrix<n, 1, DType>& Q_,
                                   const Matrix<n, 1, DType>& Qf_,
                                   const Matrix<m, 1, DType>& R_);                                   


        Matrix<nahx, 1, DType> xref; 
        
        // State space matrices 
        Matrix<n,n,DType> A;
        Matrix<n,m,DType> B;

        Matrix<mhx, mhx, DType> P;
        Matrix<mhx,   1, DType> q;
        QP<mhx,0,0,DType>* QuadProg;
        bool alloc;

};                        

}

#include "controller.tpp"