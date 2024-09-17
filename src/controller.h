
#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"
#include "ALQPSolver.h"

using namespace BLA;
using namespace ALQPS;

namespace LMPC {

template<int m, int n, int k, int p, typename DType = float>
void insert_at(Matrix<m, n, DType>& M, Matrix<k, p, DType> SubM, int row, int col);
 
template <int n, typename DType = float>
struct DiagonalMatrix : public MatrixBase<DiagonalMatrix<n>, n, n, DType> {
    Matrix<n, 1, DType> diagonal;

    DiagonalMatrix();

    DiagonalMatrix(const Matrix<n, 1, DType>& diag);

    BLA::Matrix<n, n, DType> toMatrix() const; 

    DType operator()(int row, int col) const;  
};

template<int n, int m, int hx, typename DType = float>
class ControllerWeights {
    public:

        ControllerWeights(const Matrix<n, 1, DType>& Q_, 
                          const Matrix<n, 1, DType>& Qf_, 
                          const Matrix<m, 1, DType>& R_);
        
        void print_weights();

        static constexpr int nahx = static_cast<int>((n+m)*hx);
        static constexpr int mhx = static_cast<int>(m*hx);

        Matrix<nahx, nahx, DType> get_Q();
        Matrix<mhx,   mhx, DType> get_R();

    private:
        // Tracking and control weights   
        Matrix<nahx, nahx, DType> Q;
        Matrix<mhx, mhx, DType> R;
};

template<int n, int m, int hx, typename DType = float>
class PredictionMatrices {
    public:

        PredictionMatrices(const Matrix<n,n, DType>& A_, 
                           const Matrix<n,m, DType>& B_);

        static constexpr int nhx = static_cast<int>(n*hx);
        static constexpr int mhx = static_cast<int>(m*hx);

        Matrix<nhx, n, DType> get_predmat_A();
        Matrix<nhx, mhx, DType> get_predmat_B();

    private:
        // Prediction matrices    
        Matrix<nhx,   n, DType> Ap;
        Matrix<nhx, mhx, DType> Bp;
};

template<int n, int m, int hx, typename DType = float>
class CondensedMPC {
    public: 
        static constexpr int na = static_cast<int>(n+m);
        static constexpr int nahx = static_cast<int>(na*hx);
        static constexpr int mhx = static_cast<int>(m*hx);
 
        CondensedMPC(const LTIModel<n, m>& model_,
                     const Matrix<n, 1, DType>& Q_,
                     const Matrix<n, 1, DType>& Qf_,
                     const Matrix<m, 1, DType>& R_,
                     const Matrix<n, 1, DType>& xref_); 

        // destructor
        ~CondensedMPC(){
            if(alloc) {
                delete weights;
                delete PredMat;
                delete QuadProg;
            }
        }

        Matrix<mhx, mhx, DType> get_P();
        Matrix<mhx,   1, DType> get_q();

    private: 
        void update_QP_matrices(const Matrix<static_cast<int>(n+m), 1, DType>& x);                             

        Matrix<nahx, 1, DType> xref; 
        ControllerWeights<n, m, hx, DType>* weights; 
        PredictionMatrices<na, m, hx, DType>* PredMat;
        Matrix<mhx, mhx, DType> P;
        Matrix<mhx,   1, DType> q;
        QP<mhx,0,0>* QuadProg;
        LTIModel<n, m> model;        
        bool alloc;

};                        

}

#include "controller.tpp"