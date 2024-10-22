#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"
#include "ALQPSolver.h"


using namespace ALQPS;

const QPparams params_default;

namespace LMPC {

// linear MPC class
template<int n, int m, int hx, int pu = 0, int px = 0, typename DType = float>
class lMPC {
    public: 
        // dimension used across all functions 
        static constexpr int na = n+m;
        static constexpr int nahx = na*hx;
        static constexpr int nhx = n*hx;        
        static constexpr int mhx = m*hx;
        static constexpr int phx = (px+pu)*hx;

        // constrained MPC 
        lMPC(const LTIModel<n, m, DType>& model_,
            const Matrix<n, 1, DType>& Q_,
            const Matrix<n, 1, DType>& Qf_,
            const Matrix<m, 1, DType>& R_,
            const Matrix<n, 1, DType>& xref_,
            const Constraints<n, m, pu, px, DType>& con,
            const QPparams params = params_default);

        // unconstrained MPC 
        lMPC(const LTIModel<n, m, DType>& model_,
            const Matrix<n, 1, DType>& Q_,
            const Matrix<n, 1, DType>& Qf_,
            const Matrix<m, 1, DType>& R_,
            const Matrix<n, 1, DType>& xref_,
            const QPparams params = params_default);

        Matrix<m, 1, DType> get_control(const Matrix<n,1, DType>& x); 

    private:
        // adapt state and control matrix to augmented state 
        void state_augmentation(Matrix<n,n,DType>& A, const Matrix<n,m,DType>& B);

        // build prediction matrices for future augmented state prediction
        void build_prediction_matrices(Matrix<nahx,  na, DType>& Ap, 
                                       BlockLowerTriangularToeplitzMatrix<na, m, hx, DType>& Bp, 
                                       const Matrix<na, na, DType>& Aaug, 
                                       const Matrix<na,  m, DType>& Baug);

        // build prediction matrices for future state prediction
        void build_prediction_matrices(Matrix<nhx,  n, DType>& Ap, 
                                       BlockLowerTriangularToeplitzMatrix<n, m, hx, DType>& Bp, 
                                       const Matrix<n, n, DType>& A, 
                                       const Matrix<n, m, DType>& B);
        
        // pass quadratic cost matrix to solver
        void initialize_objective();

        // pass linear cost term to solver
        void update_objective(const Matrix<n, 1, DType>& x);

        // pass quadratic inequality constraint matrix to solver
        void initialize_constraints();

        // pass linear inequality constraint term to solver
        void update_constraints(const Matrix<n,1,DType>& x);

        void initialize_state_constraints(int& abs_index, 
                                          int rel_index_start, 
                                          int rel_index_end, 
                                          const Matrix<n+n,1>& h_const_part,
                                          const Matrix<nhx, mhx, DType>& G_const_part,
                                          Matrix<phx,mhx,DType>& G);

        void update_state_constraints(int& abs_index, 
                                      int rel_index_start, 
                                      int rel_index_end, 
                                      const Matrix<nhx,1>& h_offset,
                                      Matrix<phx,1,DType>& h);

        void initialize_input_constraints(int& abs_index, 
                                          int rel_index_start, 
                                          int rel_index_end, 
                                          const Matrix<m+m,1>& h_const_part,
                                          const Matrix<mhx, mhx, DType>& G_const_part,
                                          Matrix<phx,mhx,DType>& G);

        void update_input_constraints(int& abs_index, 
                                      int rel_index_start, 
                                      int rel_index_end, 
                                      const Matrix<m,1,DType>& h_offset,
                                      Matrix<phx,1,DType>& h);

        // enable printing free memory 
        bool debugging; 

        // horizon wide reference (setpoint)
        Matrix<nahx, 1, DType> xref; 
        
        // State space matrices 
        Matrix<n,n,DType> A;
        Matrix<n,m,DType> B;

        // Augmented state space matrices
        Matrix<na, na, DType> Aaug; 
        Matrix<na,  m, DType> Baug;

        // Controller weights
        Matrix<n, 1, DType> Q;
        Matrix<n, 1, DType> Qf;
        Matrix<m, 1, DType> R;

        // controls at current operational point 
        Matrix<m,1,DType> u_op;

        // limits on states and controls
        Matrix<m+m,1> _limits_u;
        Matrix<n+n,1> _limits_x;        
        
        // prediction matrices for unaugmented state (used for constraint assembly)
        Matrix<nhx,  n, DType> Ap_c;
        BlockLowerTriangularToeplitzMatrix<n, m, hx, DType> Bp_c;
        Matrix<phx,1,DType> h_const;

        // prediction matrices for augmented state (used for cost function assembly)
        Matrix<nahx,  na, DType> Ap;
        Matrix<mhx,nahx, DType> Bp_scaled;

        QP<mhx,0,phx,DType> QuadProg;
    
};                        

}

#include "controller.tpp"