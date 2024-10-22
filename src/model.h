#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace LMPC {

template<int n, int m, int pu, int px, typename DType = float>
class Constraints {
    public: 
        Matrix<m+m,1> _limits_u;
        Matrix<n+n,1> _limits_x;

        int num_state_constraints(){
            int num_constraints = 0;
            for (int i = 0; i < n+n; i++){
                if (!isnan(_limits_x(i,0))) 
                    num_constraints += 1;
            }
            return num_constraints;            
        }
        int num_input_constraints(){
            int num_constraints = 0;
            for (int i = 0; i < m+m; i++){
                if (!isnan(_limits_u(i,0))) 
                    num_constraints += 1;
            }
            return num_constraints;             
        }                                
        Constraints(const Matrix<m+m,1>& limits_u, 
                    const Matrix<n+n,1>& limits_x)
                    : _limits_u(limits_u),
                      _limits_x(limits_x){
                        // constexpr int nc = num_state_constraints()+num_input_constraints();
                        // static_assert(nc!=p,"Please provide the correct number of constraints!");
        }; 
        Constraints(const Matrix<m+m,1>& limits_u)
                    : _limits_u(limits_u){
                      _limits_x.Fill(0.0/0.0);
                        // constexpr int nc = num_input_constraints();                                                
                        // static_assert(nc!=p, "Please provide the correct number of constraints!");
        };                      
        Constraints(const Matrix<n+n,1>& limits_x)
                    : _limits_x(limits_x){        
                        _limits_u.Fill(0.0/0.0);                        
                        // constexpr int nc = num_input_constraints();                                                
                        // static_assert(nc!=p, "Please provide the correct number of constraints!");
        };                       
};

template<int n, int m, typename DType = float>
class LTIModel {
    public:
        LTIModel(const Matrix<n,n,DType>& A_, 
                 const Matrix<n,m,DType>& B_, 
                 const DType& h_)
                    : A(A_),
                      B(B_),
                      h(h_) {};

        Matrix<n,n> get_A(){
            return A;
        }
        Matrix<n,m> get_B(){
            return B;
        }        
        Matrix<n,1,DType> simulate_time_step(const Matrix<n,1,DType>& x, const Matrix<m,1,DType>& u){
            return A*x+B*u; 
        }

    private:
        // State space matrices 
        Matrix<n,n,DType> A;
        Matrix<n,m,DType> B;
        // discretization step size 
        DType h;
};

}