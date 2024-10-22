#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace LMPC {

// constraint class which takes in limits on controls and state  
template<int n, int m, int pu, int px, typename DType = float>
class Constraints {
    public: 
        // numerical values are constraints while NaN entries represent inactive limits
        // vectors first list upper bounds in ascending order, then lower bounds in ascending order
        // e.g. m = 2; 
        // limits: 
        // u1_max = 1, u1_min = -5, u2_max = no_constraint, u2_min = -3  
        // _limits_u = {u1_max, u2_max, u1_min, u2_min} => {1, 0.0/0.0, -5, -3}; 
        Matrix<n+n,1> _limits_x;    // state constraints 
        Matrix<m+m,1> _limits_u;    // input constraints 

        // Constructor that initializes object with input and state constraints
        Constraints(const Matrix<m+m,1>& limits_u, 
                    const Matrix<n+n,1>& limits_x)
                    : _limits_u(limits_u),
                      _limits_x(limits_x){
                    test_num_input_constraints();
                    test_num_state_constraints();
        }; 

        // Constructor that initializes object with only input constraints
        Constraints(const Matrix<m+m,1>& limits_u)
                    : _limits_u(limits_u){
                      _limits_x.Fill(0.0/0.0);
                    test_num_input_constraints();
        };                      

        // Constructor that initializes object with only state constraints
        Constraints(const Matrix<n+n,1>& limits_x)
                    : _limits_x(limits_x){        
                      _limits_u.Fill(0.0/0.0);
                    test_num_state_constraints();                        
        }; 

    private:
        // since compiler avr compiler doesnt support exceptions, use while loop   
        void throw_exception(String exception){
            Serial.println(exception);
            while (1){
                // do nothing
            }
        }
        
        // test if template parameter px matches provided number of state constraints    
        void test_num_state_constraints(){
            int num_constraints = 0;
            for (int i = 0; i < n+n; i++){
                if (!isnan(_limits_x(i,0))) 
                    num_constraints += 1;
            }
            if (num_constraints!=px)
                throw_exception(F("Please provide the correct number of state constraints!"));
                        
        }
        
        // test if template parameter pu matches provided number of input constraints    
        void test_num_input_constraints(){
            int num_constraints = 0;
            for (int i = 0; i < m+m; i++){
                if (!isnan(_limits_u(i,0))) 
                    num_constraints += 1;
            }
            if (num_constraints!=pu)
                throw_exception(F("Please provide the correct number of input constraints!"));
        }                                

};

// Class which stores state equation matrices of a LTI State space model (full state access is assumed) 
template<int n, int m, typename DType = float>
class LTIModel {
    public:
        // constructor which takes in state and input matrix 
        LTIModel(const Matrix<n,n,DType>& A_, 
                 const Matrix<n,m,DType>& B_)
                    : A(A_),
                      B(B_){};
        
        // acess state matrix                  
        Matrix<n,n> get_A(){
            return A;
        }
        
        // acess input matrix                  
        Matrix<n,m> get_B(){
            return B;
        }        
        
        // iterate over state equation onece          
        Matrix<n,1,DType> simulate_time_step(const Matrix<n,1,DType>& x, const Matrix<m,1,DType>& u){
            return A*x+B*u; 
        }

    private:
        // State space matrices 
        Matrix<n,n,DType> A;
        Matrix<n,m,DType> B;
};

}