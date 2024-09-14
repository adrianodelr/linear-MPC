#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace LMPC {

template<int n, int m>
class LTIModel {
    public:
        LTIModel(const Matrix<n,n>& A_, const Matrix<n,m>& B_, const double& h_);
    private:
        // State space matrices 
        Matrix<n,n> A;
        Matrix<n,m> B;
        // discretization step size 
        double h;
}

}