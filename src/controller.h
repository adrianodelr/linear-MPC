
#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"

using namespace BLA;

namespace LMPC {

template<int n, int m, int hx>
class ControllerWeights {
    public:
        ControllerWeights(const Matrix<n,1>& Q_, 
                          const Matrix<n,1>& Qf_, 
                          const Matrix<m,1>& R_);
    private:
        static constexpr size_t nhx = (size_t)n*hx;
        static constexpr size_t mhx = (size_t)n*hx; 
        // Tracking and control weights   
        Matrix<nhx,nhx> Q;
        Matrix<mhx,mhx> R;
};

}