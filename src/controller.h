
#include <Arduino.h> 
#include <BasicLinearAlgebra.h>
#include "model.h"

namespace LMPC {

template<int n, int m, int hx>
class ControllerWeights {
    public:
        ControllerWeights(const Matrix<n,1>& Q_, 
                          const Matrix<n,1>& Qf_, 
                          const Matrix<m,1>& R_, 
                          const size_t& hx);
    private:
        const size_t nhx = (size_t)n*hx;
        const size_t mhx = (size_t)n*hx; 
        // Tracking and control weights   
        Matrix<nhx,nhx> Q;
        Matrix<mhx,mhx> R;
}

}