#include "controller.h"

namespace LMPC {

template<int n, int m, int hx>
ControllerWeights<n,m,hx>::ControllerWeights(const Matrix<n,1>& Q_, 
                                             const Matrix<n,1>& Qf_, 
                                             const Matrix<m,1>& R_){
        Q.Fill(0.0);
        R.Fill(0.0);
        
        size_t c = 0;   
        // Reference tracking weights 
        for (size_t i = 0; i < nhx-n; i++){
            Q(i,i) = Q_(c,1);
            if (c==n-1) c=0;  
        }
        // Reference tracking final weights 
        for (size_t i = nhx-n; i < nhx; i++){
            Q(i,i) = Qf_(c,1);
            if (c==n-1) c=0;  
        }
        // Control effort weights 
        for (size_t i = 0; i < mhx; i++){
            R(i,i) = R_(c,1);
            if (c==m-1) c=0;  
        }
    }

}