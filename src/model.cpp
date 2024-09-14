
#include "model.h"

namespace LMPC {

template<int n, int m>
LTIModel<n,m>::LTIModel(const Matrix<n,n>& A_, 
                        const Matrix<n,m>& B_, 
                        const double& h_)
                        : A(A_),
                          B(B_),
                          h(h_) {}; 

}
