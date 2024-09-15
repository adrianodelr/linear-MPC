#pragma once

using namespace BLA;

namespace LMPC {

template<int n, int m, int hx>
ControllerWeights<n,m,hx>::ControllerWeights(const Matrix<n>& Q_, 
                                             const Matrix<n>& Qf_, 
                                             const Matrix<m>& R_){
    Q.Fill(0.0);
    R.Fill(0.0);

    size_t c = 0;   
    // Reference tracking weights 
    for (size_t i = 0; i < nhx-n; i++){
        Q(i,i) = Q_(c);
        if (c < n-1) {
            c+=1;
        }
        else {
            c=0;
        }
    }

    // Reference tracking final weights 
    for (size_t i = nhx-n; i < nhx; i++){
        Q(i,i) = Qf_(c);
        if (c < n-1) {
            c+=1;
        }
        else {
            c=0;
        }
    }

    // Control effort weights 
    for (size_t i = 0; i < mhx; i++){
        R(i,i) = R_(c);
        if (c < m-1) {
            c+=1;
        }
        else {
            c=0;
        }
    }
};

template<int n, int m, int hx>
void ControllerWeights<n,m,hx>::print_weights(){
    Serial.println("Reference tracking: ");
    Serial.println(Q);
    Serial.print("Applied efforts: ");
    Serial.println(R);
}

template<int n, int m, int hx>
PredictionMatrices<n, m, hx>::PredictionMatrices(const Matrix<n,n>& A_, 
                                                 const Matrix<n,m>& B_){
    Ap.Fill(0);
    Bp.Fill(0);

    Matrix<n,n> PowA = A_;
    Matrix<n,n> PowAB;
    PowAB.Fill(0.0);

    for (int i = 0; i < hx; i++){
        
        // Ap.Submatrix<n,n>(static_cast<int>(i*n),0) = PowA;

        RefMatrix<BLA::Matrix<nhx, n>, 4, 4> submatrix(Ap.Submatrix<4, 4>(0, 0));

        PowA *= A_; 
        
        for (size_t j = 0; j < hx; j++){
            if(i>=j){                
                
                PowAB = A_;
                for (int k = 1; k < i-j; i++){
                    PowAB *= A_;
                }
                // Bp.Submatrix<n, m>(static_cast<int>(i*n), static_cast<int>(j*m)) = PowAB*B_;
            }
        }
    }                                                         

}; 

// template<int n, int m, int hx>
// Matrix<nhx,n> PredictionMatrices<n, m, hx>::get_predmat_A(){
//     return Ap;
// }
// template<int n, int m, int hx>
// Matrix<nhx,mhx> PredictionMatrices<n, m, hx>::get_predmat_B(){
//     return Bp;
// }
}
