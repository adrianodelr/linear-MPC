#pragma once

using namespace BLA;

namespace LMPC {

template<int m, int n, int k, int p, typename DType>
void insert_at(Matrix<m,n,DType>& M, Matrix<k,p,DType> SubM, int row, int col){
    for (int i = 0; i < k; i++){
    
        for (int j = 0; j < p; j++){
            M(row+i,col+j) = SubM(i,j); 
    
        }
    }
};

template <int n, typename DType>
DiagonalMatrix<n, DType>::DiagonalMatrix() {
    diagonal.Fill(0.0f);
}

template <int n, typename DType>
DiagonalMatrix<n, DType>::DiagonalMatrix(const Matrix<n, 1, DType>& diag) : diagonal(diag) {}

template <int n, typename DType>
DType DiagonalMatrix<n, DType>::operator()(int row, int col) const {
    // If it's on the diagonal and it's not larger than the matrix dimensions then return the element
    if (row == col)
        return diagonal(row);
    else
        return static_cast<DType>(0.0f);
}

template<int n, int m, int hx, typename DType>
ControllerWeights<n, m, hx, DType>::ControllerWeights(const Matrix<n, 1, DType>& Q_, 
                                                      const Matrix<n, 1, DType>& Qf_, 
                                                      const Matrix<m, 1, DType>& R_){
    Q.Fill(0.0);
    R.Fill(0.0);

    size_t c = 0;   
    // Reference tracking weights 
    for (size_t i = 0; i < nhx-n; i++){
        Q(i,i) = Q_(c);
        if (c < n-1)
            c+=1;
        else
            c=0;
    };

    // Reference tracking final weights 
    for (size_t i = nhx-n; i < nhx; i++){
        Q(i,i) = Qf_(c);
        if (c < n-1)
            c+=1;
        else
            c=0;
    };

    // Control effort weights 
    for (size_t i = 0; i < mhx; i++){
        R(i,i) = R_(c);
        if (c < m-1)
            c+=1;
        else
            c=0;
    };
};

template<int n, int m, int hx, typename DType>
void ControllerWeights<n,m,hx,DType>::print_weights(){
    Serial.println("Reference tracking: ");
    Serial.println(Q);
    Serial.print("Applied efforts: ");
    Serial.println(R);
}

template<int n, int m, int hx, typename DType>
PredictionMatrices<n, m, hx, DType>::PredictionMatrices(const Matrix<n, n, DType>& A_, 
                                                        const Matrix<n, m, DType>& B_){
    Ap.Fill(0);
    Bp.Fill(0);

    Matrix<n, n, DType> PowA = A_;
    Matrix<n, n, DType> PowAB;
    PowAB.Fill(0.0);

    DiagonalMatrix<n, DType> Id; 
    Id.diagonal.Fill(1.0);

    for (int i = 0; i < hx; i++){

        insert_at(Ap, PowA, static_cast<int>(i*n), 0);
        // Ap.Submatrix<n,n>(static_cast<int>(i*n),0) = PowA;

        PowA *= A_; 
        
        for (size_t j = 0; j < hx; j++){
            if(i>=j){                
                
                PowAB = Id;

                for (int k = 0; k < i-j; k++){
                    PowAB *= A_;
                }
                insert_at(Bp, PowAB*B_, static_cast<int>(i*n), static_cast<int>(j*m));
                // Bp.Submatrix<n, m>(static_cast<int>(i*n), static_cast<int>(j*m)) = PowAB*B_;
            }
        }
    }
}; 

template<int n, int m, int hx, typename DType>
Matrix<static_cast<int>(n*hx),n, DType> PredictionMatrices<n, m, hx, DType>::get_predmat_A(){
    return Ap;
}

template<int n, int m, int hx, typename DType>
Matrix<static_cast<int>(n*hx),static_cast<int>(m*hx), DType> PredictionMatrices<n, m, hx, DType>::get_predmat_B(){
    return Bp;
}
}
