#pragma once

extern unsigned int __heap_start;
extern void *__brkval;

using namespace BLA;

namespace LMPC {


int freeMemory() {
  int free_memory;
  if ((int)__brkval == 0) {
    free_memory = ((int)&free_memory) - ((int)&__heap_start);
  } else {
    free_memory = ((int)&free_memory) - ((int)__brkval);
  }
  return free_memory;
}

template<int n, int m, typename DType>
void printMatrix(const Matrix<n, m, DType>& mat) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            Serial.print(mat(i, j));
            Serial.print(" ");
        }
        Serial.println();
    }
}

template<int m, int n, int k, int p, typename DType>
void insert_at(Matrix<m, n, DType>& M, const Matrix<k, p, DType>& SubM, int row, int col){
    for (int i = 0; i < k; i++){
        for (int j = 0; j < p; j++){
            M(row+i,col+j) = SubM(i,j); 
        }
    }
};

template<int m, int n, int k, typename DType>
void insert_at(Matrix<m, n, DType>& M, const DiagonalMatrix<k, DType>& SubM, int row, int col){
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k; j++){
            M(row+i,col+j) = SubM(i,j); 
        }
    }
};

template <int n, typename DType>
DiagonalMatrix<n, DType>::DiagonalMatrix() {
    diagonal.Fill(static_cast<DType>(0.0f));
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

template <int n, typename DType>
Matrix<n, n, DType> DiagonalMatrix<n, DType>::toMatrix() const {
    
    Matrix<n, n, DType> fullMatrix;
    fullMatrix.Fill(static_cast<DType>(0.0)); 
    
    for (int i = 0; i < n; i++) {
        fullMatrix(i, i) = diagonal(i); 
    }
    return fullMatrix;
}

template<int n, int m, int hx, typename DType>
void CondensedMPC<n, m, hx, DType>::update_QP_matrices(const Matrix<CondensedMPC<n, m, hx, DType>::na, 1, DType>& x,
                                                       const Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::na, DType>& Ap,
                                                       const Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::mhx, DType>& Bp,
                                                       const Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::nahx, DType>& Q, 
                                                       const Matrix<CondensedMPC<n, m, hx, DType>::mhx,CondensedMPC<n, m, hx, DType>::mhx, DType>& R){
    
    P = ~Bp*Q*Bp+R; 
    q = ~Bp*Q*(Ap*x - xref);    

    Serial.println("P");
    printMatrix(P);

    Serial.println("q");
    printMatrix(q);

    // pass quadratic and linear coefficient matrix to solver 
    QuadProg->update(P,q,{},{},{},{}); 
}

template<int n, int m, int hx, typename DType>
void CondensedMPC<n, m, hx, DType>::build_prediction_matrices(Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::na, DType>& Ap, 
                                                              Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::mhx, DType>& Bp, 
                                                              const Matrix<CondensedMPC<n, m, hx, DType>::na, CondensedMPC<n, m, hx, DType>::na, DType>& Aaug, 
                                                              const Matrix<CondensedMPC<n, m, hx, DType>::na, m, DType>& Baug){
    Ap.Fill(static_cast<DType>(0.0));
    Bp.Fill(static_cast<DType>(0.0));

    Matrix<na, na, DType> PowA = Aaug;
    Matrix<na, na, DType> PowAB;

    DiagonalMatrix<na, DType> Id; 
    Id.diagonal.Fill(static_cast<DType>(1.0));

    for (int i = 0; i < hx; i++){
        insert_at(Ap, PowA, static_cast<int>(i*na), 0);
        // Ap.Submatrix<n,n>(static_cast<int>(i*n),0) = PowA;

        PowA *= Aaug; 

        for (size_t j = 0; j < hx; j++){
            if(i>=j){                
                
                PowAB = Id.toMatrix();

                for (int k = 0; k < i-j; k++){
                    PowAB *= Aaug;

                }
                insert_at(Bp, PowAB*Baug, static_cast<int>(i*na), static_cast<int>(j*m));
                // Bp.Submatrix<n, m>(static_cast<int>(i*n), static_cast<int>(j*m)) = PowAB*B_;
            }
        }
    }
}

template<int n, int m, int hx, typename DType>
void CondensedMPC<n, m, hx, DType>::build_weight_matrices(Matrix<CondensedMPC<n, m, hx, DType>::nahx,CondensedMPC<n, m, hx, DType>::nahx, DType>& Q, 
                                                          Matrix<CondensedMPC<n, m, hx, DType>::mhx,CondensedMPC<n, m, hx, DType>::mhx, DType>& R, 
                                                          const Matrix<n, 1, DType>& Q_, 
                                                          const Matrix<n, 1, DType>& Qf_, 
                                                          const Matrix<m, 1, DType>& R_){
    Q.Fill(static_cast<DType>(0));
    R.Fill(static_cast<DType>(0));

    size_t c = 0;   
    size_t na = n+m;
    // Reference tracking weights 
    for (size_t i = 0; i < nahx-na; i++){
        if (c < n)
            Q(i,i) = Q_(c);
        else 
            Q(i,i) = static_cast<DType>(0.0);
        if (c < na-1)
            c+=1;
        else
            c=0;
    };

    // Reference tracking final weights 
    for (size_t i = nahx-na; i < nahx; i++){
        if (c < n)
            Q(i,i) = Qf_(c);
        else 
            Q(i,i) = static_cast<DType>(0.0);
        if (c < na-1)
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
CondensedMPC<n, m, hx, DType>::CondensedMPC(const LTIModel<n, m, DType>& model_,
                                            const Matrix<n, 1, DType>& Q_,
                                            const Matrix<n, 1, DType>& Qf_,
                                            const Matrix<m, 1, DType>& R_,
                                            const Matrix<n, 1, DType>& xref_)
                                             : A(model_.get_A()),
                                               B(model_.get_B()),
                                               QuadProg(new QP<mhx,0,0,DType>()), 
                                               alloc(true){
    Serial.print("free memory after when entering constructor: ");
    Serial.println(freeMemory());

    // horizon wide reference 
    xref.Fill(static_cast<DType>(0.0));

    size_t c = 0; 
    for (size_t i = 0; i < nahx; i++){
        if (c < n){
            xref(i,0) = xref_(c);
        } 
        else {
            xref(i,0) = static_cast<DType>(0.0);
        }
        if (c < na-1)
            c+=1;
        else
            c=0; 
    };

    Serial.print("free memory after filling xref: ");
    Serial.println(freeMemory());

    // Tracking and control weights   
    Matrix<nahx, nahx, DType> Q;
    Matrix<mhx, mhx, DType> R;

    build_weight_matrices(Q, R, Q_, Qf_, R_);

    Serial.print("free memory after building controller weights: ");
    Serial.println(freeMemory());

    // // QP in terms of relative controls is defined via state augmentation
    Matrix<na, na, DType> Aaug; 
    Matrix<na,  m, DType> Baug;

    DiagonalMatrix<m, DType> Id; 
    Id.diagonal.Fill(static_cast<DType>(1.0));

    Aaug.Fill(static_cast<DType>(0.0)); 
    Baug.Fill(static_cast<DType>(0.0)); 
    
    // augment state transition matrix 
    insert_at(Aaug, A, 0, 0); 
    insert_at(Aaug, B, 0, n);
    insert_at(Aaug, Id, n, n);

    // augment control matrix 
    insert_at(Baug, B, 0, 0);
    insert_at(Baug, Id, n, 0);

    Serial.print("free memory after building augmented matrices: ");
    Serial.println(freeMemory());

    Matrix<nahx,  na, DType> Ap;
    Matrix<nahx, mhx, DType> Bp;

    // // finally the 'augmented' prediction matrices     
    build_prediction_matrices(Ap, Bp, Aaug, Baug);

    Serial.print("free memory after building prediction matrices: ");
    Serial.println(freeMemory());

    // // initialize the QP at initial state [0.0, 0.0, ....]
    Matrix<na, 1, DType> xinit;
    xinit.Fill(static_cast<DType>(0.0));  

    update_QP_matrices(xinit, Ap, Bp, Q, R);
    Serial.print("free memory after updating qp: ");
    Serial.println(freeMemory());

}

template<int n, int m, int hx, typename DType>
Matrix<CondensedMPC<n, m, hx, DType>::mhx,
       CondensedMPC<n, m, hx, DType>::mhx, 
       DType> 
CondensedMPC<n, m, hx, DType>::get_P(){
    return P;
}

template<int n, int m, int hx, typename DType>
Matrix<CondensedMPC<n, m, hx, DType>::mhx, 
       1, 
       DType> 
CondensedMPC<n, m, hx, DType>::get_q(){
    return q;
}


}
