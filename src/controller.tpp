#pragma once

using namespace BLA;

namespace LMPC {

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
void insert_at(Matrix<m, n, DType>& M, Matrix<k, p, DType> SubM, int row, int col){
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

template <int n, typename DType>
BLA::Matrix<n, n, DType> DiagonalMatrix<n, DType>::toMatrix() const {
    
    BLA::Matrix<n, n, DType> fullMatrix;
    fullMatrix.Fill(0.0); 
    
    for (int i = 0; i < n; i++) {
        fullMatrix(i, i) = diagonal(i); 
    }
    return fullMatrix;
}

template<int n, int m, int hx, typename DType>
ControllerWeights<n, m, hx, DType>::ControllerWeights(const Matrix<n, 1, DType>& Q_, 
                                                      const Matrix<n, 1, DType>& Qf_, 
                                                      const Matrix<m, 1, DType>& R_){
    Q.Fill(0.0);
    R.Fill(0.0);

    size_t c = 0;   
    size_t na = n+m;
    Serial.print("nahx");
    Serial.println(nahx);

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
void ControllerWeights<n,m,hx, DType>::print_weights(){
    Serial.println("Reference tracking: ");
    Serial.println(Q);
    Serial.print("Applied efforts: ");
    Serial.println(R);
}

template<int n, int m, int hx, typename DType>
Matrix<ControllerWeights<n, m, hx, DType>::nahx, 
       ControllerWeights<n, m, hx, DType>::nahx, 
       DType> 
ControllerWeights<n, m, hx, DType>::get_Q(){
    return Q;
}

template<int n, int m, int hx, typename DType>
Matrix<ControllerWeights<n, m, hx, DType>::mhx, 
       ControllerWeights<n, m, hx, DType>::mhx, 
       DType>  
ControllerWeights<n, m, hx, DType>::get_R(){
    return R;
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
Matrix<PredictionMatrices<n, m, hx, DType>::nhx,
       n, 
       DType> 
PredictionMatrices<n, m, hx, DType>::get_predmat_A(){
    return Ap;
}

template<int n, int m, int hx, typename DType>
Matrix<PredictionMatrices<n, m, hx, DType>::nhx, 
       PredictionMatrices<n, m, hx, DType>::mhx, 
       DType> 
PredictionMatrices<n, m, hx, DType>::get_predmat_B(){
    return Bp;
}

template<int n, int m, int hx, typename DType>
void CondensedMPC<n, m, hx, DType>::update_QP_matrices(const Matrix<static_cast<int>(n+m), 1, DType>& x){
    const auto Ap = PredMat -> get_predmat_A();
    const auto Bp = PredMat -> get_predmat_B(); 
    const auto Q = weights -> get_Q();
    const auto R = weights -> get_R();

    Matrix<mhx, mhx, DType> BQBR = ~Bp*Q*Bp+R; 
    Matrix<mhx, mhx, DType> R_ = R;  
    

    Serial.print("~Bp*Q*Bp");    
    printMatrix(~Bp*Q*Bp); 

    Serial.print("R");    
    printMatrix(R); 

    Serial.print("~Bp*Q*Bp(0,0)+R(0,0): ");    
    Serial.println((~Bp*Q*Bp)(0,0)+R(0,0));    

    Serial.print("~Bp*Q*Bp+R: ");    
    printMatrix((~Bp*Q*Bp)+R);    


    // Matrix<mhx,   1, DType> q_ = ~Bp*Q*(Ap*x - xref);
    // P = ~Bp*Q*Bp + R; 
    // q = ~Bp*Q*(Ap*x - xref);

    // Serial.print("R: ");    
    // Serial.println(R);

    // Serial.print("Bp: ");    
    // Serial.println(Bp);

    // Serial.print("Q: ");    
    // Serial.println(Q);

    // Serial.print("Ap: ");    
    // Serial.println(Ap);

    // Serial.print("xref: ");    
    // Serial.println(xref);

    // Serial.print("x: ");    
    // Serial.println(x);

    // Serial.print("P_");    
    // printMatrix(P_); 




    // Serial.print("~Bp*Q*Bp: ");    
    // Serial.println(~Bp*Q*Bp);
    
    // Serial.print("P_: ");    
    // Serial.println(P_);
    // Serial.print("q: ");    
    // Serial.println(q);

    // pass quadratic and linear coefficient matrix to solver 
    // QuadProg->update(P,q,{},{},{},{}); 
}

template<int n, int m, int hx, typename DType>
CondensedMPC<n, m, hx, DType>::CondensedMPC(const LTIModel<n, m>& model_,
                                            const Matrix<n, 1, DType>& Q_,
                                            const Matrix<n, 1, DType>& Qf_,
                                            const Matrix<m, 1, DType>& R_,
                                            const Matrix<n, 1, DType>& xref_)
                                             : model(model_),
                                               QuadProg(new QP<static_cast<int>(m*hx),0,0>()), 
                                               alloc(true){
    
    weights = new ControllerWeights<n, m, hx>(Q_, Qf_, R_);

    // horizon wide reference 
    xref.Fill(0.0);
    size_t c = 0; 
    Serial.println("im here0");
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
    Serial.println("im here1");
    // QP in terms of relative controls is defined via state augmentation
    Matrix<na, na, DType> Aaug; 
    Matrix<na,  m, DType> Baug;

    DiagonalMatrix<m, DType> Id; 
    Id.diagonal.Fill(1.0);

    Aaug.Fill(0.0); 
    Baug.Fill(0.0); 
    
    // augment state transition matrix 
    insert_at(Aaug, model.get_A(), 0, 0); 
    insert_at(Aaug, model.get_B(), 0, n);
    insert_at(Aaug, Id.toMatrix(), m, n);

    // augment control matrix 
    insert_at(Baug, model.get_B(), 0, 0);
    insert_at(Baug, Id.toMatrix(), m, 0);
    Serial.println("im here2");

    // finally the 'augmented' prediction matrices 
    PredMat = new PredictionMatrices<na, m, hx, DType>(Aaug,Baug);
    Serial.println("im here3");

    // initialize the QP at initial state [0.0, 0.0, ....]
    Matrix<na, 1, DType> xinit;
    xinit.Fill(0.0);  
    Serial.println("im here4");

    update_QP_matrices(xinit);
    
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
