#pragma once

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
        Serial.println(c);  
        Q(i,i) = Q_(c);
        if (c < n-1) c+=1;
        else c=0;
    }
    // Reference tracking final weights 
    for (size_t i = nhx-n; i < nhx; i++){
        Q(i,i) = Qf_(c);
        if (c < n-1) c+=1;
        else c=0;
    }
    // Control effort weights 
    for (size_t i = 0; i < mhx; i++){
        R(i,i) = R_(c);
        if (c < m-1) c+=1;
        else c=0;
    }
};

template<int n, int m, int hx>
void ControllerWeights<n,m,hx>::print_weights(){
    Serial.println("Reference tracking: ");
    Serial.println(Q);
    Serial.print("Applied efforts: ");
    Serial.println(R);
}

}
