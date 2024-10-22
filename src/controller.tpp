#pragma once

namespace LMPC {

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::state_augmentation(Matrix<n,n,DType>& A, const Matrix<n,m,DType>& B){
    // QP in terms of relative controls is defined via state augmentation
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
}

template<int n, int m, int hx, int pu, int px, typename DType>
Matrix<m, 1, DType> lMPC<n, m, hx, pu, px, DType>::get_control(const Matrix<n,1, DType>& x){

    if (debugging){
        Serial.print("Free RAM when calling the controller: ");
        Serial.print(freeMemory());            
        Serial.println(" kB");    
    }
    // update constraint matrices
    if (pu+px!=0){
        update_constraints(x);
    }

    // pass linear coefficient matrix to solver 
    update_objective(x);
    
    auto sol_obj = QuadProg.solve();     

    Matrix<m, 1, DType> Delta_u; 
    for (int i = 0; i < m; i++){
        Delta_u(i,0) = sol_obj.x(i,0);
    } 
    u_op += Delta_u;
    return u_op;
} 

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::update_objective(const Matrix<n, 1, DType>& x){
    // build the augmented prediction matrices     
    BlockLowerTriangularToeplitzMatrix<na, m, hx, DType> Bp;
    build_prediction_matrices(Ap, Bp, Aaug, Baug);
    // pass q to the QP solver 
    // q = Bp_scaled*(Ap*(x && u_op) - xref) replaces  q = (~Bp*Q)*(Ap*(x && u_op) - xref);    
    QuadProg.update_q(Bp_scaled*(Ap*(x && u_op) - xref));        
}


template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::initialize_objective(){
    // build the augmented prediction matrices     
    BlockLowerTriangularToeplitzMatrix<na, m, hx, DType> Bp;
    build_prediction_matrices(Ap, Bp, Aaug, Baug);
    
    size_t c = 0;   
    size_t na = n+m;
    // replaces 
    // Matrix<mhx, mhx, DType> P = (~Bp*Q)*Bp+R;
    // : 
    Bp_scaled = ~Bp; 
    for (int i = 0; i < nahx; i++){
        DType qii;
        if (c<n){
            if (i < nahx-na){
                qii = Q(c);
            }
            else{
                qii = Qf(c);
            }
        }
        else {
            qii = 0.0;
        }
        if (c < na-1){
            c+=1;
        }
        else{
            c=0;        
        }
        for (int j = 0; j < mhx; j++){
            Bp_scaled(j,i) *= qii;
        }
    }
    Matrix<mhx, mhx, DType> P = Bp_scaled*Bp;
    for (int j = 0; j < mhx; j++){
        P(j,j) += R(c);
        if (c < m-1)
            c+=1;
        else
            c=0;        
    }
    QuadProg.update_Q(P);        
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::build_prediction_matrices(Matrix<nahx, na, DType>& Ap, 
                                                                      BlockLowerTriangularToeplitzMatrix<na, m, hx, DType>& Bp, 
                                                                      const Matrix<na, na, DType>& Aaug, 
                                                                      const Matrix<na, m, DType>& Baug){
    
    Ap.Fill(static_cast<DType>(0.0));
    
    DiagonalMatrix<na, DType> Id; 
    Id.diagonal.Fill(static_cast<DType>(1.0));
    
    Matrix<na, na, DType> PowA = Id;

    for (int i = 0; i < hx; i++){
        // diagonalBlocks_Bp[i] = PowA*Baug; 
        Bp.blocks[i] = PowA*Baug;
        PowA *= Aaug;         
        insert_at(Ap, PowA, i*na, 0);
    } 
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::build_prediction_matrices(Matrix<nhx, n, DType>& Ap, 
                                                                      BlockLowerTriangularToeplitzMatrix<n, m, hx, DType>& Bp, 
                                                                      const Matrix<n, n, DType>& A, 
                                                                      const Matrix<n, m, DType>& B){
    
    Ap.Fill(static_cast<DType>(0.0));
    
    DiagonalMatrix<n, DType> Id; 
    Id.diagonal.Fill(static_cast<DType>(1.0));
    
    Matrix<n, n, DType> PowA = Id;

    for (int i = 0; i < hx; i++){
        // diagonalBlocks_Bp[i] = PowA*Baug; 
        Bp.blocks[i] = PowA*B;
        PowA *= A;         
        insert_at(Ap, PowA, i*n, 0);
    } 
}


template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::update_constraints(const Matrix<n,1,DType>& x){

    Matrix<phx,1,DType> h = h_const;

    int rcount_total = 0;
    if (pu!=0){
        // upper limits
        update_input_constraints(rcount_total, 0, m-1, -u_op, h); 
        // lower limits
        update_input_constraints(rcount_total, m, m+m-1, u_op, h); 
    }
    if (px!=0){
        // state constraints 
        int count = 0;
        Matrix<mhx, 1, DType> L_uk;
        for (int i=0; i < mhx; i++){
            L_uk(i,0) = u_op(count,0);
            if (count < m-1)
                count+=1;
            else 
                count =m;
        }
        Matrix<nhx, 1, DType> limits_x_const_offset = -Ap_c*x-Bp_c*L_uk;
        // upper limits
        update_state_constraints(rcount_total, 0, n-1, limits_x_const_offset, h);         
        // lower limits
        update_state_constraints(rcount_total, n, n+n-1, limits_x_const_offset, h);                 
    }
    QuadProg.update_h(h);
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::update_state_constraints(int& abs_index, 
                                                                    int rel_index_start, 
                                                                    int rel_index_end, 
                                                                    const Matrix<nhx,1>& h_offset,
                                                                    Matrix<phx,1,DType>& h){
    int ncount = rel_index_start;
    for (int i=0; i < nhx; i++){
        if (!isnan(_limits_x(ncount,0))){
            h(abs_index,0) += h_offset(i,0);                       
            abs_index+=1;
        }
        if (ncount < rel_index_end)
            ncount+=1;
        else 
            ncount =rel_index_start;
    }    
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::update_input_constraints(int& abs_index, 
                                                                    int rel_index_start, 
                                                                    int rel_index_end, 
                                                                    const Matrix<m,1,DType>& h_offset,
                                                                    Matrix<phx,1,DType>& h){          
    int mcount = rel_index_start;
    for (int i=0; i < mhx; i++){
        if (!isnan(_limits_u(mcount,0))){
            h(abs_index,0) += h_offset(mcount,0);                       
            abs_index+=1;
        }
        if (mcount < rel_index_end)
            mcount+=1;
        else 
            mcount =rel_index_start;
    }
}


template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::initialize_state_constraints(int& abs_index, 
                                                                        int rel_index_start, 
                                                                        int rel_index_end, 
                                                                        const Matrix<n+n,1>& h_const_part,
                                                                        const Matrix<nhx, mhx, DType>& G_const_part,
                                                                        Matrix<phx,mhx,DType>& G){
    int ncount = rel_index_start;
    for (int i=0; i < nhx; i++){
        if (!isnan(_limits_x(ncount,0))){
            h_const(abs_index,0) = h_const_part(ncount,0);                       
            for (int j = 0; j <= mhx; j++){
                G(abs_index,j) = G_const_part(i,j);
            }
            abs_index+=1;
        }
        if (ncount < rel_index_end)
            ncount+=1;
        else 
            ncount =rel_index_start;
    }    
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::initialize_input_constraints(int& abs_index, 
                                                                        int rel_index_start, 
                                                                        int rel_index_end, 
                                                                        const Matrix<m+m,1>& h_const_part,
                                                                        const Matrix<mhx, mhx, DType>& G_const_part,
                                                                        Matrix<phx,mhx,DType>& G){
                                    
    int mcount = rel_index_start;
    for (int i=0; i < mhx; i++){
        if (!isnan(_limits_u(mcount,0))){
            h_const(abs_index,0) = h_const_part(mcount,0);                       
            for (int j = 0; j <= mhx; j++){
                G(abs_index,j) = G_const_part(i,j);
            }
            abs_index+=1;
        }
        if (mcount < rel_index_end)
            mcount+=1;
        else 
            mcount =rel_index_start;
    }
}

template<int n, int m, int hx, int pu, int px, typename DType>
void lMPC<n, m, hx, pu, px, DType>::initialize_constraints(){
    Matrix<phx,mhx,DType> G;

    DiagonalMatrix<m, DType> Id; 
    Id.diagonal.Fill(static_cast<DType>(1.0));

    BlockLowerTriangularToeplitzMatrix<m, m, hx, DType> CDelta;
    for (int i = 0; i < hx; i++){
        CDelta.blocks[i] = Id;
    }    
    int rcount_total = 0;

    if (pu!=0){
        // upper limits 
        initialize_input_constraints(rcount_total, 0, m-1,  _limits_u, -(-CDelta), G);        
        // lower limits
        initialize_input_constraints(rcount_total, m, m+m-1, -_limits_u, -CDelta, G);        
    }
    if (px!=0){
        build_prediction_matrices(Ap_c, Bp_c, A, B);
        Matrix<nhx, mhx, DType> Bp_CDelta = Bp_c*CDelta;

        // state constraints upper limit        
        initialize_state_constraints(rcount_total, 0, n-1, _limits_x, Bp_CDelta, G);
        // state constraints lower limit        
        initialize_state_constraints(rcount_total, n, n+n-1, -_limits_x, -Bp_CDelta, G);
    }

    // update static constraint matrix 
    QuadProg.update_G(G);
}

template<int n, int m, int hx, int pu, int px, typename DType>
lMPC<n, m, hx, pu, px, DType>::lMPC(const LTIModel<n, m, DType>& model_,
                                               const Matrix<n, 1, DType>& Q_,
                                               const Matrix<n, 1, DType>& Qf_,
                                               const Matrix<m, 1, DType>& R_,
                                               const Matrix<n, 1, DType>& xref_,
                                               const Constraints<n, m, pu, px, DType>& con,
                                               const QPparams params)
                                             : A(model_.get_A()),
                                               B(model_.get_B()),
                                               Q(Q_),
                                               Qf(Qf_),
                                               R(R_),
                                               _limits_u(con._limits_u),
                                               _limits_x(con._limits_x),
                                               debugging(params.debugging){
    
    if (debugging){
        Serial.print("Free RAM when building the controller: ");
        Serial.print(freeMemory());            
        Serial.println(" kB");    
    }

    QuadProg = QP<mhx,0,phx,DType>(params);
    
    // horizon wide reference 
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
    state_augmentation(A, B);

    u_op.Fill(static_cast<DType>(0.0));

    initialize_constraints();

    initialize_objective();

}


template<int n, int m, int hx, int pu, int px, typename DType>
lMPC<n, m, hx, pu, px, DType>::lMPC(const LTIModel<n, m, DType>& model_,
                                               const Matrix<n, 1, DType>& Q_,
                                               const Matrix<n, 1, DType>& Qf_,
                                               const Matrix<m, 1, DType>& R_,
                                               const Matrix<n, 1, DType>& xref_,
                                               const QPparams params)
                                             : A(model_.get_A()),
                                               B(model_.get_B()),
                                               Q(Q_),
                                               Qf(Qf_),
                                               R(R_),
                                               debugging(params.debugging){
    if (debugging){
        Serial.print("Free RAM when building the MPC: ");
        Serial.print(freeMemory());            
        Serial.println(" kB");    
    }

    QuadProg = QP<mhx,0,phx,DType>(params);
    
    // horizon wide reference 
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
    state_augmentation(A, B);

    u_op.Fill(static_cast<DType>(0.0));

    initialize_objective();    
}

}
