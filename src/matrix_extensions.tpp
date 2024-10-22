#pragma once


namespace LMPC {

// build diagonal matrix with zeros on main diagonal
template <int n, typename DType>
DiagonalMatrix<n, DType>::DiagonalMatrix() {
    diagonal.Fill(static_cast<DType>(0.0f));
}

// build diagonal matrix with vector diag on main diagonal
template <int n, typename DType>
DiagonalMatrix<n, DType>::DiagonalMatrix(const Matrix<n, 1, DType>& diag) : diagonal(diag) {}

// overloading of function operator
template <int n, typename DType>
DType DiagonalMatrix<n, DType>::operator()(int row, int col) const {
    // If it's on the diagonal and it's not larger than the matrix dimensions then return the element
    if (row == col)
        return diagonal(row);
    else
        return static_cast<DType>(0.0f);
}

// initialize block lower triangular toeplitz matrix with given diagonal blocks
template<int nB, int mB, int numB, typename DType>
BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType>::BlockLowerTriangularToeplitzMatrix(const Matrix<nB, mB, DType> diagonalBlocks[numB]) {
    for (int i = 0; i < numB; ++i) {
        blocks[i] = diagonalBlocks[i];  
    };
};

// initialize block lower triangular toeplitz matrix with zero diagonal blocks
template<int nB, int mB, int numB, typename DType>
BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType>::BlockLowerTriangularToeplitzMatrix() {
    for (int i = 0; i < numB; ++i) {
        blocks[i].Fill(static_cast<DType>(0));  
    };
};

// convert block lower triangular toeplitz matrix to BLA matrix
template<int nB, int mB, int numB, typename DType>
Matrix<nB*numB,mB*numB, DType> BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType>::toMatrix(){
    Matrix<nB*numB,mB*numB, DType> Mat;
    Mat.Fill(static_cast<DType>(0));

    for (int i = 0; i < nB*numB; i++) {
        for (int j = 0; j <= i; j++){
            int nStart = i*nB;
            int mStart = j*mB; 

            for (int r = 0; r < nB; ++r) {
                for (int c = 0; c < mB; ++c) {
                    Mat(nStart + r, mStart + c) = blocks[i - j](r, c);
                }
            }
        }  
    };
    return Mat; 
};

// overloading of function operator
template<int nB, int mB, int numB, typename DType>
DType BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType>::operator()(int n, int m) const {
    int blockRow = n / nB;  // Determine the block index (row)
    int blockCol = m / mB;  // Determine the block index (column)

    // Ensure lower triangular structure: i >= j
    if (blockRow < blockCol) {
        return static_cast<DType>(0);  
    }

    // Compute the relative position inside the block
    int intraBlockRow = n % nB;
    int intraBlockCol = m % mB;

    return blocks[blockRow - blockCol](intraBlockRow, intraBlockCol);
}


// print a BLA matrix
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

// print a BlockLowerTriangularToeplitzMatrix
template<int nB, int mB, int numB, typename DType>
void printMatrix(const BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType> & mat) {
    for (int i = 0; i < nB * numB; ++i) {
        for (int j = 0; j < mB * numB; ++j) {
            Serial.print(mat(i, j));
            Serial.print(" ");
        }
        Serial.println();
    }
}

// insert submatrix at given row, col in matrix 
template<int m, int n, int k, int p, typename DType>
void insert_at(Matrix<m, n, DType>& M, const Matrix<k, p, DType>& SubM, int row, int col){
    for (int i = 0; i < k; i++){
        for (int j = 0; j < p; j++){
            M(row+i,col+j) = SubM(i,j); 
        }
    }
};

// insert diagonal submatrix at given row, col in matrix 
template<int m, int n, int k, typename DType>
void insert_at(Matrix<m, n, DType>& M, const DiagonalMatrix<k, DType>& SubM, int row, int col){
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k; j++){
            M(row+i,col+j) = SubM(i,j); 
        }
    }
};


}