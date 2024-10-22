#pragma once

#include <Arduino.h> 
#include <BasicLinearAlgebra.h>

namespace LMPC {

// economic diagonal matrices
template <int n, typename DType = float>
struct DiagonalMatrix : public MatrixBase<DiagonalMatrix<n>, n, n, DType> {
    Matrix<n, 1, DType> diagonal;
    
    // build diagonal matrix with zeros on main diagonal
    DiagonalMatrix();

    // build diagonal matrix with vector diag on main diagonal
    DiagonalMatrix(const Matrix<n, 1, DType>& diag);

    // overloading of function operator
    DType operator()(int row, int col) const;  
};

// economic lower triangular toeplitz matrices in block form
template<int nB, int mB, int numB, typename DType = float>
class BlockLowerTriangularToeplitzMatrix : public BLA::MatrixBase<BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType>, nB*numB, mB*numB, DType> {
    public:
        // initialize diagonal blocks
        BlockLowerTriangularToeplitzMatrix(const BLA::Matrix<nB, mB, DType> diagonalBlocks[numB]);

        // initialize diagonal blocks with zeros
        BlockLowerTriangularToeplitzMatrix();

        // overloading of function operator
        DType operator()(int n, int m) const; 
        // allocate list of matrices, each repeated along a block diagonal band 
        // first block is on the main diagonal 
        Matrix<nB, mB, DType> blocks[numB];

        // returns a BLA matrix 
        Matrix<nB*numB,mB*numB, DType> toMatrix();  
};

// print a BLA matrix
template<int n, int m, typename DType>
void printMatrix(const Matrix<n, m, DType>& mat);

// print a BlockLowerTriangularToeplitzMatrix
template<int nB, int mB, int numB, typename DType>
void printMatrix(const BlockLowerTriangularToeplitzMatrix<nB, mB, numB, DType> & mat);

// insert submatrix at given row, col in matrix 
template<int m, int n, int k, int p, typename DType = float>
void insert_at(Matrix<m, n, DType>& M, const Matrix<k, p, DType>& SubM, int row, int col);

// insert diagonal submatrix at given row, col in matrix 
template<int m, int n, int k, typename DType>
void insert_at(Matrix<m, n, DType>& M, const DiagonalMatrix<k, DType>& SubM, int row, int col);

}

#include "matrix_extensions.tpp"