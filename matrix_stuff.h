/*
 * matrix_stuff.h
 *
 *  Created on: Nov 8, 2022
 *      Author: norms
 */

// Forward declaration to use State_Getter class
class State_Getter;



#ifndef MATRIX_STUFF_H_
#define MATRIX_STUFF_H_

// Matrix block, represents only non zero values and saves time when matrices have a lot of zeros
class Matrix_Block
{
	public:
	// row and column number of origin. this is the offset
	int row;
	int col;
	// row and column lengths.
	// jlength is always the number of dimensions the particle is in. in this case, 3. ilength is however many constraints that are of the same type we have. Number of rows.
	int rows;
	int cols;
	// float array containing the matrix. matrix will be ilength * jlength. Number of columns of matrix. 1d array.
	float *matrix;
};

// Returns the transpose of a symmetric matrix
void matrix_transpose_symmetric(float *A, int rows, int cols, float *transposed_matrix);

// Returns the transpose of a nonsymmetric matrix
void matrix_transpose_nonsymmetric(float *A, int rows, int cols, float *matrix_transposed);

// Matrix multiplication
void matrix_multiplication_ROW_BY_ROW(float *A, float *B, int *A_dim, int *B_dim, float *result);

// Normal matrix multiplication
void matrix_multiplication(float *A, float *B, int *A_dim, int *B_dim, float *result);

// Take inverse of a diagonal matrix, represented as a vector
void inverse_diagonal_matrix(float *matrix, int length, float *result);

// Jacobian times the diagonal inverse mass matrix. Edits a list of resulting matrix blocks as the resulting m x 3n matrix is sparse.
void matrix_blocks_times_diagonal(Matrix_Block *block_list, int num_blocks, float *inverse_mass_matrix, Matrix_Block *result_blocks);

void diagonal_times_vector(float *diagonal, float *vector, int n, float *result);

// Multiplies a list of matrix blocks times a vector. Results in a vector
void matrix_blocks_times_vector(Matrix_Block *block_list, int num_blocks, float *vector, float *result);

// Multiplies the transpose of a list of matrix blocks times a vector. Results in a vector.
void matrix_blocks_transpose_times_vector(Matrix_Block *block_list, int num_blocks, float *vector, float *result);

// Function that either returns the magnitude of the vector, or the magnitude of the largest component.
// Used by linbcg, itol is the protocol
float calculate_relevant_norm(float *v, int n, int itol);

// Function that solves the Ax = B
void solve_for_x(int n, float *A_diagonal, float *b, float *x);

// Transpose matrix blocks
void matrix_blocks_transpose(Matrix_Block *block_list, int num_blocks, Matrix_Block *result_blocks);

// Function that calculates A*x: A is always symmetric in our problem. A= JM^-1JT
void A_times_x(State_Getter *state, float *x, float *result);

// zeros all elements in array
void zero_vector(float  *vector, int length);

// Fills array with ones
void ones_vector(float *vector, int length);

#endif /* MATRIX_STUFF_H_ */
