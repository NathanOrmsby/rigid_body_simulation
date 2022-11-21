#include "vectors.h"
#include "rigid_bodies.h"
#include <vector>
#include <iostream>
#include "matrix_stuff.h"
#include <math.h>

#include "get_state.h"

// THIS IS IT, THE COOLEST SHIT I'VE EVER CODED. LINEAR ALGEBRA IN C++

// Define a matrix block structure. These will be relevant in jacobian and jacobian time derivative matrices
// Each matrix block will pertain to a type of constraint
// The idea is to utilize these in sparse matrices so time is not wasted on zero elements

void matrix_transpose_symmetric(double *A, int rows, int cols, double *matrix_transposed)
{
	// Rows is number of rows, cols is number of cols. Gives dimension of matrix
	// Swap the row indexes with the col indexes
	// First and last elements always the same
	matrix_transposed[0] = A[0];
	matrix_transposed[rows * cols - 1] = A[rows * cols - 1];
	for (int i = 1; i < rows * cols - 1; i++)
	{
		// Find which row im in
		int row_num = i / cols;
		// Find column im in. Current position minus first position in the row
		int col_num = i - row_num * cols;

		// Find swapped position and put it in new transposed matrix
		matrix_transposed[col_num * cols + row_num] = A[i];

	}
}

void matrix_transpose_nonsymmetric(double *A, int rows, int cols, double *matrix_transposed)
{
	// Assymmetric transpose has opposite dimensions
	int new_cols = rows;
	int new_rows = cols;
	// First and last are the same
	matrix_transposed[0] = A[0];
	matrix_transposed[new_rows * new_cols - 1] = A[rows * cols - 1];

	for (int i = 1; i < rows * cols - 1; i++)
	{
		// Find which row im in
		int row_num = i / cols;
		// Find column im in. Current position minus first position in the row
		int col_num = i - row_num * cols;

		// Find swapped position and put it in new transposed matrix, which has a different amount of columns
		matrix_transposed[col_num * new_cols + row_num] = A[i];
	}
}

// Given two matrices and dimensions, perform matrix multiplication using 1d arrays.
// THIS IS ROW TO ROW MULTIPLICATION, SO TAKING A TRANSPOSE OF THE SECOND MATRIX IS ESSENTIAL
void matrix_multiplication_ROW_BY_ROW(double *A, double *B, int *A_dim, int *B_dim, double *result)
{
	// Check if matrix multiplication is possible.
	if (A_dim[1] != B_dim[0])
	{
		std::cout << "Columns in A: " << A_dim[1] << " DO NOT MATCH ROWS IN B: " << B_dim[0] << "! MATRIX MULTIPLICATION IS NOT POSSIBLE!" << std::endl;
		return;
	}

	// Get the transpose of B. makes life easier. This isnt actually matrix multiplication, so this is fine.
	double B_T[B_dim[0] * B_dim[1]];
	int B_T_dim[2];
	// If the matrix is symmetric
	if (B_dim[0] == B_dim[1])
	{
		matrix_transpose_symmetric(B, B_dim[0], B_dim[1], B_T);
		B_T_dim[0] = B_dim[0];
		B_T_dim[1] = B_dim[1];
	}
	else
	{
		matrix_transpose_nonsymmetric(B, B_dim[0], B_dim[1], B_T);
		B_T_dim[0] = B_dim[1];
		B_T_dim[1] = B_dim[0];
	}

	// Perform row to row multiplication
	// Iterate through number of rows in matrix A
	int count = 0;
	for (int row_A = 0; row_A < A_dim[0]; row_A++)
	{
		// Array of row elements in this row.
		double row_A_elements[A_dim[1]];
		// Fill up row elements
		for (int i = 0; i < A_dim[1]; i++)
		{
			// Same as formula used in transpose functions
			row_A_elements[i] = A[row_A * A_dim[1] + i];
		}
		for (int row_B_T = 0; row_B_T < B_T_dim[0]; row_B_T++)
		{
			// Define the row elements to be multiplied by
			// Array of row elements in B
			double row_B_T_elements[B_T_dim[1]];
			for (int i = 0; i < B_T_dim[1]; i++)
			{
				row_B_T_elements[i] = B_T[row_B_T * B_T_dim[1] + i];
			}

			// element by element multiplication. A_cols and B_transpose cols are the same
			double sum = 0;
			for (int i = 0; i < A_dim[1]; i++)
			{
				sum += row_A_elements[i] * row_B_T_elements[i];
			}

			// Add the sum to the new matrix
			result[count] = sum;
			count++;
		}
	}
}

// Normal matrix multiplication
void matrix_multiplication(double *A, double *B, int *A_dim, int *B_dim, double *result)
{
	// Check if matrix multiplication is possible.
	if (A_dim[1] != B_dim[0])
	{
		std::cout << "Columns in A: " << A_dim[1] << " DO NOT MATCH ROWS IN B: " << B_dim[0] << "! MATRIX MULTIPLICATION IS NOT POSSIBLE!" << std::endl;
		return;
	}

	// No transpose taken

	// Perform row to row multiplication
	// Iterate through number of rows in matrix A
	int count = 0;
	for (int row_A = 0; row_A < A_dim[0]; row_A++)
	{
		// Array of row elements in this row.
		double row_A_elements[A_dim[1]];
		// Fill up row elements
		for (int i = 0; i < A_dim[1]; i++)
		{
			// Same as formula used in transpose functions
			row_A_elements[i] = A[row_A * A_dim[1] + i];
		}
		for (int col_B = 0; col_B < B_dim[1]; col_B++)
		{
			// Define the column elements to be multiplied by
			// Array of column elements in B. As many as there are rows
			double col_B_elements[B_dim[0]];
			for (int i = 0; i < B_dim[0]; i++)
			{
				col_B_elements[i] = B[i * B_dim[1] + col_B];
			}

			// element by element multiplication. A_cols and B_transpose cols are the same
			double sum = 0;
			for (int i = 0; i < A_dim[1]; i++)
			{
				sum += row_A_elements[i] * col_B_elements[i];
			}

			// Add the sum to the new matrix
			result[count] = sum;
			count++;
		}
	}
}

// Take inverse of a diagonal matrix, represented as a vector
void inverse_diagonal_matrix(double *matrix, int length, double *result)
{
	// Inverse of a diagonal matrix is just the reciprocal of every element
	for (int i = 0; i < length; i++)
	{
		result[i] = 1 / matrix[i];
	}
}

// Jacobian times the diagonal inverse mass matrix. Edits a list of resulting matrix blocks as the resulting m x 3n matrix is sparse.
void matrix_blocks_times_diagonal(Matrix_Block *block_list, int num_blocks, double *diagonal, Matrix_Block *result_blocks)
{
	// Result blocks have same offset as input blocks
	// Loop through input blocks
	for (int i = 0; i < num_blocks; i++)
	{
		// Resulting block matrix is same length as input block matrix
		result_blocks[i].row = block_list[i].row;
		result_blocks[i].col = block_list[i].col;
		result_blocks[i].rows = block_list[i].rows;
		result_blocks[i].cols = block_list[i].cols;


		// Perform the multiplication. Multiply element in input block matrix times element in the mass matrix offset by the column of the element in the block matrix
		// The column of the block element in the entire matrix is its origin column plus its relative column
		for (int j = 0; j < block_list[i].rows * block_list[i].cols; j++)
		{
			// Relative row number
			int row_num = j / block_list[i].cols;
			// Find column im in IN THIS BLOCK. Current position minus first position in the row IN THIS BLOCK.
			int col_num = j - row_num * block_list[i].cols;

			result_blocks[i].matrix[j] = block_list[i].matrix[j] * diagonal[block_list[i].col + col_num];
		}
	}
}

// Diagonal matrix times a vector. Just element by element multiplication, diagonal matrix represented as a vector.
// n is the size
void diagonal_times_vector(double *diagonal, double *vector, int n, double *result)
{
	for (int i = 0; i < n; i++)
	{
		result[i] = diagonal[i] * vector[i];
	}
}

// Multiplies a list of matrix blocks times a vector. Results in a vector
void matrix_blocks_times_vector(Matrix_Block *block_list, int num_blocks, double *vector, double *result)
{
	// Loop through input blocks
	// Loop through input blocks
	for (int i = 0; i < num_blocks; i++)
	{
		// Assume vector is already zeroed
		// Each matrix block represents a certain number of rows and columns in the matrix
		for (int j = 0; j < block_list[i].rows * block_list[i].cols; j++)
		{

			// Relative row number
			// Find column im in IN THIS BLOCK. Current position minus first position in the row IN THIS BLOCK.
			// The indice of the resulting vector that gets edited is the same as the row of the element in the matrix
			result[block_list[i].row + (j / block_list[i].cols)] += block_list[i].matrix[j] * vector[block_list[i].col + (j - (j / block_list[i].cols) * block_list[i].cols)];
		}
	}
}

// Multiplies the transpose of a list of matrix blocks times a vector. Results in a vector.
void matrix_blocks_transpose_times_vector(Matrix_Block *block_list, int num_blocks, double *vector, double *result)
{
	// Loop through input blocks
	for (int i = 0; i < num_blocks; i++)
	{
		// Perform the multiplication. Multiply element in input block matrix times element in the mass matrix offset by the column of the element in the block matrix
		// Assume vector is already zeroed
		// The column of the block element in the entire matrix is its origin column plus its relative column
		for (int j = 0; j < block_list[i].rows * block_list[i].cols; j++)
		{
			// This is the transpose, reverse the row and column numbers
			// Also reverse the row and col origin in result
			// Relative row number.
			// Find column im in IN THIS BLOCK. Current position minus first position in the row IN THIS BLOCK.

			// Instead of editing the current element in a block matrix, we edit the current position in the vector
			result[block_list[i].col + (j - (j / block_list[i].cols)* block_list[i].cols)] += block_list[i].matrix[j] * vector[block_list[i].row + (j / block_list[i].cols)];
		}
	}
}

// NOT USED
// Transpose matrix blocks
void matrix_blocks_transpose(Matrix_Block *block_list, int num_blocks, Matrix_Block *result_blocks)
{
	for (int i = 0; i < num_blocks; i++)
	{
		// Swap the row and col origin
		int tmp = result_blocks[i].row;
		result_blocks[i].row = block_list[i].col;
		result_blocks[i].col = tmp;

		// Swap the dimensions of the block matrix
		tmp = result_blocks[i].rows;
		result_blocks[i].rows = block_list[i].cols;
		result_blocks[i].cols = tmp;
	}
}

// NOT USED
//Example function for row reduction
void row_reduction(double matrix[][4], int rows, int cols)
{
	// lead position in the row, this is what gets set to 1
	int lead = 0;

	// Assume matrix is the augmented matrix form. Ax = b so the last column is b.
	while (lead < cols - 1)
	{
		// Loop through the rows
		for (int row = 0; row < rows; row++)
		{
			// Calculate divisor and multiplier
			double d = matrix[lead][lead];
			double m = matrix[row][lead];

			// Loop through the columns
			for (int col = 0; col < cols; col++)
			{
				// If we're at a row with an identity, then just divide to make the identity 1
				if (row == lead)
				{
					matrix[row][col] /= d;
				}
				else
				{
					matrix[row][col] -= (m / d) * matrix[lead][col];
				}
			}
			lead++;
		}
	}
}

// Function that either returns the magnitude of the vector, or the magnitude of the largest component.
// Used by linbcg, itol is the protocol
double calculate_relevant_norm(double *v, int n, int itol)
{
	// If itol is 1, just calculate the magnitude of the vector
	if (itol == 1)
	{
		double magnitude = 0.0;
		for (int i = 0; i < n; i++)
		{
			magnitude += (v[i] * v[i]);
		}
		return sqrt(magnitude);
	}
	// Otherwise return the magnitude of the largest component
	else
	{
		int max_index = 0;
		for (int i = 0; i < n; i++)
		{
			if (fabs(v[i]) > fabs(v[max_index]))
			{
				max_index = i;
			}
		}
		return fabs(v[max_index]);
	}
}

// Function that calculates A*x: A is always symmetric in our problem. A= J* M^-1 *J_T
// Symmetric means A_T = A
void A_times_x(State_Getter *state, double *x, double *result)
{
	// Zero the vector
	zero_vector(result, state->num_constraints);
	// First multiply the jacobian transpose times x
	double first_res[state->num_bodies * 3];
	zero_vector(first_res, state->num_bodies * 3);
	matrix_blocks_transpose_times_vector(state->jacobian, state->num_constraints, x, first_res);

	// Then multiply mass inverse times first_res
	double second_res[state->num_bodies * 3];
	zero_vector(second_res, state->num_bodies * 3);
	diagonal_times_vector(state->inverse_mass_matrix, first_res, 3 * state->num_bodies, second_res);

	// Then multiply the jacobian times the second result, and thats it
	matrix_blocks_times_vector(state->jacobian, state->num_constraints, second_res, result);
}

// Function that solves the A_diagonal * x = B.
void solve_for_x(int n, double *A_diagonal, double *b, double *x)
{
	// n is the length of the arrays
	// A_diagonal is the diagonal part of A. We just plug in identity here because we dont calc A.
	// Plugging in identity just means dividing by 1, so its just assigning shit right now.
	// If the diagonal is equal to zero, then the x component is equal to the b component
	for (int i = 0; i < n + 1; i++)
	{
		if (A_diagonal[i] == 0.0)
		{
			x[i] = b[i];
		}
		else
		{
			x[i] = b[i] / A_diagonal[i];
		}
	}
}

// zeros all elements in array
void zero_vector(double *vector, int length)
{
	for (int i = 0; i < length; i++)
	{
		vector[i] = 0;
	}
}

void ones_vector(double *vector, int length)
{
	for (int i = 0; i < length; i++)
	{
		vector[i] = 1;
	}
}


