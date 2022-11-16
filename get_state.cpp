#include <iostream>
#include "matrix_stuff.h"
#include "constraint_bodies.h"
#include "rigid_bodies.h"
#include "vectors.h"
#include "utils.h"
#include "get_state.h"
#include "bc_gradient.h"
#include "springs.h"

// Functions involved with calculating the state of the system

void State_Getter::get_current_state(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, Spring_2 *spring2s)
{
	// Results in a net force vector for the current state assigned to the class.
	// Net force vector is the sum of the external force vector and the constraint force vector.

	// *** Creating the matrices to start ***
	// Create state vector and state vector derivative
	create_state_vector(mass_list);
	create_state_vector_derivative(mass_list);
	// Create inverse mass matrix
	create_inverse_mass_matrix(mass_list);
	// Create external force vector
	create_force_ext_vector(mass_list, spring2s);
//	std::cout << "External force vector " << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << force_ext_vector[i] << " ";
//	}
//	std::cout << std::endl;
	// Create jacobian
	create_jacobian(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
	// Create jacobian derivative
	create_jacobian_derivative(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// For spring and damping coefficients to remove numerical drift
	// Create constraint vector
	create_constraint_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
	// Create constraint derivative vector
	create_constraint_derivative_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// DEBUGGING
//	std::cout << "Printing out state_vector" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << state_vector[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing out state_vector derivative" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << state_vector_derivative[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing out inverse mass matrix" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << inverse_mass_matrix[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing out external force vector" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << force_ext_vector[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing out jacobian" << std::endl;
//	print_matrix_block_list(jacobian, num_constraints);
//
//	std::cout << "Printing out jacobian time derivative" << std::endl;
//	print_matrix_block_list(jacobian_derivative, num_constraints);

	// Now calculate b.
	float b[num_constraints];
	zero_vector(b, num_constraints);
	calculate_b(b);

//	std::cout << "Printing b vector" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << b[i] << " ";
//	}
//	std::cout << std::endl;

	// Solve for the lambda vector using biconjugate gradient method. try itol 1 to start
	float x[num_constraints];
	zero_vector(x, num_constraints);

	// Try putting the max at 2 * N
	// Pointer to iteration count, and error. Dynamically allocated
	int *iter = (int *)malloc(sizeof(int));
	float *err = (float *)malloc(sizeof(float));

	biconjugate_gradient(this, b, x, 1, TOL, 10 * num_constraints, iter, err);

	// Free dynamically allocated iter and err.
	free(iter);
	free(err);

	// DEBUGGING
//	std::cout << "Total iterations was: " << *iter << std::endl;
//	std::cout << "Error in Ax - B is: " << *err << std::endl;

//	std::cout << "Printing the lambda vector" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << x[i] << " ";
//	}
//	std::cout << std::endl;

	// Now we have x, which is our lambda_vector.
	// Multiply by jacobian transpose to get constraint forces. F = J_T * lambda
	// Allocate space
	constraint_force_vector = (float *)malloc(3 * num_bodies * sizeof(float));
	zero_vector(constraint_force_vector, 3 * num_bodies);

	matrix_blocks_transpose_times_vector(jacobian, num_constraints, x, constraint_force_vector);

//	std::cout << "Printing the constraint force vector" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << constraint_force_vector[i] << " ";
//	}
//	std::cout << std::endl;
	// Find the net force vector
	calculate_net_force_vector();
}

// Create the state_vector
void State_Getter::create_state_vector(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	state_vector = (float *)malloc(3 * num_bodies * sizeof(float));
	// Loop through the rigid bodies
	int count = 0;
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		state_vector[count] = mass_list[i].pos.x;
		// y position
		state_vector[count + 1] = mass_list[i].pos.y;
		// angle
		state_vector[count + 2] = mass_list[i].angle;
		// Increment count
		count += 3;
	}
}

// Creating the matrices and vectors
void State_Getter::create_state_vector_derivative(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	state_vector_derivative = (float *)malloc(3 * num_bodies * sizeof(float));
	// Loop through the rigid bodies
	int count = 0;
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		state_vector_derivative[count] = mass_list[i].linear_vel.x;
		// y position
		state_vector_derivative[count + 1] = mass_list[i].linear_vel.y;
		// angle
		state_vector_derivative[count + 2] = mass_list[i].angular_vel;
		count += 3;
	}
}

void State_Getter::create_inverse_mass_matrix(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	inverse_mass_matrix = (float *)malloc(3 * num_bodies * sizeof(float));
	// Loop through the rigid bodies
	int count = 0;
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		inverse_mass_matrix[count] = 1.0 / mass_list[i].mass;
		// y position
		inverse_mass_matrix[count + 1] = 1.0 / mass_list[i].mass;
		// angle
		inverse_mass_matrix[count + 2] = 1.0 / mass_list[i].moi;
		count += 3;
	}
}

void State_Getter::create_force_ext_vector(Circular_Rigid_Body *mass_list, Spring_2 *spring2s)
{
	// Allocate space for the external force vector
	force_ext_vector = (float *)malloc(3 * num_bodies * sizeof(float));
	// Loop through masses and assign force vectors
	int count = 0;
	for (int i = 0; i < num_bodies; i++)
	{
		// Determine the external force on each mass
		mass_list[i].set_initial_force_ext();
		// Yoink the numbers from each mass: fx, fy, torque
		force_ext_vector[count] = mass_list[i].force_ext.x;
		force_ext_vector[count + 1] = mass_list[i].force_ext.y;
		force_ext_vector[count + 2] = mass_list[i].torque;
		count += 3;
	}
//	std::cout << "Number of spring2s: " << num_spring2s << std::endl;
	// Apply spring force to relevant masses
	for (int i = 0; i < num_spring2s; i++)
	{
		spring2s[i].determine_spring_force(mass_list);
//		std::cout << "Spring force of spring " << i << std::endl;
//		std::cout << "On mass 0: x: " << spring2s[i].sf0.x << " y: " << spring2s[i].sf0.y << std::endl;
//		std::cout << "On mass 1: x: " << spring2s[i].sf1.x << " y: " << spring2s[i].sf1.y << std::endl;
		spring2s[i].apply_spring_force(mass_list, force_ext_vector);
	}
}

void State_Getter::create_constraint_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Vector of scalar constraint values for each constraint
	// Allocate space
	constraint_vector = (float *)malloc(num_constraints * sizeof(float));

	int count = 0;
	// Loop through lists of each type of constraint, these are all the constraint functions
	// Rigid_Bar_1s.
	for (int i = 0; i < num_bar1s; i++)
	{
		// Each rigid_bar_1 applies a constraint to the vector
		constraint_vector[count] = bar1s[i].constraint(mass_list);
		// Increment count
		count++;
	}
	// Rigid_Bar_2s
	for (int i = 0; i < num_bar2s; i++)
	{
		// Each rigid_bar_2 applies a constraint to the vector
		constraint_vector[count] = bar2s[i].constraint(mass_list);
		// Increment count
		count++;
	}

}

void State_Getter::create_constraint_derivative_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Vector of constraint time derivative values for each constraint
	// Allocate space
	constraint_derivative_vector = (float *)malloc(num_constraints * sizeof(float));

	int count = 0;
	// Loop through lists of each type of constraint, these are all the constraint functions
	// Rigid_Bar_1s.
	for (int i = 0; i < num_bar1s; i++)
	{
		// Each rigid_bar_1 applies a constraint to the vector
		constraint_derivative_vector[count] = bar1s[i].constraint_time_derivative(mass_list);
		// Increment count
		count++;
	}
	// Rigid_Bar_2s
	for (int i = 0; i < num_bar2s; i++)
	{
		//std::cout << "count is : " << count << std::endl;
		//std::cout << "i is: " << i << std::endl;
		// Each rigid_bar_2 applies a constraint to the vector
		constraint_derivative_vector[count] = bar2s[i].constraint_time_derivative(mass_list);
		// Increment count
		count++;
	}

}

void State_Getter::create_jacobian(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Jacobian is made up of a list of matrix blocks representing nonzero values. There will be as many matrix blocks as there are constraints
	// Dynamically allocate space for the matrix block list m x 3n
	jacobian = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));

	// Now I need to create my matrix blocks
	// Loop through each list of constraints
	// Count to represent how many matrix_blocks have been created
	int count = 0;
	// Loop through lists of each type of constraint, these are all the constraint functions
	// Rigid_Bar_1s.
	for (int i = 0; i < num_bar1s; i++)
	{
		// Each rigid_bar_1 only applies a partial derivative to the x and y of 1 attached mass. this contributes two nonzero terms to the jacobian
		jacobian[count].rows = 1;
		jacobian[count].cols = 2;
		// Row offset: 1 row per constraint object
		jacobian[count].row = count;
		// Column offset: column of origin of matrix block starts at the mass attached to the rigid rod.
		// 3 derivatives per mass so multiplied by 3.
		jacobian[count].col = 3 * bar1s[i].attached_mass;

		// Allocate the space needed for the matrix in the matrix block.
		// Rigid_Bar_1 constraints contribute 2 non zero terms for the 1 particle affected. partial x and partial y
		jacobian[count].matrix = (float *)malloc(2 * sizeof(float));

		// Fill up the matrix
		//jacobian[count].matrix[0] = bar1s[i].partial_x(mass_list);
		jacobian[count].matrix[0] = bar1s[i].jacobian_entry_x(mass_list);
		jacobian[count].matrix[1] = bar1s[i].jacobian_entry_y(mass_list);

		// Increment count
		count++;
	}
	// Rigid_Bar_2s
	for (int i = 0; i < num_bar2s; i++)
	{
		// Each rigid_bar_2 applies a partial derivative to the x and y of both attached masses. this contributes four nonzero terms to the jacobian
		jacobian[count].rows = 1;
		jacobian[count].cols = 5;
		// Row offset: 1 row per constraint object
		jacobian[count].row = count;
		// Column offset: column of origin of matrix block starts at the first mass the rigid rod is attached to.
		// 3 derivatives per mass so multiplied by 3.
		// Assume minimum attached mass is first in the attached mass list
		jacobian[count].col = 3 * bar2s[i].attached_masses[0];
		// Allocate the space needed for the matrix in the matrix block.
		// Rigid_Bar_2 constraints contribute 4 non zero terms for the 2 particles affected. Need 5 total terms. x1, y1, theta1, x2, y2.
		// partial theta1 is just 0.
		jacobian[count].matrix = (float *)malloc(5 * sizeof(float));

		// Fill up the matrix
		jacobian[count].matrix[0] = bar2s[i].jacobian_entry_x1(mass_list);
		jacobian[count].matrix[1] = bar2s[i].jacobian_entry_y1(mass_list);
		jacobian[count].matrix[2] = 0.0;
		jacobian[count].matrix[3] = bar2s[i].jacobian_entry_x2(mass_list);
		jacobian[count].matrix[4] = bar2s[i].jacobian_entry_y2(mass_list);

		// Increment count
		count++;
	}
}

void State_Getter::create_jacobian_transpose(void)
{
	// Dynamically allocate the space for it
	jacobian_transpose = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));

	// Create it
	matrix_blocks_transpose(jacobian, num_constraints, jacobian_transpose);
}

void State_Getter::create_jacobian_derivative(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Jacobian is made up of a list of matrix blocks representing nonzero values. There will be as many matrix blocks as there are constraints
		// Dynamically allocate space for the matrix block list m x 3n
		jacobian_derivative = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));

		// Now I need to create my matrix blocks
		// Loop through each list of constraints
		// Count to represent how many matrix_blocks have been created
		int count = 0;
		// Loop through lists of each type of constraint, these are all the constraint functions
		// Rigid_Bar_1s.
		for (int i = 0; i < num_bar1s; i++)
		{
			// Each rigid_bar_1 only applies a double partial derivative multiplied by velocity term for both x and y of 1 attached mass.
			// this contributes two nonzero terms to the jacobian
			jacobian_derivative[count].rows = 1;
			jacobian_derivative[count].cols = 2;
			// Row offset: 1 row per constraint object
			jacobian_derivative[count].row = count;
			// Column offset: column of origin of matrix block starts at the mass attached to the rigid rod.
			// 3 derivatives per mass so multiplied by 3.
			jacobian_derivative[count].col = 3 * bar1s[i].attached_mass;

			// Allocate the space needed for the matrix in the matrix block.
			// Rigid_Bar_1 constraints contribute 2 non zero terms for the 1 particle affected. partial x and partial y
			jacobian_derivative[count].matrix = (float *)malloc(2 * sizeof(float));

			// Fill up the matrix
			//jacobian[count].matrix[0] = bar1s[i].partial_x(mass_list);
			jacobian_derivative[count].matrix[0] = bar1s[i].jacobian_derivative_entry_x(mass_list);
			jacobian_derivative[count].matrix[1] = bar1s[i].jacobian_derivative_entry_y(mass_list);

			// Increment count
			count++;
		}
		// Rigid_Bar_2s
		for (int i = 0; i < num_bar2s; i++)
		{
			// Each rigid_bar_2 applies a partial derivative to the x and y of both attached masses. this contributes four nonzero terms to the jacobian
			jacobian_derivative[count].rows = 1;
			jacobian_derivative[count].cols = 5;
			// Row offset: 1 row per constraint object
			jacobian_derivative[count].row = count;
			// Column offset: column of origin of matrix block starts at the first mass the rigid rod is attached to.
			// 3 derivatives per mass so multiplied by 3.
			// Assume minimum attached mass is first in the attached mass list
			jacobian_derivative[count].col = 3 * bar2s[i].attached_masses[0];
			// Allocate the space needed for the matrix in the matrix block.
			// Rigid_Bar_2 constraints contribute 4 non zero terms for the 2 particles affected. Need 5 total terms. x1, y1, theta1, x2, y2.
			// partial theta1 is just 0.
			jacobian_derivative[count].matrix = (float *)malloc(5 * sizeof(float));

			// Fill up the matrix
			jacobian_derivative[count].matrix[0] = bar2s[i].jacobian_derivative_entry_x1(mass_list);
			jacobian_derivative[count].matrix[1] = bar2s[i].jacobian_derivative_entry_y1(mass_list);
			jacobian_derivative[count].matrix[2] = 0.0;
			jacobian_derivative[count].matrix[3] = bar2s[i].jacobian_derivative_entry_x2(mass_list);
			jacobian_derivative[count].matrix[4] = bar2s[i].jacobian_derivative_entry_y2(mass_list);

			// Increment count
			count++;
		}
}

// Calculate the B vector evaluating the right side: -J_dot * q_dot - J * M^-1 * F_ext
void State_Getter::calculate_b(float *b)
{
	// First calculate the left side. J_dot * q_dot
	float lhs[num_constraints];
	zero_vector(lhs, num_constraints);

	matrix_blocks_times_vector(jacobian_derivative, num_constraints, state_vector_derivative, lhs);

//	std::cout << "Printing lhs. jdot * qdot" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << lhs[i] << " ";
//	}
//	std::cout << std::endl;

	// Then multiply it by negative 1.
	for (int i = 0; i < num_constraints; i++)
	{
		lhs[i] = -1 * lhs[i];
	}

	// Now do the right hand side.
	// First inverse mass times F_ext
	float first_res[3 * num_bodies];
	zero_vector(first_res, 3 * num_bodies);
	diagonal_times_vector(inverse_mass_matrix, force_ext_vector, 3 * num_bodies, first_res);

//	std::cout << "Printing rhs. " << std::endl;
//	std::cout << "Printing M-1 * Fext" << std::endl;
//	for (int i = 0; i < 3 * num_bodies; i++)
//	{
//		std::cout << first_res[i] << " ";
//	}
//	std::cout << std::endl;

	// Then J times that
	float rhs[num_constraints];
	zero_vector(rhs, num_constraints);
	matrix_blocks_times_vector(jacobian, num_constraints, first_res, rhs);

	// b = lhs - rhs
	for (int i = 0; i < num_constraints; i++)
	{
		b[i] = lhs[i] - rhs[i];
	}

	// Add the spring and damping coefficient vectors

	// Safety check
	if (ks == 0.0 || kd == 0.0)
	{
		ks = 0.5;
		kd = 0.5;
	}

	float ks_times_constraint[num_constraints];
	float kd_times_constraint_derivative[num_constraints];
	for (int i = 0; i < num_constraints; i++)
	{
		ks_times_constraint[i] = ks * constraint_vector[i];
	}
	for (int i = 0; i < num_constraints; i++)
	{
		kd_times_constraint_derivative[i] = kd * constraint_derivative_vector[i];
	}

//	std::cout << "Printing constraint and constraint derivative" << std::endl;
//	std::cout << "Constraint vector" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << constraint_vector[i] << " ";
//	}
//	std::cout << std::endl;
//	std::cout << "Constraint vector derivative" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << constraint_derivative_vector[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing ksc and kdcd" << std::endl;
//	std::cout << "ksc" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << ks_times_constraint[i] << " ";
//	}
//	std::cout << std::endl;
//	std::cout << "kdcd" << std::endl;
//	for (int i = 0; i < num_constraints; i++)
//	{
//		std::cout << kd_times_constraint_derivative[i] << " ";
//	}
//	std::cout << std::endl;

	// Subtract the sum from b
	for (int i = 0; i < num_constraints; i++)
	{
		b[i] = b[i] - (ks_times_constraint[i] + kd_times_constraint_derivative[i]);
	}
}

void State_Getter::calculate_net_force_vector(void)
{
	// Allocate space
	net_force_vector = (float *)malloc(3 * num_bodies * sizeof(float));
	zero_vector(net_force_vector, 3 * num_bodies);

	for (int i = 0; i < 3 * num_bodies; i++)
	{
		net_force_vector[i] = force_ext_vector[i] + constraint_force_vector[i];
	}
}

void State_Getter::free_all(void)
{
	// Free the vectors
	free(state_vector);
	free(state_vector_derivative);
	free(inverse_mass_matrix);
	free(force_ext_vector);
	free(constraint_vector);
	free(constraint_derivative_vector);
	free(constraint_force_vector);
	free(net_force_vector);

	// Free the matrix block structs
	for (int i = 0; i < num_constraints; i++)
	{
		// Free the matrix arrays in each block
		free(jacobian[i].matrix);
	}
	free(jacobian);

	for (int i = 0; i < num_constraints; i++)
	{
		// Free the matrix arrays in each block
		free(jacobian_derivative[i].matrix);
	}
	free(jacobian_derivative);
}
