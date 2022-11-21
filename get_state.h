/*
 * get_state.h
 *
 *  Created on: Nov 12, 2022
 *      Author: norms
 */

#ifndef GET_STATE_H_BKUP_
#define GET_STATE_H_BKUP_

#include <iostream>
#include "matrix_stuff.h"
#include "constraint_bodies.h"
#include "rigid_bodies.h"
#include "vectors.h"
#include "utils.h"

class Spring_2;

// Functions involved with calculating the state of the system
class State_Getter
{
	public:

	// Number of stuff
	int num_bodies;
	int num_constraints;
	int num_spring2s;

	// Different vectors and matrices
	double *state_vector;
	double *state_vector_derivative;
	double *inverse_mass_matrix;
	double *force_ext_vector;
	double *constraint_vector;
	double *constraint_derivative_vector;
	Matrix_Block *jacobian;
	Matrix_Block *jacobian_transpose;
	Matrix_Block *jacobian_derivative;
	double *constraint_force_vector;
	double *net_force_vector;

	// Spring and damping coefficents
	double ks = 0.0;
	double kd = 0.0;

	void get_current_state(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, Spring_2 *spring2s);

	// Create the matrices and vectors
	void create_state_vector(Circular_Rigid_Body *mass_list);

	void create_state_vector_derivative(Circular_Rigid_Body *mass_list);

	void create_inverse_mass_matrix(Circular_Rigid_Body *mass_list);

	void create_force_ext_vector(Circular_Rigid_Body *mass_list, Spring_2 *spring2s);

	void create_constraint_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s);

	void create_constraint_derivative_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s);

	void create_jacobian(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s);

	void create_jacobian_transpose(void);

	void create_jacobian_derivative(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s);

	void calculate_b(double *b);

	void calculate_net_force_vector(void);

	void free_all(void);
};


#endif /* GET_STATE_H_BKUP_ */
