/*
 * rk4.cpp
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */

#include "euler.h"
#include "rigid_bodies.h"
#include "vectors.h"
#include "rk4.h"

#include "get_state.h"
#include "utils.h"
#include "my_timer.h"
// Runge Kutta 4th Order to solve differential equations

void rk4(State_Getter *state, Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, Spring_2 *spring2s, double dt, int loop_count)
{

	// k1 step
	// Get net force at t0
	// Calculate the net force at the starting position
	state->get_current_state(mass_list, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Copy net force list into k1
	double k1_force[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k1_force[i] = state->net_force_vector[i];
	}

	// DEBUGGING
//	if (loop_count == 1)
//	{
//		std::cout << "DEBUGGING" << std::endl;
//		std::cout << "\nPrinting mass properties" << std::endl;
//		for (int i = 0; i < state->num_bodies; i++)
//		{
//			std::cout << "Mass " << i << std::endl;
//			std::cout << "Pos: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//			std::cout << "Vel: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		}
//		std::cout << "\n\nPrinting jacobian" << std::endl;
//		print_matrix_block_list(state->jacobian, state->num_constraints);
//		std::cout << "\n\nPrinting jacobian derivative" << std::endl;
//		print_matrix_block_list(state->jacobian_derivative, state->num_constraints);
//
//	}

	state->free_all();
	// K1 velocity
	double k1_vel[3 * state->num_bodies];
	for (int i = 0; i < state->num_bodies; i++)
	{
		k1_vel[3 * i] = mass_list[i].linear_vel.x;
		k1_vel[3 * i + 1] = mass_list[i].linear_vel.y;
		k1_vel[3 * i + 2] = mass_list[i].angular_vel;
	}

//	if (loop_count == 1)
//	{
//		std::cout << "Printing k1 force" << std::endl;
//		for (int i = 0; i < 3 * state->num_bodies; i++)
//		{
//			std::cout << k1_force[i] << " ";
//		}
//		std::cout << "\n\nPrinting k1 velocity" << std::endl;
//		for (int i = 0; i < 3 * state->num_bodies; i++)
//		{
//			std::cout << k1_vel[i] << " ";
//		}
//		std::cout << std::endl;
//	}

	// k2 step
	// Make new mass_list for k2 state.
	Circular_Rigid_Body mass_list_k2[state->num_bodies];

	// Copy properties of mass_list

	copy_mass_list(mass_list, mass_list_k2, state->num_bodies);

	// Push k2 state forward dt / 2 using k1 as net force
	for (int i = 0; i < state->num_bodies; i++)
	{
		// Update velocities: THIS USES k1 / 2
		mass_list_k2[i].linear_vel.x += (k1_force[3 * i] / mass_list_k2[i].mass) * 0.5 * dt;
		mass_list_k2[i].linear_vel.y += (k1_force[3 * i + 1] / mass_list_k2[i].mass * 0.5) * dt;
		mass_list_k2[i].angular_vel += (k1_force[3 * i + 2] / mass_list_k2[i].moi * 0.5) * dt;

		// Update positions: Using k1 velocity. Initial position velocity
		mass_list_k2[i].pos.x += k1_vel[3 * i] * 0.5 * dt;
		mass_list_k2[i].pos.y += k1_vel[3 * i + 1] * 0.5 * dt;
		mass_list_k2[i].angle += k1_vel[3 * i + 2] * 0.5 * dt;
	}

	// New state_getter for k2
	State_Getter k2_state;
	k2_state.num_bodies = state->num_bodies;
	k2_state.num_constraints = num_bar1s + num_bar2s;
	k2_state.num_spring2s = state->num_spring2s;

	// Get state for k2
	k2_state.get_current_state(mass_list_k2, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k2 net force
	double k2_force[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k2_force[i] = k2_state.net_force_vector[i];
	}

	// Get k2 velocity vector
	double k2_vel[3 * state->num_bodies];
	for (int i = 0; i < state->num_bodies; i++)
	{

		k2_vel[3 * i] = mass_list_k2[i].linear_vel.x;
		k2_vel[3 * i + 1] = mass_list_k2[i].linear_vel.y;
		k2_vel[3 * i + 2] = mass_list_k2[i].angular_vel;
	}

	// Free all matrices in k2_state
	k2_state.free_all();

	// k3 step
	// Make new mass_list for k3 state
	Circular_Rigid_Body mass_list_k3[state->num_bodies];
	copy_mass_list(mass_list, mass_list_k3, state->num_bodies);

	// Push the mass_list_k3 state forward dt / 2 using k2 as net force
	for (int i = 0; i < state->num_bodies; i++)
	{
		// Update velocities
		mass_list_k3[i].linear_vel.x += (k2_force[3 * i] / mass_list_k3[i].mass) * 0.5 * dt;
		mass_list_k3[i].linear_vel.y += (k2_force[3 * i + 1] / mass_list_k3[i].mass) * 0.5 * dt;
		mass_list_k3[i].angular_vel += (k2_force[3 * i + 2] / mass_list_k3[i].moi) * 0.5 * dt;

		// Update positions
		mass_list_k3[i].pos.x += k2_vel[3 * i] * 0.5 * dt;
		mass_list_k3[i].pos.y += k2_vel[3 * i + 1] * 0.5 * dt;
		mass_list_k3[i].angle += k2_vel[3 * i + 2] * 0.5 * dt;
	}

	// New state_getter for k3
	State_Getter k3_state;
	k3_state.num_bodies = state->num_bodies;
	k3_state.num_constraints = num_bar1s + num_bar2s;
	k3_state.num_spring2s = state->num_spring2s;

	// Get state for k3
	k3_state.get_current_state(mass_list_k3, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k3 net force
	double k3_force[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k3_force[i] = k3_state.net_force_vector[i];
	}

	// Get k3 velocity vector
	double k3_vel[3 * state->num_bodies];
	for (int i = 0; i < state->num_bodies; i++)
	{

		k3_vel[3 * i] = mass_list_k3[i].linear_vel.x;
		k3_vel[3 * i + 1] = mass_list_k3[i].linear_vel.y;
		k3_vel[3 * i + 2] = mass_list_k3[i].angular_vel;
	}

	// Free all matrices in k3_state
	k3_state.free_all();

	// k4 step
	// Make new mass_list for k4 state
	Circular_Rigid_Body mass_list_k4[state->num_bodies];
	copy_mass_list(mass_list, mass_list_k4, state->num_bodies);

	// Push the k4 state forward dt using k3 as net force and k3 vel
	for (int i = 0; i < state->num_bodies; i++)
	{
		// Update velocities
		mass_list_k4[i].linear_vel.x += (k3_force[3 * i] / mass_list_k4[i].mass) * dt;
		mass_list_k4[i].linear_vel.y += (k3_force[3 * i + 1] / mass_list_k4[i].mass) * dt;
		mass_list_k4[i].angular_vel += (k3_force[3 * i + 2] / mass_list_k4[i].moi) * dt;

		// Update positions
		mass_list_k4[i].pos.x += k3_vel[3 * i] * dt;
		mass_list_k4[i].pos.y += k3_vel[3 * i + 1] * dt;
		mass_list_k4[i].angle += k3_vel[3 * i + 2] * dt;
	}

	// New state_getter for k4
	State_Getter k4_state;
	k4_state.num_bodies = state->num_bodies;
	k4_state.num_constraints = num_bar1s + num_bar2s;
	k4_state.num_spring2s = state->num_spring2s;

	// Get state for k4
	k4_state.get_current_state(mass_list_k4, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k4 net force
	double k4_force[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k4_force[i] = k4_state.net_force_vector[i];
	}

	// Get k4 velocity vector
	double k4_vel[3 * state->num_bodies];
	for (int i = 0; i < state->num_bodies; i++)
	{

		k4_vel[3 * i] = mass_list_k4[i].linear_vel.x;
		k4_vel[3 * i + 1] = mass_list_k4[i].linear_vel.y;
		k4_vel[3 * i + 2] = mass_list_k4[i].angular_vel;
	}

	// Free all matrices in k4_state
	k4_state.free_all();

	// Calculate the weighted rk4 net force (k1_force + 2*k2_force + 2*k3_force + k4_force)
	double rk4_force[3 * state->num_bodies];
	double rk4_vel[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		rk4_force[i] = k1_force[i] + 2 * k2_force[i] + 2 * k3_force[i] + k4_force[i];
	}
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

	// Calculate the weighted rk4 velocity (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

	// Push the initial state forward dt using the weighted average rk4 force as net force
	for (int i = 0; i < state->num_bodies; i++)
	{
		// Update velocities
		mass_list[i].linear_vel.x += ((rk4_force[3 * i] / mass_list[i].mass) / 6.0) * dt;
		mass_list[i].linear_vel.y += ((rk4_force[3 * i + 1] / mass_list[i].mass) / 6.0) * dt;
		mass_list[i].angular_vel += ((rk4_force[3 * i + 2] / mass_list[i].moi) / 6.0) * dt;

		// Update positions
		mass_list[i].pos.x += (rk4_vel[3 * i] / 6.0)  * dt;
		mass_list[i].pos.y += (rk4_vel[3 * i + 1] / 6.0)  * dt;
		mass_list[i].angle += (rk4_vel[3 * i + 2] / 6.0)  * dt;
	}
}

// Fills a mass list for another state
void copy_mass_list(Circular_Rigid_Body *old_list, Circular_Rigid_Body *new_list, int num_bodies)
{
	for (int i = 0; i < num_bodies; i++)
	{
		new_list[i].pos.x = old_list[i].pos.x;
		new_list[i].pos.y = old_list[i].pos.y;
		new_list[i].angle = old_list[i].angle;
		new_list[i].linear_vel.x = old_list[i].linear_vel.x;
		new_list[i].linear_vel.y = old_list[i].linear_vel.y;
		new_list[i].angular_vel = old_list[i].angular_vel;
		new_list[i].force_ext.x = old_list[i].force_ext.x;
		new_list[i].force_ext.y = old_list[i].force_ext.y;
		new_list[i].torque = old_list[i].torque;
		new_list[i].mass = old_list[i].mass;
		new_list[i].radius = old_list[i].radius;
		new_list[i].moi = old_list[i].moi;
		new_list[i].color = old_list[i].color;
	}
}
