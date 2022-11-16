/*
 * rk4.cpp
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */

#include "euler.h"
#include "rigid_bodies.h"
#include "get_state.h"
#include "vectors.h"
#include "rk4.h"
// Runge Kutta 4th Order to solve differential equations

void rk4(State_Getter *state, Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, Spring_2 *spring2s, float dt)
{

	// k1 step
	// Get net force at t0
	// Calculate the net force at the starting position
	state->get_current_state(mass_list, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Copy net force list into k1
	float k1[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k1[i] = state->net_force_vector[i];
	}
	state->free_all();

//	std::cout << "Printing mass properties BEFORE k1 calculation" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

	// k2 step
	// Make new mass_list for k2 state
	Circular_Rigid_Body mass_list_k2[state->num_bodies];
	copy_mass_list(mass_list, mass_list_k2, state->num_bodies);

	// Update position for k2_masses using k1 as net force with dt / 2 time step
	float k2_step = (float)(dt / 2.0);
	euler_method(mass_list_k2, state->num_bodies, k1, k2_step);

//	std::cout << "Printing mass properties AFTER updating mass_list_k2" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

	// New state_getter for k2
	State_Getter k2_state;
	k2_state.num_bodies = state->num_bodies;
	k2_state.num_constraints = num_bar1s + num_bar2s;
	k2_state.ks = 0.5;
	k2_state.kd = 0.5;
	k2_state.num_spring2s = state->num_spring2s;

	// Get state for k2
	k2_state.get_current_state(mass_list_k2, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k2 net force
	float k2[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k2[i] = k2_state.net_force_vector[i];
	}
	// Free all matrices in k2_state
	k2_state.free_all();

	// k3 step
	// Make new mass_list for k3 state
	Circular_Rigid_Body mass_list_k3[state->num_bodies];
	copy_mass_list(mass_list, mass_list_k3, state->num_bodies);

	// Update position for k3 using k2 as net force
	float k3_step = (float)(dt / 2.0);
	euler_method(mass_list_k3, state->num_bodies, k2, k3_step);

//	std::cout << "Printing mass properties AFTER updating mass_list_k3" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

	// New state_getter for k3
	State_Getter k3_state;
	k3_state.num_bodies = state->num_bodies;
	k3_state.num_constraints = num_bar1s + num_bar2s;
	k3_state.num_spring2s = state->num_spring2s;
	k3_state.ks = 0.5;
	k3_state.kd = 0.5;

	// Get state for k3
	k3_state.get_current_state(mass_list_k3, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k3 net force
	float k3[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k3[i] = k3_state.net_force_vector[i];
	}
	// Free all matrices in k3_state
	k3_state.free_all();

	// k4 step
	// Make new mass_list for k4 state
	Circular_Rigid_Body mass_list_k4[state->num_bodies];
	copy_mass_list(mass_list, mass_list_k4, state->num_bodies);

	// Update position for k4_masses using k3 as net force
	euler_method(mass_list_k4, state->num_bodies, k3, dt);

//	std::cout << "Printing mass properties AFTER updating mass_list_k4" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

	// New state_getter for k4
	State_Getter k4_state;
	k4_state.num_bodies = state->num_bodies;
	k4_state.num_constraints = num_bar1s + num_bar2s;
	k4_state.num_spring2s = state->num_spring2s;
	k4_state.ks = 0.5;
	k4_state.kd = 0.5;

	// Get state for k4
	k4_state.get_current_state(mass_list_k4, bar1s, num_bar1s, bar2s, num_bar2s, spring2s);

	// Fill up k4 net force
	float k4[3 * state->num_bodies];
	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		k4[i] = k4_state.net_force_vector[i];
	}
	// Free all matrices in k4_state
	k4_state.free_all();

	// Calculate the weighted rk4 net force (k1 + 2*k2 + 2*k3 + k4)
	float rk4_force[3 * state->num_bodies];

	for (int i = 0; i < 3 * state->num_bodies; i++)
	{
		rk4_force[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
	}

	// Apply weighted force using Euler
	euler_method(mass_list, state->num_bodies, rk4_force, dt);

//	std::cout << "Printing mass properties AFTER applying weighted force" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

//	std::cout << "Printing mass properties of rk4" << std::endl;
//	for (int i = 0; i < state->num_bodies; i++)
//	{
//		std::cout << "Mass " << i << std::endl;
//		std::cout << "Position is: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
//		std::cout << "Velocity is: vx: " << mass_list[i].linear_vel.x << " vy: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Acceleration is: ax: " << mass_list[i].linear_accel.x << " ay: " << mass_list[i].linear_accel.y << std::endl;
//	}

//	std::cout << "Printing all calculated net force vectors" << std::endl;
//	std::cout << "Printing k1:" << std::endl;
//	for (int i = 0; i < 3 * state->num_bodies; i++)
//	{
//		std::cout << k1[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing k2:" << std::endl;
//	for (int i = 0; i < 3 * state->num_bodies; i++)
//	{
//		std::cout << k2[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing k3:" << std::endl;
//	for (int i = 0; i < 3 * state->num_bodies; i++)
//	{
//		std::cout << k3[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing k4:" << std::endl;
//	for (int i = 0; i < 3 * state->num_bodies; i++)
//	{
//		std::cout << k4[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing weighted:" << std::endl;
//	for (int i = 0; i < 3 * state->num_bodies; i++)
//	{
//		std::cout << rk4_force[i] / 6 << " ";
//	}
//	std::cout << std::endl;

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
		new_list[i].linear_accel.x = old_list[i].linear_accel.x;
		new_list[i].linear_accel.y = old_list[i].linear_accel.y;
		new_list[i].angular_accel = old_list[i].angular_accel;
		new_list[i].force_ext.x = old_list[i].force_ext.x;
		new_list[i].force_ext.y = old_list[i].force_ext.y;
		new_list[i].torque = old_list[i].torque;
		new_list[i].mass = old_list[i].mass;
		new_list[i].radius = old_list[i].radius;
		new_list[i].moi = old_list[i].moi;
		new_list[i].color = old_list[i].color;

		//		float x = old_list[i].pos.x;
		//		float y = old_list[i].pos.y;
		//		float angle = old_list[i].angle;
		//		float vx = old_list[i].linear_vel.x;
		//		float vy = old_list[i].linear_vel.y;
		//		float angular_vel = old_list[i].angular_vel;
		//		float ax = old_list[i].linear_accel.x;
		//		float ay = old_list[i].linear_accel.y;
		//		float angular_accel = old_list[i].angular_accel;
		//		float force_ext_x = old_list[i].force_ext.x;
		//		float force_ext_y = old_list[i].force_ext.y;
		//		float torque = old_list[i].torque;
		//		float mass = old_list[i].mass;
		//		float radius = old_list[i].radius;
		//		float moi = old_list[i].moi;
		//		unsigned int color = old_list[i].color;
		//		new_list[i].pos.x = x;
		//		new_list[i].pos.y = y;
		//		new_list[i].angle = angle;
		//		new_list[i].linear_vel.x = vx;
		//		new_list[i].linear_vel.y = vy;
		//		new_list[i].angular_vel = angular_vel;
		//		new_list[i].linear_accel.x = ax;
		//		new_list[i].linear_accel.y = ay;
		//		new_list[i].angular_accel = angular_accel;
		//		new_list[i].force_ext.x = force_ext_x;
		//		new_list[i].force_ext.y = force_ext_y;
		//		new_list[i].torque = torque;
		//		new_list[i].mass = mass;
		//		new_list[i].radius = radius;
		//		new_list[i].moi = moi;
		//		new_list[i].color = color;
	}

}
