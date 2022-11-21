#include <iostream>

#include "get_state.h"
#include "rigid_bodies.h"
#include "vectors.h"

// Euler's Method
void euler_method(Circular_Rigid_Body *mass_list, int num_bodies, double *net_force_vector, double dt)
{
	// Loop through masses
	for (int i = 0; i < num_bodies; i++)
	{
		// Torque
		mass_list[i].torque = net_force_vector[3 * i + 2] / mass_list[i].moi;

		// Update the linear and angular velocities
		mass_list[i].linear_vel.x += (net_force_vector[3 * i] / mass_list[i].mass) * dt;
		mass_list[i].linear_vel.y += (net_force_vector[3 * i + 1] / mass_list[i].mass) * dt;
		mass_list[i].angular_vel += (net_force_vector[3 * i + 2] / mass_list[i].moi) * dt;

		// Update the position and orientation
		mass_list[i].pos.x += mass_list[i].linear_vel.x * dt;
		mass_list[i].pos.y += mass_list[i].linear_vel.y * dt;
		mass_list[i].angle += mass_list[i].angular_vel * dt;
	}
}


