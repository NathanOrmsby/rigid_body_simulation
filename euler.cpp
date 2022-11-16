#include <iostream>
#include "rigid_bodies.h"
#include "get_state.h"
#include "vectors.h"

// Start with Euler's Method

void euler_method(Circular_Rigid_Body *mass_list, int num_bodies, float *net_force_vector, float dt)
{
	// Get the current state
//	std::cout << "Made it to euler method" << std::endl;
	//std::cout << "Time step is" << dt << std::endl;
	for (int i = 0; i < num_bodies; i++)
	{
//		std::cout << "i is: " << i << std::endl;
		// Update the acceleration for each mass
		// a = Fnet / m, ac = torque / moi
		// ax
		mass_list[i].linear_accel.x = net_force_vector[3 * i] / mass_list[i].mass;
//		std::cout << "made it here 1" << std::endl;
		mass_list[i].linear_accel.y = net_force_vector[3 * i + 1] / mass_list[i].mass;

		// Torque
		mass_list[i].torque = net_force_vector[3 * i + 2] / mass_list[i].moi;

		// Update the linear and angular velocities
		// Debugging
		mass_list[i].linear_vel.x += mass_list[i].linear_accel.x * dt;
		mass_list[i].linear_vel.y += mass_list[i].linear_accel.y * dt;
		mass_list[i].angular_vel += mass_list[i].angular_accel * dt;

		// Update the position and orientation
		mass_list[i].pos.x += mass_list[i].linear_vel.x * dt;
		mass_list[i].pos.y += mass_list[i].linear_vel.y * dt;
		mass_list[i].angle += mass_list[i].angular_vel * dt;
	}
}


