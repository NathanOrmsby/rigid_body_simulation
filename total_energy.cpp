// Calculates the kinetic plus potential energies of masses in the scene

#include "rigid_bodies.h"
#include "total_energy.h"
double total_energy(Circular_Rigid_Body *mass_list, int num_bodies)
{
	double tot_kinetic = 0.0;
	double tot_gp = 0.0;
	// Calculate the kinetic and potential energies of each mass
	for (int i = 0; i < num_bodies; i++)
	{
		double vx = mass_list[i].linear_vel.x;
		double vy = mass_list[i].linear_vel.y;
		tot_kinetic += 0.5 * mass_list[i].mass * ((vx * vx) + (vy * vy));
		tot_gp += mass_list[i].mass * 9.8 * mass_list[i].pos.y;
	}
	return tot_kinetic + tot_gp;
}


