#ifndef SPRING_H
#define SPRING_H

// Forward declarations
class Vector;
class Circular_Rigid_Body;

// Spring_2 connects and applies forces to two masses
class Spring_2
{
	// Attaches two masses with a straight line. Has equilibrium length.
	// Emits forces inwards and outwards on the two masses in the direction the spring is oriented depending on length change
	public:

	// Equilibrium length
	double eq_length;
	// Spring constant
	double spring_const;
	// Point 1 and 2. Where spring attaches
	int attached_masses[2];

	// Spring forces on each mass
	Vector sf0;
	Vector sf1;

	// Color
	unsigned int color;

	void determine_spring_force(Circular_Rigid_Body *mass_list);

	void apply_spring_force(Circular_Rigid_Body *mass_list, double *force_ext_vector);

	void draw_spring(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width);
};

#endif
