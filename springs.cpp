#include "renderer.h"
#include "vectors.h"
#include "rigid_bodies.h"
#include "springs.h"

// Spring 2 connects two masses.
void Spring_2::determine_spring_force(Circular_Rigid_Body *mass_list)
{
	Vector point0_to_1 = {mass_list[attached_masses[1]].pos.x - mass_list[attached_masses[0]].pos.x, mass_list[attached_masses[1]].pos.y - mass_list[attached_masses[0]].pos.y};
	// Find distance of length or contraction
	double dist = vector_magnitude(point0_to_1);

	// Find center point
	Vector center = {mass_list[attached_masses[0]].pos.x + point0_to_1.x / 2, mass_list[attached_masses[0]].pos.y + point0_to_1.y / 2};
	// Calculate length change
	double length_change = dist - eq_length;

	// Determine magnitude
	double sf_mag = -1 * spring_const * length_change;

	// Determine directions to each point from center
	Vector center_to_point0 = {mass_list[attached_masses[0]].pos.x - center.x, mass_list[attached_masses[0]].pos.y - center.y};
	Vector center_to_point1 = {mass_list[attached_masses[1]].pos.x - center.x, mass_list[attached_masses[1]].pos.y - center.y};

	Vector dir0 = vector_direction(center_to_point0);
	Vector dir1 = vector_direction(center_to_point1);

	// Calculate spring forces on each point
	sf0 = {sf_mag * dir0.x, sf_mag * dir0.y};
	sf1 = {sf_mag * dir1.x, sf_mag * dir1.y};
}

void Spring_2::apply_spring_force(Circular_Rigid_Body *mass_list, double *force_ext_vector)
{
	// Applies spring force onto masses it is attached to
	// Mass 0
	force_ext_vector[3 * attached_masses[0]] += sf0.x;
	force_ext_vector[3 * attached_masses[0] + 1] += sf0.y;
	// Mass 1
	force_ext_vector[3 * attached_masses[1]] += sf1.x;
	force_ext_vector[3 * attached_masses[1] + 1] += sf1.y;
}
void Spring_2::draw_spring(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width)
{
	// Draw the spring
	draw_line(to_screen(mass_list[attached_masses[0]].pos.x), to_screen(mass_list[attached_masses[0]].pos.y), to_screen(mass_list[attached_masses[1]].pos.x), to_screen(mass_list[attached_masses[1]].pos.y), color, first_pixel, pixel_buffer_width);
}

