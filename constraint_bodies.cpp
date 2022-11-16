#include <cmath>
#include <iostream>
#include "rigid_bodies.h"
#include "constraint_bodies.h"
#include "vectors.h"
#include "renderer.h"

// The rigid bar classes. NEVER CHANGES LENGTH BETWEEN POINTS
// Rigid bar 1 methods
// This rigid bar attaches at a pivot at point0, and a mass at point1

float Rigid_Bar_1::constraint(Circular_Rigid_Body *mass_list)
{
	// Given a state vector q, return a scalar value of this constraint function.
	// For single attachments, where one end is fixed. Connects a maximum of two masses
	// STATE VECTOR HAS 3 elements per mass, x, y, and theta

	// First constraint equation
	// Set the pivot to the origin
	// Shift x1 and y1 so they are relative to the pivot
	float x1 = mass_list[attached_mass].pos.x;
	float y1 = mass_list[attached_mass].pos.y;

	// Calculate the constraint
	float length = vector_magnitude({initial_point.x - pivot.x, initial_point.y - pivot.y});
	return ((x1 - pivot.x) * (x1 - pivot.x)) + ((y1 - pivot.y) * (y1 - pivot.y)) - (length * length);
}

float Rigid_Bar_1::constraint_time_derivative(Circular_Rigid_Body *mass_list)
{
	// Returns a value for the constraint time derivative
	float x1 = mass_list[attached_mass].pos.x;
	float y1 = mass_list[attached_mass].pos.y;
	float x1_dot = mass_list[attached_mass].linear_vel.x;
	float y1_dot = mass_list[attached_mass].linear_vel.y;
	return 2.0 * x1_dot * (x1 - pivot.x) + 2.0 * y1_dot * (y1 - pivot.y);

}

// Partial of constraint function with respect to x
float Rigid_Bar_1::jacobian_entry_x(Circular_Rigid_Body* mass_list)
{
	return 2.0 * (mass_list[attached_mass].pos.x - pivot.x);
}

// Partial of constraint function with respect to y
float Rigid_Bar_1::jacobian_entry_y(Circular_Rigid_Body *mass_list)
{
	return 2.0 * (mass_list[attached_mass].pos.y - pivot.y);
}
// partial of constraint time derivative wrt x and y
float Rigid_Bar_1::jacobian_derivative_entry_x(Circular_Rigid_Body *mass_list)
{
	return 2.0 * mass_list[attached_mass].linear_vel.x;
}

float Rigid_Bar_1::jacobian_derivative_entry_y(Circular_Rigid_Body *mass_list)
{
	return 2.0 * mass_list[attached_mass].linear_vel.y;
}

// Assign initial point given attached mass
void Rigid_Bar_1::determine_initial_point(Circular_Rigid_Body *mass_list)
{
	initial_point = mass_list[attached_mass].pos;
}

// Draw to screen: Either rod or line
void Rigid_Bar_1::draw_to_screen(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width)
{
	if (rod == 0)
	{
		draw_line(to_screen(pivot.x), to_screen(pivot.y), to_screen(mass_list[attached_mass].pos.x), to_screen(mass_list[attached_mass].pos.y), color, first_pixel, pixel_buffer_width);
	}
	else
	{
		draw_rod(to_screen(pivot.x), to_screen(pivot.y), to_screen(mass_list[attached_mass].pos.x), to_screen(mass_list[attached_mass].pos.y), thickness, color, first_pixel, pixel_buffer_width);
	}
}

// Rigid bar 2 functions
// Attaches to two moving mass objects

float Rigid_Bar_2::constraint(Circular_Rigid_Body *mass_list)
{

		float x1 = mass_list[attached_masses[0]].pos.x;
		float y1 = mass_list[attached_masses[0]].pos.y;
		float x2 = mass_list[attached_masses[1]].pos.x;
		float y2 = mass_list[attached_masses[1]].pos.y;

		// Calculate the constraint
		float length = vector_magnitude({initial_point1.x - initial_point0.x, initial_point1.y - initial_point0.y});
		return ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) - (length * length);
}

float Rigid_Bar_2::constraint_time_derivative(Circular_Rigid_Body *mass_list)
{
	float x1 = mass_list[attached_masses[0]].pos.x;
	float y1 = mass_list[attached_masses[0]].pos.y;
	float x2 = mass_list[attached_masses[1]].pos.x;
	float y2 = mass_list[attached_masses[1]].pos.y;

	float x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	float y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	float x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	float y2_dot = mass_list[attached_masses[1]].linear_vel.y;

	return 2.0 * x1_dot * (x1 - x2) + 2.0 * y1_dot * (y1 - y2) + 2.0 * x2_dot * (x2 - x1) + 2.0 * y2_dot * (y2 - y1);
}

float Rigid_Bar_2::jacobian_entry_x1(Circular_Rigid_Body *mass_list)
{

	float x1 = mass_list[attached_masses[0]].pos.x;
	float x2 = mass_list[attached_masses[1]].pos.x;
	return -2.0 * (x2 - x1);

}

float Rigid_Bar_2::jacobian_entry_y1(Circular_Rigid_Body *mass_list)
{

	float y1 = mass_list[attached_masses[0]].pos.y;
	float y2 = mass_list[attached_masses[1]].pos.y;
	return -2.0 * (y2 - y1);
}

float Rigid_Bar_2::jacobian_entry_x2(Circular_Rigid_Body *mass_list)
{

	float x1 = mass_list[attached_masses[0]].pos.x;
	float x2 = mass_list[attached_masses[1]].pos.x;
	return 2.0 * (x2 - x1);
}

float Rigid_Bar_2::jacobian_entry_y2(Circular_Rigid_Body *mass_list)
{

	float y1 = mass_list[attached_masses[0]].pos.y;
	float y2 = mass_list[attached_masses[1]].pos.y;
	return 2.0 * (y2 - y1);
}

// partial of constraint time derivative wrt x1, y2, x2, y2
float Rigid_Bar_2::jacobian_derivative_entry_x1(Circular_Rigid_Body *mass_list)
{
	float x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	float x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	return 2.0 * (x1_dot - x2_dot);
}

float Rigid_Bar_2::jacobian_derivative_entry_y1(Circular_Rigid_Body *mass_list)
{
	float y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	float y2_dot = mass_list[attached_masses[1]].linear_vel.y;
	return 2.0 * (y1_dot - y2_dot);
}

float Rigid_Bar_2::jacobian_derivative_entry_x2(Circular_Rigid_Body *mass_list)
{
	float x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	float x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	return 2.0 * (x2_dot - x1_dot);
}

float Rigid_Bar_2::jacobian_derivative_entry_y2(Circular_Rigid_Body *mass_list)
{
	float y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	float y2_dot = mass_list[attached_masses[1]].linear_vel.y;
	return 2.0 * (y2_dot - y1_dot);
}

void Rigid_Bar_2::determine_initial_points(Circular_Rigid_Body *mass_list)
{
	initial_point0 = mass_list[attached_masses[0]].pos;
	initial_point1 = mass_list[attached_masses[1]].pos;
}

void Rigid_Bar_2::draw_to_screen(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width)
{
	if (rod == 0)
	{
		draw_line(to_screen(mass_list[attached_masses[0]].pos.x), to_screen(mass_list[attached_masses[0]].pos.y), to_screen(mass_list[attached_masses[1]].pos.x), to_screen(mass_list[attached_masses[1]].pos.y), color, first_pixel, pixel_buffer_width);
	}
	else
	{
		draw_rod(to_screen(mass_list[attached_masses[0]].pos.x), to_screen(mass_list[attached_masses[0]].pos.y), to_screen(mass_list[attached_masses[1]].pos.x), to_screen(mass_list[attached_masses[1]].pos.y), thickness, color, first_pixel, pixel_buffer_width);
	}
}

