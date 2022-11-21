/*
 * rigid_bar.h
 *
 *  Created on: Nov 7, 2022
 *      Author: norms
 */
#include "vectors.h"
#include "rigid_bodies.h"
#ifndef RIGID_BAR_H_
#define RIGID_BAR_H_

// The rigid bar classes. NEVER CHANGES LENGTH BETWEEN POINTS

// This rigid bar attaches at a pivot at point0, and a mass at point1
class Rigid_Bar_1
{
	public:

	Vector pivot;
	Vector initial_point;
	int attached_mass;
	unsigned int color;
	// Draw a rod or not
	int rod = 0;
	int thickness;

	// Constraint equation of rigid rod
	double constraint(Circular_Rigid_Body *mass_list);
	double constraint_time_derivative(Circular_Rigid_Body *mass_list);
	// Constraint derivatives
	// Partial derivatives wrt x and y
	double jacobian_entry_x(Circular_Rigid_Body* mass_list);
	double jacobian_entry_y(Circular_Rigid_Body *mass_list);

	// Partial of constraint time derivative wrt x and y.
	double jacobian_derivative_entry_x(Circular_Rigid_Body *mass_list);
	double jacobian_derivative_entry_y(Circular_Rigid_Body *mass_list);

	void determine_initial_point(Circular_Rigid_Body *mass_list);

	// Draw to screen
	void draw_to_screen(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width);


};

// Attaches to two moving mass objects
class Rigid_Bar_2
{
	public:

	Vector initial_point0;
	Vector initial_point1;
	// Tells which masses the rigid rod is connected to
	int attached_masses[2];

	// Color
	unsigned int color;

	// Draw a rod or not
	int rod = 0;
	int thickness;

	// Constraint equation of rigid rod connecting two moving masses
	double constraint(Circular_Rigid_Body *mass_list);
	double constraint_time_derivative(Circular_Rigid_Body *mass_list);

	// Partial of constraint function wrt x1, y1, x2, y2
	double jacobian_entry_x1(Circular_Rigid_Body *mass_list);
	double jacobian_entry_y1(Circular_Rigid_Body *mass_list);
	double jacobian_entry_x2(Circular_Rigid_Body *mass_list);
	double jacobian_entry_y2(Circular_Rigid_Body *mass_list);

	// partial of constraint time derivative wrt x1, y1, x2, y2
	double jacobian_derivative_entry_x1(Circular_Rigid_Body *mass_list);
	double jacobian_derivative_entry_y1(Circular_Rigid_Body *mass_list);
	double jacobian_derivative_entry_x2(Circular_Rigid_Body *mass_list);
	double jacobian_derivative_entry_y2(Circular_Rigid_Body *mass_list);

	void determine_initial_points(Circular_Rigid_Body *mass_list);

	// Draw to screen
	void draw_to_screen(Circular_Rigid_Body *mass_list, void *first_pixel, int pixel_buffer_width);

};



#endif /* RIGID_BAR_H_ */
