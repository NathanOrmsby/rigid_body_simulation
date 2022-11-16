#include "vectors.h"
#include "renderer.h"
#include "rigid_bodies.h"
// Class for a mass that acts like a point mass, and is drawn in the shape of a circle


void Circular_Rigid_Body::draw_to_screen(void *first_pixel, int pixel_buffer_width, int pixel_buffer_height)
{
	draw_circle(to_screen(pos.x), to_screen(pos.y), to_screen(radius), color, first_pixel, pixel_buffer_width, pixel_buffer_height);
}

void Circular_Rigid_Body::determine_moi(void)
{
	moi = (mass * (radius * radius)) / 2.0;
}

void Circular_Rigid_Body::set_initial_force_ext(void)
{
	// Initial external force is gravitational
	// Gravitational
	force_ext.y = mass * -9.8;
	force_ext.x = 0.0;
	torque = 0.0;

	// Then apply spring force if a spring is attached using the spring class method
}





