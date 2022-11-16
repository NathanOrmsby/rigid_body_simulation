
// Include header of Vector class because there are actual variables being defined
#include "vectors.h"

#ifndef MASS_H_
#define MASS_H_
// Class for a circular rigid body that is drawn in the shape of a circle

class Circular_Rigid_Body
{
	public:
	Vector pos;
	float angle;
	Vector linear_vel;
	float angular_vel;
	Vector linear_accel;
	float angular_accel;
	Vector force_ext;
	float torque;
	float mass;
	float radius;
	float moi;
	unsigned int color;

	void draw_to_screen(void *first_pixel, int pixel_buffer_width, int pixel_buffer_height);

	void determine_moi(void);

	void set_initial_force_ext(void);
};



#endif /* MASS_H_ */
