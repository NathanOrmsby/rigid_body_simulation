/*
 * euler.h
 *
 *  Created on: Nov 14, 2022
 *      Author: norms
 */

#include "rigid_bodies.h"
#include "get_state.h"
#include "vectors.h"

#ifndef EULER_H_
#define EULER_H_

// Eulers method
void euler_method(Circular_Rigid_Body *mass_list, int num_bodies, float *net_force_vector, float dt);

#endif /* EULER_H_ */
