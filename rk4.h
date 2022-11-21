/*
 * rk4.h
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */

#ifndef RK4_H_
#define RK4_H_

class State_Getter;
class Circular_Rigid_Body;
class Rigid_Bar_1;
class Rigid_Bar_2;
class Spring_2;

// Runge Kutta 4th Order to solve differential equations

void rk4(State_Getter *state, Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, Spring_2 *spring2s, double dt, int loop_count);

// Fills a mass list for another state
void copy_mass_list(Circular_Rigid_Body *old_list, Circular_Rigid_Body *new_list, int num_bodies);

#endif /* RK4_H_ */
