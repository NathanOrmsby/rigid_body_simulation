/*
 * bc_gradient.h
 *
 *  Created on: Nov 14, 2022
 *      Author: norms
 */
// Forward declaration because we only use a class pointer
class State_Getter;

#ifndef BC_GRADIENT_H_
#define BC_GRADIENT_H_

#define EPS 1.0e-14

// tol is convergence tolerance, default value of 0.001 normally
#define TOL 0.001

void biconjugate_gradient(State_Getter *state, double *b, double *x, int itol, double tol, int itmax, int *iter, double *err);



#endif /* BC_GRADIENT_H_ */
