/*
 * vectors.h
 *
 *  Created on: Nov 7, 2022
 *      Author: norms
 */

#ifndef VECTORS_H_
#define VECTORS_H_

// Vector class
class Vector
{
	public:

	// Coordinates for x and y directions
	double x;
	double y;
};

class Screen_Vector
{
	public:

	// Coordinates for x and y directions
	int x;
	int y;
};

// List of vector functions

// Magnitude
double vector_magnitude(Vector);

// Direction, returns normalized vector
Vector vector_direction(Vector);

// Addition
Vector vector_addition(Vector, Vector);

// Dot product: Returns a scaler
double vector_dp(Vector, Vector);

// Angle between two vectors: Dot product method
double angle_between_vectors(Vector, Vector);

#endif /* VECTORS_H_ */
