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
	float x;
	float y;
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
float vector_magnitude(Vector);

// Direction, returns normalized vector
Vector vector_direction(Vector);

// Addition
Vector vector_addition(Vector, Vector);

// Dot product: Returns a scaler
float vector_dp(Vector, Vector);

// Angle between two vectors: Dot product method
float angle_between_vectors(Vector, Vector);

#endif /* VECTORS_H_ */
