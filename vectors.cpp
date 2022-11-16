#include <cmath>
#include "vectors.h"
// Vector class


// List of vector functions

// Magnitude
float vector_magnitude(Vector a)
{
	return sqrt((a.x * a.x) + (a.y * a.y));
}

// Direction, returns normalized vector
Vector vector_direction(Vector a)
{
	Vector unit;
	float mag = vector_magnitude(a);

	unit.x = a.x / mag;
	unit.y = a.y / mag;

	return unit;
}

// Addition
Vector vector_addition(Vector a, Vector b)
{
	Vector Result;
	Result.x = a.x + b.x;
	Result.y = a.y + b.y;

	return Result;
}

// Dot product: Returns a scaler
float vector_dp(Vector a, Vector b)
{
	return a.x * b.x + a.y * b.y;
}

// Angle between two vectors: Dot product method
float angle_between_vectors(Vector a, Vector b)
{
	float angle = acos(vector_dp(a, b) / (vector_magnitude(a) * vector_magnitude(b)));
	return angle;
}
