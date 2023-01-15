#include <cmath>
#include "vectors.h"
// Vector class


// List of vector functions

// Minimum and Maximum
Vector minimum_vector(Vector a, Vector b)
{
	Vector result;
	// X coordinate
	if (a.x < b.x)
	{
		result.x = a.x;
	}
	else
	{
		result.x = b.x;
	}
	// Y coordinate
	if (a.y < b.y)
	{
		result.y = a.y;
	}
	else
	{
		result.y = b.y;
	}
	return result;
}

Vector maximum_vector(Vector a, Vector b)
{
	Vector result;
	// X coordinate
	if (a.x < b.x)
	{
		result.x = b.x;
	}
	else
	{
		result.x = a.x;
	}
	// Y coordinate
	if (a.y < b.y)
	{
		result.y = b.y;
	}
	else
	{
		result.y = a.y;
	}
	return result;
}

// Magnitude
double vector_magnitude(Vector a)
{
	return sqrt((a.x * a.x) + (a.y * a.y));
}

// Direction, returns normalized vector
Vector vector_direction(Vector a)
{
	Vector unit;
	double mag = vector_magnitude(a);

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
double vector_dp(Vector a, Vector b)
{
	return a.x * b.x + a.y * b.y;
}

// Angle between two vectors: Dot product method
double angle_between_vectors(Vector a, Vector b)
{
	double angle = (double)acos(vector_dp(a, b) / (vector_magnitude(a) * vector_magnitude(b)));
	return angle;
}
