
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "vectors.h"
#include "utils.h"
#include "renderer.h"

// FUNCTIONS THAT HAVE TO DO WITH DRAWING TO THE SCREEN
// Move from world to screen coordinates
// SCALE FROM WORLD TO SCREEN
float SCALE = 100.0;
int to_screen(double coordinate)
{
	return (int)(SCALE * coordinate);
}

double to_world(int coordinate)
{
	return (double)(coordinate / SCALE);
}

// Clearing the screen
void clear_screen(unsigned int color, void *first_pixel, int pixel_buffer_width, int pixel_buffer_height)
{
	unsigned int* pixel = (unsigned int*)first_pixel;
	for (int y = 0; y < pixel_buffer_height; y++)
	{
		for (int x = 0; x < pixel_buffer_width; x++)
		{
			*pixel++ = color;
		}
	}
}

// Circle
inline int clamp(int min, int val, int max)
{
	if (val < min) return min;
	if (val > max) return max;
	return val;
}

void draw_circle(int x, int y, int radius, unsigned int color, void *first_pixel, int pixel_buffer_width, int pixel_buffer_height)
{
	y = clamp(radius, y, pixel_buffer_height - radius);
	x = clamp(radius, x, pixel_buffer_width);
	// loop over top right quarter of square
	for (int row = 0; row < radius + 1; row++)
	{
		// Assign current pixel address
		unsigned int* pixel = (unsigned int*)first_pixel + x + (y+row)*pixel_buffer_width;
		for (int col = 0; col < radius + 1; col++)
		{
			// Use distance formula. reduces to below form if considering xprime = x + col, yprime = y + row
			float distance = sqrt((float)(col * col) + (float)(row * row));

			if (distance <= radius)
			{
				// Draw the top quarter circle
				*pixel = color;
				// Flip over the horizontal: same x value
				unsigned int *hfpixel = pixel - 2*row*pixel_buffer_width;
				*hfpixel = color;
				// Flip over the vertical same y value
				unsigned int *vfpixel = pixel - 2 * col;
				*vfpixel = color;
				// Flip over both
				unsigned int *vhfpixel = pixel - 2 * col - 2*row*pixel_buffer_width;
				*vhfpixel = color;
			}
			// Increment pixel
			pixel++;
		}
	}
}

// Calculate pixel numbers of a line to screen giving world coordinates

void draw_line(int x0, int y0, int x1, int y1, unsigned int color, void *first_pixel, int pixel_buffer_width)
{
	// Check for vertical line so no division by zero
	float m;
	if (x0 == x1)
	{
		m = 9999;
	}
	else
	{
		// Find slope of the line given the two points
		m = (float)(y1 - y0) / (float)(x1 - x0);
	}

	// Find the y intercept of the line
	float b = y0 - m * x0;

	// Protocol for horizontalish lines, 1 y value per x value.
	if (fabs(m) <= 1.5)
	{
		// Check if the line is positive sloped
		// Calculate array of y values for each x coordinate within range
		float *y_array;
		if (x0 < x1)
		{
			y_array = (float*)malloc(((x1 + 1) - x0) * sizeof(float));
		}
		else
		{
			y_array = (float*)malloc(((x0 + 1) - x1) * sizeof(float));
		}

		// Fill the y_array using line equation solving for y
		int count = 0;
		if (x0 < x1)
		{
			for (int x = x0; x < x1 + 1; x++)
			{
				y_array[count] = m * x + b;
				count++;
			}
		}
		else
		{
			for (int x = x1; x < x0 + 1; x++)
			{
				y_array[count] = m * x + b;
				count++;
			}
		}

		// Loop through pixels
		count = 0;
		int point_count = 0;
		if (m < 0)
		{
			// Line goes right
			if (x0 < x1)
			{

				for (int row = y1; row < y0 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					// Reset count to zero to count with columns
					count = 0;
					for (int col = x0; col < x1 + 1; col++)
					{
						// Calculate vertical distance from point
						float vert_dist = (float)fabs(y_array[count] - row);

						// If distance is smaller than 1, color the pixel
						if (vert_dist < 1)
						{
							point_count++;
							*pixel = color;
						}
						// Increment the pixel address count
						pixel++;
						// Increment count
						count++;
					}
				}
			}
			// Line goes left
			else
			{
				for (int row = y0; row < y1 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x1 + row*pixel_buffer_width;

					// Reset count to zero to count with columns
					count = 0;
					for (int col = x1; col < x0 + 1; col++)
					{
						// Calculate vertical distance from point
						float vert_dist = (float)fabs(y_array[count] - row);

						// If distance is smaller than 1, color the pixel
						if (vert_dist < 1)
						{
							point_count++;
							*pixel = color;
						}
						// Increment the pixel address count
						pixel++;
						// Increment count
						count++;
					}
				}
			}
		}
		// Positive slope
		else
		{
			// Line goes right
			if (x0 < x1)
			{
				for (int row = y0; row < y1 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					// Reset count to zero to count with columns
					count = 0;
					for (int col = x0; col < x1 + 1; col++)
					{
						// Calculate vertical distance from point
						float vert_dist = (float)fabs(y_array[count] - row);

						// If distance is smaller than 1, color the pixel
						if (vert_dist < 1)
						{
							point_count++;
							*pixel = color;
						}
						// Increment the pixel address count
						pixel++;
						// Increment count
						count++;
					}
				}
			}
			// Line goes left
			else
			{
				for (int row = y1; row < y0 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x1 + row*pixel_buffer_width;

					// Reset count to zero to count with columns
					count = 0;
					for (int col = x1; col < x0 + 1; col++)
					{
						// Calculate vertical distance from point
						float vert_dist = (float)fabs(y_array[count] - row);

						// If distance is smaller than 1, color the pixel
						if (vert_dist < 1)
						{
							point_count++;
							*pixel = color;
						}
						// Increment the pixel address count
						pixel++;
						// Increment count
						count++;
					}
				}
			}
		}
		// Free the y_array
		free(y_array);
	}
	// Protocol for verticalish lines
	else
	{
		// Calculate array of y values for each x coordinate within range
		float *x_array;
		if (y0 < y1)
		{
			x_array = (float*)malloc(((y1 + 1) - y0) * sizeof(float));
		}
		else
		{
			x_array = (float*)malloc(((y0 + 1) - y1) * sizeof(float));
		}

		// Fill the x_array using line equation solving for x
		int count = 0;
		// Vertical line
		if (m == 9999)
		{
			// Line going up
			if (y0 < y1)
			{
				for (int y = y0; y < y1 + 1; y++)
				{
					x_array[count] = x0;
					count++;
				}
			}
			// Line going down
			else
			{
				for (int y = y1; y < y0 + 1; y++)
				{
					x_array[count] = x0;
					count++;
				}
			}
		}
		// positive slope
		else if (m > 0)
		{
			// Line going right
			if (y0 < y1)
			{
				for (int y = y0; y < y1 + 1; y++)
				{
					x_array[count] = (y - b) / m;
					count++;
				}
			}
			// Line going left
			else
			{
				for (int y = y1; y < y0 + 1; y++)
				{
					x_array[count] = (y - b) / m;
					count++;
				}
			}
		}
		// negative slope
		else
		{
			// Line going right
			if (y0 > y1)
			{
				for (int y = y1; y < y0 + 1; y++)
				{
					x_array[count] = (y - b) / m;
					count++;
				}
			}
			// Line going left
			else
			{
				for (int y = y0; y < y1 + 1; y++)
				{
					x_array[count] = (y - b) / m;
					count++;
				}
			}
		}

		// Loop through pixels
		count = 0;

		// For vertical line.
		if (m == 9999)
		{
			// Line from bottom to top
			if (y0 < y1)
			{
				for (int row = y0; row < y1 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					for (int col = x0; col < x1 + 1; col++)
					{
						// If line is vertical, all points are on same column
						if (col == x0)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}
				}
			}
			// Line from top to bottom
			else
			{
				for (int row = y1; row < y0 + 1; row++)
				{
					// x0 == x1 in this case
					// Assign current pixel address: x0 == x1 in this case
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					for (int col = x0; col < x1 + 1; col++)
					{
						// If distance is smaller than 1, color the pixel
						if (col == x0)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}
				}
			}
		}
		else if (m > 0)
		{
			// Line going from left to right
			if (y0 < y1)
			{
				// Count with rows for vertical lines
				count = 0;
				//printf("y0: %d, y1: %d\n", y0, y1);
				for (int row = y0; row < y1 + 1; row++)
				{
					//printf("Count is: %d", count);
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					for (int col = x0; col < x1 + 1; col++)
					{
						// Calculate horizontal distance from point
						float horizontal_dist = (float)fabs(x_array[count] - col);

						// If distance is smaller than 1, color the pixel
						if (horizontal_dist < 1)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}

					// Increment the count
					count++;
				}
			}
			// Line going from right to left
			else
			{
				// Count with rows for vertical lines
				count = 0;
				//printf("y0: %d, y1: %d\n", y0, y1);
				for (int row = y1; row < y0 + 1; row++)
				{
					//printf("Count is: %d", count);
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x1 + row*pixel_buffer_width;

					for (int col = x1; col < x0 + 1; col++)
					{
						// Calculate horizontal distance from point
						float horizontal_dist = (float)fabs(x_array[count] - col);

						// If distance is smaller than 1, color the pixel
						if (horizontal_dist < 1)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}

					// Increment the count
					count++;
				}
			}
		}
		// slope is negative
		else
		{
			// Line going from left to right
			// Count with rows for vertical lines
			count = 0;
			if (y0 > y1)
			{
				for (int row = y1; row < y0 + 1; row++)
				{
					// Assign current pixel address
					unsigned int* pixel = (unsigned int*)first_pixel + x0 + row*pixel_buffer_width;

					for (int col = x0; col < x1 + 1; col++)
					{
						// Calculate horizontal distance from point
						float horizontal_dist = (float)fabs(x_array[count] - col);

						// If distance is smaller than 1, color the pixel
						if (horizontal_dist < 1)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}

					// Increment count
					count++;
				}
			}
			// Line going from right to left
			else
			{
				// Count with rows
				count = 0;
				for (int row = y0; row < y1 + 1; row++)
				{
					// Assign current pixel address
					unsigned int *pixel = (unsigned int*)first_pixel + x1 + row*pixel_buffer_width;

					for (int col = x1; col < x0 + 1; col++)
					{
						// Calculate horizontal distance from point
						float horizontal_dist = (float)fabs(x_array[count] - col);

						// If distance is smaller than 1, color the pixel
						if (horizontal_dist < 1)
						{
							*pixel = color;
						}

						// Increment the pixel address count
						pixel++;
					}
					// Increment count
					count++;
				}
			}
		}

		// Free the x_array
		free(x_array);
	}
}

// Line for pendulum. First point will always be on top. Draws a rectangle around a center line acting as a thick rod or string. GIVE IN SCREEN COORDINATES
void draw_rod(float x0, float y0, float x1, float y1, float thickness, unsigned int color, void *first_pixel, int pixel_buffer_width)
{

	// Calculate the slope of the line
	// Check for vertical line so no division by zero
	float m;
	if (x0 == x1)
	{
		m = 9999;
	}
	else
	{
		// Find slope of the line given the two points
		m = (y1 - y0) / (x1 - x0);
	}

	// Find the four points of the rectangle

	// Find perpendicular slope
	float m_perp;
	// Check for horizontal line
	if (m == 0)
	{
		m_perp = 9999;
	}
	// Check for vertical line
	else if (m == 9999)
	{
		m_perp = 0;
	}
	else
	{
		m_perp = -1 / m;
	}

	// Find angle between thickness vector and horizontal
	float angle;
	// Convert thickness slope to vector with same slope
	if (m_perp == 9999)
	{
		// We have a 90 degree angle
		angle = 1.5707;
	}
	else if (m_perp == 0)
	{
		angle = 0;
	}
	else
	{
		Vector perp;
		perp = {1, 1 * m_perp};
		Vector unit_x;
		unit_x = {1, 0};
		angle = angle_between_vectors(perp, unit_x);
	}

	// Find x and y coordinates of new points
	float rx0, ry0, rx1, ry1, rx2, ry2, rx3, ry3;
	// Top right
	rx0 = x0 + (thickness / 2) * cos(angle);
	ry0 = y0 + (thickness / 2) * sin(angle);
	// Bottom right
	rx1 = x1 + (thickness / 2) * cos(angle);
	ry1 = y1 + (thickness / 2) * sin(angle);
	// Bottom left
	rx2 = x1 - (thickness / 2) * cos(angle);
	ry2 = y1 - (thickness / 2) * sin(angle);
	// Top left
	rx3 = x0 - (thickness / 2) * cos(angle);
	ry3 = y0 - (thickness / 2) * sin(angle);

	// Draw the borders of the rectangle

	draw_line(rx0, ry0, rx1, ry1, color, first_pixel, pixel_buffer_width);
	draw_line(rx1, ry1, rx2, ry2, color, first_pixel, pixel_buffer_width);
	draw_line(rx2, ry2, rx3, ry3, color, first_pixel, pixel_buffer_width);
	draw_line(rx3, ry3, rx0, ry0, color, first_pixel, pixel_buffer_width);

	// Fill in the borders. NOT WORKING

	// *** TODO **** THIS ONLY WORKS IF LINE GOES LEFT TO RIGHT AND TOP TO BOTTOM. FIX THIS
//	int min_y = minimum(ry1, ry2);
//	int max_y = maximum(ry0, ry3);
//	int min_x = minimum(rx2, rx3);
//	int max_x = maximum(rx0, rx1);
//
//	if (min_y == 9999) min_y = ry1;
//	if (max_y == 9999) max_y = ry0;
//	if (min_x == 9999) min_x = rx2;
//	if (max_x == 9999) max_x = rx0;
//
//	//printf("ymin: %d, ymax: %d, xmin: %d, xmax: %d\n", min_y, max_y, min_x, max_x);
//
//	for (int row = min_y + 1; row < max_y; row++)
//	{
//		// Assign current pixel address
//		unsigned int *pixel = (unsigned int*)first_pixel + min_x + row*pixel_buffer_width;
//
//		bool fill = false;
//		for (int col = min_x; col < max_x + 1; col++)
//		{
//			if (*pixel == color && *(pixel - 1) != color)
//			{
//				//printf("First pixel!\n");
//				fill = true;
//			}
//			else if (*pixel == color && *(pixel + 1) != color)
//			{
//				//printf("Last pixel!\n");
//				fill = false;
//			}
//			else if (fill && *pixel != color)
//			{
//				*pixel = color;
//			}
//
//			// Increment pixel
//			pixel++;
//		}
//	}
}

//void draw_spring_coils(float x0, float y0, float x1, float y1, int num_coils, unsigned int color)
//{
//	// Draws a spring between two points with defined number of coils
//
//	// Calculate slope between two points
//	// Check for vertical line so no division by zero
//	float m;
//	if (x0 == x1)
//	{
//		m = 9999;
//	}
//	else
//	{
//		// Find slope of the line given the two points
//		m = (float)(y1 - y0) / (float)(x1 - x0);
//	}
//
//	// Find points given coils
//
//}


