/*
 * renderer.h
 *
 *  Created on: Nov 7, 2022
 *      Author: norms
 */

#ifndef RENDERER_H_
#define RENDERER_H_

// FUNCTIONS THAT HAVE TO DO WITH DRAWING TO THE SCREEN
// Move from world to screen coordinates and vice versa
int to_screen(double coordinate);

double to_world(int coordinate);

// Clearing the screen
void clear_screen(unsigned int color, void *first_pixel, int pixel_buffer_width, int pixel_buffer_height);

// Clamp function, ensures value is between the minimum and the maximum. Prevents drawing outside the screen
inline int clamp(int min, int val, int max);

// Draws a circle around a given position. Has safety check so that circle may not be drawn outside the screen preventing crashes
void draw_circle(int x, int y, int radius, unsigned int color, void *first_pixel, int pixel_buffer_width, int pixel_buffer_height);

// Draws a line between two integer point values that can be any orientation.
void draw_line(int x0, int y0, int x1, int y1, unsigned int color, void *first_pixel, int pixel_buffer_width);

// Draws a rectangle around a line to imitate a rod of certain thickness
void draw_rod(float x0, float y0, float x1, float y1, float thickness, unsigned int color, void *first_pixel, int pixel_buffer_width);

#endif /* RENDERER_H_ */
