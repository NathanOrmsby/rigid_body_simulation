#include <windows.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <chrono>
#include <time.h>
#include "math.h"
#include <cmath>
#include <stdint.h>

// Custom header files
#include "matrix_stuff.h"
#include "rigid_bodies.h"
#include "constraint_bodies.h"
#include "utils.h"
#include "renderer.h"
#include "euler.h"
#include "get_state.h"
#include "rk4.h"
#include "my_timer.h"
#include "springs.h"
#include "total_energy.h"


using namespace std;

// Globals
bool running = true;

struct Render_State
{
	int height, width;
	void *memory;

	BITMAPINFO bitmap_info;
};

Render_State render_state;

LRESULT CALLBACK window_callback(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {

	LRESULT result;

	// If user interacts with the window e.g. clicks the close button
	switch (uMsg) {
		case WM_CLOSE:
		case WM_DESTROY: {
			running = false;
		} break;

		case WM_SIZE: {
			RECT rect;
			GetClientRect(hwnd, &rect);
			render_state.width = rect.right - rect.left;
			render_state.height = rect.bottom - rect.top;

			int buffer_size = render_state.width * render_state.height * sizeof(unsigned int);

			// Check if the pointer already has memory to be freed before reassigning
			if (render_state.memory) {
				VirtualFree(render_state.memory, 0, MEM_RELEASE);
			}

			render_state.memory = VirtualAlloc(0, buffer_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);

			render_state.bitmap_info.bmiHeader.biSize = sizeof(render_state.bitmap_info.bmiHeader);
			render_state.bitmap_info.bmiHeader.biWidth = render_state.width;
			render_state.bitmap_info.bmiHeader.biHeight = render_state.height;
			render_state.bitmap_info.bmiHeader.biPlanes = 1;
			render_state.bitmap_info.bmiHeader.biBitCount = 32;
			render_state.bitmap_info.bmiHeader.biCompression = BI_RGB;
		} break;

		// For other cases do the window default procedure
		default: {
			result = DefWindowProc(hwnd, uMsg, wParam, lParam);
		}
	}
	return result;
}

int WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd) {

	// Create window class
	WNDCLASS window_class = {};
	window_class.style = CS_HREDRAW | CS_VREDRAW;
	window_class.lpszClassName = "Game Window Class";
	window_class.lpfnWndProc = window_callback;

	// Register Class
	RegisterClass(&window_class);

	// Create window
	HWND window = CreateWindow(window_class.lpszClassName, "My First Game!", WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, 1280, 720, 0, 0, hInstance, 0);
	HDC hdc = GetDC(window);

	// ********* INITIAL CONDITIONS *************

	// Where the initial conditions are set.

	// MASSES
	int num_masses = 2;

	Circular_Rigid_Body masses[num_masses];

	masses[0].pos = {to_world(render_state.width / 2 + 100), to_world(render_state.height / 2 + 100)};
	masses[1].pos = {masses[0].pos.x, masses[0].pos.y + to_world(100)};



	// Fill in the rest of the parameters for the masses
	for (int i = 0; i < num_masses; i++)
	{
		masses[i].angle = 0.0;
		masses[i].linear_vel = {0.0, 0.0};
		masses[i].angular_vel = 0.0;
		masses[i].force_ext = {0.0, 0.0};
		masses[i].torque = 0.0;
		masses[i].mass = 1.0;
		masses[i].radius = to_world(10);
		masses[i].determine_moi();
		masses[i].color = 0x0000FF;
	}
	// Constraints lists
	int num_bar1s = 1;
	int num_bar2s = 1;

	Rigid_Bar_1 bar1s[num_bar1s];
	bar1s[0].pivot = {to_world((render_state.width / 2)), to_world(render_state.height / 2)};
	bar1s[0].attached_mass = 0;
	bar1s[0].determine_initial_point(masses);
	bar1s[0].rod = 0;
	bar1s[0].thickness = to_world(5);
	bar1s[0].color = 0xFFFFFF;


	Rigid_Bar_2 bar2s[num_bar2s];


	bar2s[0].attached_masses[0] = 0;
	bar2s[0].attached_masses[1] = 1;
	bar2s[0].determine_initial_points(masses);
	bar2s[0].rod = 0;
	bar2s[0].thickness = to_world(5);
	bar2s[0].color = 0xFFFFFF;

	// Springs lists
	//int num_spring1s = 0;
	int num_spring2s = 0;

	Spring_2 spring2s[num_spring2s];
//	spring2s[0].attached_masses[0] = 1;
//	spring2s[0].attached_masses[1] = 4;

	for (int i = 0; i < num_spring2s; i++)
	{
		spring2s[i].eq_length = to_world(75);
		spring2s[i].spring_const = 5.0;
		spring2s[i].color = 0xFFFF00;
	}



	// Create the State_Getter class for current state
	State_Getter current_state;
	// Number of masses and constraints
	current_state.num_bodies = num_masses;
	current_state.num_constraints = num_bar1s + num_bar2s;
	current_state.num_spring2s = num_spring2s;

	// Timer stuff. Limiting the framerate, setting dt. SUPER AWESOME, to the microsecond
	LARGE_INTEGER ticks_per_second, start, current;
	int64_t ticks_per_loop;
	unsigned int loop_count = 0;

	// Find number of ticks per second
	QueryPerformanceFrequency(&ticks_per_second);

	// Find number of ticks per loop. I want 300 loops per second. Set dt
	ticks_per_loop = ticks_per_second.QuadPart / 1000;
	double dt = 1.0 / 1000.0;

	// Start the counter
	QueryPerformanceCounter(&start);

	// time the rk4 process. took 14 microseconds
//	my_timer_init();
//
//	my_timer_t timer = my_timer_start();

	rk4(&current_state, masses, bar1s, num_bar1s, bar2s, num_bar2s, spring2s, dt, loop_count);

//	int64_t time = my_timer_end(&timer);
//
//	std::cout << "rk4 took " << time << " microseconds." << std::endl;

//  return 0;

	// Game loop
	while (running) {
		// Input
		MSG message;
		while (PeekMessage(&message, window, 0, 0, PM_REMOVE)) {
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		loop_count++;

		// rk4 method
		rk4(&current_state, masses, bar1s, num_bar1s, bar2s, num_bar2s, spring2s, dt, loop_count);

		 //Draw a frame if time is right
		if (loop_count >= 17)
		{

			// Clear screen
			clear_screen(0xff0000, render_state.memory, render_state.width, render_state.height);

			// Draw all masses
			for (int i = 0; i < num_masses; i++)
			{
				masses[i].draw_to_screen(render_state.memory, render_state.width, render_state.height);
			}

			// Draw constraint objects
			for (int i = 0; i < num_bar1s; i++)
			{
				bar1s[i].draw_to_screen(masses, render_state.memory, render_state.width);
			}

			for (int i = 0; i < num_bar2s; i++)
			{
				bar2s[i].draw_to_screen(masses, render_state.memory, render_state.width);
			}

			for (int i = 0; i < num_spring2s; i++)
			{
				spring2s[i].draw_spring(masses, render_state.memory, render_state.width);
			}

			// Reset loop count to prepare for next frame
			loop_count = 0;

//			std::cout << "Constraints" << std::endl;
			std::cout << "c1: " << bar1s[0].constraint(masses) << std::endl;
			std::cout << "Total energy: " << total_energy(masses, num_masses) << std::endl;
		}

		// Render
		StretchDIBits(hdc, 0, 0, render_state.width, render_state.height, 0, 0, render_state.width, render_state.height, render_state.memory, &render_state.bitmap_info, DIB_RGB_COLORS, SRCCOPY);
		// Timer do while loop until next cycle
		do
		{
			QueryPerformanceCounter(&current);
		}
		while(current.QuadPart - start.QuadPart < ticks_per_loop);

		start = current;
	}

	// Success
	return 0;

}
