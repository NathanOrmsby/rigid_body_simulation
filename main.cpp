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
#include "get_state.h"
#include "utils.h"
#include "renderer.h"
#include "euler.h"
#include "rk4.h"
#include "my_timer.h"
#include "springs.h"


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
	HWND window = CreateWindow(window_class.lpszClassName, "My First Game!", WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, 1920, 1080, 0, 0, hInstance, 0);
	HDC hdc = GetDC(window);

	// Timer stuff
	//TimerMsec my_timer;
	//my_timer.Reset();
	//float time = 0.0;
	//float prev_time;
	//float dt;
	// Current time for managing the rendering of frames

//

	// ********* INITIAL CONDITIONS *************

	// Where the initial conditions are set.

	// MASSES
	int num_masses = 8;

	Circular_Rigid_Body masses[num_masses];
	masses[0].pos = {to_world((render_state.width / 2) - (render_state.width / 4)), to_world(render_state.height - 350)};

	masses[1].pos = {masses[0].pos.x - to_world(30), to_world(render_state.height - 400)};

	masses[2].pos = {masses[1].pos.x - to_world(70), to_world(render_state.height - 500)};

	masses[3].pos = {to_world(render_state.width / 2), to_world(render_state.height - 500)};

	masses[4].pos = {to_world(render_state.width / 2), to_world(render_state.height - 400)};

	masses[7].pos = {to_world((render_state.width / 2) + (render_state.width / 4)), to_world(render_state.height - 350)};

	masses[6].pos = {masses[7].pos.x + to_world(30), to_world(render_state.height - 400)};

	masses[5].pos = {masses[6].pos.x + to_world(70), to_world(render_state.height - 500)};

	// Fill in the rest of the parameters for the masses
	for (int i = 0; i < num_masses; i++)
	{
		masses[i].angle = 0.0;
		masses[i].linear_vel = {0.0, 0.0};
		masses[i].angular_vel = 0.0;
		masses[i].linear_accel = {0.0, 0.0};
		masses[i].angular_accel = 0.0;
		masses[i].force_ext = {0.0, 0.0};
		masses[i].torque = 0.0;
		masses[i].mass = 0.8;
		masses[i].radius = to_world(10);
		masses[i].determine_moi();
		masses[i].color = 0x0000FF;
	}
	// Constraints lists
	int num_bar1s = 2;
	int num_bar2s = 5;
	Rigid_Bar_1 bar1s[num_bar1s];
	bar1s[0].pivot = {to_world((render_state.width / 2) - (render_state.width / 4)), to_world(render_state.height - 250)};
	bar1s[0].attached_mass = 0;
	bar1s[0].determine_initial_point(masses);
	bar1s[0].rod = 1;
	bar1s[0].thickness = to_world(5);
	bar1s[0].color = 0xFFFFFF;

	bar1s[1].pivot = {to_world((render_state.width / 2) + (render_state.width / 4)), to_world(render_state.height - 250)};
	bar1s[1].attached_mass = 7;
	bar1s[1].determine_initial_point(masses);
	bar1s[1].rod = 1;
	bar1s[1].thickness = to_world(5);
	bar1s[1].color = 0xFFFFFF;

	Rigid_Bar_2 bar2s[num_bar2s];
	bar2s[0].attached_masses[0] = 0;
	bar2s[0].attached_masses[1] = 1;
	bar2s[0].determine_initial_points(masses);
	bar2s[0].rod = 1;
	bar2s[0].thickness = to_world(5);
	bar2s[0].color = 0xFFFFFF;

	bar2s[1].attached_masses[0] = 1;
	bar2s[1].attached_masses[1] = 2;
	bar2s[1].determine_initial_points(masses);
	bar2s[1].rod = 1;
	bar2s[1].thickness = to_world(5);
	bar2s[1].color = 0xFFFFFF;

	bar2s[2].attached_masses[0] = 3;
	bar2s[2].attached_masses[1] = 4;
	bar2s[2].determine_initial_points(masses);
	bar2s[2].rod = 1;
	bar2s[2].thickness = to_world(5);
	bar2s[2].color = 0xFFFFFF;

	bar2s[3].attached_masses[0] = 5;
	bar2s[3].attached_masses[1] = 6;
	bar2s[3].determine_initial_points(masses);
	bar2s[3].rod = 1;
	bar2s[3].thickness = to_world(5);
	bar2s[3].color = 0xFFFFFF;

	bar2s[4].attached_masses[0] = 6;
	bar2s[4].attached_masses[1] = 7;
	bar2s[4].determine_initial_points(masses);
	bar2s[4].rod = 1;
	bar2s[4].thickness = to_world(5);
	bar2s[4].color = 0xFFFFFF;

	// Springs lists
	int num_spring1s = 0;
	int num_spring2s = 4;

	Spring_2 spring2s[num_spring2s];
	spring2s[0].attached_masses[0] = 1;
	spring2s[0].attached_masses[1] = 4;

	spring2s[1].attached_masses[0] = 2;
	spring2s[1].attached_masses[1] = 3;

	spring2s[2].attached_masses[0] = 3;
	spring2s[2].attached_masses[1] = 5;

	spring2s[3].attached_masses[0] = 4;
	spring2s[3].attached_masses[1] = 6;

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
	current_state.ks = 0.5;
	current_state.kd = 0.5;


	// Timer stuff. Limiting the framerate, setting dt. SUPER AWESOME, to the microsecond
	LARGE_INTEGER ticks_per_second, start, current;
	int64_t ticks_per_loop;
	unsigned int loop_count = 0;

	// Find number of ticks per second
	QueryPerformanceFrequency(&ticks_per_second);

	// Find number of ticks per loop. I want 300 loops per second. Set dt
	ticks_per_loop = ticks_per_second.QuadPart / 300;
	float dt = 1.0 / 300;

	// Start the counter
	QueryPerformanceCounter(&start);

	// time the rk4 process. took 14 microseconds
//	my_timer_init();
//
//	my_timer_t timer = my_timer_start();
//
//	rk4(&current_state, masses, bar1s, num_bar1s, bar2s, num_bar2s, dt);
//
//	int64_t time = my_timer_end(&timer);
//
//	std::cout << "rk4 took " << time << " microseconds." << std::endl;

//	return 0;

	// Game loop
	while (running) {
		// Input
		MSG message;
		while (PeekMessage(&message, window, 0, 0, PM_REMOVE)) {
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		// Euler method

		// Checking masses
//		std::cout << "Printing mass properties" << std::endl;
//		for (int i = 0; i < current_state.num_bodies; i++)
//		{
//			std::cout << "Mass " << i << std::endl;
//			std::cout << "Position is: x: " << masses[i].pos.x << " y: " << masses[i].pos.y << std::endl;
//			std::cout << "Velocity is: vx: " << masses[i].linear_vel.x << " vy: " << masses[i].linear_vel.y << std::endl;
//			std::cout << "Acceleration is: ax: " << masses[i].linear_accel.x << " ay: " << masses[i].linear_accel.y << std::endl;
//		}
		// Calculate net forces
//		euler_state.get_current_state(euler_masses, bar1s, num_bar1s, bar2s, num_bar2s);
		// Plug net forces into ode solver and update values
//		euler_method(euler_masses, current_state.num_bodies, (&euler_state)->net_force_vector, dt);

//		std::cout << "Printing euler calculated net force: " << std::endl;
//		for (int i = 0; i < 3 * euler_state.num_bodies; i++)
//		{
//			std::cout << euler_state.net_force_vector[i] << " ";
//		}
//		std::cout << std::endl;

//		euler_state.free_all();
		// rk4 method
		loop_count++;
		// Do rk4
		// Try rk4
		rk4(&current_state, masses, bar1s, num_bar1s, bar2s, num_bar2s, spring2s, dt);



		//std::cout << "loop count is: " << loop_count << std::endl;

		 //Draw a frame if time is right
		if (loop_count >= 5)
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
//			std::cout << "c1: " << bar1s[0].constraint(masses) << std::endl;
//			std::cout << "c2: " << bar2s[0].constraint(masses) << std::endl;

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

		// Free all the shit
		//current_state.free_all();

//		std::cout << "Printing net force vector" << std::endl;
//		for (int i = 0; i < 3 * current_state.num_bodies; i++)
//		{
//			std::cout << current_state.net_force_vector[i] << " ";
//		}
//		std::cout << std::endl;
//
//		std::cout << "Printing out constraint function" << std::endl;
//		std::cout << "Constraint function: c = " << bar1s[0].constraint(masses) << std::endl;
//		std::cout << "Printing out state vector" << std::endl;
//		for (int i = 0; i < 3 * current_state.num_bodies; i++)
//		{
//			std::cout << current_state.state_vector[i] << " ";
//		}
//		std::cout << std::endl;
//
//		std::cout << "Printing position of mass" << std::endl;
//		std::cout << "x: " << masses[0].pos.x << " y: " << masses[0].pos.y << std::endl;
//		std::cout << "Printing velocity of mass" << std::endl;
//		std::cout << "vx: " << masses[0].linear_vel.x << " vy: " << masses[0].linear_vel.y << std::endl;
//		std::cout << "Printing acceleration of mass" << std::endl;
//		std::cout << "ax: " << masses[0].linear_accel.x << " ay: " << masses[0].linear_accel.y << std::endl;




		// Print out the external force vector
//		if (count == 0)
//		{
//			std::cout << "Printing external force vector" << std::endl;
//			for (int i = 0; i < 3 * current_state.num_bodies; i++)
//			{
//				std::cout << current_state.force_ext_vector[i] << " ";
//			}
//			std::cout << std::endl;
//
//			std::cout << "Printing constraint force vector" << std::endl;
//			for (int i = 0; i < 3 * current_state.num_bodies; i++)
//			{
//				std::cout << current_state.constraint_force_vector[i] << " ";
//			}
//			std::cout << std::endl;
//
//			std::cout << "Printing net force vector" << std::endl;
//			for (int i = 0; i < 3 * current_state.num_bodies; i++)
//			{
//				std::cout << current_state.net_force_vector[i] << " ";
//			}
//			std::cout << std::endl;
//
//			std::cout << "Printing out mass properties" << std::endl;
//			for (int i = 0; i < current_state.num_bodies; i++)
//			{
//				std::cout << "ax: " << masses[i].linear_accel.x << " ay: " << masses[i].linear_accel.y << std::endl;
//				std::cout << "vx: " << masses[i].linear_vel.x << " vy: " << masses[i].linear_vel.y << std::endl;
//				std::cout << "x: " << masses[i].pos.x << " y: " << masses[i].pos.y << std::endl;
//			}
//
//			std::cout << "Printing out constraint function" << std::endl;
//			std::cout << "Constraint function: c = " << bar1s[0].constraint_equation(masses) << std::endl;
//			std::cout << "Printing out state vector" << std::endl;
//			for (int i = 0; i < 3 * current_state.num_bodies; i++)
//			{
//				std::cout << current_state.state_vector[i] << " ";
//			}
//			std::cout << std::endl;
//		}




		// Print out the jacobian and jacobian time derivative
//		std::cout << "Printing the jacobian" << std::endl;
//		print_matrix_block_list((&current_state)->jacobian, 1);

//		current_state.free_all();
//
//		std::cout << "made it here 3" << std::endl;
//
//		time += dt;
//		count++;


		// Clear screen

		//draw_line(0, 0, 500, 500, 0x0000FF, render_state.memory, render_state.width);

	}

	// Success
	return 0;

}
