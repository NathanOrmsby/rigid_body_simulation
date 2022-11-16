#include <stdio.h>
#include <math.h>
#include "matrix_stuff.h"
#include "get_state.h"
#include "bc_gradient.h"


// Linear biconjugate gradient method: Solves systems of linear equations aka our matrix problem to solve for constraint forces
// Uses doubles for more precision
// m is num_constraints
void biconjugate_gradient(State_Getter *state, float *b, float *x, int itol, float tol, int itmax, int *iter, float *err)
{
	// This algorithm MINIMIZES Gradient F = Ax - b to converge within specified value
	// Generates a succession of search directions and improved minimizers xk.
	// After N iterations of an NxN matrix you arrive at the minimizer over the entire space, aka the solution
	// Double numbers for stuff
	float ak,akden,bk,bkden,bknum,b_norm,dx_norm,x_norm,zm1_norm,z_norm;

	// Set num constraint so i dont have to write it everywhere
	int m = state->num_constraints;

//	std::cout << "Printing elements of x before biconjugate gradient" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << x[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Printing b vector" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << b[i] << " ";
//	}
//	std::cout << std::endl;

	// All the vectors
	float p[m];
	float p_bar[m];
	float r[m];
	float r_bar[m];
	float z[m];
	float z_bar[m];

	// zero all the vectors
	zero_vector(p, m);
	zero_vector(p_bar, m);
	zero_vector(r, m);
	zero_vector(r_bar, m);
	zero_vector(z, m);
	zero_vector(z_bar, m);


	// Start iteration count
	*iter = 0;

	// Sets A * x = r
	A_times_x(state, x, r);

//	std::cout << "Printing elements of x after A_times_x(x, r) biconjugate gradient" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << x[i] << " ";
//	}
//	std::cout << std::endl;

//	std::cout << "Printing elements of r after r = A*x" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << r[i] << " ";
//	}
//	std::cout << std::endl;

	// residual vector "r" is calculated by subtracting r from b.
	for (int i = 0; i < m; i++)
	{
		r[i] = b[i] - r[i];
		//r_bar[i] = r[i];
	}

//	std::cout << "Printing initial residual vector" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << r[i] << " ";
//	}
//	std::cout << std::endl;

	// Minimum residual version of algorithm: sets r_bar = A * r at start
//	std::cout << "Printing initial residual bar vector" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << r_bar[i] << " ";
//	}
//	std::cout << std::endl;
	A_times_x(state, r, r_bar);

//	std::cout << "Printing after residual bar vector" << std::endl;
//	for (int i = 0; i < m; i++)
//	{
//		std::cout << r_bar[i] << " ";
//	}
//	std::cout << std::endl;

	// Plug in an identity matrix for preconditioned M, we don't calculate A
	float identity[m];
	ones_vector(identity, m);

	// Different itol protocols
	if (itol == 1)
	{
		// Calculates magnitude of b
		b_norm = calculate_relevant_norm(b, m, itol);
		// Sets z = r
		solve_for_x(m, identity, r, z);

//		std::cout << "Printing b_norm: " << std::endl;
//		std::cout << "b_norm is: " << b_norm << std::endl;
//		std::cout << "Printing z" << std::endl;
//		for (int i = 0; i < m; i++)
//		{
//			std::cout << z[i] << " ";
//		}
//		std::cout << std::endl;
	}
	else if (itol == 2)
	{
		// Sets z = b
		solve_for_x(m, identity, b, z);
		// Sets b_norm equal to magnitude of largest component of z
		b_norm = calculate_relevant_norm(z, m, itol);
		// Sets z equal to r
		solve_for_x(m, identity, r, z);
	}
	else if (itol == 3 || itol == 4)
	{
	solve_for_x(m, identity, b, z);
	b_norm = calculate_relevant_norm(z, m, itol);
	solve_for_x(m, identity, r, z);
	z_norm = calculate_relevant_norm(z, m, itol);
	}
	else
	{
		std::cout << "Itol protocols only 1-4, exiting..." << std::endl;
		exit(1);
	}

	// Main loop
	while (*iter <= itmax)
	{
		// Increment the count
		(*iter)++;

		//std::cout << "iter is: " << *iter << std::endl;

		// Sets z_bar = r_bar
		solve_for_x(m, identity, r_bar, z_bar);

//		std::cout << "Printing z_bar" << std::endl;
//		for (int i = 0; i < m; i++)
//		{
//			std::cout << z_bar[i] << " ";
//		}
//		std::cout << std::endl;

		// Calculate bk_numerator
		bknum = 0.0;
		for (int i = 0; i < m; i++)
		{
			bknum += z[i] * r_bar[i];
		}
//		std::cout << "Printing bk numerator" << std::endl;
//		std::cout << "bk_num: " << bknum << std::endl;

		// Calculate coefficient bk, and direction vectors p and pp.
		// First iteration count, bk is not calculated, so bknum is now k+1 and bkden is k
		if (*iter == 1)
		{
			for (int i = 0; i < m; i++)
			{
				p[i] = z[i];
				p_bar[i] = z_bar[i];
			}

//			std::cout << "Printing p" << std::endl;
//			for (int i = 0; i < m; i++)
//			{
//				std::cout << p[i] << " ";
//			}
//			std::cout << std::endl;
//			std::cout << "Printing p bar" << std::endl;
//			for (int i = 0; i < m; i++)
//			{
//				std::cout << p_bar[i] << " ";
//			}
//			std::cout << std::endl;
		}
		else
		{
			bk = bknum / bkden;
			for (int i = 0; i < m; i++)
			{
				p[i] = bk * p[i] + z[i];
				p_bar[i] = bk * p_bar[i] + z_bar[i];
			}
		}
		// Set bk_denominator = bk_numerator
		bkden = bknum;

//		std::cout << "Printing bk denominator" << std::endl;
//		std::cout << "bk_den: " << bkden << std::endl;

		// Set z = A * p
		A_times_x(state, p, z);

//		std::cout << "Printing z" << std::endl;
//		for (int i = 0; i < m; i++)
//		{
//			std::cout << z[i] << " ";
//		}
//		std::cout << std::endl;

		akden = 0.0;
		for (int i = 0; i < m; i++)
		{
			akden += z[i] * p_bar[i];
		}
//		std::cout << "Printing ak_denominator" << std::endl;
//		std::cout << "akden: " << akden << std::endl;
		// Calculate ak
		ak = bknum / akden;
//		std::cout << "Printing ak" << std::endl;
//		std::cout << "ak: " << ak << std::endl;
//
//		std::cout << "Printing p bar" << std::endl;
//		for (int i = 0; i < m; i++)
//		{
//			std::cout << p_bar[i] << " ";
//		}
//		std::cout << std::endl;

		// Set z_bar = A_T * p_bar

		A_times_x(state, p_bar, z_bar);
//		std::cout << "Printing z bar" << std::endl;
//		for (int i = 0; i < m; i++)
//		{
//			std::cout << z_bar[i] << " ";
//		}
//		std::cout << std::endl;

		// Calculate improved minimizer x_k+1 = x_k + ak * pk and new residuals r and r_bar
		// Eq 2.7.32 in numerical_recipes_in_c
		for (int i = 0; i < m; i++)
		{
			x[i] += ak * p[i];
			r[i] -= ak * z[i];
			r_bar[i] -= ak * z_bar[i];
		}

		// Set z = r
		solve_for_x(m, identity, r, z);

		// Calculate error depending on protocol chosen
		// Checks if |Ax - b| / |b| < tol
		if (itol == 1)
		{
			*err = calculate_relevant_norm(r, m, itol) / b_norm;
		}
		// Other itols have different checks. See page 86 in numerical_recipes_in_c
		else if (itol == 2)
		{
			*err = calculate_relevant_norm(z, m, itol) / b_norm;
		}
		else if (itol == 3 || itol == 4)
		{
			zm1_norm = z_norm;
			z_norm = calculate_relevant_norm(z, m, itol);
			if (fabs(zm1_norm - z_norm) > EPS * z_norm)
			{
				dx_norm = fabs(ak) * calculate_relevant_norm(p, m, itol);
				*err = z_norm / fabs(zm1_norm - z_norm) * dx_norm;
			}
			// Error may not be accurate, loop again
			else
			{
				*err = z_norm / b_norm;
				continue;
			}
			x_norm = calculate_relevant_norm(x, m, itol);
			if (*err <= 0.5 * x_norm)
			{
				*err /= x_norm;
			}
			// Error may not be accurate, loop again
			else
			{
				*err = z_norm / b_norm;
				continue;
			}
		}

		// Print out relevant info
		//std::cout << "iter_count: " << *iter << " error: " << *err << std::endl;

		// If the error is less than tol, then break
		if (*err <= tol)
		{
			//std::cout << "Reached error margin" << std::endl;
			break;
		}
	}
}







