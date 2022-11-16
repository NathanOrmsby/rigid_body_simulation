/*
 * utils.h
 *
 *  Created on: Nov 7, 2022
 *      Author: norms
 */
#include "matrix_stuff.h"
#ifndef UTILS_H_
#define UTILS_H_

int minimum(int a, int b);
int maximum(int a, int b);
void free_matrix_block_list(Matrix_Block *block_list, int num_blocks);
void print_matrix_block_list(Matrix_Block *block_list, int num_blocks);

#endif /* UTILS_H_ */
