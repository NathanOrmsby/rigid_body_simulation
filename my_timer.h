/*
 * my_timer.h
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */

#ifndef MY_TIMER_H_
#define MY_TIMER_H_

#include <windows.h>
#include <stdint.h>

typedef struct
{
	LARGE_INTEGER start, end;
} my_timer_t;

static LARGE_INTEGER ticks_per_second = {.QuadPart = 1};
static int64_t ticks_per_microsecond = 1;

// You must call my_timer_init() to use it. otherwise tick rate is set to 1 which is bad.
static inline void my_timer_init()
{
	QueryPerformanceFrequency(&ticks_per_second);
	ticks_per_microsecond = ticks_per_second.QuadPart / 1000000;
}

static inline my_timer_t my_timer_start()
{
	my_timer_t timer;
	QueryPerformanceCounter(&timer.start);
	return timer;
}

// Returns time in microseconds
static inline int64_t my_timer_end(my_timer_t *timer)
{
	QueryPerformanceCounter(&timer->end);
	return (timer->end.QuadPart - timer->start.QuadPart) / ticks_per_microsecond;
}
#endif /* MY_TIMER_H_ */
