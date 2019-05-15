/* timerc.h
 * this files contains 
 *				timer header functions IMP
 * 	that can be compiled with nvcc without -std=c++11 flag
 */

		/* function prototype */

#ifndef TIMERC_P_H
#define TIMERC_P_H
#include <stdio.h>

inline void cstart();
inline void cend(float * cputime);
inline void gstart();
inline void gend(float * gputime);

#endif

#ifndef TIMERC_H
#define TIMERC_H

		/* function implementation */
//--------------------------------------------------------
//	CPU
//--------------------------------------------------------
clock_t cpu_start, cpu_end;

//cstart
inline void cstart() { cpu_start = clock(); }

//cend
inline void cend(float * cputime) 
{ 
	cpu_end = clock(); 
	*cputime = (float) (cpu_end - cpu_start) / CLOCKS_PER_SEC * 1000.0;
}

//------------------------------------------------------------------------
//	GPU
//------------------------------------------------------------------------


static void _init();
static cudaEvent_t gpu_start, gpu_end;

//gstart
inline void gstart(){
	_init();
	cudaEventRecord(gpu_start, 0);
}

//gend
inline void gend(float * gputime){
	cudaEventRecord(gpu_end, 0);
	cudaEventSynchronize(gpu_end);
	cudaEventElapsedTime(gputime, gpu_start, gpu_end);
}

static int init = 0;
static void _init()
{
	if(!init){
		cudaEventCreate(&gpu_start);
		cudaEventCreate(&gpu_end);
		init = 1;
	}
}

#endif
