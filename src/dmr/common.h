#include <sys/time.h>
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#define GCC_EXTENSION
//#define OPENMP_3_1

double rtclock() {
	struct timezone Tzp;
	struct timeval Tp;
	int stat;
	stat = gettimeofday (&Tp, &Tzp);
	if (stat != 0) printf("Error return from gettimeofday: %d",stat);
	return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

template <class T>
inline T my_fetch_add(T *ptr, T val) {
#ifdef ENABLE_OPENMP
#ifdef GCC_EXTENSION
	return __sync_fetch_and_add(ptr,val);
#endif
#ifdef OPENMP_3_1
	T old;
	#pragma omp atomic capture
	{old = *ptr; *ptr += val;}
	return old;
#endif
#else
	T old; old = *ptr; *ptr += val;
	return old;
#endif
}

template <class T>
inline T my_fetch_sub(T *ptr, T val) {
#ifdef ENABLE_OPENMP
#ifdef GCC_EXTENSION
	return __sync_fetch_and_sub(ptr,val);
#endif
#ifdef OPENMP_3_1
	T old;
	#pragma omp atomic capture
	{old = *ptr; *ptr -= val;}
	return old;
#endif
#else
	T old; old = *ptr; *ptr -= val;
	return old;
#endif
}
;

template <class T>
inline T my_compare_swap(T *ptr, T old_val, T new_val) {
#ifdef ENABLE_OPENMP
#ifdef GCC_EXTENSION
	return __sync_val_compare_and_swap(ptr,old_val,new_val);
#endif
#ifdef OPENMP_3_1
	T old = *ptr;
	#pragma omp critical
	{
	if(*ptr == old_val) {
		*ptr = new_val;
	}
	}
	return old;
#endif
#else
	T old = *ptr;
	if(*ptr == old_val) *ptr = new_val;
	return old;
#endif
}
;

#ifdef ENABLE_OPENMP
void __syncthreads() {
#ifdef GCC_EXTENSION
	#pragma omp barrier
	//__sync_synchronize();
#else
	#pragma omp barrier
#endif
}
#endif
