//Definition of DoubleBuffer
//Author: Xuhao Chen <cxh@illinois.edu>
#ifndef DB_H
#define DB_H
template <typename T>
struct DoubleBuffer {
	// Pair of device buffer pointers
	T *d_buffers[2];

	///  Selector into \p d_buffers (i.e., the active/valid buffer)
	int selector;

	/// \brief Constructor
	DoubleBuffer() {
		selector = 0;
		d_buffers[0] = NULL;
		d_buffers[1] = NULL;
	}
	
	/// \brief Constructor
	DoubleBuffer(
			T *d_current,         ///< The currently valid buffer
			T *d_alternate)       ///< Alternate storage buffer of the same size as \p d_current
	{
		selector = 0;
		d_buffers[0] = d_current;
		d_buffers[1] = d_alternate;
	}

	/// \brief Return pointer to the currently valid buffer
	T* Current() { return d_buffers[selector]; }
};
#endif
