/**
 * \brief Sorts key-value pairs into ascending order.
 *
 * \par
 * - The sorting operation requires a pair of key buffers and a pair of value
 *   buffers.  Each pair is wrapped in a DoubleBuffer structure whose member
 *   DoubleBuffer::Current() references the active buffer.  The currently-active
 *   buffer may be changed by the sorting operation.
 * \par
 * The code snippet below illustrates the sorting of a device vector of \p int keys
 * with associated vector of \p int values.
 * \par
 * \code
 * // Declare, allocate, and initialize device pointers for sorting data
 * int  num_items;          // e.g., 7
 * int  *d_key_buf;         // e.g., [8, 6, 7, 5, 3, 0, 9]
 * int  *d_key_alt_buf;     // e.g., [        ...        ]
 * int  *d_value_buf;       // e.g., [0, 1, 2, 3, 4, 5, 6]
 * int  *d_value_alt_buf;   // e.g., [        ...        ]
 * ...
 * // Create a set of DoubleBuffers to wrap pairs of device pointers
 * DoubleBuffer<int> d_keys(d_key_buf, d_key_alt_buf);
 * DoubleBuffer<int> d_values(d_value_buf, d_value_alt_buf);
 * // Determine temporary device storage requirements
 * SortPairs(d_temp_storage, temp_storage_bytes, d_keys, num_items);
 * // d_keys.Current()      <-- [0, 3, 5, 6, 7, 8, 9]
 * // d_values.Current()    <-- [5, 4, 3, 1, 2, 0, 6]
 * \endcode
 *
 * \tparam Key      <b>[inferred]</b> Key type
 * \tparam Value    <b>[inferred]</b> Value type
 *
 * Author: Xuhao Chen <cxh@illinois.edu>
 */

#include <map>
#include "db.h"
template <typename Key, typename Value>
void sortPairs(DoubleBuffer<Key>keys, DoubleBuffer<Value>values, int num_items) {
//	printf("Sorting %d elements\n", num_items);
	int i;
	std::multimap<Key,Value> my_map;
	typename std::multimap<Key,Value>::iterator iter;
	for(i=0; i<num_items; i++) {
		my_map.insert(std::pair<Key,Value>(keys.Current()[i], values.Current()[i]));
	}
	i = 0;
	for(iter=my_map.begin(); iter!=my_map.end(); ++iter) {
		keys.Current()[i] = iter->first;
		values.Current()[i] = iter->second;
		i++;
	}
	return;
}
