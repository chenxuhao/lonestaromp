#pragma once

struct Worklist {
	int *wl;
	int length;
	int index;

	Worklist(size_t nsize) {
		index = 0;
		length = (int) nsize; 
		wl = (int *) calloc(nsize, sizeof(int));
	}
	~Worklist() { }

	void display_items() {
		printf("WL: ");
		for(int i = 0; i < length; i++)
			printf("%d %d, ", i, wl[i]);
		printf("\n");
		return;
	}

	void reset() {
		index = 0;
	}

	int nitems(){
		return index;
	}

	int push(int item) {
		int old_index = my_fetch_add<int>(&index, 1);
		if(old_index >= length)
			return 0;
		wl[old_index] = item;
		return 1;
	}

	int pop(int &item) {
		int old_index = my_fetch_sub<int>(&index, 1);
		if(old_index <= 0) {
			index = 0;
			return 0;
		}
		item = wl[old_index - 1];
		return 1;
	}
};

struct Worklist2: public Worklist {
	Worklist2(int nsize) : Worklist(nsize) {}
	template <typename T>
	int push_1item(int nitem, int item) {
		typename T::TempStorage temp_storage;
		int queue_index;
		int total_items = 0;

		T(temp_storage).ExclusiveSum(nitem, nitem, total_items);
		//printf("t: %d\n", total_items);
		queue_index = my_fetch_add<int>(&index, total_items);
		//printf("queueindex: %d %d %d\n", queue_index, nitem+n_items, total_items);
		if(nitem == 1) {
			if(queue_index + nitem >= length) {
				printf("exceeded length: %d %d %d\n", queue_index, nitem, length);
				return 0;
			}
			wl[queue_index + nitem] = item;
		}
		return total_items;
	}

	template <typename T>
	int push_nitems(int n_items, int *items, int threads_per_block) {
		typename T::TempStorage temp_storage;
		int queue_index;
		int total_items;
		T(temp_storage).ExclusiveSum(n_items, n_items, total_items);
		queue_index = my_fetch_add<int>(&index, total_items);
		//printf("queueindex: %d %d %d\n", queue_index, nitem+n_items, total_items);
		for(int i = 0; i < n_items; i++) {
			//printf("pushing %d to %d\n", items[i], queue_index + n_item + i);
			if(queue_index + n_items + i >= length) {
				printf("exceeded length: %d %d %d %d\n", queue_index, n_items, i, length);
				return 0;
			}
			wl[queue_index + n_items + i] = items[i];
		}
		return total_items;
	}

	int pop_id(int id, int &item) {
		if(id < index) {
			item = wl[id];
			return 1;
		}
		return 0;
	} 
};
