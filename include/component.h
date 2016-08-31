struct ComponentSpace {
	ComponentSpace(unsigned nelements);

	unsigned numberOfElements();
	unsigned numberOfComponents();
	bool isBoss(unsigned element);
	unsigned find(unsigned lelement, bool compresspath = true);
	bool unify(unsigned one, unsigned two);
	void print1x1();
	void print();
	void copy(ComponentSpace &two);
	void dump_to_file(const char *F);
	void allocate();
	void deallocate();
	void init();
//	unsigned numberOfComponentsHost();

	unsigned nelements;
	unsigned ncomponents,	// number of components.
			 *complen, 		// lengths of components.
			 *ele2comp;		// components of elements.
};

ComponentSpace::ComponentSpace(unsigned nelements) {
	this->nelements = nelements;
	allocate();
	init();
}

void ComponentSpace::dump_to_file(const char *F) {
	static FILE *f;
	static unsigned *mem;

	if(!f) {
		f = fopen(F, "w");
	}
	int i;
	for(i = 0; i < nelements; i++) {
		fprintf(f, "%d %d\n", i, ele2comp[i]);
	}
	fprintf(f, "\n");
}

void ComponentSpace::print1x1() {
	printf("\t\t-----------------\n");
	for (unsigned ii = 0; ii < nelements; ++ii) {
		printf("\t\t%d -> %d\n", ii, ele2comp[ii]);
	}	
	printf("\t\t-----------------\n");
}

void print1x1(ComponentSpace cs) {
	cs.print1x1();
}

void ComponentSpace::print() {
}

unsigned ComponentSpace::numberOfElements() {
	return nelements;
}

unsigned ComponentSpace::numberOfComponents() {
	return ncomponents;
}

void ComponentSpace::allocate() {
	complen = (unsigned *)malloc(nelements * sizeof(unsigned));
	ele2comp = (unsigned *)malloc(nelements * sizeof(unsigned));
}

void ComponentSpace::deallocate() {
	free(complen);
	free(ele2comp);
}

void initcs(unsigned nelements, unsigned *complen, unsigned *ele2comp) {
	for (int id = 0; id < nelements; id++) {
		complen[id]	= 1;
		ele2comp[id] = id;
	}
}

void ComponentSpace::init() {
	// init the elements.
	initcs(nelements, complen, ele2comp);
	// init number of components.
	ncomponents = nelements;
}

bool ComponentSpace::isBoss(unsigned element) {
	return ele2comp[element] == element;
}

unsigned ComponentSpace::find(unsigned lelement, bool compresspath/*= true*/) {
	// do we need to worry about concurrency in this function?
	// for other finds, no synchronization necessary as the data-structure is a tree.
	// for other unifys, synchornization is not required considering that unify is going to affect only bosses, while find is going to affect only non-bosses.
	unsigned element = lelement;
	while (isBoss(element) == false) {
		element = ele2comp[element];
	}
	if (compresspath) ele2comp[lelement] = element;	// path compression.
	return element;
}

bool ComponentSpace::unify(unsigned one, unsigned two) {
	// if the client makes sure that one component is going to get unified as a source with another destination only once, then synchronization is unnecessary.
	// while this is true for MST, due to load-balancing in if-block below, a node may be source multiple times.
	// if a component is source in one thread and destination is another, then it is okay for MST.
	do {
		if(!isBoss(one)) return false;
		if(!isBoss(two)) return false;

		unsigned onecomp = one;
		unsigned twocomp = two;

		if (onecomp == twocomp) return false; // "duplicate" edges due to symmetry

		unsigned boss = twocomp;
		unsigned subordinate = onecomp;
		if (boss < subordinate) {			// break cycles by id.
			boss = onecomp;
			subordinate = twocomp;
		}
		// merge subordinate into the boss.
		unsigned oldboss = my_compare_swap(&ele2comp[subordinate],subordinate,boss);
		//unsigned oldboss = atomicCAS(&ele2comp[subordinate], subordinate, boss);
		if (oldboss != subordinate) {	// someone else updated the boss.
			// we need not restore the ele2comp[subordinate], as union-find 
			// ensures correctness and complen of subordinate doesn't matter.
			one = oldboss;
			two = boss;
			return false;
		} else {
//			printf("\t\tunifying %d -> %d (%d)\n", subordinate, boss);
//#ifdef ENABLE_OPENMP
			//atomicAdd(&complen[boss], complen[subordinate]);
//			my_fetch_add<unsigned>(&complen[boss], complen[subordinate]);
//#else
			complen[boss] += complen[subordinate];
//#endif
			// complen[subordinate] doesn't matter now, since find() will find its boss.

//#ifdef ENABLE_OPENMP
			//unsigned ncomp = atomicSub(ncomponents, 1);
//			my_fetch_sub<unsigned>(&ncomponents, 1);
//#else
			// a component has reduced.
			ncomponents --;
//#endif
			return true;
		}
	} while (true);
}
