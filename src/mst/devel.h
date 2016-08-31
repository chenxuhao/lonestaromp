#pragma once

/* Development routines */

void print_comp_mins(ComponentSpace cs, Graph graph, foru *minwtcomponent, unsigned *goaheadnodeofcomponent, unsigned *partners, bool *pin)
{
	foru *cminwt;
	unsigned *cgah, *cpart;
	unsigned *ele2comp;
	bool *cpin;

	ele2comp = (unsigned *) calloc(cs.nelements, sizeof(unsigned));
	cgah = (unsigned *) calloc(cs.nelements, sizeof(unsigned));
	cpart = (unsigned *) calloc(cs.nelements, sizeof(unsigned));
	cminwt = (foru *) calloc(cs.nelements, sizeof(unsigned));
	cpin = (bool *) calloc(cs.nelements, sizeof(bool));

	for(int i = 0; i < cs.nelements; i++)
	{
		if(ele2comp[i] == i && cminwt[i] != MYINFINITY && cpin[cgah[i]])
			printf("CM %d %d %d %d\n",  i, cminwt[i], cgah[i], cpart[i]);
	}

	free(ele2comp);
	free(cgah);
	free(cminwt);
}

