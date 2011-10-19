#ifndef _T_WORLD_H_

#include "globals.h"
#include "viech.h"
#include "chemistry.h"
#include <vector>

using namespace std;

class t_world{
public:
	int inhibitor;
	double max_fitness;
	double ave_fitness;
	t_chemistry *chemistry;
    vector<t_viech*> population;
    vector<int> order;
    int currentTime;
	double fittness_for_compound[608];
	double compound_avail[608];
	t_world();
	~t_world();
	double limits;
	int32_t _viecher_n;
	double max_energy;
	unsigned char feld[mx+1][my+1];
	void setup(void);
	void load_setup(int hm,int number, char *path);
	void run_all_viecher(void);
	void test_run_viecher(int32_t n,int iterations);
	void test_run_activity_viecher(int32_t n,int iterations,int inactive_gene,double test_activity, int overexpressed_gene, double overexpression);
	void depopulate(void);
};

#define _T_WORLD_H_
#endif