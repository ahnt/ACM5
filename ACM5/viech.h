#ifndef _T_VIECH_H_

#include "chemistry.h"
#include <stdint.h>
#include <vector>
#include <string>

using namespace std;

struct t_gene{
	double n;
	unsigned char ident;
	unsigned char direction;
	int domain[5];
	int A,B,wo_A,wo_B,AP,BP;
	double affinity[5];
	signed char delta_e;
};

class t_viech{
public:
    t_viech *ancestor;
    int nrPointingAtMe;
    int born;
    string changeLog;
	int hmvg[16];
	int speed_hits;
	int32_t age;
	unsigned char heading;
	double speed;
	double energy;
	double xpos,ypos;
	double total_n;
	double n_comp[608];
	int n_comp_hits[608];
	t_gene *genes;
	double *activity;
	int32_t hmg;
	long nr;
	vector<vector<unsigned char> > chromosome;
	t_viech *next,*last;
	t_viech();
	~t_viech();
	void show(t_chemistry *chemistry);
	void fill_chr_rand(int32_t length_A,int32_t length_B);
	void inherit(t_viech *from,double mutationRate,double duplicationRate,double deletionRate);
	void setup_proteom(t_chemistry *chemistry);
	void setup_total_n(double a);
	int32_t how_many_genes(void);
	int32_t how_many_working_genes(void);
	void step(t_chemistry *chemistry);
	void calc_nucleotide_needs(void);
	void save_proteom_network(int name,char *path, t_chemistry *chemistry);
    void showProteomNetwork(t_chemistry *chemistry);
	void save_genome(int name);
	void load_genome(char *filename);
	int ready_to_divide(double level,double *fitness_for_compound);
	double get_total_n(void);
    void saveLOD(FILE *lod,FILE *comment,FILE *data);
    t_viech* getLMRCA(void);
};
#define _T_VIECH_H_
#endif