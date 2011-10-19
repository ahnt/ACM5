#ifndef _T_VIECH_H_

#include "chemistry.h"
#include <stdint.h>

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
	unsigned char *chr_A,*chr_B;
	int32_t length_chr_A,length_chr_B;
	t_viech *next,*last;
	t_viech();
	~t_viech();
	void show(t_chemistry *chemistry);
	void fill_chr_rand(int32_t length_A,int32_t length_B);
	void fill_chr_klone_of(t_viech *from);
	void make_perfekt_klone_of(t_viech *ancestor);
	void mutate(int32_t how_many);
	void insertion(int32_t how_many);
	void deletion(int32_t how_many);
	void duplication(int32_t how_many);
	void setup_proteom(t_chemistry *chemistry);
	void setup_total_n(double a);
	int32_t how_many_genes(void);
	int32_t how_many_working_genes(void);
	void step(t_chemistry *chemistry);
	void calc_nucleotide_needs(void);
	void save_proteom_network(int name,char *path, t_chemistry *chemistry);
	void save_genome(int name);
	void load_genome(int name,char *path);
	int ready_to_divide(double level,double *fitness_for_compound);
	double get_total_n(void);
};
#define _T_VIECH_H_
#endif