#ifndef _T_CHEMISTRY_H_

#include "globals.h"
#include <stdint.h>

struct t_pcr{
	int16_t A,B;
	int8_t wo_A,wo_B;
};

class t_chemistry{
public:
	int32_t compound[highest_compound+1];
	unsigned char compound_content[highest_compound+1][4];
	//double n[mx+1][my+1][highest_compound+1];
	double xp[highest_compound+1];
	double yp[highest_compound+1];
	double np[highest_compound+1];
	//short q_aff[0x1000000];
	t_pcr *q_pcr;
	t_chemistry();
	~t_chemistry();
	void setup(void);
	int is_valid_compound(int32_t to_test);
	double affinity(int32_t a, int32_t b);
	unsigned char contains_site(int32_t A, int32_t B);
	unsigned char split(int32_t O1,int32_t O2,int32_t O3,int32_t O4,int32_t &O5,int32_t &O6);
	void show_str_compound(int32_t i);
	unsigned char delta_e(int32_t O1,int32_t O2,int32_t O3,int32_t O4);
	unsigned char length(int32_t O);
	double c(double x, double y,int32_t comp);
	void show_concentrations(void);
	unsigned char quick_split(int A, int B, int wo_A, int wo_B, int &AP, int &BP);
	int inttolinked(int v);
};

#define _T_CHEMISTRY_H_
#endif