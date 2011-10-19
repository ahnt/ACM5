#ifndef _T_GLOBALS_H_
       
#define gene_code 0
#define max_nr_of_viecher 1000 //5000 //20000
#define size_of_clean_sweep 1000
#define highest_compound 608
#define feeding_level 53
#define max_pcr 5020279
#define start_chr_length_A 3000 //3000 //1500  //400
#define start_chr_length_B 3000 //3000	//1500 //400
#define mx 31
#define my 31
#define viecher_per_spot 16
//#define start_viecher mx*mx*(viecher_per_spot/2)
#define start_viecher 1000
#define start_compounds 100.0//0.10//0.0003//0.00027//0.000275//50//0.5//5//0.011
//#define synth_rate_enhancement 10000000
//#define synth_rate_enhancement 1000000
//#define synth_rate_enhancement 10000000000
// #define synth_rate_enhancement 1000000000

#define speed_enhancement   10000000000//100000000
#define turning_enhancement 100000000000
#define protein_synth_rate_enhancement 1
#define max_age 200

#define n_test_loops 30
#define repression_rate 0.0
#define over_expression_rate 2.0

#define sudden_death_val 63
#define limit 0.0000001


const unsigned char id_active[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
const unsigned char id_what[16]  ={0,1,2,2,0,1,2,2,3,4,3,4,0,1,2,2}; 

#define _T_GLOBALS_H_
#endif
