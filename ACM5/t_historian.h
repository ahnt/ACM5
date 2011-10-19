#ifndef _T_HISTORIAN_H_

class t_history{
public:
	t_history *ancestor;
	char was;
	__int64 when;
	int wo;
	unsigned char *vorher;
	int size;
	int n_points;
};

class t_historian{
public:
	int main_loop_counter;
	t_historian();
	void born_kloned(t_history *&ancestor,t_history *&offspring);
	void born_rand(t_history *&h);
	void born_loaded(t_history *&h,int ident);
	void add_event(t_history *&h,char was, int wo, unsigned char *vorher, int sizeofvorher);
	void remove_h(t_history *&h);
	void save_history(char *name,t_history *h);
};

#define _T_HISTORIAN_H_
#endif
