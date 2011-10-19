#include "t_historian.h"
#include "t_globals.h"
#include <stdlib.h>
#include <stdio.h>

#define historian_on

t_historian::t_historian()
{
	this->main_loop_counter=0;
}


void t_historian::born_rand(t_history *&h)
{
#ifdef historian_on
	if (h==NULL)
	{
		h=new t_history;
		h->ancestor=NULL;
		h->was='R';
		h->when=main_loop_counter;
		h->wo=255;
		h->n_points=1; //1
		h->vorher=NULL;
		h->size=0;
	}
#endif
}

void t_historian::born_loaded(t_history *&h,int ident)
{
#ifdef historian_on	
	if (h==NULL)
	{
		h=new t_history;
		h->ancestor=NULL;
		h->was='L';
		h->when=main_loop_counter;
		h->wo=ident;
		h->n_points=1; //1
		h->vorher=NULL;
		h->size=0;
	}
#endif
}

void t_historian::add_event(t_history *&h, char was, int wo, unsigned char *vorher, int sizeofvorher)
{
#ifdef historian_on
	t_history *nh;
	//h->n_points++;
	nh=new t_history;
	nh->ancestor=h;
	nh->was=was;
	nh->when=main_loop_counter;
	nh->wo=wo;
	//nh->vorher=new unsigned char[sizeofvorher];
	nh->vorher=vorher;
	nh->size=sizeofvorher;
	nh->n_points=1; //1
	h=nh;
#else
	delete vorher;
	vorher=NULL;
#endif
}

void t_historian::remove_h(t_history *&h)
{
#ifdef historian_on
	h->n_points--;
	if (h->n_points<1)
	{
		if(h->ancestor!=NULL)
			remove_h(h->ancestor);
		//if (h->vorher!=NULL)
		delete h->vorher;
		delete h;
		h=NULL;
	}
#endif
}

void t_historian::save_history(char *name, t_history *h)
{
	#ifdef historian_on
	FILE *f;
	t_history *r;
	int z;
	int i;
	f=fopen(name,"w+t");
	r=h;
	while (r!=NULL)
	{
		fprintf(f,"%c	%i	%i	%i	",r->was,(int)r->when,(int)r->wo,(int)r->size);
		if((r->size!=0)&&(r->vorher!=NULL))
			for(z=0;z<r->size;z++) 
				fprintf(f,"%i	",r->vorher[z]);
		fprintf(f,"\n");
		r=r->ancestor;
	}
	fclose(f);
#endif
}

void t_historian::born_kloned(t_history *&ancestor, t_history *&offspring)
{
#ifdef historian_on
	offspring=ancestor;
	ancestor->n_points++;
#endif
}