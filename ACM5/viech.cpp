#include "viech.h"
#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

//int synth_rate_enhancement=1000000000;
double MOL=10000000000000000000;

t_viech::t_viech()
{
	this->chr_A=NULL;
	this->chr_B=NULL;
	this->genes=NULL;
	this->activity=NULL;
    this->nrPointingAtMe=1;
    this->ancestor=NULL;
    this->energy=0.0;
}

t_viech::~t_viech()
{
	if(chr_A!=NULL)
	{
		delete chr_A;
		chr_A=NULL;
	}
	if(chr_B!=NULL)
	{
		delete chr_B;
		chr_B=NULL;
	}
	if(genes!=NULL)
	{
		delete genes;
		genes=NULL;
	}
	if (activity!=NULL)
	{
		delete activity;
		activity=NULL;
	}
    if(ancestor!=NULL){
        ancestor->nrPointingAtMe--;
        if(ancestor->nrPointingAtMe==0)
            delete ancestor;
    }
}

void t_viech::show(t_chemistry *chemistry)
{
	int z;
	for(z=0;z<608;z++)
		if(this->n_comp[z]!=0) 
		{
			printf("%f	",this->n_comp[z]);
			chemistry->show_str_compound(chemistry->compound[z]);
			printf("\n");
		}
}

void t_viech::fill_chr_rand(int length_A, int length_B)
{
	int32_t z;
	this->length_chr_A=length_A;
	this->length_chr_B=length_B;
	chr_A=new unsigned char[this->length_chr_A];
	chr_B=new unsigned char[this->length_chr_B];
	for(z=0;z<this->length_chr_A;z++)
		this->chr_A[z]=rand()&255;
	for(z=0;z<this->length_chr_B;z++)
		this->chr_B[z]=rand()&255;
}
int32_t t_viech::how_many_genes(void)
{
	int z;
	int32_t i=0;
	for(z=0;z<length_chr_A;z++)
		if ((chr_A[z])<=gene_code) //0 
		{
			i++;
			//z=(z+16)&262128;
		}
	for(z=0;z<length_chr_B;z++)
		if ((chr_B[z])<=gene_code) //0
		{
			i++;
			//z=(z+16)&262128;
		}
	//return i;
		return i;
}
int32_t t_viech::how_many_working_genes(void)
{
	int z;
	return this->hmg;
}

void t_viech::setup_proteom(t_chemistry *chemistry)
{
	int z;
	int g;
	int i;
	int RAP,RBP;
	this->hmg=this->how_many_genes();
	this->genes=new t_gene[hmg];
	this->age=0;
	for(z=0;z<608;z++)
	{
		this->n_comp_hits[z]=0;
		this->n_comp[z]=(double)0.0;
	}
	g=0;
	for(z=0;z<this->length_chr_A;z++)
		if(this->chr_A[z]<=gene_code)
		{
			for(i=0;i<5;i++)
			genes[g].domain[i]=	((int)this->chr_A[(z+3+(i*3))%length_chr_A]<<16)+
								((int)this->chr_A[(z+4+(i*3))%length_chr_A]<<8)+
								((int)this->chr_A[(z+5+(i*3))%length_chr_A]);
			genes[g].n=((double)chr_A[(z+1)%length_chr_A]+(double) 1.0)/(double)256000.0;
			genes[g].ident=id_what[chr_A[(z+2)%length_chr_A]&15];
			g++;
		}
	for(z=0;z<this->length_chr_B;z++)
		if(this->chr_B[z]<=gene_code)
		{
			for(i=0;i<5;i++)
			genes[g].domain[i]=	((int)this->chr_B[(z+3+(i*3))%length_chr_B]<<16)+
								((int)this->chr_B[(z+4+(i*3))%length_chr_B]<<8)+
								((int)this->chr_B[(z+5+(i*3))%length_chr_B]);
			genes[g].n=((double)chr_B[(z+1)%length_chr_B]+(double) 1.0)/(double)256000.0;
			genes[g].ident=id_what[chr_B[(z+2)%length_chr_B]&15];
			g++;
		}
	for(z=0;z<hmg;z++)
	{
		switch(this->genes[z].ident){
			case 0: //import protein
				genes[z].A=genes[z].domain[0]%608;
				this->genes[z].affinity[0]=(chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[1])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[2])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[3])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[4]))/(double) 4.0;
				break;
			case 1: //export protein
				genes[z].A=genes[z].domain[0]%608;
				this->genes[z].affinity[0]=(chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[1])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[2])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[3])+
											chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[4]))/(double) 4.0;
				this->n_comp_hits[genes[z].A]++;
				break;
			case 2: //mol reaktion
				i=genes[z].domain[0]%max_pcr;
				genes[z].A=chemistry->q_pcr[i].A;
				genes[z].B=chemistry->q_pcr[i].B;
				genes[z].wo_A=chemistry->q_pcr[i].wo_A;
				genes[z].wo_B=chemistry->q_pcr[i].wo_B;
				chemistry->quick_split(genes[z].A,genes[z].B,genes[z].wo_A,genes[z].wo_B,RAP,RBP);
				genes[z].AP=chemistry->inttolinked(RAP);
				genes[z].BP=chemistry->inttolinked(RBP);
				this->genes[z].affinity[0]=(chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[1])+
											chemistry->affinity(chemistry->compound[genes[z].B],genes[z].domain[2])+
											chemistry->affinity(chemistry->compound[genes[z].AP],genes[z].domain[3])+
											chemistry->affinity(chemistry->compound[genes[z].BP],genes[z].domain[4]))/(double) 4.0;
				this->n_comp_hits[genes[z].A]++;
				this->n_comp_hits[genes[z].B]++;
				break;
			case 3://flagellum
				i=genes[z].domain[0]%max_pcr;
				genes[z].A=chemistry->q_pcr[i].A;
				genes[z].B=chemistry->q_pcr[i].B;
				genes[z].wo_A=chemistry->q_pcr[i].wo_A;
				genes[z].wo_B=chemistry->q_pcr[i].wo_B;
				chemistry->quick_split(genes[z].A,genes[z].B,genes[z].wo_A,genes[z].wo_B,RAP,RBP);
				genes[z].AP=chemistry->inttolinked(RAP);
				genes[z].BP=chemistry->inttolinked(RBP);
				this->genes[z].affinity[0]=(chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[1])+
											chemistry->affinity(chemistry->compound[genes[z].B],genes[z].domain[2])+
											chemistry->affinity(chemistry->compound[genes[z].AP],genes[z].domain[3])+
											chemistry->affinity(chemistry->compound[genes[z].BP],genes[z].domain[4]))/(double) 4.0;
				this->n_comp_hits[genes[z].A]++;
				this->n_comp_hits[genes[z].B]++;
				break;
			case 4://cillium
				i=genes[z].domain[0]%max_pcr;
				genes[z].A=chemistry->q_pcr[i].A;
				genes[z].B=chemistry->q_pcr[i].B;
				genes[z].wo_A=chemistry->q_pcr[i].wo_A;
				genes[z].wo_B=chemistry->q_pcr[i].wo_B;
				chemistry->quick_split(genes[z].A,genes[z].B,genes[z].wo_A,genes[z].wo_B,RAP,RBP);
				genes[z].AP=chemistry->inttolinked(RAP);
				genes[z].BP=chemistry->inttolinked(RBP);
				this->genes[z].affinity[0]=(chemistry->affinity(chemistry->compound[genes[z].A],genes[z].domain[1])+
											chemistry->affinity(chemistry->compound[genes[z].B],genes[z].domain[2])+
											chemistry->affinity(chemistry->compound[genes[z].AP],genes[z].domain[3])+
											chemistry->affinity(chemistry->compound[genes[z].BP],genes[z].domain[4]))/(double) 4.0;
				this->n_comp_hits[genes[z].A]++;
				this->n_comp_hits[genes[z].B]++;
				break;
	}
	//	printf("%f\n",genes[z].n);
	}
	total_n=this->get_total_n();
	energy=(double)0.0;
	for(z=0;z<608;z++) if(this->n_comp_hits[z]==0) this->n_comp_hits[z]=1;
	heading=rand()&255;
}

void t_viech::setup_total_n(double a)
{
	int z;
	printf("function missing setuptotal_n\n");
}

void t_viech::step(t_chemistry *chemistry)
{
	double n_add[608],n_sub[608];
	int z,d;
	double alteration,l;
	double tumble_force=0.0;
	double tt;
	int speed_hits=1;
	this->age++;
	for(z=0;z<608;z++) 
	{
		n_add[z]=(double)0.0;
		n_sub[z]=(double)0.0;
	}
	speed=0;
	for(z=0;z<hmg;z++)
		switch(genes[z].ident){
			case 0: //import protein
				if(chemistry->np[genes[z].A]>0.0)
				{
					alteration=chemistry->c(this->xpos,this->ypos,this->genes[z].A)*(genes[z].n/total_n)*(genes[z].affinity[0]);
					n_add[genes[z].A]+=alteration;
					chemistry->np[genes[z].A]-=alteration;
					chemistry->np[608]-=alteration;
				}
				break;
			case 1: //export protein
				if(n_comp[genes[z].A]>0.0)
				{
					alteration=((n_comp[genes[z].A]/total_n)/(double)n_comp_hits[genes[z].A])*(genes[z].n/total_n)*(genes[z].affinity[0]);
					n_sub[genes[z].A]+=alteration;
					chemistry->np[genes[z].A]+=alteration;
					chemistry->np[608]+=alteration;
				}
				break;
			case 2: //molecular reaktion
				if((n_comp[genes[z].A]>0.0)&&(n_comp[genes[z].B]>0.0))
				{
					alteration=	((n_comp[genes[z].A]/total_n)/(double)n_comp_hits[genes[z].A])*
								((n_comp[genes[z].B]/total_n)/(double)n_comp_hits[genes[z].B])*
								(genes[z].n/total_n)*(genes[z].affinity[0]);
					n_sub[genes[z].A]+=alteration;
    				n_sub[genes[z].B]+=alteration;
					n_add[genes[z].AP]+=alteration;
					n_add[genes[z].BP]+=alteration;
				}
				break;
			case 3: //flagellum
				if((n_comp[genes[z].A]>0.0)&&(n_comp[genes[z].B]>0.0))
				{
					alteration=	((n_comp[genes[z].A]/total_n)/(double)n_comp_hits[genes[z].A])*
								((n_comp[genes[z].B]/total_n)/(double)n_comp_hits[genes[z].B])*
								(genes[z].n/total_n)*(genes[z].affinity[0]);
					n_sub[genes[z].A]+=alteration;
    				n_sub[genes[z].B]+=alteration;
					chemistry->np[genes[z].AP]+=alteration;
					chemistry->np[608]+=alteration;
					chemistry->np[genes[z].BP]+=alteration;
					chemistry->np[608]+=alteration;
					speed+=alteration;
					speed_hits+=n_comp_hits[genes[z].A]+n_comp_hits[genes[z].B];
				}
				break;
			case 4: //cillium
				if((n_comp[genes[z].A]>0.0)&&(n_comp[genes[z].B]>0.0))
				{
					alteration=	((n_comp[genes[z].A]/total_n)/(double)n_comp_hits[genes[z].A])*
								((n_comp[genes[z].B]/total_n)/(double)n_comp_hits[genes[z].B])*
								(genes[z].n/total_n)*(genes[z].affinity[0]);
					n_sub[genes[z].A]+=alteration;
    				n_sub[genes[z].B]+=alteration;
					chemistry->np[genes[z].AP]+=alteration;
					chemistry->np[608]+=alteration;
					chemistry->np[genes[z].BP]+=alteration;
					chemistry->np[608]+=alteration;
					tumble_force+=alteration;
				}
				break;
	}
	z=0;
	for(z=0;z<608;z++)
	{
		if((n_add[z]-n_sub[z])<(double)0.0)
			z=z;
		this->n_comp[z]=n_comp[z]+n_add[z]-n_sub[z];
		if(n_comp[z]<(double)0.0)
			n_comp[z]=(double) 0.0;
		//	printf("%f	",(float)n_comp[z]);
	}
	//printf("\n");
	total_n=this->get_total_n();
	tt=tumble_force+speed;
	if(tt!=(double)0.0)
	{
		if((rand()&1)==1)
		heading=(heading+(rand()&15)+(int)((rand()&255)*(tumble_force/tt)))&255;
		else
		heading=(heading-(rand()&15)+(int)((rand()&255)*(tumble_force/tt)))&255;
		if(log(speed/tt)>=-(double)1.0) 
			speed=(double)1.0;
		else speed=(double)1.0/(-(double)1.0*log(speed/tt));
		xpos=xpos+sin(heading*M_PI/127)*speed;
		ypos=ypos+cos(heading*M_PI/127)*speed;
	}

}
void t_viech::calc_nucleotide_needs(void)
{
}

void t_viech::fill_chr_klone_of(t_viech *from)
{
	int z;
	if(chr_A!=NULL)
	{
		delete chr_A;
		chr_A=NULL;
	}
	if(chr_B!=NULL)
	{
		delete chr_B;
		chr_B=NULL;
	}
	this->length_chr_A=from->length_chr_A;
	this->length_chr_B=from->length_chr_B;
	this->chr_A=new unsigned char [this->length_chr_A];
	this->chr_B=new unsigned char [this->length_chr_B];
	for(z=0;z<this->length_chr_A;z++)
		this->chr_A[z]=from->chr_A[z];
	for(z=0;z<this->length_chr_B;z++)
		this->chr_B[z]=from->chr_B[z];
    ancestor=from;
    from->nrPointingAtMe++;
}

void t_viech::mutate(int32_t how_many)
{
	int z;
	int wo;
	unsigned char *old;
	for(z=0;z<how_many;z++)
		if((rand()&1)==0)
		{
			wo=rand()%length_chr_A;
			old=new unsigned char[1];
			old[0]=chr_A[wo];
			chr_A[wo]=rand()&255;
		}
		else
		{
			wo=rand()%length_chr_B;
			old=new unsigned char[1];
			old[0]=chr_B[wo];
			chr_B[wo]=rand()&255;
		}
	
}


void t_viech::insertion(int how_many)
{
	
}

void t_viech::deletion(int how_many)
{
	int z;
	int d_pos;
	unsigned char *reserve;
	unsigned char *old;
	if((rand()&1)==0)
	{
		if (this->length_chr_A>how_many)
		{
			d_pos=(rand()%(this->length_chr_A-how_many));
			reserve=new unsigned char[length_chr_A];
			for(z=0;z<this->length_chr_A;z++)
				reserve[z]=this->chr_A[z];
			delete chr_A; chr_A=NULL;
			chr_A=new unsigned char[this->length_chr_A-how_many];
			for(z=0;z<d_pos;z++)
				chr_A[z]=reserve[z];
			for(z=d_pos;z<length_chr_A-how_many;z++)
				chr_A[z]=reserve[z+how_many];
			length_chr_A-=how_many;
			
			old=new unsigned char[how_many];	
			for(z=d_pos;z<d_pos+how_many;z++)
				old[z-d_pos]=reserve[z];
			delete reserve;
		}
	}
	else
	{
		if (this->length_chr_B>how_many)
		{
			d_pos=(rand()%(this->length_chr_B-how_many));
			reserve=new unsigned char[length_chr_B];
			for(z=0;z<this->length_chr_B;z++)
				reserve[z]=this->chr_B[z];
			delete chr_B; chr_B=NULL;
			chr_B=new unsigned char[this->length_chr_B-how_many];
			for(z=0;z<d_pos;z++)
				chr_B[z]=reserve[z];
			for(z=d_pos;z<length_chr_B-how_many;z++)
				chr_B[z]=reserve[z+how_many];
			length_chr_B-=how_many;
			old=new unsigned char[how_many];	
			for(z=d_pos;z<d_pos+how_many;z++)
				old[z-d_pos]=reserve[z];
			delete reserve;
		}
	}
}

void t_viech::duplication(int32_t how_many)
{
	int d_pos;
	int z;
	unsigned char *reserve;
	
	if((rand()&1)==0)
	{
		if (this->length_chr_A>how_many)
		{
			d_pos=rand()%(this->length_chr_A-how_many);
			reserve=new unsigned char[this->length_chr_A];
			for(z=0;z<this->length_chr_A;z++)
				reserve[z]=this->chr_A[z];
			delete chr_A; chr_A=NULL;
			chr_A=new unsigned char[length_chr_A+how_many];
			for(z=0;z<d_pos+how_many;z++)
				chr_A[z]=reserve[z];
			for(z=d_pos+how_many;z<length_chr_A+how_many;z++)
				chr_A[z]=reserve[z-how_many];
			this->length_chr_A+=how_many;
			delete reserve;
		}
	}
	else
	{
		if (this->length_chr_B>how_many)
		{
			d_pos=rand()%(this->length_chr_B-how_many);
			reserve=new unsigned char[this->length_chr_B];
			for(z=0;z<this->length_chr_B;z++)
				reserve[z]=this->chr_B[z];
			delete chr_B; chr_B=NULL;
			chr_B=new unsigned char[length_chr_B+how_many];
			for(z=0;z<d_pos+how_many;z++)
				chr_B[z]=reserve[z];
			for(z=d_pos+how_many;z<length_chr_B+how_many;z++)
				chr_B[z]=reserve[z-how_many];
			this->length_chr_B+=how_many;
			delete reserve;
		}
	}
}

void t_viech::save_proteom_network(int name,char *path, t_chemistry *chemistry)
{
    printf("WARING PROTEOM NOT SAVED!\n");
    /*
	char n[200],n2[100];
	FILE *f;
	int z;
	int i,j;
	strcpy(n,path);
	_itoa_s(name,n2,100,10);
	strcat(n,n2);
	strcat(n,"_PNW.txt");
	f=fopen(n,"w+t");
	for(z=0;z<hmg;z++)
	{
		switch (genes[z].ident)
		{
			case 0:
				fprintf(f,"IMPORT:	");
				for(i=0;i<12;i++)
					fprintf(f,"%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				fprintf(f,"\n");
				break;
			case 1:
				fprintf(f,"EXPORT:	");
				for(i=0;i<12;i++)
					fprintf(f,"%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				fprintf(f,"\n");
				break;
			case 2:
				fprintf(f,"MOLREACT:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_B) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].B]>>(i*2))&3);

				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
					fprintf(f,"%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				fprintf(f,"	");
				for(i=0;i<12;i++)
					fprintf(f,"%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				fprintf(f,"\n");
				break;
			case 3:
				fprintf(f,"FLAGELLUM:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].B]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				}
				fprintf(f,"\n");
				break;
			case 4:
				fprintf(f,"CILIA:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].B]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				}
				fprintf(f,"	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) fprintf(f,"|");
					fprintf(f,"%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				}
				fprintf(f,"\n");
				break;
		}
	}
	for(z=0;z<hmg;z++)
		switch(genes[z].ident)
	{
		case 0:fprintf(f,"IMPORT:	%i\n",genes[z].A); break;
		case 1:fprintf(f,"EXPORT:	%i\n",genes[z].A); break;
		case 2:fprintf(f,"REACTION:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
		case 3:fprintf(f,"FLAGELLUM:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
		case 4:fprintf(f,"CILIA:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
	}

	fclose(f);
     */
}

void t_viech::save_genome(int name)
{
    printf("WARING NO save_gemome FUNCTION!\n");
/*	
	char n[100];
	FILE *f;
	int z;
	int i;
	_itoa_s(name,n,100,10);
	strcat(n,".txt");
	f=fopen(n,"w+t");
	fprintf(f,"%i	%i\n",this->length_chr_A,this->length_chr_B);
	for(z=0;z<this->length_chr_A;z++)
		fprintf(f,"%i\n",this->chr_A[z]);
	for(z=0;z<this->length_chr_B;z++)
		fprintf(f,"%i\n",this->chr_B[z]);
	fclose(f);
 */
}

void t_viech::load_genome(int name, char *path)
{
    printf("GENOME NOT LOADED!\n");
    /*
	char n[200],nu[200];
	FILE *f;
	int z;
	int i;
	unsigned char c;
	_itoa_s(name,nu,200,10);
	strcpy(n,path);
	strcat(n,nu);
	strcat(n,".txt");
	f=fopen(n,"r+t");
	fscanf(f,"%i	%i\n",&length_chr_A,&length_chr_B);
	chr_A=new unsigned char[length_chr_A];
	chr_B=new unsigned char[length_chr_B];
	for(z=0;z<this->length_chr_A;z++)
	{
		fscanf(f,"%i",&i);
		this->chr_A[z]=i&255;
	}
	for(z=0;z<this->length_chr_B;z++)
	{
		fscanf(f,"%i",&i);
		this->chr_B[z]=i&255;
	}
	fclose(f);
     */
}

void t_viech::make_perfekt_klone_of(t_viech *ancestor)
{
/*	int z,i,t;
	this->chr_A=ancestor->chr_A;
	this->chr_B=ancestor->chr_B;
	this->length_chr_A=ancestor->length_chr_A;
	this->length_chr_B=ancestor->length_chr_B;
	historian.born_kloned(ancestor->history,this->history);

	this->age=0;
	this->energy=0;
	this->heading=rand()&255;
	this->speed=0;
	for(z=0;z<5;z++)
		this->protein_synth[z]=(double) 1.0;
	if (genes!=NULL)
	{
		delete genes;
		genes=NULL;
	}
	hmg=ancestor->hmg;
	calc_nucleotide_needs();
	genes=new t_gene[hmg];
	activity=new double[hmg];
	for(z=0;z<608;z++)
		this->n_comp[z]=(double)0.0;
	for(z=0;z<hmg;z++)
		activity[z]=(double)1.0;

	for(i=0;i<hmg;i++)
	{
	}*/
	printf("important function missing\n");
}

int t_viech::ready_to_divide(double level,double *fitness_for_compound)
{
	double r;
	int z;
	r=(double)1.0;
	for(z=0;z<608;z++)
	{
		//r+=this->n_comp[z]*fitness_for_compound[z];
		if((n_comp[z]>(double)0.0)&&(fitness_for_compound[z]>(double)0.0)) 
			r*=(double)1.1+(this->n_comp[z]*fitness_for_compound[z]);
	}
	energy=r;
	if(r>level)
		return 1;
	else 
		return 0;
}


double t_viech::get_total_n(void)
{
	int z;
	double t;
	t=(double) 1.0;
	for(z=0;z<hmg;z++)
		t+=genes[z].n;
	for(z=0;z<608;z++)
		t+=n_comp[z];
	return t;
}