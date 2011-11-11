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
	this->chromosome.clear();
	this->genes=NULL;
	this->activity=NULL;
    this->nrPointingAtMe=1;
    this->ancestor=NULL;
    this->energy=0.0;
    changeLog.clear();
}

t_viech::~t_viech()
{
    chromosome.clear();
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
void t_viech::inherit(t_viech *from,double mutationRate,double duplicationRate,double deletionRate){
    from->nrPointingAtMe++;
    ancestor=from;
    chromosome=from->chromosome;
    for(int i=0;i<chromosome.size();i++)
        for(int j=0;j<chromosome[i].size();j++)
                if(((double)rand()/(double)RAND_MAX)<mutationRate){
                    chromosome[i][j]=rand()&255;//(chromosome[i][j]&(255-(3<<k)))+((rand()&3)<<k);
                    char n[100];
                    sprintf(n,"M %i %i\n",i,j);
                    changeLog.append(n);
                }
    if(((double)rand()/(double)RAND_MAX)<duplicationRate){
        int width,fromChr,toChr,startFromC,startToC;
        vector<unsigned char> thePiece;
        fromChr=(int)(rand()%chromosome.size());
        toChr=(int)(rand()%chromosome.size());
        if(chromosome[toChr].size()<20000){
            width=64+(rand()&511);
            startFromC=(int)(rand()%(chromosome[fromChr].size()-width));
            startToC=(int)(rand()%(chromosome[toChr].size()-width));
            thePiece.assign(chromosome[fromChr].begin()+startFromC,chromosome[fromChr].begin()+startFromC+width);
            chromosome[toChr].insert(chromosome[toChr].begin()+startToC,thePiece.begin(),thePiece.end());
            char n[1000];
            sprintf(n,"D %i %i %i %i %i\n",fromChr,toChr,width,startFromC,startToC);
            changeLog.append(n);
        }
    }
    if(((double)rand()/(double)RAND_MAX)<deletionRate){
        int width,from,chr;
        chr=rand()%(int)chromosome.size();
        if(chromosome[chr].size()>1000){
            width=rand()&255;
            from=1+(rand()%((int)chromosome[chr].size()-width));
            chromosome[chr].erase(chromosome[chr].begin(),chromosome[chr].begin()+width);
            char n[100];
            sprintf(n, "E %i %i %i\n",chr,from,width);
        }
    }
}


void t_viech::fill_chr_rand(int length_A, int length_B)
{
    chromosome.resize(2);
    chromosome[0].resize(length_A);
    chromosome[1].resize(length_B);
    for(int i=0;i<chromosome.size();i++)
        for(int j=0;j<chromosome[i].size();j++)
		this->chromosome[i][j]=(unsigned char)(rand()&255);
}
int32_t t_viech::how_many_genes(void)
{
	int32_t k=0;
    for(int i=0;i<chromosome.size();i++)
        for(int j=0;j<chromosome[i].size();j++)
            if (chromosome[i][j]<=gene_code)
                k++;
		return k;
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
    for(int i=0;i<chromosome.size();i++)
        for(int j=0;j<chromosome[i].size();j++)
            if(this->chromosome[i][j]<=gene_code)
            {
                for(int k=0;k<5;k++)
                    genes[g].domain[k]=	((int)this->chromosome[i][(j+3+(k*3))%this->chromosome[i].size()]<<16)+
                                        ((int)this->chromosome[i][(j+4+(k*3))%this->chromosome[i].size()]<<8)+
                                        ((int)this->chromosome[i][(j+5+(k*3))%this->chromosome[i].size()]);
                genes[g].n=((double)this->chromosome[i][(j+1)%this->chromosome[i].size()]+(double) 1.0)/(double)256000.0;
                genes[g].ident=id_what[chromosome[i][(j+2)%chromosome[i].size()]&15];
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


void t_viech::save_proteom_network(char *filename, t_chemistry *chemistry)
{
	char n[200],n2[100];
	FILE *f;
	int z;
	int i,j;
	f=fopen(filename,"w+t");
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
}

void t_viech::save_proteom_toDot(char *filename, t_chemistry *chemistry){
    char n[200],n2[100];
	FILE *f;
	int z;
	int i,j;
	f=fopen(filename,"w+t");
	fprintf(f,"digraph brain {\n");
	for(z=0;z<hmg;z++)
		switch(genes[z].ident)
	{
		case 0:
            fprintf(f,"out%i    ->  I%i;\n",genes[z].A,z); 
            fprintf(f,"I%i    ->  %i;\n",z,genes[z].A); 
            break;
		case 1:
            fprintf(f,"O%i   ->  out%i;\n",z,genes[z].A); 
            fprintf(f,"%i   ->  O%i;\n",genes[z].A,z); 
            break;
		case 2:
            fprintf(f,"%i       ->  G%i;\n",genes[z].A,z);
            fprintf(f,"%i       ->  G%i;\n",genes[z].B,z);
            fprintf(f,"G%i       ->  %i;\n",z,genes[z].AP);
            fprintf(f,"G%i       ->  %i;\n",z,genes[z].BP);
            break;
		case 3:
            fprintf(f,"%i       ->  F%i;\n",genes[z].A,z);
            fprintf(f,"%i       ->  F%i;\n",genes[z].B,z);
            fprintf(f,"F%i       ->  %i;\n",z,genes[z].AP);
            fprintf(f,"F%i       ->  %i;\n",z,genes[z].BP);
            break;
		case 4:
            fprintf(f,"%i       ->  C%i;\n",genes[z].A,z);
            fprintf(f,"%i       ->  C%i;\n",genes[z].B,z);
            fprintf(f,"C%i       ->  %i;\n",z,genes[z].AP);
            fprintf(f,"C%i       ->  %i;\n",z,genes[z].BP);
            break;
	}
    
	fprintf(f,"}\n");
	fclose(f);
}


void t_viech::showProteomNetwork(t_chemistry *chemistry){
	int z;
	int i,j;
	for(z=0;z<hmg;z++)
	{
		switch (genes[z].ident)
		{
			case 0:
				printf("IMPORT:	");
				for(i=0;i<12;i++)
					printf("%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				printf("\n");
				break;
			case 1:
				printf("EXPORT:	");
				for(i=0;i<12;i++)
					printf("%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				printf("\n");
				break;
			case 2:
				printf("MOLREACT:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_B) printf("|");
					printf("%i",(chemistry->compound[genes[z].B]>>(i*2))&3);
                    
				}
				printf("	");
				for(i=0;i<12;i++)
					printf("%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				printf("	");
				for(i=0;i<12;i++)
					printf("%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				printf("\n");
				break;
			case 3:
				printf("FLAGELLUM:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].B]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				}
				printf("\n");
				break;
			case 4:
				printf("CILIA:	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].A]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].B]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].AP]>>(i*2))&3);
				}
				printf("	");
				for(i=0;i<12;i++)
				{
					if(i==genes[z].wo_A) printf("|");
					printf("%i",(chemistry->compound[genes[z].BP]>>(i*2))&3);
				}
				printf("\n");
				break;
		}
	}
	for(z=0;z<hmg;z++)
		switch(genes[z].ident)
	{
		case 0:printf("IMPORT:	%i\n",genes[z].A); break;
		case 1:printf("EXPORT:	%i\n",genes[z].A); break;
		case 2:printf("REACTION:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
		case 3:printf("FLAGELLUM:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
		case 4:printf("CILIA:	%i	%i	%i	%i\n",genes[z].A,genes[z].B,genes[z].AP,genes[z].BP); break;
	}
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

void t_viech::load_genome(char *filename)
{
	FILE *f;
	int z;
	int i;
    int length_chr_A,length_chr_B;
	f=fopen(filename,"r+t");
	fscanf(f,"%i	%i\n",&length_chr_A,&length_chr_B);
    chromosome.resize(2);
    chromosome[0].resize(length_chr_A);
    chromosome[1].resize(length_chr_B);
	for(z=0;z<length_chr_A;z++)
	{
		fscanf(f,"%i",&i);
		this->chromosome[0][z]=i&255;
	}
	for(z=0;z<length_chr_B;z++)
	{
		fscanf(f,"%i",&i);
		this->chromosome[1][z]=i&255;
	}
	fclose(f);
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
			r*=(double)1.01+(this->n_comp[z]*fitness_for_compound[z]);
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

void t_viech::saveLOD(FILE *lod,FILE *comment,FILE *data){
    if(ancestor!=NULL)
        ancestor->saveLOD(lod,comment, data);
    fprintf(comment,">%i\n%s\n",born,changeLog.c_str());
    fprintf(data, "%i   %i\n",(int)chromosome[0].size(),(int)chromosome[1].size());
    for(int i=0;i<chromosome.size();i++){
        for(int j=0;j<chromosome[i].size();j++)
            fprintf(data,"%i\n",chromosome[i][j]);
        fprintf(data,"\n");
    }
    fprintf(lod,"%i %f\n",born,energy);
}

t_viech* t_viech::getLMRCA(void){
        
}


