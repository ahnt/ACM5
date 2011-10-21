#include "world.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

t_world::t_world()
{
	int x,y;
	max_fitness=(double) 0.0;
	ave_fitness=(double) 0.0;
    population.clear();
	chemistry = new t_chemistry;
	chemistry->setup();
	for(x=0;x<mx+1;x++)
		for(y=0;y<my+1;y++)
			feld[x][y]=0;
	limits=(double)limit;
	for(x=0;x<608;x++)
	{
		if(x<feeding_level)
		{
			this->compound_avail[x]=(double)1-((double)x/(double)608);
			this->fittness_for_compound[x]=(double) 0.0;
		}
		else 
		{
			this->compound_avail[x]=(double)0.0;
			this->fittness_for_compound[x]=(((double)x*(double)x)/((double)608*(double)608));
			/*if (x<141) this->fittness_for_compound[x]=(double)0.0;
			if ((x>230)&&(x<375)) this->fittness_for_compound[x]=(double)0.0;*/
		}
		//printf("%i	%f\n",x,fittness_for_compound[x]);
	}
	inhibitor=5000;
    currentTime=0;
}

t_world::~t_world()
{
}

void t_world::setup()
{
	int z;
	printf("	* SETUP VIECHER RANDOM...");
	_viecher_n=start_viecher;
    order.clear();
	for(z=0;z<start_viecher;z++)
	{
		t_viech *dummy = new t_viech;
        /*
        if(z==0){
            dummy->load_genome((char*)"startGenome.txt");
        }
        else{
            dummy->inherit(population[0], 0.1, 0.0, 0.1);
        }//*/
        dummy->fill_chr_rand(start_chr_length_A,start_chr_length_B);
		dummy->setup_proteom(chemistry);
        //if(z==0)
        //    dummy->showProteomNetwork(chemistry);
		dummy->xpos=0;
		dummy->ypos=0;
        dummy->born=currentTime;
        population.push_back(dummy);
        order.push_back(z);
	}
	printf("	...COMPLETED *\n");
	chemistry->np[608]=(double) 0.0;
	for(z=0;z<608;z++)
	{
		chemistry->np[z]=(double) start_compounds*(double) compound_avail[z];
		chemistry->np[608]+=(double) start_compounds*(double) compound_avail[z];
		chemistry->xp[z]=16-(rand()%33);
		chemistry->yp[z]=16-(rand()%33);
		while ((chemistry->xp[z]==0)||(chemistry->yp[z]==0))
		{
			chemistry->xp[z]=64-(rand()%129);
			chemistry->yp[z]=64-(rand()%129);
		}
	}

}

void t_world::run_all_viecher(bool replace,double mutationRate,double duplicationRate,double deletionRate)
{
	int32_t k;
	int u;
	//double max_fitness;
	int16_t l,m;
	int born_break=0;
	int nx,ny;
	int z=start_viecher;
    int who;
	int pos_in_line=0;
	long loop_counter=0;
	for(k=0;k<population.size();k++){ 
        int a=rand()%population.size();
        int b=rand()%population.size();
        int c=order[a];
        order[a]=order[b];
        order[b]=c;
    }
    currentTime++;
	max_fitness=(double)0.0;
	ave_fitness=(double)0.0;
	for(k=0;k<population.size();k++)
	{
		population[order[k]]->step(chemistry);
		population[order[k]]->ready_to_divide(this->limits,this->fittness_for_compound);
		if(population[order[k]]->energy>max_fitness)
			max_fitness=population[order[k]]->energy;
		ave_fitness+=population[order[k]]->energy;
	}
	ave_fitness=ave_fitness/(double)population.size();
	for(k=0;k<population.size();k++)
		if ((((rand()&sudden_death_val)==sudden_death_val)||(population[k]->age>max_age))){
			for(l=0;l<608;l++)
			{
				chemistry->np[l]+=population[k]->n_comp[l];
				chemistry->np[608]+=population[k]->n_comp[l];
			}
            population[k]->nrPointingAtMe--;
            if(population[k]->nrPointingAtMe==0)
                delete population[k];
            population[k]=NULL;
		}
    if(replace)
        for(k=0;k<population.size();k++)
            if(population[k]==NULL){
                do{
                    who=rand()%population.size();
                }while((population[who]==NULL)||(population[who]->born==currentTime)||(((double)rand()/(double)RAND_MAX)>(population[who]->energy/max_fitness)));
                t_viech *dummy=new t_viech;
                dummy->inherit(population[who],mutationRate,duplicationRate,deletionRate);
                dummy->born=currentTime;
                dummy->setup_proteom(chemistry);	
                for(m=0;m<608;m++)
                {
                    population[who]->n_comp[m]=population[who]->n_comp[m]/(double)4.0;
                    chemistry->np[m]+=population[who]->n_comp[m]*(double)2.0;
                    chemistry->np[608]+=population[who]->n_comp[m]*(double)2.0;
                    if(m<feeding_level)
                    {
                        dummy->n_comp[m]=chemistry->c(population[who]->xpos,population[who]->ypos,m)/(double)1000000.0;
                        chemistry->np[608]-=dummy->n_comp[m];
                    }
                }
                population[who]->get_total_n();
                dummy->get_total_n();
                dummy->xpos=population[who]->xpos+(double)(((rand()&255)-(rand()&255))/(double)2560);
                dummy->ypos=population[who]->ypos+(double)(((rand()&255)-(rand()&255))/(double)2560);
                population[k]=dummy;
            }
}

void t_world::depopulate(void)
{
	int i;
    for(i=0;i<population.size();i++){
        population[i]->nrPointingAtMe--;
        if(population[i]->nrPointingAtMe==0)
            delete population[i];
    }
    population.clear();
}



