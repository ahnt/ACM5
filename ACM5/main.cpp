//
//  main.cpp
//  ACM5
//
//  Created by Arend on 10/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <time.h>
#include "globals.h"
#include "chemistry.h"
#include "world.h"

#define updates 100000
int main (int argc, const char * argv[])
{
    t_world *myWorld=new t_world;
    srand((unsigned int)time(NULL));
    myWorld->setup();
    FILE *lod,*data,*comment;
    int update;
    for(update=0;update<updates;update++){
        myWorld->run_all_viecher(true,0.001,0.01,0.01);
        printf("%i  %f  %i  %i  %i  %i\n",
               myWorld->currentTime,
               myWorld->max_fitness,
               (int)myWorld->population.size(),
               myWorld->population[0]->how_many_genes(),
               (int)myWorld->population[0]->chromosome[0].size(),
               myWorld->population[0]->hmg);
    }
    lod=fopen("LOD.txt","w+t");
    data=fopen("DATA.txt","w+t");
    comment=fopen("COMMENT.txt","w+t");
    myWorld->population[0]->saveLOD(lod, comment, data);
    return 0;
}

