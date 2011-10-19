//
//  main.cpp
//  ACM5
//
//  Created by Arend on 10/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "globals.h"
#include "chemistry.h"
#include "world.h"

int main (int argc, const char * argv[])
{
    t_world *myWorld=new t_world;
    myWorld->setup();
    while(true){
        myWorld->run_all_viecher();
        printf("%i  %f  %i  %i\n",
               myWorld->currentTime,
               myWorld->max_fitness,
               (int)myWorld->population.size(),
               myWorld->population[0]->how_many_genes());
    }
    return 0;
}

