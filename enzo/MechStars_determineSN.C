/*
    Probabilistically determines supernova based on analytic
    starburst99 simulations.  fits taken from Hopkins 2017

    07/2019: Azton Wells
 */

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
#include "phys_constants.h"

int determineSN(float age, int* nSNII, int* nSNIA){

    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return SUCCESS;
    }
}