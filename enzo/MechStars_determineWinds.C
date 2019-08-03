/*
    Determines wind feedback parameters according to fits in Hopkins 2017:
    These fits are known to be erroneous, need to re-run and fit using SB99 sims.

    07/2019: Azton Wells
 */

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"

int determineWinds(float age, float* eWinds, float* mWinds, float* zWinds ){
    *eWinds = *mWinds = *zWinds = 0.0;
    return SUCCESS;
}