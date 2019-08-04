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

int determineWinds(float age, float* eWinds, float* mWinds, float* zWinds,
                        float massMsun, float zZsun, float TimeUnits, float dtFixed){

    bool oldEnough = (age < 0.01)?(false):(true);
    float windE = 0,  windM = 0, windZ = 0.0;
    float wind_factor = 0.0;
    float e_factor = 0.0;
 // I dont want to deal with new particles
    printf("Computing Winds for age = %f, Msun = %e\n", age, massMsun);
    if (StellarWinds && oldEnough){

        if (0.01 < age < 1.0){
            wind_factor =4.763
                         * min((0.01 + zZsun), 1.0) * age;
        }
        if (1 <= age < 3.5){
            wind_factor = 4.763*min(0.01+zZsun, 1.0)* 
                pow(age, 1.45+0.08*min(log(zZsun), 1.0));
        }
        if (3.5 <= age <= 100){
            wind_factor = 29.4*pow(age/3.5, 3.25)*0.0042;
        
        }
        if (100 < age){
            wind_factor = 0.42*pow(age/1000, -1.1)/(19.81/log(age));
            e_factor = 4.83;
        }
        if (age < 100){
            e_factor = 5.94e4/pow((1-age/2.5), 1.4)
                + pow(age/50.0, 5.0)+4.83;
        }
        windM = massMsun * wind_factor*
                        dtFixed*TimeUnits/3.1557e16; //Msun/Gyr
        printf("First winds mass = %e\n", windM);
        if (windM > massMsun){
            printf("Winds too large Mw = %e, Mp = %e",
                windM, massMsun);
            windM = 0.125*massMsun; // limit loss to huge if necessary.
        }
        windZ = 0.0278+ 0.0041* min(max(zZsun, 1.65), 5.0)*windM;
        windE = e_factor * 1e12 * windM; // what the actual fuck are the units here???
        *mWinds = windM;
        *zWinds = windZ;
        *eWinds = windE;
        printf("Winds Mass = %e\n",*mWinds);
        printf("Winds Energy = %e\n", *eWinds);
        printf("wind_factor = %f\n", wind_factor);
        printf("energy_factor = %f\n", e_factor);
        printf("Metallicity = %e\n", zZsun);
    }
    
    return SUCCESS;
}