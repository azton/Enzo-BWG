/*
    Routine to determine feedback quantities and couple them to the 
    Grid.  Coupling follows the implementation of Hopkins 2017 with
    modification for Enzo's fixed grid

    07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"

    int determineSN(float age, int* nSNII, int* nSNIA, float massMsun, 
                    float TimeUnits, float dtFixed);
    int determineWinds(float age, float* eWinds, float* zWinds, float* mWinds,
                        float massMsun, float zZsun, float TimeUnits, float dtFixed);
    int checkCreationCriteria(float* Density, float* Metals,
                        float* Temperature,float* DMField,
                        float* Vel1, float* Vel2, float* Vel3, 
                        float* CoolingTime, int GridDim,
                        float* shieldedFraction, float* freeFallTime, 
                        float* dynamicalTime, int i, int j, int k, 
                        float Time, float* RefinementField, float CellWidth);
    int FindField(int field, int farray[], int numfields);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);





int grid::MechStars_FeedbackRoutine(int level, float* mu_field)
{

    printf("IN FEEDBACK ROUTINE\n");
    if (NumberOfParticles == 0) {
        printf("NO PARITCLES FOUND\n");
        return SUCCESS;
    }
    float stretchFactor = 1.0;//1.5/sin(M_PI/10.0);
    /* Get units to use */

    int dim, i, j, k, index, size, field, GhostZones = DEFAULT_GHOST_ZONES;
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    
    /* Compute size (in floats) of the current grid. */
    
    size = 1;
    for (dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                        HMNum, H2INum, H2IINum, DINum, DIINum, HDINum;
    /* Find fields: density, total energy, velocity1-3. */

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                        Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    /* Set the units */
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
            TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;    
    } 
    double dx = CellWidth[0][0];    
    MassUnits = DensityUnits*pow(LengthUnits*dx, 3)/SolarMass;

    /* 
        get metallicity field and set flag; assumed true thoughout feedback
        since so many quantities are metallicity dependent
     */
    int MetallicityField = FALSE, MetalNum;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
        != -1)
        MetallicityField = TRUE;
    else
        MetalNum = 0;
    
    int numSN = 0;

    /* Begin Iteration of all particles */
    printf("\nIterating all particles\n");
    for (int pIndex=0; pIndex < NumberOfParticles; pIndex++){
    
        /* Selection criteria */

        if (ParticleType[pIndex] != 2.0
                || ParticleMass[pIndex] < 0.0
                ||ParticleAttribute[0][pIndex] < 0.0){
            continue;
        }
        // if (StarMakerAgeCutoff)
        //     if ((Time-ParticleAttribute[0][pIndex])
        //         *TimeUnits/(150*3.1557e7) > 150)
        //         continue;

        /* get index of cell hosting particle */
        float xp = ParticlePosition[0][pIndex];
        float yp = ParticlePosition[1][pIndex];
        float zp = ParticlePosition[2][pIndex];

        int ip = (xp-CellLeftEdge[0][0]-0.5*dx)/dx;
        int jp = (yp-CellLeftEdge[1][0]-0.5*dx)/dx;
        int kp = (zp-CellLeftEdge[2][0]-0.5*dx)/dx;



        /* error check particle position; Cant be on the border or outside grid
            If on border, reposition to within grid for CIC deposit */
        float age = (Time-ParticleAttribute[0][pIndex])*TimeUnits/3.1557e13;// Myr

        float gridDx = GridDimension[0]*dx;
        float gridDy = GridDimension[1]*dx;
        float gridDz = GridDimension[2]*dx;
       /* Keep particle 2 cells from edge since we cant transfer to 
            neighboring grids */
        float borderDx = (stretchFactor+1)*dx;
        if (xp > CellLeftEdge[0][0]+gridDx 
            || xp < CellLeftEdge[0][0]
            || yp > CellLeftEdge[1][0]+gridDy
            || yp < CellLeftEdge[1][0]
            || zp > CellLeftEdge[2][0]+gridDz
            || zp < CellLeftEdge[2][0]){
            fprintf(stderr, "Particle %d out of grid!\nage: %d, pos: %f, %f, %f GridEdge: %f %f %f", pIndex,
                age, xp, yp, zp, CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]
                );
            exit(137);
            }
        int shifted = 0;

        if (xp < CellLeftEdge[0][0]+borderDx){
            xp = CellLeftEdge[0][0]+borderDx+0.5*dx;
            shifted ++;
        }
        if (xp > CellLeftEdge[0][0]+gridDx-borderDx){
            xp -= CellLeftEdge[0][0]+gridDx-borderDx-0.5*dx;
            shifted = 1;
        }
        if (yp < CellLeftEdge[1][0]+borderDx){
            yp += CellLeftEdge[1][0]+borderDx+0.5*dx;
            shifted = 1;
        }
        if (yp > CellLeftEdge[1][0]+gridDx-borderDx){
            yp -= CellLeftEdge[1][0]+gridDx-borderDx-0.5*dx;
            shifted = 1;
        }
        if (zp < CellLeftEdge[2][0]+borderDx){
            zp += CellLeftEdge[2][0]+borderDx+0.5*dx;
            shifted = 1;
        }
        if (zp > CellLeftEdge[2][0]+gridDx-borderDx){
            zp -= CellLeftEdge[2][0]+gridDx-borderDx-0.5*dx;
            shifted = 1;
        }
        if (shifted > 0){
            //if (debug)
                fprintf(stderr, "Particle position shifted away from grid: %f %f %f\n", xp, yp, zp);
            ip = int((xp - CellLeftEdge[0][0])/dx + 0.5);
            jp = int((yp - CellLeftEdge[1][0])/dx + 0.5);
            kp = int((zp - CellLeftEdge[2][0])/dx + 0.5);
        }
        /* REMOVED: Check for continual formation, i guess. Only really done because
            Hopkins did it.  We can just make more stars next timestep I guess. 
            On the other hand, the function is already written... */
            /* create some stuff.  This is a lot of overhead for something
            optional... */

        /* Start actual feedback: Supernova calculations */
        index = ip+jp*GridDimension[0]+kp*GridDimension[0]*GridDimension[1];
        int nSNII = 0;
        int nSNIA = 0;
        int SNMassEjected = 0;
        
        /* determine how many supernova events */
        if (SingleSN)
        {
            determineSN(age, &nSNII, &nSNIA, ParticleMass[pIndex]*MassUnits,
                         TimeUnits, dtFixed);
            numSN += nSNII+nSNIA;
            if (numSN > 0)
                printf("\n\nSUPERNOVAE!!!! %d %d\n\n", nSNII, nSNIA);
            if (nSNII > 0 || nSNIA > 0){
                /* set feedback qtys based on number and types of events */
                    /* 1e51 erg per sn */
                float energySN = (nSNII + nSNIA)*1e51;
            
                    /*10.5 Msun ejecta for type II and IA*/
                SNMassEjected = (nSNII+nSNIA)*10.5;
                    /* Metal yeilds from starburst 99 */
                float zZsun = BaryonField[MetalNum][index]/BaryonField[DensNum][index]/0.02;
                float SNMetalEjected = nSNII*(1.91+0.0479*max(zZsun, 1.65));
                SNMetalEjected += nSNIA*(1.4);
                MechStars_DepositFeedback(energySN, SNMassEjected, SNMetalEjected,
                            &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                            &ParticlePosition[0][pIndex], &ParticlePosition[1][pIndex], &ParticlePosition[2][pIndex],
                            ip, jp, kp, size, mu_field);
            }
        }
        
        float windEnergy=0, windMass=0, windMetals=0;
        /* Do the same for winds. Cooling Radius is very small,
        So almost no energy is coupled, but some mass may be. */
        
        
        if (StellarWinds)
        {
            float zZsun = min(BaryonField[MetalNum][index]/BaryonField[DensNum][index]/0.02, 1e-8);
            determineWinds(age, &windEnergy, &windMass, &windMetals, 
                            ParticleMass[pIndex]*MassUnits, zZsun,
                            TimeUnits, dtFixed);
            if (windEnergy > 0){
            MechStars_DepositFeedback(windEnergy, windMass, windMetals,
                        &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                        &ParticlePosition[0][pIndex], &ParticlePosition[1][pIndex], &ParticlePosition[2][pIndex],
                        ip, jp, kp, size, mu_field);
            }
        }
        printf("Subtracting off mass %e\n",(windMass+SNMassEjected));
        ParticleMass[pIndex] -= (windMass+SNMassEjected)/MassUnits;
        printf("Post-feedback MP = %e\n", ParticleMass[pIndex]*MassUnits);
    }// end iteration over particles

    delete [] mu_field;
    delete [] BaryonField[NumberOfBaryonFields];   // refinement flag field
    BaryonField[NumberOfBaryonFields] = NULL;
    return SUCCESS;
}