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

    extern "C"  void FORTRAN_NAME(cic_deposit)(float* xPosition, float* yPosition,
        float* zPosition, int* gridRank, int* nParticles, 
        float* DepositQuantity, float* FieldToDepositTo,
        float* leftEdge, int* xGridDim, int* yGridDim, 
        int* zGridDim, float* gridDx);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
    int transformComovingWithStar(float* Density, float* Metals, 
        float* Vel1, float* Vel2, float* Vel3, 
        float up, float vp, float wp,
        int sizeX, int sizeY, int sizeZ, int direction);
    int FindField(int field, int farray[], int numfields);


int grid::MechStars_DepositFeedback(float ejectaEnergy, 
                        float ejectaMass, float ejectaMetal,
                        float* up, float* vp, float* wp,
                        float* xp, float* yp, float* zp,
                        int ip, int jp, int kp,
                        int size, float* muField){
    
    /*
     This routine will create an isocahedron of coupling particles, where we determine
        the feedback quantities.  The vertices of the isocahedron are ~coupled particles
        and all have radius dx from the source particle. 
        Each vertex particle will then be CIC deposited to the grid!
    */
    int index = ip+jp*GridDimension[0]+kp*GridDimension[0]*GridDimension[1];
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    
    /* Compute size (in floats) of the current grid. */
    float stretchFactor =1.25;//1.5/sin(M_PI/10.0);  // How far should cloud particles be from their host
                                // in units of dx
    /* Cloud geometry:
    1 = direct couple to neighbor cells along axes
    2 = icosahedron coupling
    3 = dual icosahedtron
    4 = dodecahedtron 
    5 = dodecahedron with isocahedron*/
    int geometry = 5; 
    size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    float DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
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
            &TimeUnits, &VelocityUnits, 
            &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;    
    } 
    double dx = CellWidth[0][0];        

    /*debug */
    fprintf(stdout, "depositing quantities: Energy %e, Mass %e, Metals %e\n",
            ejectaEnergy, ejectaMass, ejectaMetal);

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

    /* set other units that we need */
    MassUnits = DensityUnits*pow(LengthUnits*dx, 3)/SolarMass;
    float EnergyUnits = DensityUnits*pow(LengthUnits*dx, 3) 
                    * VelocityUnits*VelocityUnits;//[g cm^2/s^2] -> code_energy
    float MomentaUnits = MassUnits*SolarMass*VelocityUnits;  

    /* Make small 5^3 copys of fields to work with. These will conatin the added deposition
        of all quantities and be coupled to the grid after the cic deposition. */
    float *density = new float[size];
    float *metals = new float[size];
    float *u = new float[size];
    float *v = new float[size];
    float *w = new float[size];
    float *totalEnergy = new float[size];
    for (int i=0; i<size; i++){
        density[i] = 0.0;
        metals[i] = 0.0;
        u[i] = 0.0;
        v[i] =0.0;
        w[i] = 0.0;
        totalEnergy[i] = 0.0;
    }
    /* Transform coordinates so that metals is fraction (rho metal/rho baryon)
        u, v, w -> respective momenta.  Use -1 to reverse transform after.*/
    /* these temp arrays are implicitly comoving with the star! */
    /* Array of coordinates of the isocahedron vertices scaled by r=dx */
    
    float phi = (1.0+sqrt(5))/2.0; //Golden Ratio
    float iphi = 1.0/phi; // inverse GR
    int nCouple = 0; // number of coupled neighbors
    float A = 0.0; //scaling to transform vectors into position space
    /* Particle Vectors */
    float *CloudParticleVectorX;
    float *CloudParticleVectorY;
    float *CloudParticleVectorZ;
    
/* Direct neighbor coupling */
    if (geometry == 1){

    
        A = stretchFactor*dx; //for isocahedron
        nCouple = 6;

/* coupling direct to neighbor cells */
        CloudParticleVectorX = new float[nCouple]{0, 0, 1, -1, 0, 0};
        CloudParticleVectorY = new float[nCouple]{1, -1, 0, 0, 0, 0};
        CloudParticleVectorZ= new float[nCouple]{0, 0, 0, 0, 1, -1};
    }
/* ISOCOHEDRON COUPLING */
    if (geometry == 2){
        nCouple = 12;
        A = stretchFactor*dx/(2.0*sin(2*M_PI/5));
        CloudParticleVectorX = new float[nCouple]{0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi};
        CloudParticleVectorY = new float[nCouple]{1, 1, -1, -1, phi, -phi, -phi, phi, 0, 0, 0, 0};
        CloudParticleVectorZ = new float[nCouple]{phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1};        
    }
// /* dual Isocahedron coupling */
    if (geometry == 3){
        nCouple = 24;
        A = stretchFactor*dx/(2.0*sin(2*M_PI/5));
        CloudParticleVectorX = new float[nCouple]{0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi,
                                        0, 0, 0, 0, phi, phi, -phi, -phi, 1, 1, -1, -1};
        CloudParticleVectorY = new float[nCouple]{1, 1, -1, -1, phi, -phi, -phi, phi, 0, 0, 0, 0,
                                        phi, phi, -phi, -phi, 1, -1, -1, 1, 0, 0, 0, 0};
        CloudParticleVectorZ = new float[nCouple]{phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1,
                                        1, -1, -1, 1, 0, 0, 0, 0, phi, -phi, -phi, phi};
    }
/* DODECAHEDRON CLOUD */
    if (geometry == 4){
            nCouple = 20;
            A = stretchFactor*dx/(sqrt(3)*sin(2*M_PI/5));

            CloudParticleVectorX = new float[nCouple]{1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
                                            iphi, iphi, -iphi, -iphi, phi, phi, -phi, -phi};
            CloudParticleVectorY = new float[nCouple]{1,1,-1,-1, 1, 1, -1, -1, iphi, iphi, -iphi, -iphi,
                                            phi, -phi, phi,-phi, 0, 0, 0, 0};
            CloudParticleVectorZ = new float[nCouple]{1,-1, 1,-1, 1,-1, 1,-1, phi,-phi, phi,-phi, 
                                            0, 0, 0, 0, iphi, -iphi, iphi, -iphi};
    }
    /* A DODECAHEDRON+ISOCAHEDRON */
    if (geometry == 5){
        nCouple = 32;
        A = stretchFactor*dx;

        /* Dodec points followed by isoca points */
        CloudParticleVectorX = new float[nCouple]{1, 1, 1, 1, -1, -1, -1, -1,
                                             0, 0, 0, 0,
                                            iphi, iphi, -iphi, -iphi, phi, phi, -phi, -phi, 0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi};
        CloudParticleVectorY = new float[nCouple]{1,1,-1,-1, 1, 1, -1, -1, iphi, iphi, -iphi, -iphi,
                                            phi, -phi, phi,-phi, 0, 0, 0, 0,1, 1, -1, -1, phi, -phi, -phi, phi, 0, 0, 0, 0};
        CloudParticleVectorZ = new float[nCouple]{1,-1, 1,-1, 1,-1, 1,-1, 
                                                phi,-phi, phi,-phi, 
                                            0, 0, 0, 0, iphi, -iphi, iphi, -iphi, phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1};
                            
    }   
            /* Set position of feedback cloud particles */

            float CloudParticlePositionX [nCouple];
            float CloudParticlePositionY [nCouple];
            float CloudParticlePositionZ [nCouple];

                /*all possible values of x,y,z with origin particle at x=y=z=0.0 */
            for (int cpInd = 0; cpInd < nCouple; cpInd++){
                float norm = sqrt(CloudParticleVectorX[cpInd]*CloudParticleVectorX[cpInd]
                    + CloudParticleVectorY[cpInd]*CloudParticleVectorY[cpInd]
                    + CloudParticleVectorZ[cpInd]*CloudParticleVectorZ[cpInd]);
                /* in this cloud, take the coupling particle position as 0.5, 0.5, 0.5 */
                CloudParticlePositionX[cpInd] = *xp+CloudParticleVectorX[cpInd]/norm*
                            A-CellLeftEdge[0][0];
                CloudParticlePositionY[cpInd] = *yp+CloudParticleVectorY[cpInd]/norm*
                            A-CellLeftEdge[1][0];
                CloudParticlePositionZ[cpInd] = *zp+CloudParticleVectorZ[cpInd]/norm*
                            A-CellLeftEdge[2][0];
                /* turn the vectors back into unit-vectors */
                CloudParticleVectorZ[cpInd] *= norm;
                CloudParticleVectorY[cpInd] *= norm;
                CloudParticleVectorX[cpInd] *= norm;  
                
            }
    /* Each particle gets 1/12 of energy, momenta, mass, and metal.  There
    are no varying vector / scalar weights to worry about.  The momenta coupled
    is simply \hat(r_ba) p/12 for r_ba the vector from source to coupled
    particle.  */
    /* cooling radius of the host cell  */
    float zZsun = BaryonField[MetalNum][index]/BaryonField[DensNum][index]/0.02;
    float fz = (zZsun < 0.01)? (2.0): (pow(zZsun, -0.14));

    /* conversions */
    float mhCode = mh/MassUnits/SolarMass; // hydrogen mass in code
    float densityFactor = max(BaryonField[DensNum][index]*DensityUnits/mh/0.6, 1e-8);
    float CoolingRadius = 28.4 *
        pow(densityFactor, -3.0/7.0)
        *pow(ejectaEnergy/1.0e51, 2.0/7.0)* fz;
    printf("cooling radius [pc] = %f\n %f %f %f\n", 
            CoolingRadius, densityFactor, ejectaEnergy/1e51, fz);
    /* Calculate coupled energy scaled by reduction to account for unresolved
    cooling, then use that energy to calculate momenta*/
    float coupledEnergy = ejectaEnergy;
    float nBaryons = BaryonField[DensNum][index]
                    *DensityUnits/mh/muField[index];
    printf("Dx [pc] = %f\n", dx*LengthUnits/pc_cm);
    float dxRatio = stretchFactor*dx*LengthUnits/pc_cm/CoolingRadius;
    int usePt = 0;
    float coupledMomenta = 0.0;
    /* Hopkins uses ratio of masses to determine how to couple.
        Radius here is well-known and fixed, so we use that instead */
    if (dxRatio > 1.0){ 
        if (ejectaEnergy < 1e10 || dxRatio > 100){
            coupledEnergy = 0.0;
            coupledMomenta = 0.0;
        }else{
            coupledEnergy = ejectaEnergy*pow(dxRatio, -6.5);
            usePt = 1;
            
        /* Determine coupled momenta if rc < dx use scaled energy
        else couple p_ej*(1+dx/r_cool) */
            printf("Using P_t with Nb = %f, E= %e\n",nBaryons, coupledEnergy/1e51);
            float Efactor = coupledEnergy/1e51;
            coupledMomenta = 4.8e5*pow(nBaryons, -1.0/7.0)
                * pow(Efactor, 13.0/14.0) * fz; //Msun*km/s
            coupledMomenta *= (SolarMass*1e5); //g*cm/s
        }
    } else {
        printf("Directly calculating momenta using energy = %e and mass = %e\n", 
                    ejectaEnergy, ejectaMass*SolarMass);
        coupledMomenta = pow(2.0*ejectaEnergy*(ejectaMass*SolarMass), 0.5) 
                                * pow(1.0+dxRatio*dxRatio, 0.5); //g*cm/s
        

    }
    float shellMass = 0.0, shellVelocity = 0.0;
    /* If resolution is in a range compared to Rcool and
        Analytic SNR shell mass is on, adjust the shell mass */
    if (dxRatio <= 100 & dxRatio >= 0.1 & coupledEnergy > 0
        & AnalyticSNRShellMass){
            shellVelocity = 413 *pow(nBaryons, 1.0/7.0)
                *pow(zZsun, 3.0/14.0)*pow(coupledEnergy/EnergyUnits/1e51, 1.0/14.0)
                *pow(dxRatio, 7.0/3.0);//km/s
            shellVelocity *= 1e5; //cm/s
            shellMass = coupledMomenta/shellVelocity/SolarMass;
            /* cant let host cell evacuate completely! */
            if (shellMass >= 0.125*BaryonField[DensNum][index]*MassUnits)
                shellMass = 0.125*BaryonField[DensNum][index]*MassUnits;
    }
    printf("Ejecta Mass = %f\n", ejectaMass);
    float coupledMass = shellMass+ejectaMass;
    /* rescale momentum for new shell */
    float shellMetals = zZsun*0.02 * shellMass;
    float coupledMetals = ejectaMetal + shellMetals;

    fprintf(stdout, "Coupled Momentum: %e\n", coupledMomenta/2e37/float(nCouple));
    fprintf(stdout, "Coupled shell mass: %e\n", shellMass);
    /* Reduce coupled quantities to per-particle quantity and convert to 
        code quantities.
        Hopkins has complicated weights due to complicated geometry. 
            This is more simple since our coupled particles are 
            spherically symmetric about the feedback particle*/
    
    coupledEnergy /= float(nCouple);
    coupledMass /= float(nCouple);
    printf("Coupled Mass = %e\n", coupledMass);
    coupledMetals /= float(nCouple);
    printf("Coupled metals = %e\n", coupledMetals);
    coupledMomenta /= float(nCouple);
    /* Transform coupled quantities to code units */
    printf("PreUnits energy = %e\n", coupledEnergy);
    coupledEnergy /= EnergyUnits;
    printf("PostUnits energy = %e\n", coupledEnergy);
    coupledMass /= MassUnits;
    coupledMetals /= MassUnits;
    printf("Pre units momenta = %e\n", coupledMomenta);
    coupledMomenta /= MomentaUnits;
    printf("Post units Momenta = %e\n", coupledMomenta);
    /* CIC deposit the particles with their respective quantities */
    for (int n = 0; n < nCouple; n++){
        /* set vector qtys for this particle */
        // MechStars_CICdeposit(coupledEnergy, coupledMass, coupledMetals,
        //         coupledMomenta, CloudParticlePositionX[n], 
        //         CloudParticlePositionY[n],CloudParticlePositionZ[n], 
        //         CloudParticleVectorX[n], CloudParticleVectorY[n],
        //         CloudParticleVectorZ[n], density, metals, totalEnergy, u,
        //         v, w);
        float pX = coupledMomenta*CloudParticleVectorX[n];
        float pY = coupledMomenta*CloudParticleVectorY[n];
        float pZ = coupledMomenta*CloudParticleVectorZ[n];
        int np = 1;
        FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
            &CloudParticlePositionZ[n], &GridRank, &np,&coupledMass, density, GridLeftEdge, 
            &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
        FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
            &CloudParticlePositionZ[n], &GridRank, &np,&coupledMetals, metals, GridLeftEdge, 
            &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
        if (pX != 0.0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank, &np,&pX, u, GridLeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
        if (pY != 0.0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&pY, v, GridLeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
        if (pZ != 0.0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&pZ, w, GridLeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);


        FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
            &CloudParticlePositionZ[n], &GridRank,&np,&coupledEnergy, totalEnergy, GridLeftEdge, 
            &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);

    }
    /* Deposit one negative mass particle centered on star to account for 
        shell mass leaving host cells */
    int np = 1;
    shellMass *= -1/MassUnits;
    FORTRAN_NAME(cic_deposit)(xp, yp, zp, &GridRank,&np,&shellMass, density, GridLeftEdge, 
        &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
    shellMetals *= -1/MassUnits;   
    FORTRAN_NAME(cic_deposit)(xp, yp, zp, &GridRank,&np,&shellMetals, metals, GridLeftEdge, 
        &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx);
    /* transform the grid to comoving with star ; wouldnt recommend this on root grid...*/
    transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum], 
                        BaryonField[Vel1Num],BaryonField[Vel2Num],BaryonField[Vel3Num],
                        *up, *vp, *wp, GridDimension[0], GridDimension[1],
                        GridDimension[2], 1);
    for (int index = 0; index < size; index++){ 
                
                float deltaMass = (density[index]+BaryonField[DensNum][index])
                                    /BaryonField[DensNum][index];
                /* Couple placeholder fields to the grid, account 
                    for grids that got initialized to -0.0*/
                BaryonField[DensNum][index] += density[index];
                // if (density[index] < 0.0){
                //         printf("NEGATIVE DEPOSITION DENSITY = %f\n",
                //                         density[index]);
                //         printf("%d\n", index);
                //         exit(253);
                //     }

            /*  */


                //Metals transformed back to density in transform routine
                
                BaryonField[MetalNum][index] += metals[index]; 
                BaryonField[TENum][index] += totalEnergy[index]
                            /(BaryonField[DensNum][index]*BaryonField[DensNum][index]);
                BaryonField[Vel1Num][index] += u[index];
                BaryonField[Vel2Num][index] += v[index];
                BaryonField[Vel3Num][index] += w[index];
            }


    /* Sum of feedback quantities: */

    float dsum = 0.0, zsum=0.0, psum=0.0, psqsum =0.0, esum=0.0;
    for (int i = 0; i<size; i++){
        dsum += density[i];
        zsum += metals[i];
        psum += u[i]+v[i]+w[i];
        esum += totalEnergy[i];
        psqsum += u[i]*u[i]+v[i]*v[i]+w[i]*w[i];

    }
    if (debug){
        printf("Sum Mass Deposited = %e\n", dsum*MassUnits);
        printf("Sum Metals Deposited = %e\n", zsum*MassUnits);
        printf("Sum momenta magnitude = %e\n", sqrt(psqsum)*MomentaUnits/2e37);

        printf("Sum momenta error deposited = %e\n", psum*MomentaUnits/2e37);
        printf("Sum energy = %e\n", esum * EnergyUnits);
    }
    /* Transform the grid back */

        transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum], 
                        BaryonField[Vel1Num],BaryonField[Vel2Num],
                        BaryonField[Vel3Num],*up, *vp, *wp, 
                        GridDimension[0], GridDimension[1],
                        GridDimension[2], -1);
    delete [] CloudParticleVectorX;
    delete [] CloudParticleVectorY;
    delete [] CloudParticleVectorZ;
    return SUCCESS;
}