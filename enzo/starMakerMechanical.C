/*
    Star maker paired with star_feedback_mechanical.src.  
    Uses all the criteria from Cen & Ostriker, but determines
    gas being converted to stars by analytic estimation of the 
    shielded fraction of gas from Krumholz & Gnedin 2011.
*/
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "Grid.h"
#include "StarParticleData.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int starCreationMechanical(int level, grid *tg, int NumberOfNewParticleSoFar,int MaxNewStars, float coolingTime,)
    {
        if (level < StarMakerMinimumRefinementLevel) return 0;
        // Get all our units sorted out
        float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
                TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
        if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	            &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
                fprintf(stderr, "Error in GetUnits.\n");
            return FAIL;    
        }          
        int ii = NumberOfNewParticlesSoFar;
        double G = GravConst;
        double msolar = SolarMass;
        double max_mass = StarMakerMaximumFormationMass / DensityUnits / pow(LengthUnits*CellWidth[0][0], 3.0) * msolar;
        int i,j,k;
        double sndspdC = 1.3095e8;  
        int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
        // find the useful fields well need later
        if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
            fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
            return FAIL;
        }
        float *temperature = new float[size];
        this->ComputeTemperatureField(temperature);
          int MetallicityField = FALSE, MetalNum;
        if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
            != -1)
            MetallicityField = TRUE;
        else
            MetalNum = 0;
        // analytic formula depend on metallicity of the gas
        if (MetallicityField == FALSE){
            fprintf(stderr, "StarMakerMechanical only functions with metallicity field active!!\n");
            ERROR_MESSAGE;
        }
        // Limit formation to the level specified by parameter
        float *cooling_time = new float[size];
        this->ComputeCoolingTime(cooling_time);
        /* Loop over all cells in this grid and check each one 
        for star formation criteria:
        1. is this the finest level of at this point in space?
        2. is the density > critical density (odthresh)
        3. is the flow converging?
        4. is the cooling time < dynamical time?
        5. is the M_jeans <= M_critical?
        6. does the gas have a non-zero shielded fraction?*/
        
        for ( i = 0+NumberOfGhostZones; i <= GridDimension[0]-NumberOfGhostzones; i++){
            for ( j = 0+NumberOfGhostZones; j <= GridDimension[1]-NumberOfGhostZones; j++){
                for ( k = 0+NumberOfGhostZones; j <= GridDimension[2]-NumberOfGhostZones; k++){
                    if (ii > MaxNewStars) break;
                // 1. Is this the finest level of refinement?
                    if (BaryonField[NumberOfBaryonFields] != 0.0) continue;

                // 2. Is density > odthresh
                    fprintf(stderr,"Density = %f\n",BaryonField[DensNum][i,j,k]);
                    if (BaryonField[DensNum][i,j,k] < StarMakerOverDensityThreshold) continue;
                // 3. Is flow converging?
                    double diver = 0.0;
                    if (imethod == 2){
                        diver = BaryonField[Vel1Num][i+1, j, k] - BaryonField[Vel1Num][i,j,k]
                            + BaryonField[Vel2Num][i,j+1,k]-BaryonField[Vel2Num][i,j,k]
                            + BaryonField[Vel3Num][i,j,k+1]-BaryonField[Vel3Num][i,j,k];
                    }
                    else{
                        diver = BaryonField[Vel1Num][i+1, j, k] - BaryonField[Vel1Num][i-1,j,k]
                            + BaryonField[Vel2Num][i,j+1,k]-BaryonField[Vel2Num][i,j-1,k]
                            + BaryonField[Vel3Num][i,j,k+1]-BaryonField[Vel3Num][i,j,k-1];
                    };
                    if (diver >= 0.0) continue;
                // 4. is cooling time < dynamical time or t < 1.1e4


                    float dtot = (BaryonField[DensNum][i,j,k] + dm[i,j,k])* DensityUnits;
                    float tdyn = sqrt(3.0 * pi / 32.0/G/dtot)/TimeUnits;
                    if ((tdyn < cooltime[i,j,k]) && temperature[i,j,k] > 1.1e4) 
                        continue;
                    
                //5. Is m_jeans < m_critical
                    float bmass = GridDensity[DensNum][i,j,k] * DensityUnits * pow(LengthUnits*CellWidth[0][0],3.0)/msolar;
                    float isosndsp2 = sndspdC * (temperature[i,j,k]);
                    float m_j = pi / (6.0*sqrt((BaryonField[DensNum][i,j,k]*DensityUnits)))
                                * pow((pi * isosndsp2 / G),1.5) / msolar;
                    if ((bmass > 1e3) && (m_j < bmass)) continue;
                    if ((bmass < 1e3) && (m_j < 1e3)) continue;
                //6. is there a non-zero shielded fraction of gas? 
                // estimated according to Krumholz & Gnedin 2011
                    float absGradRho = sqrt(pow(BaryonField[DensNum][i+1,j,k] - BaryonField[DensNum][i-1,j,k],2.0)
                            + pow(BaryonField[DensNum][i,j+1,k]-BaryonField[DensNum][i,j-1,k],2.0)
                            + pow(BaryonField[DensNum][i,j,k+1] - BaryonField[DensNum][i,j,k-1],2.0));
                    float tau = 434.8 * BaryonField[DensNum][i,j,k] * DensityUnits / absGradRho;
                    float phi = 0.756 * pow(1.0+3.1*(BaryonField[MetalNum][i,j,k])/0.02,0.365);
                    float psi = (0.6 * tau * 
                                (0.01 + BaryonField[MetalNum][i,j,k]/ 0.02))
                                / log(1.0 + 0.06 * phi * 0.01 * pow(phi, 2));
                    float shieldFrac = 1.0 - 3.0/(1.0 + 4.0 * psi);
                    if (shieldFrac < 0.0) continue;
                    if (shieldFrac > 1.0) shieldFrac = 1.0;
                    float t_free = sqrt(3.0 * pi / (32.0 * G * BaryonField[DensNum][i,j,k] * DensityUnits))/ TimeUnits;
                // If we've gotten here, its time to create a particle
                    ii = ii + 1;
                    fprintf(stderr, "Creating Particle\n");
                    float mp == 0.001 * BaryonField[DensNum][i,j,k];
                    if (shieldFrac * BaryonField[DensNum][i,j,k] /t_free > mp) 
                        mp = shieldFrac * BaryonField[DensNum][i,j,k] /t_free;
                    if (max_mass < mp) 
                        mp = max_mass;
                    if (mp > BaryonField[DensNum][i,j,k]) 
                        fprintf(stderr,"starMakerMechanical: Star mass > M_cell!");
                        ERROR_MESSAGE;
                    tg->ParticleMass[ii] = mp;
                    tg->ParticleAttribute[0][ii] = Time;
                    tg->ParticleAttribute[1][ii] = tdyn; // attr 1
                    if (StarMakerTrackEnergy == 1) tg->ParticleAttribute[1][ii] = 0.0;
                    tg->ParticlePosition[0][ii] = CellLeftEdge[0] + (float(i)-0.5) * CellWidth[0][0];
                    tg->ParticlePosition[1][ii] = CellLeftEdge[1] + (float(j)-0.5) * CellWidth[0][0];
                    tg->ParticlePosition[2][ii] = CellLeftEdge[2] + (float(k)-0.5) * CellWidth[0][0];

                    tg->ParticleVelocity[0][ii] = BaryonField[Vel1Num][i,j,k];
                    tg->ParticleVelocity[1][ii] =BaryonField[Vel2Num][i,j,k];
                    tg->ParticleVelocity[2][ii] = BaryonField[Vel3Num][i,j,k];

                    if (MetallicityField) tg->ParticleAttribute[2][ii] = BaryonField[MetalNum][i,j,k];
                    else BaryonField[MetalNum][ii] = 0.0; // attr 2

                    BaryonField[DensNum][i,j,k] -= mp;

                    if (ii == MaxNewStars) break;

                }
            }
        }
        if (ii > MaxNewStars)
            fprintf(stderr, "starMakerMech: Reached max count");
        delete [] cooling_time;
        delete [] temperature;
        return ii;
    }