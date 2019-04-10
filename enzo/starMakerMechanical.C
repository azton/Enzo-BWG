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
#include "StarParticleData.h"
#include "phys_constants.h"
int starMakerMechanical(int *nx, int *ny, int *nz,
            float *d, float *dm, float *temp, 
            float *u, float *v, float *w,
            float *cooltime, float *dt, 
            float *r, float *metal, float *zfield1, 
            float *zfield2, float *dx, FLOAT *t, 
            float *z,int *procnum,float *dunits, 
            float *x1, float *v1, float *t1,
            int *nmax, FLOAT *xstart, FLOAT *ystart, 
            FLOAT *zstart,int *ibuff,int *imetal, hydro_method *imethod, 
            float *mintdyn,float *odthresh, float *massff, 
            float *smthrest, int *level, int *np, 
            FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	        float *mp, float *tdp, float *tcp, float *metalf,
	        int *imetalSNIa, float *metalSNIa, float *metalfSNIa, float *exptime,
            float *max_form_mass)
    {
        int ii = *np;
        double G = GravConst;
        double msolar = SolarMass;
        double max_mass = *max_form_mass / *dunits / pow(*dx * *x1, 3.0) * msolar;
        int i,j,k;
        double sndspdC = 1.3095e8;
        /* Loop over all cells in this grid and check each one 
        for star formation criteria:
        1. is this the finest level of at this point in space?
        2. is the density > critical density (odthresh)
        3. is the flow converging?
        4. is the cooling time < dynamical time?
        5. is the M_jeans <= M_critical?
        6. does the gas have a non-zero shielded fraction?*/
        for ( i = 0+*ibuff; i <= *nx-*ibuff; i++){
            for ( j = 0+*ibuff; j <= *ny-*ibuff; j++){
                for ( k = 0+*ibuff; j <= *nz-*ibuff; j++){
                    if (ii > *nmax) break;
                // 1. Is this the finest level of refinement?
                    
                    if (r[i,j,k] != 0.0) continue;

                // 2. Is density > odthresh
                    if (d[i,j,k] < *odthresh) continue;
                // 3. Is flow converging?
                    double diver = 0.0;
                    if (*imethod == 2){
                        diver = u[i+1, j, k] - u[i,j,k]
                            + v[i,j+1,k]-v[i,j,k]
                            + w[i,j,k+1]-w[i,j,k];
                    }
                    else{
                        diver = u[i+1, j, k] - u[i-1,j,k]
                            + v[i,j+1,k]-v[i,j-1,k]
                            + w[i,j,k+1]-w[i,j,k-1];
                    };
                    if (diver >= 0.0) continue;
                // 4. is cooling time < dynamical time or t < 1.1e4
                    double dtot = (d[i,j,k] + dm[i,j,k])* (*dunits);
                    double tdyn = sqrt(3.0 * pi / 32.0/G/dtot)/(*t1);
                    if ((tdyn < cooltime[i,j,k]) && temp[i,j,k] > 1.1e4) 
                        continue;
                    
                //5. Is m_jeans < m_critical
                    double bmass = d[i,j,k] * (*dunits)* pow(((*x1)*(*dx)),3.0)/msolar;
                    double isosndsp2 = sndspdC * (temp[i,j,k]);
                    double m_j = pi / (6.0*sqrt((d[i,j,k]*(*dunits))))
                                * pow((pi * isosndsp2 / G),1.5) / msolar;
                    if ((bmass > 1e3) && (m_j < bmass)) continue;
                    if ((bmass < 1e3) && (m_j < 1e3)) continue;
                //6. is there a non-zero shielded fraction of gas? 
                // estimated according to Krumholz & Gnedin 2011
                    double absGradRho = sqrt(pow(d[i+1,j,k] - d[i-1,j,k],2.0)
                            + pow(d[i,j+1,k]-d[i,j-1,k],2.0)
                            + pow(d[i,j,k+1] - d[i,j,k-1],2.0));
                    double tau = 434.8 * d[i,j,k] * *dunits / absGradRho;
                    double phi = 0.756 * pow(1.0+3.1*(metal[i,j,k])/0.02,0.365);
                    double psi = (0.6 * tau * 
                                (0.01 + metal[i,j,k]/ 0.02))
                                / log(1.0 + 0.06 * phi * 0.01 * pow(phi, 2));
                    double shieldFrac = 1.0 - 3.0/(1.0 + 4.0 * psi);
                    if (shieldFrac < 0.0) continue;
                    if (shieldFrac > 1.0) shieldFrac = 1.0;
                    double t_free = sqrt(3.0 * pi / (32.0 * G * d[i,j,k] * *dunits))/ *t1;
                // If we've gotten here, its time to create a particle
                    ii = ii + 1;
                    mp[ii] == 0.001 * d[i,j,k];
                    if (shieldFrac * d[i,j,k] /t_free > mp[ii]) mp[ii] = shieldFrac * d[i,j,k] /t_free;
                    if (max_mass < mp[ii]) mp[ii] = max_mass;
                    //mp[ii] = max(0.001 * d[i,j,k],
                      //  min(shieldFrac * d[i,j,k]/t_free, max_mass));
                    if (mp[ii] > d[i,j,k]) 
                        fprintf(stderr,"starMakerMechanical: Star mass > M_cell!");
                    tcp[ii] = *t;
                    tdp[ii] = tdyn;
                    xp[ii] = *xstart + (double(i)-0.5) * *dx;
                    yp[ii] = *ystart + (double(j)-0.5) * *dx;
                    zp[ii] = *zstart + (double(k)-0.5) * *dx;

                    up[ii] = u[i,j,k];
                    vp[ii] = v[i,j,k];
                    wp[ii] = w[i,j,k];

                    if (*imetal==1) metalf[ii] = metal[i,j,k];
                    else metalf[ii] = 0.0;

                    d[i,j,k] -= mp[ii];

                    if (ii == *nmax) break;

                }
            }
        }
        *np = ii;
        if (ii >= *nmax)
            fprintf(stderr, "starMakerMech: Reached max count");
        return 0;
    }