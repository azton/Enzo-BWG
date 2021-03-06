/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE STAR PARTICLES
/
/  written by: Greg Bryan
/  date:       February, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

#define STAR_PARTICLE_NUMBER_START 1000000000

/* Number of Star particles. */

SPEXTERN int NumberOfStarParticles;

/* Star particle parameters. */

SPEXTERN float StarMakerOverDensityThreshold;
SPEXTERN float StarMakerMassEfficiency;
SPEXTERN float StarMakerMinimumMass;
SPEXTERN float StarMakerMinimumDynamicalTime;
SPEXTERN float StarMassEjectionFraction;
SPEXTERN float StarMetalYield;
SPEXTERN float StarEnergyToThermalFeedback;
SPEXTERN float StarFeedbackKineticFraction;
SPEXTERN float StarMakerExplosionDelayTime;
SPEXTERN float StarEnergyToStellarUV;
SPEXTERN float StarEnergyToQuasarUV;
/* adding distributed starmaker and feedback */
SPEXTERN int StarFeedbackDistRadius;
SPEXTERN int StarFeedbackDistCellStep;
SPEXTERN int StarFeedbackDistTotalCells;
/* mechanical feedback */
SPEXTERN int StellarWinds;
SPEXTERN int SingleSN;
SPEXTERN float StarMakerMaximumFormationMass;
SPEXTERN float StarMakerMaximumMass;
SPEXTERN int DepositUnresolvedEnergyAsThermal;
SPEXTERN int StarMakeLevel;
SPEXTERN int NEvents;
SPEXTERN int AnalyticSNRShellMass;
SPEXTERN int UnrestrictedSN;