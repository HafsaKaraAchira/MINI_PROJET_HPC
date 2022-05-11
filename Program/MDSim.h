#ifndef MDSIM_H
#define MDSIM_H

#ifndef VEC_H
#include "vec.h"
#endif

#ifndef MAT_H
#include "mat.h"
#endif

#ifndef MDPARTICLELIST_H
#include "MDParticleList.h"
#include "Gaz.h"
#endif


typedef struct graphData_t{
	double eKin, ePot;
	graphData_t* next;
} graphData_t;

class MDSim{

public:
	MDSim(double,int,int);
	~MDSim() {};

	void initSim(bool, vec, int, double);
	void velocityVerletStep(bool);
	void updateGraphs();
	void resetGraphs();
	void resetRadialDistribution();
	vec getRadialDistribution();
	int getHistogramResolution();
	void resetDirectionalDistribution();
	mat getDirectionalDistribution();
	double getDt();
	double getT();
	int getDim();
	double getEPotMin();
	Gaz* getGaz();

	//  initialize positions on simple cubic lattice, also calls function to initialize velocities
	void initialize();  
	//  update positions and velocities using Velocity Verlet algorithm 
	//  print particle coordinates to file for rendering via VMD or other animation software
	//  return 'instantaneous pressure'
	double VelocityVerlet(double dt, int iter);  
	//  Compute Force using F = -dV/dr
	//  solve F = ma for use in Velocity Verlet
	void computeAccelerations();
	//  Numerical Recipes function for generation gaussian distribution
	double gaussdist();
	//  Initialize velocities according to user-supplied initial Temperature (Tinit)
	void initializeVelocities();
	//  Compute total potential energy from particle coordinates
	double Potential();
	//  Compute mean squared velocity from particle velocities
	double MeanSquaredVelocity();
	//  Compute total kinetic energy from particle mass and velocities
	double Kinetic();

    MDParticleList* particles;
	bool pause;
	int numSteps;
	graphData_t *graphDataFirst, *graphDataLast; // graphDataFirst contains max values and leads to values for t=0

private:
	double refreshVerletLists(bool, bool);
	double velocityVerletForce();

	Gaz* gaz;
	double (*pot)(vec), dt, t, r_inter, r_verlet, eKin, ePot, ePotMin, histogramLength;
	vec (*f)(vec), simBox, radial;
    mat directional;
	mat r,v,a,F;
	int histogramResolution;
	bool periodic;
};

#endif
