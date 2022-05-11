#include <cmath>
#include <iostream>
#include "vec.h"
#include "mat.h"
#include "MDParticle.h"
#include "MDParticleList.h"
#include "MDParticleListEntry.h"
//#include "AndersonThermostat.h"
#include "Gaz.h"
#include "MDSim.h"
using namespace std;

MDParticle::MDParticle(int dim) : r(0, dim), v(0, dim), a(0, dim), verletList(new MDParticleList()) {}

MDParticle::MDParticle(vec r_, vec v_) : dim(r_.get_dim()), r(r_), v(v_), verletList(new MDParticleList()) {}

MDParticle::~MDParticle(){

}

// Deleting a list also deletes its entries
MDParticleList::~MDParticleList(){
	this->clear();
}

int MDParticleList::getLength(){
	return length;
}

MDParticleListEntry* MDParticleList::getFirst(){
	return first;
}

MDParticleListEntry* MDParticleList::getLast(){
	return last;
}

void MDParticleList::addParticle(MDParticle* particle){
	if (this->last) this->last = (this->last->next = new MDParticleListEntry(particle, this, this->last));
	else this->last = (this->first = new MDParticleListEntry(particle, this, this->last));
	this->length++;
}

void MDParticleList::clear(){
	while (this->first) delete this->first; // ParticleListEntry destructor handles everything (see below).
}

// Deleting a list entry keeps list intact.
MDParticleListEntry::~MDParticleListEntry(){
	if (this->prior) this->prior->next = this->next;
	if (this->next) this->next->prior = this->prior;
	if (!this->prior) list->first = this->next;
	if (!this->next) list->last = this->prior;
	this->list->length--;
}

MDParticleListEntry* MDParticleListEntry::getNext(){
	return next;
}

MDParticleListEntry* MDParticleListEntry::getPrior(){
	return prior;
}

MDParticle* MDParticleListEntry::getThis(){
	return thisParticle;
}

bool MDParticleListEntry::isFirst(){
	return !prior;
}

bool MDParticleListEntry::isLast(){
	return !next;
}

MDSim::MDSim(double dt_,int N,int NumSteps) :
numSteps(NumSteps),dt(dt_), t(0),particles(new MDParticleList()),gaz(0)
,pause(false),radial(0, 201), histogramResolution(201),
directional(201, 201),graphDataFirst(0),graphDataLast(0), periodic(false), simBox(0, 2)
{
	this->gaz = new Gaz(this,N);
	this->r = mat(this->gaz->N,2);
	this->v = mat(this->gaz->N,2);
	this->a = mat(this->gaz->N,2);
	this->F = mat(this->gaz->N,2);
	//this->thermo = new AndersonThermostat(this, 0, 10);
};

void MDSim::initSim(bool periodic_, vec simBox_, int histogramResolution_ = 201, double histogramLength_ = 0.5){
	// restore starting conditions
	// **************************************************
	int i, j, k;
	this->histogramResolution = histogramResolution_;
	this->histogramLength = histogramLength_;
	this->simBox = simBox_;

	MDParticleListEntry* entry;
	
	this->periodic = periodic_;
	this->t = 0;
	this->radial = vec(0, this->histogramResolution);
	this->directional = mat(this->histogramResolution, this->histogramResolution);
	//init pot energy
	this->ePot = (this->ePotMin = this->refreshVerletLists(true, true));
	// calculate kinetic energy and total momentum
	// vec p_ges(0, this->dim);
	this->eKin = 0;
	for (entry = this->particles->getFirst(); entry; entry = entry->getNext()){
		// p_ges += entry->getThis()->v;
		this->eKin += entry->getThis()->v * entry->getThis()->v;
	}
	this->eKin *= 0.5;

	this->resetGraphs();
	this->updateGraphs();
	this->pause = true;
}

void MDSim::velocityVerletStep(bool countRadial){
	if (!this->pause) {
	}
}

double MDSim::refreshVerletLists(bool calc, bool countRadial){
	int i;
	double r_abs, pot_ = 0;
	vec r_, f_;
	MDParticleListEntry* entry;
	if (calc) for (entry = this->particles->getFirst(); entry; entry = entry->getNext()) entry->getThis()->a *= 0;
	for (entry = this->particles->getFirst(); entry; entry = entry->getNext()){
		entry->getThis()->verletList->clear();
		MDParticleListEntry* verlet;
		for (verlet = entry->getNext(); verlet; verlet = verlet->getNext()) {
			r_ = entry->getThis()->r - verlet->getThis()->r;
			// nearest image convention:
			if (this->periodic) for (i = 1; i <= 2; i++) if (abs(r_[i]) > 0.5 * this->simBox[i]) r_[i] -= copysign(this->simBox[i], r_[i]);
			r_abs = r_.v_abs();
			if (countRadial) {
				double hr = this->histogramResolution;
				double hl = this->histogramLength;
				if (r_abs < hl) this->radial[(int)(r_abs/hl * hr) + 1]++;
				if ((abs(r_[1]) < hl) && (abs(r_[2]) < hl)){
					this->directional
						[(int)((r_[1] / (2.0 * hl) + 0.5)*hr) + 1]
						[(int)((r_[2] / (2.0 * hl) + 0.5)*hr) + 1]++;
					this->directional
						[(int)((r_[1] / (-2.0 * hl) + 0.5)*hr) + 1]
						[(int)((r_[2] / (-2.0 * hl) + 0.5)*hr) + 1]++;
				}
			}
			if (this->r_verlet > 0) if (r_abs > this->r_verlet) continue;
			entry->getThis()->verletList->addParticle(verlet->getThis());
			if (calc) if ((r_abs < this->r_inter) || r_inter == 0){
				f_ = this->f(r_);
				pot_ += this->pot(r_);
				entry->getThis()->a += f_; verlet->getThis()->a -= f_;
			}
		}
	}
	return pot_;
}

double MDSim::velocityVerletForce(){
	int i;
	double pot_ = 0, r_inter_2 = this->r_inter*this->r_inter;
	vec f_, r_;
	MDParticleListEntry* entry;
	for (entry = this->particles->getFirst(); entry; entry = entry->getNext()) entry->getThis()->a *= 0; // reset force
	// address every pair of particles only once:
	for (entry = this->particles->getFirst(); entry; entry = entry->getNext()) {
		MDParticleListEntry* verlet;
		for (verlet = entry->getThis()->verletList->getFirst(); verlet; verlet = verlet->getNext()){
			r_ = entry->getThis()->r - verlet->getThis()->r;
			if (this->periodic) for (i = 1; i <= 2; i++) if (abs(r_[i]) > 0.5 * this->simBox[i]) r_[i] -= copysign(this->simBox[i], r_[i]);
			if (r_inter_2 > 0) if ((r_*r_) > r_inter_2) continue;
			f_ = this->f(r_);
			pot_ += this->pot(r_);
			entry->getThis()->a += f_; verlet->getThis()->a -= f_;
		}
	}
	return pot_;
}


// has to be called after the velocityVerletStep function to extend the graphs by the last time step
void MDSim::updateGraphs(){
	if (this->pause) return;
	this->graphDataLast = (this->graphDataLast->next = new graphData_t);
  
	/*
	for (i=0; i<sim->numSteps+1; i++) {

    // This updates the positions and velocities using Newton's Laws
    // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
    // which is a Kinetic Theory of gasses concept of Pressure
    	Press = sim->VelocityVerlet(dt, i+1, sim->getGaz()->tfp);
    	Press *= PressFac;

    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Now we would like to calculate somethings about the system:
    //  Instantaneous mean velocity squared, Temperature, Pressure
    //  Potential, and Kinetic Energy
    //  We would also like to use the IGL to try to see if we can extract the gas constant
    	mvs = MeanSquaredVelocity();
    	KE = Kinetic();
    	PE = Potential();

    	// Temperature from Kinetic Theory
    	Temp = m*mvs/(3*kB) * TempFac;

    	// Instantaneous gas constant and compressibility - not well defined because
    	// pressure may be zero in some instances because there will be zero wall collisions,
    	// pressure may be very high in some instances because there will be a number of collisions
    	gc = NA*Press*(Vol*VolFac)/(N*Temp);
    	Z  = Press*(Vol*VolFac)/(N*kBSI*Temp);

    	Tavg += Temp;
    	Pavg += Press;

    	fprintf(ofp,"  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
  	}

  	// Because we have calculated the instantaneous temperature and pressure, 
  	// we can take the average over the whole simulation here
  	sim->getGaz()->Pavg /= sim->numSteps;
	sim->getGaz()->Tavg /= sim->numSteps;
	sim->getGaz()->Z = sim->getGaz()->Pavg*(sim->getGaz()->Vol*sim->getGaz()->VolFac)/(sim->getGaz()->N*sim->getGaz()->kBSI*(sim->getGaz()->Tavg));
  	sim->getGaz()->gc = sim->getGaz()->NA*sim->getGaz()->Pavg*(sim->getGaz()->Vol*sim->getGaz()->VolFac)/(sim->getGaz()->N*sim->getGaz()->Tavg);
	*/

	// The first list entry (graphDataFirst) contains the maximum values of the graphs (for normalization).
	// Assignment of new values and query for wether they are greater than last maximum can be combined:
	if ((this->graphDataLast->ePot = this->ePot) > this->graphDataFirst->ePot) this->graphDataFirst->ePot = this->graphDataLast->ePot;
	if ((this->graphDataLast->eKin = this->eKin) > this->graphDataFirst->eKin) this->graphDataFirst->eKin = this->graphDataLast->eKin;
	this->graphDataLast->next = 0;
	if (this->ePot < this->ePotMin) this->ePotMin = this->ePot;
}

void MDSim::resetGraphs(){
	graphData_t* graphData = this->graphDataFirst;
	if (graphData){
		graphData->ePot = (graphData->eKin = 0);
		graphData = graphData->next;
		while (graphData){
			graphData_t* temp = graphData;
			graphData = graphData->next;
			delete temp;
		}
		this->graphDataLast = this->graphDataFirst;
	}
	else{
		this->graphDataFirst = (this->graphDataLast = new graphData_t);
		this->graphDataFirst->eKin = (this->graphDataFirst->ePot = 0);
		this->graphDataFirst->next = 0;
	}
}

void MDSim::resetRadialDistribution(){
	this->radial *= 0;
}

vec MDSim::getRadialDistribution(){
	vec result = this->radial;
	result.normalize();
	return result;
}

int MDSim::getHistogramResolution(){
	return this->histogramResolution;
}

void MDSim::resetDirectionalDistribution(){
	this->directional *= 0;
}

mat MDSim::getDirectionalDistribution(){
	return this->directional / this->directional.get_max();
}

double MDSim::getDt(){
	return this->dt;
}

double MDSim::getT(){
	return this->t;
}

double MDSim::getEPotMin(){
	return this->ePotMin;
}

Gaz* MDSim::getGaz(){
	return this->gaz;
}

void MDSim::initialize() {
   int n, p, i, j, k;
   double pos;

   // Number of atoms in each direction
   n = int(ceil(pow(this->gaz->N, 1.0/3)));

   //  spacing between atoms along a given direction
   pos = this->gaz->L / n;
   
   //  index for number of particles assigned positions
   p = 0;
   //  initialize positions
   for (i=0; i<n; i++) {
     for (j=0; j<n; j++) {
       for (k=0; k<n; k++) {
         if (p<this->gaz->N) {

           r[p][0] = (i + 0.5)*pos;
           r[p][1] = (j + 0.5)*pos;
           r[p][2] = (k + 0.5)*pos;
         }
         p++;
       }
     }
   }

   // Call function to initialize velocities
   initializeVelocities();

   /***********************************************
 *   Uncomment if you want to see what the initial positions and velocities are
   printf("  Printing initial positions!\n");
   for (i=0; i<this->gaz->N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
   }

   printf("  Printing initial velocities!\n");
   for (i=0; i<this->gaz->N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
   }
   */
 


}   


//  Function to calculate the averaged velocity squared
double MDSim::MeanSquaredVelocity() { 

  double vx2 = 0;
  double vy2 = 0;
  double vz2 = 0;
  double v2;

  for (int i=0; i<this->gaz->N; i++) {
  
    vx2 = vx2 + v[i][0]*v[i][0];
    vy2 = vy2 + v[i][1]*v[i][1];
    vz2 = vz2 + v[i][2]*v[i][2];

  }
  v2 = (vx2+vy2+vz2)/this->gaz->N;


  //printf("  Average of x-component of velocity squared is %f\n",v2);
  return v2;
}

//  Function to calculate the kinetic energy of the system
double MDSim::Kinetic() { //Write Function here!  

  double v2, kin;

  kin =0.;
  for (int i=0; i<this->gaz->N; i++) {

    v2 = 0.;  
    for (int j=0; j<3; j++) {

      v2 += v[i][j]*v[i][j];
    
    }
    kin += this->gaz->m*v2/2.;

  }

  //printf("  Total Kinetic Energy is %f\n",this->gaz->N*mvs*m/2.);
  return kin;

}


// Function to calculate the potential energy of the system
double MDSim::Potential() {
  double quot, r2, rnorm, term1, term2, Pot;
  int i, j, k;

  Pot=0.;
  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<this->gaz->N; j++) {

      if (j!=i) {
        r2=0.;
        for (k=0; k<3; k++) {
          r2 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
        }
        rnorm=sqrt(r2);
        quot=this->gaz->sigma/rnorm;
        term1 = pow(quot,12.);
        term2 = pow(quot,6.);

        Pot += 4*this->gaz->epsilon*(term1 - term2);

      }
    }
  }

  return Pot;
}



//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
void MDSim::computeAccelerations() {
  int i, j, k;
  double f, rSqd;
  double rij[3]; // position of i relative to j


  for (i = 0; i < this->gaz->N; i++) {  // set all accelerations to zero
    for (k = 0; k < 3; k++) {
      a[i][k] = 0;
    }
  } 
  for (i = 0; i < this->gaz->N-1; i++) {   // loop over all distinct pairs i,j
    for (j = i+1; j < this->gaz->N; j++) {
      // initialize r^2 to zero
      rSqd = 0;

      for (k = 0; k < 3; k++) {
        //  component-by-componenent position of i relative to j
        rij[k] = r[i][k] - r[j][k];
        //  sum of squares of the components
        rSqd += rij[k] * rij[k];
      }

      //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
      f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
      for (k = 0; k < 3; k++) {
        //  from F = ma, where m = 1 in natural units!
        a[i][k] += rij[k] * f;
        a[j][k] -= rij[k] * f;
      }
    }
  }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double MDSim::VelocityVerlet(double dt, int iter) {
  int i, j, k;
 
  double psum = 0.;

  //  Compute accelerations from forces at current position
  computeAccelerations();
  //  Update positions and velocity with current velocity and acceleration
  //printf("  Updated Positions!\n");
  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<3; j++) {
      r[i][j] += v[i][j]*dt + 0.5*a[i][j]*dt*dt;

      v[i][j] += 0.5*a[i][j]*dt;
    }
    //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
  }
  //  Update accellerations from updated positions
  computeAccelerations();
  //  Update velocity with updated acceleration
  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<3; j++) {
      v[i][j] += 0.5*a[i][j]*dt;
    }
  }

  // Elastic walls
  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<3; j++) {
       if (r[i][j]<0.) {
         v[i][j] *=-1.; //- elastic walls
         psum += 2*this->gaz->m*fabs(v[i][j])/dt;  // contribution to pressure from "left" walls
       }
       if (r[i][j]>=this->gaz->L) {
         v[i][j]*=-1.;  //- elastic walls
         psum += 2*this->gaz->m*fabs(v[i][j])/dt;  // contribution to pressure from "right" walls
       }
    }
  }

  for (i=0; i<this->gaz->N; i++) {
    printf("%s",this->gaz->atype);
    for (j=0; j<3; j++) {
      printf("  %12.10e ",r[i][j]);
    }
    printf("\n");
  }
  //fprintf(fp,"\n \n");

  return psum/(6*this->gaz->L*this->gaz->L);
}


void MDSim::initializeVelocities() {

  int i, j;

  for (i=0; i<this->gaz->N; i++) {

    for (j=0; j<3; j++) {
      //  Pull a number from a Gaussian Distribution
      v[i][j] = gaussdist();

    }
  }

  // Vcm = sum_i^this->gaz->N  m*v_i/  sum_i^this->gaz->N  M
  // Compute center-of-mas velocity according to the formula above
  double vCM[3] = {0, 0, 0};

  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<3; j++) {

      vCM[j] += this->gaz->m*this->v[i][j];

    }
  }

 
  for (i=0; i<3; i++) vCM[i] /= this->gaz->N*this->gaz->m;

  //  Subtract out the center-of-mass velocity from the
  //  velocity of each particle... effectively set the
  //  center of mass velocity to zero so that the system does
  //  not drift in space!
  for (i=0; i<this->gaz->N; i++) {
    for (j=0; j<3; j++) {

      v[i][j] -= vCM[j];

    }
  }

 //  Now we want to scale the average velocity of the system 
 //  by a factor which is consistent with our initial temperature, Tinit
 double vSqdSum, lambda;
 vSqdSum=0.;
 for (i=0; i<this->gaz->N; i++) {
   for (j=0; j<3; j++) {

     vSqdSum += v[i][j]*v[i][j];

   }
 }

 lambda = sqrt( 3*(this->gaz->N-1)*this->gaz->TempInit/vSqdSum);

 for (i=0; i<this->gaz->N; i++) {
   for (j=0; j<3; j++) {

     v[i][j] *= lambda;

   }
 }
}


//  Numerical recipes Gaussian distribution number generator
double MDSim::gaussdist() {
  static bool available = false;
  static double gset;
  double fac, rsq, v1, v2;
  if (!available) {
  do {
    v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
    v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
    rsq = v1 * v1 + v2 * v2;
  } while (rsq >= 1.0 || rsq == 0.0);

  fac = sqrt(-2.0 * log(rsq) / rsq);
  gset = v1 * fac;
  available = true;

  return v2*fac;
  } else {

  available = false;
  return gset;

}
}
