#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "vec.h"
#include "MDSim.h"
#include "Gaz.h"
#include <cmath>
#include <cstdint>
#include <iostream>
//using namespace std;

class MDParticleList;
class MDParticleListEntry;

Gaz::Gaz(MDSim* sim_,int N_) : sim(sim_),N(N_){};

void Gaz::setT(double temp_){
	this->temp = (-2.0)*temp_;
}

double Gaz::getT(){
	return (-0.5)*this->temp;
}

void Gaz::execute(){

	if (this->sim->pause) return;

	MDParticleListEntry* entry;
	for (entry = this->sim->particles->getFirst(); entry; entry = entry->getNext()) {
		
	}
}