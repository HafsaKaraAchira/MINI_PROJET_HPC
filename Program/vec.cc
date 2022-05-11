#include <cmath>
#include <iostream>
#include "vec.h"
using namespace std;


//vecteur null
vec::vec(){
	this->comp = 0;
	this->dim = 0;
}

//vecteur de dimension n sans initialisation
vec::vec(int n){ 
	this->comp = new double[n];
	this->dim = n;
}


//vecteur de dimension n avec initialisation
vec::vec(double* v, int n){
	int i;
	this->dim = n;
	this->comp = new double[n];
        if (!v) for (i = 0; i < n; i++) this->comp[i] = 0.0;
	else for (i = 0; i < n; i++) this->comp[i] = v[i];
}

//init avec un autre vecteur
vec::vec(const vec &other){
	int i;
	this->dim = other.dim;
	this->comp = new double[other.dim];
	for (i = 0; i < this->dim; i++) this->comp[i] = other.comp[i];
}



vec::~vec(){
	if (this->comp) delete[] this->comp;
}

//= operator
vec& vec::operator=(const vec &right){
	if (this == &right) return (*this);
	if (this->dim != right.dim){
		if (this->comp) delete[] this->comp;
		this->comp = new double[right.dim];
	}
	int i;
	for (i = 0; i < right.dim; i++) this->comp[i] = right.comp[i];
	this->dim = right.dim;
	return (*this);
}

//+= opertaor
vec& vec::operator+=(const vec &right){
	if (this->dim != right.dim){
		cerr << "vector addition failed: vectors of different dimension" << endl;
		return (*this);
	}
	int i;
	for (i = 0; i < right.dim; i++) this->comp[i] += right.comp[i];
	return (*this);
}

// - operator
vec& vec::operator-=(const vec &right){
	(*this) += (-1)*right;
	return (*this);
}

// * operator
vec& vec::operator*=(const double right){
	int i;
	for (i = 0; i < this->dim; i++) this->comp[i] *= right;
	return (*this);
}


// / operator
vec& vec::operator/=(const double right){
	(*this) *= (1 / right);
	return (*this);
}

// access vec[i]
double& vec::operator[](const int i){
    if ((this->dim < i)||(i<1)) cerr << "vector component index higher than vector dimension" << endl;
	return this->comp[i-1];
}

// get dim of vec
int vec::get_dim(){
	return this->dim;
}

// vec2
double vec::v_abs(){
	return sqrt((*this)*(*this));
}

// normalization
void vec::normalize(){
        (*this) /= this->v_abs();
}
