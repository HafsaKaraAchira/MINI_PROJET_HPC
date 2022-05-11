#ifndef GAZ_H
#define GAZ_H

class MDParticleList;
class MDSim;

class Gaz{

public:
	Gaz(MDSim*,int);
	~Gaz();
    void setT(double);
	double getT();
    void execute();

    // Number of particles
    int N;
    //  Lennard-Jones parameters in natural units!
    double sigma = 1.;
    double epsilon = 1.;
    double m = 1.;
    double kB = 1.;

    double NA = 6.022140857e23;
    double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

    //  Size of box, which will be specified in natural units
    double L;
    //  Initial Temperature in Natural Units
    double TempInit=41.17700;  //2;
    //density
    double rho=35000.0;
    // atom type
    char atype[10];
    // dt in natural units of time s.t. in SI it is 5 f.s.
    double dt ;
    //volume
    double Vol;
    // current temperature
    double Temp ;
    //current pressure
    double Press;
    //the average Pressure
    double Pavg ;
    //the average Temperature
    double Tavg ;
  	//unities factors
    double VolFac, TempFac, PressFac, timefac;
    //
    double KE, PE, mvs, gc, Z;

private:
    double temp;
	MDSim* sim;
};

#endif