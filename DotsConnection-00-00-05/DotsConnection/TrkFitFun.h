#ifndef TRKFITFUN_H
#define TRKFITFUN_H

#include <string>
#include <vector>

namespace TrkFitFun
{
	/* exterior product */
	void expd(double veca[3], double vecb[3], double val[3]);

	/* distance of two lines */
	int dist2Line(double sta[3], double stb[3], double veca[3], double vecb[3],double &d, double xyz_a[3], double xyz_b[3], int fgZcal = 1);

	/* doca of helix track and wire(line) */
	double docaHelixWire(double trkpar[], double wirest[], double wirev[], double &zwire, double zini);
	bool   getDoca(double trkpar[], double wpos[], double &doca,double whitPos[], double zini);
	double getPhiIni(double trkpar[], double rLayer, double pos[]);

	/* derivative of doca to track paramters */
	bool getDeriLoc(int ipar, double helix[], double &deri, double wpos[], double zini);
	
	/* number of iteration in doca calculation */
	extern int gNiter;
	
	const double CC   = 2.99792458E10;       // cm/sec, light velocity
	const double PI   = 3.141592653;
	const double PI2  = 6.283185307;
	const double HFPI = 1.570796327;
	const double BFIELD = 1.0;	/* Tesla */
	
	const int NTRKPAR = 5;	/* number of track parameters */
	
	const int gNsamLC = 100;
	const double gStepLC[5] = {0.001, 0.001, 0.00001, 0.0001, 0.0001};	/* units of dr&dz are cm */
}

#endif
