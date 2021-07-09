#include <math.h>
#include <iostream>
#include "DotsConnection/TrkFitFun.h"
#include "TSpline.h"
#include "KalFitAlg/helix/Helix.h"

using namespace std;
//using namespace CLHEP;

int TrkFitFun::gNiter;

void TrkFitFun::expd(double veca[3], double vecb[3], double val[3]){
     val[0] = veca[1]*vecb[2] - veca[2]*vecb[1];
     val[1] = veca[2]*vecb[0] - veca[0]*vecb[2];
     val[2] = veca[0]*vecb[1] - veca[1]*vecb[0];
}

int TrkFitFun::dist2Line(double sta[3], double stb[3],
			 double veca[3], double vecb[3],
			 double &d, double xyz_a[3], double xyz_b[3], int fgZcal){
     int i;
     double vst[3];   // P0 - W0
     double vp[3];    // (P * W) / |P * W|
     double modvp;
     double m2;
     double seca[3], secb[3];
     double tpa = 0.0;
     double twb = 0.0;

     for(i=0; i<3; i++) vst[i] = sta[i] - stb[i];  // vector P0-W0

     //cout<<"TrkFitFun::dist2Line() veca=("<<veca[0]<<", "<<veca[1]<<", "<<veca[2]<<")"<<endl;
     //cout<<"TrkFitFun::dist2Line() vecb=("<<vecb[0]<<", "<<vecb[1]<<", "<<vecb[2]<<")"<<endl;
     TrkFitFun::expd(veca, vecb, vp); // exterior product

     //cout<<"TrkFitFun::dist2Line() vp=("<<vp[0]<<", "<<vp[1]<<", "<<vp[2]<<")"<<endl;
     m2 = vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2];
     //cout<<"TrkFitFun::dist2Line() m2="<<m2<<endl;
     modvp = sqrt(m2);
     for(i=0; i<3; i++) vp[i] /= modvp;  // (P * W) / |P * W|

     d = 0.0;
     for(i=0; i<3; i++)  d += vst[i] * vp[i];
     //cout<<"TrkFitFun::dist2Line() d="<<d<<endl;

     if( 0 == fgZcal ) return 1;

     double veca00 = veca[0]*veca[0];
     double veca11 = veca[1]*veca[1];
     double veca22 = veca[2]*veca[2];

     double veca01 = veca[0]*veca[1];
     double veca02 = veca[0]*veca[2];
     double veca12 = veca[1]*veca[2];

     double vecb00 = vecb[0]*vecb[0];
     double vecb11 = vecb[1]*vecb[1];
     double vecb22 = vecb[2]*vecb[2];
     double vecb01 = vecb[0]*vecb[1];
     double vecb02 = vecb[0]*vecb[2];
     double vecb12 = vecb[1]*vecb[2];

     seca[0] = veca[1]*vecb01 + veca[2]*vecb02 - veca[0]*(vecb11 + vecb22);
     seca[1] = veca[0]*vecb01 + veca[2]*vecb12 - veca[1]*(vecb00 + vecb22);
     seca[2] = veca[0]*vecb02 + veca[1]*vecb12 - veca[2]*(vecb00 + vecb11);

     secb[0] = vecb[1]*veca01 + vecb[2]*veca02 - vecb[0]*(veca11 + veca22);
     secb[1] = vecb[0]*veca01 + vecb[2]*veca12 - vecb[1]*(veca00 + veca22);
     secb[2] = vecb[0]*veca02 + vecb[1]*veca12 - vecb[2]*(veca00 + veca11);

     for(i=0; i<3; i++){
	  tpa += seca[i] * (sta[i] - stb[i]);
	  twb += secb[i] * (stb[i] - sta[i]);
     }
     tpa /= m2;
     twb /= m2;
     
     for(i=0; i<3; i++) 
     {
	     xyz_a[i] = veca[i] * tpa + sta[i];
	     xyz_b[i] = vecb[i] * twb + stb[i];
     }

     return 1;
}


double TrkFitFun::docaHelixWire(double trkpar[], double wirest[], double wirev[], double &zwire, double phiIni)
{
	HepPoint3D pivot(0,0,0);
	HepVector  a(5);
	for(int i=0; i<5; i++) a(i+1)=trkpar[i];
	KalmanFit::Helix aHelix(pivot, a);
	double x0 = trkpar[0] * cos(trkpar[1]);
	double y0 = trkpar[0] * sin(trkpar[1]);
	double z0 = trkpar[3];
	double phi0 = trkpar[1] + HFPI;
	if(phi0 > PI2) phi0 -= PI2;
	double g = 1e13 / (CC * BFIELD * trkpar[2]); // alpha/kappa in cm/(GeV/c)
	double tanl = trkpar[4];

	double trkst[3];// point on track
	double trkv[3];
	//double phi = phi0 + (zini - z0) / (g * tanl);
	double phi = phiIni;
	double dphi;

	double doca;
	double ztrk;
	double xyz_trk[3];
	double xyz_wire[3];
	double phiNew;
	int iter = 0;
	//for(iter=0; iter<10; iter++){
	//    trkst[0] = x0 + g * (sin(phi) - sin(phi0));
	//    trkst[1] = y0 + g * (-cos(phi) + cos(phi0));
	//    trkst[2] = z0 + g * tanl * (phi - phi0);

	//    trkv[0] = cos(phi);
	//    trkv[1] = sin(phi);
	//    trkv[2] = tanl;

	//    Alignment::dist2Line(trkst, wirest, trkv, wirev, doca, ztrk, zwire);

	//    phiNew = phi0 + (ztrk - z0) / (g*tanl);
	//    if(fabs(phiNew - phi) < 1.0E-8) break;
	//    phi = phiNew;
	//}
	for(iter=0; iter<10; iter++){
		dphi = phi - phi0;
		if(dphi > PI) dphi -= PI2;
		if(dphi < -PI) dphi += PI2;

		trkst[0] = x0 + g * (sin(phi) - sin(phi0));
		trkst[1] = y0 + g * (-cos(phi) + cos(phi0));
		// 	  trkst[2] = z0 + g * tanl * (phi - phi0);
		trkst[2] = z0 + g * tanl * dphi;

		trkv[0] = cos(phi);
		trkv[1] = sin(phi);
		trkv[2] = tanl;
		//cout<<"TrkFitFun::docaHelixWire: trkv = "<<setw(10)<<trkv[0]<<setw(10)<<trkv[1]<<setw(10)<<trkv[2]<<endl;

		// wirest: coordinate of the point on wire according to zc
		dist2Line(trkst, wirest, trkv, wirev, doca, xyz_trk, xyz_wire);

		//ztrk=xyz_trk[2];
		//phiNew = phi0 + (ztrk - z0) / (g*tanl);
		//cout<<"phiNew = "<<phiNew<<endl;

		HepPoint3D newPivot(xyz_trk[0],xyz_trk[1],xyz_trk[2]);
		aHelix.pivot(newPivot);
		phiNew=aHelix.phi0()+HFPI;
		//cout<<"phiNew2 = "<<phiNew<<endl;
		
		double delPhi = phiNew - phi;
		while(delPhi>M_PI)  delPhi-=M_PI*2.0;
		while(delPhi<-M_PI) delPhi+=M_PI*2.0;
		if(fabs(delPhi) < 1.0E-8) break;
		phi = phiNew;
	}
	zwire=xyz_wire[2];

	gNiter = iter;

	return doca;
}

bool TrkFitFun::getDoca(double trkpar[], double wpos[], double &doca,
			double whitPos[], double zini){
     int i = 0;
     double zp;      // z of the point above in the plane of the wire
     double xyz[3];  // coordinate of the point on wire according to zc
     double dxyz[3]; // orientation of the tangent line at the point above

     double ten = wpos[6];
     double a = 9.47e-05 / (2 * ten); // a = density(g/mm)/2T(g)
     double dx = wpos[0] - wpos[3]; // the differential of xyz between the end planes
     double dy = wpos[1] - wpos[4];
     double dz = wpos[2] - wpos[5]; // 
     double length = sqrt(dx*dx + dz*dz);
     //cout<<"length = "<<length<<endl;

     double ztan = 0.0;  // z of the doca point in the tangent line
     if(whitPos[2] < 0.5*length)  ztan = whitPos[2];

     double zc=0.0;  // z of the calculated point of the wire
     //if( Alignment::gFlagMag ) zc = zini;

     // alf is the angle between z and the projection of the wire on xz
     double sinalf = dx / sqrt(dx*dx + dz*dz);
     double cosalf = dz / sqrt(dx*dx + dz*dz);
     double tanalf = dx / dz;
     //cout<<"cosalf="<<cosalf<<endl;

     /*
     double posIni[3];
     double rLayer = sqrt((wpos[3] * wpos[3]) + (wpos[4] * wpos[4]));
     double phiIni = getPhiIni(trkpar, rLayer, posIni);
     //cout<<"TrkFitFun::getDoca(): phiIni = "<<phiIni<<endl;
     */

     HepPoint3D pivot(0,0,0);
     HepVector  par(5);
     for(int i=0; i<5; i++) par(i+1)=trkpar[i];
     KalmanFit::Helix aHelix(pivot, par);
     double tmp = (ztan-wpos[5])/dz;
     pivot.setX(tmp*dx+wpos[3]);
     pivot.setY(tmp*dy+wpos[4]);
     pivot.setZ(tmp*dz+wpos[5]);
     aHelix.pivot(pivot);
     double phiIni=aHelix.phi0()+HFPI;
     //cout<<"TrkFitFun::getDoca(): phiIni2 = "<<phiIni<<endl;

     if(dz < 20){
	  std::cout << "ERROR: wire position error in getdocaLine() !!!"
		    << std::endl;
	  std::cout<<"dz = "<<dz<<std::endl;
	  return false;
     }

     while( 1 ){
	  i++;
	  if(i > 5){
	       return false;
	  }
	  zp = zc / cosalf;

	  xyz[0] = (zc - wpos[5]) * tanalf + wpos[3];
	  xyz[1] = a*zp*zp + (wpos[1] - wpos[4])*zp/length
	       + 0.5*(wpos[1] + wpos[4]) - a*length*length/4.0;
	  xyz[2] = zc;

	  dxyz[0] = sinalf;
	  dxyz[1] = 2.0 * a * zp + (wpos[1] - wpos[4]) / length;
	  dxyz[2] = cosalf;

	  //if( Alignment::gFlagMag ) doca = docaHelixWire(trkpar, xyz, dxyz, ztan, phiIni);
	  //else doca = docaLineWire(trkpar, xyz, dxyz, ztan);
	  doca = docaHelixWire(trkpar, xyz, dxyz, ztan, phiIni);
	  //cout<<"TrkFitFun::getDoca(): i, doca, ztan, zc = "<<i<<", "<<doca<<", "<<ztan<<", "<<zc<<endl;

	  if( fabs(zc-ztan) < 0.5 )  break;
	  /*else if( fabs(ztan) > (0.5*length) ){
	       doca = 99999.0;
	       break;
	  }
	  */ // --- comment out by wangll, 2020-04-27
	  zc = ztan;
     }
     whitPos[2] = ztan;
     zp = ztan / cosalf;
     whitPos[1] = a*zp*zp + (wpos[1] - wpos[4])*zp/length
	  + 0.5*(wpos[1] + wpos[4]) - a*length*length/4.0;
     whitPos[0] = (ztan - wpos[5]) * tanalf + wpos[3];

     return true;
}

double TrkFitFun::getPhiIni(double trkpar[], double rLayer, double pos[]){
     double dr = trkpar[0];
     double fi0 = trkpar[1];
     double kap = trkpar[2];
     double rw = rLayer;

     double phi0 = fi0 + HFPI;
     if(phi0 > PI2) phi0 -= PI2;
     double g = 1.0e13 / (CC * BFIELD * kap); // alpha/kappa in cm/(GeV/c)

     double aa = rw*rw - (dr-g)*(dr-g) - g*g;
     double bb = 2*g*(dr - g);
     //cout<<"TrkFitFun::getPhiIni() aa,bb="<<setw(10)<<aa<<setw(10)<<bb<<endl;
     double cc = aa/bb;
     double dd = M_PI;
     if(fabs(cc)<=1) 
	     dd=acos(cc);	// dd (0, PI)

     double phi;
     if(kap > 0) phi = phi0 + dd;
     else phi = phi0 - dd;
     //cout<<"TrkFitFun::getPhiIni() phi0,dd="<<setw(10)<<phi0<<setw(10)<<dd<<endl;

     if(phi > PI2) phi -= PI2;
     if(phi < 0) phi += PI2;

     double x0 = dr * cos(fi0);
     double y0 = dr * sin(fi0);
     pos[0] = x0 + g * (sin(phi) - sin(phi0));
     pos[1] = y0 + g * (-cos(phi) + cos(phi0));
     //      pos[2] = trkpar[3] + g * trkpar[4] * (phi - phi0);
     if(kap > 0) pos[2] = trkpar[3] + g * trkpar[4] * dd;
     else pos[2] = trkpar[3] - g * trkpar[4] * dd;

     return phi;
}


bool TrkFitFun::getDeriLoc(int ipar, double helix[], double &deri, double wpos[], double zini)
{
	int i;
	double doca;
	double sampar[NTRKPAR];
	double xxLC[gNsamLC];
	double yyLC[gNsamLC];
	double aPointOnWire[3]={0,0,0};

	for(i=0; i<NTRKPAR; i++) sampar[i] = helix[i];
	double startpar = helix[ipar] - 0.5*gStepLC[ipar]*(double)gNsamLC;

	for(i=0; i<gNsamLC; i++){
		sampar[ipar] = startpar + (double)i * gStepLC[ipar];
		xxLC[i] = sampar[ipar];
		//if(0==ipar || 3==ipar) xxLC[i] *= 10.;	// cm -> mm
		
		getDoca(sampar, wpos, doca, aPointOnWire, zini);

		if(NULL == doca){
			// 	     cout << "in getDeriLoc, doca = " << doca << endl;
			return false;
		}
		yyLC[i] = doca;
	}

	double par = helix[ipar];
	//if (0==ipar || 3==ipar) par*= 10.; // cm -> mm
	TSpline3* pSpline3 = new TSpline3("deri", xxLC, yyLC, gNsamLC);
	deri = pSpline3->Derivative(par);
	delete pSpline3;
	return true;
}

