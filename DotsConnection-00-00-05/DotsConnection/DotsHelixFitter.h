#ifndef DotsHelixFitter_H
#define DotsHelixFitter_H

#include "KalFitAlg/helix/Helix.h"
#include "MdcRawEvent/MdcDigi.h"
//#include "CgemRecEvent/RecCgemCluster.h"
//#include "CgemGeomSvc/ICgemGeomSvc.h"
//#include "CgemGeomSvc/CgemGeomSvc.h"
#include "MdcGeomSvc/MdcGeomSvc.h" 
#include "MdcCalibFunSvc/MdcCalibFunSvc.h"
//#include "CgemCalibFunSvc/CgemCalibFunSvc.h"
#include "MdcUtilitySvc/IMdcUtilitySvc.h"

//using namespace KalmanFit;

class DotsHelixFitter 
{
	public: 
		DotsHelixFitter();
		~DotsHelixFitter();
		//DotsHelixFitter(KalmanFit::Helix iniHelix, vector<const MdcDigi*> vecMdcDigi, vector<const RecCgemCluster*> vecCgemCluster);
		DotsHelixFitter(KalmanFit::Helix iniHelix, vector<const MdcDigi*> vecMdcDigi);
		void initialize();
		void setMaxIterations(int n=10) {myMaxIteration=n;};
		void setMinXhitsInCircleFit(int n=3) {myMinXhitsInCircleFit=n;};
		void setMinVhitsInHelixFit(int n=2)  {myMinVhitsInHelixFit=n;};
		void setMinHitsInHelixFit(int n=5)   {myMinHitsInHelixFit=n;};
		void setDchi2Converge(double cut=5)  {myDchi2Converge=cut;};
		void setDchi2Diverge(double cut=50)  {myDchi2Diverge=cut;};
		void setMaxNChi2Increase(int n=2)    {myMaxNChi2Increase=n;};
		void setChi2Diverge(double cut=1000000)  {myChi2Diverge=cut;};

		void setInitialHelix(KalmanFit::Helix aHelix);
		void setT0(double T0) {myEventT0=T0;};
		//void setCgemClusters(vector<const RecCgemCluster*> aVecCgemCluster);
		void setDChits(vector<const MdcDigi*> aVecMdcDigi, double T0);
		
		void fitCircleOnly(bool x=true) {myFitCircleOnly=x;};
		void fitModeCircle() {myFitCircleOnly=true;}
		void fitModeHelix()  {myFitCircleOnly=false; myUseAxialHitsOnly=false;}
		void useAxialHitsOnly(bool x=true) {myUseAxialHitsOnly=x;}

		KalmanFit::Helix getClassHelix() {return *myHelix;} // in cm, radian, GeV/c
		HepVector getHelix() {return myHelix_aVec;};
		const HepSymMatrix & getEa() {return myHelix_Ea;};// error matrix
		double    getChi2()  {return myChi2;};

		void fit();
		
		void calculateInitialHelix();
		void calcuMeasDeriv();
		int  calculateNewHelix();
		int  deactiveHits(double chi_cut=10, int nMax=1);
		int  activeHits(double chi_cut=10);
		int  deactiveHits(int layer_max=43, int nHit_max=1000);
		int  getNActiveHits() {return myNActiveHits;};
		void updateChi2();

		void loadOneDcDigi(const MdcDigi* aDcDigi);
		void calculateDocaFromTrk(const MdcDigi* aDcDigi);
		void updateDcDigiInfo(const MdcDigi* aDcDigi);
		double getFlightLength() { return myFlightLength; };
		double getFlightLength(const MdcDigi* aDcDigi) { calculateDocaFromTrk(aDcDigi); return myFlightLength; };
		int    getLayer(const MdcDigi* aDcDigi) { loadOneDcDigi(aDcDigi); return myLayer; };
		double getDocaFromTrk() {return myDocaFromTrk;};
		double getDocaFromTrk(const MdcDigi* aDcDigi) { calculateDocaFromTrk(aDcDigi); return myDocaFromTrk;};
		
		RecMdcHit makeRecMdcHit(const MdcDigi* aDcDigi);
		RecMdcHit makeRecMdcHit(const MdcDigi* aDcDigi, KalmanFit::Helix aHelix);
		vector<RecMdcHit> makeRecMdcHitVec(int sel=1);// sel=0: all, 1: only active

		vector<const MdcDigi*> getVecMdcDigi() {return myVecMdcDigi;}
		vector<double> getVecChiMdcDigi() {return myChiVecMdcDigi;}
		vector<int>    getVecMdcDigiIsAct() {return myMdcDigiIsActive;}

		//vector<const RecCgemCluster*> getVecCgemCluster() { return myVecCgemCluster;}
		//vector<int> getVecIsActiveCgemCluster() { return myCgemClusterIsActive;}
		//vector<double> getVecChiCgemCluster()   { return myChiVecCgemCluster;}
		
		double IntersectCylinder(double r);
		HepMatrix dxda_cgem(KalmanFit::Helix a, double phi);

		//double getRmidGapCgem(int i) 
		//{
		//	if(i>2) i=2;
		//	if(i<0) i=0;
		//	return myRmidDGapCgem[i];
		//}

	private:

		// --- if initialized
		bool                      myIsInitialized;
		
		// --- option
		bool                      myFitCircleOnly;
		bool                      myUseAxialHitsOnly;

		// --- criteria
		int                       myMaxIteration;
		int                       myMinXhitsInCircleFit;
		int                       myMinVhitsInHelixFit;
		int                       myMinHitsInHelixFit;
		double                    myDchi2Converge;
		double                    myDchi2Diverge;
		int                       myMaxNChi2Increase;
		double                    myChi2Diverge;

		// --- some common information at event level
		double                    myEventT0;// event start time

		// --- information of a DC digi (or CGEM cluster) under proccessing
		int    myLayer;
		int    myWire;
		double myWirePos[7];
		double myCharge;
		double myPosOnWire[3];
		double myPosOnTrk[3];
		int    myLeftRight;//left:0, right:1
		double myDocaFromTrk;
		double myEntranceAngle;
		double myFlightLength;
		double myDocaFromDigi;
		double myDriftDist[2];//left/right
		double myDriftDistErr[2];//left/right
		double myDriftTime;
		double myDcChi;


		// --- track parameters
		KalmanFit::Helix*         myHelix;// in cm, radian, GeV/c
		double                    myHelix_a[5];// in cm, radian, GeV/c
		HepVector                 myHelix_aVec;
		HepSymMatrix              myHelix_Ea;// error matrix
		double                    myChi2;
		int                       myNActiveHits;


		// --- digi collection for DC
		vector<const MdcDigi*>    myVecMdcDigi;
		vector<double>            myDelDVecMdcDigi;
		vector< vector<double> >  myDerivVecMdcDigi;
		vector<double>            myChiVecMdcDigi;
		vector<int>               myAmbiguityMdcDigi;// left: -1, right: 1, not sure: 0
		vector<int>               myMdcDigiIsActive;// active: 1, inactive: 0
		int                       myNumMdcDigiPerLayer[43];
		int                       myIdxMdcDigiNeighbour[43][2];// index
		int                       myWireIdMdcDigiNeighbour[43][2];// wire id
		
		// --- cluster collection for CGEM
		//vector<const RecCgemCluster*>   myVecCgemCluster;
		//vector<int>                     myCgemClusterIsActive;// active: 1, inactive: 0
		//vector<double>                  myChiVecCgemCluster;
		//vector<const RecCgemCluster*>   myVecCgemClusterX;
		//vector<const RecCgemCluster*>   myVecCgemClusterV;

		// --- index ordered by flight length
		map<double, int>          myMapFlylenIdx;

		
		// --- geometry for MDC
		MdcGeomSvc*               myMdcGeomSvc;
		int                       myWireFlag[43];
		int                       myNWire[43];
		double                    myRWires[43];
		vector<double>            myPhiWires[43];
		vector<double>            myTensionWires[43];
		vector<double*>           myEastPosWires[43];
		vector<double*>           myWestPosWires[43];
		vector<double*>           myPosWires[43];// position  of wires with z=0
		vector<double*>           myDirWires[43];// direction of wires
		
		// --- geometry for CGEM
		//CgemGeomSvc*              myCgemGeomSvc;
		//double                    myRmidDGapCgem[3];
		//double                    myR2midDGapCgem[3];
		//double                    myRVCgem[3];
		//double                    myRXCgem[3];
		//double                    myAngStereoCgem[3];
		//int                       myNSheets[3];

		// --- calibration svc
		MdcCalibFunSvc*           myMdcCalibFunSvc;
		//ICgemCalibFunSvc*         myCgemCalibSvc;

		// --- MDC utility
		IMdcUtilitySvc*  myMdcUtilitySvc;

};

#endif
