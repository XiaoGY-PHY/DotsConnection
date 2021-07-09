#ifndef  DotsConnection_H
#define  DotsConnection_H

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "MdcGeomSvc/MdcGeomSvc.h" 
#include "KalFitAlg/helix/Helix.h"
//#include "McTruth/CgemMcHit.h"
#include "MdcRawEvent/MdcDigi.h"
#include "RawDataProviderSvc/IRawDataProviderSvc.h"
#include "RawDataProviderSvc/RawDataProviderSvc.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "DotsConnection/DotsHelixFitter.h"
//#include "MdcUtilitySvc/IMdcUtilitySvc.h"
#include "MdcCalibFunSvc/MdcCalibFunSvc.h"

//using namespace KalmanFit;

class DotsConnection:public Algorithm {
	public:
		DotsConnection (const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();

	private:

		// --- properties ----------
		int myDebug; // flag for debug output
		int myNtProd; // if produce Ntuple
		double myDriftTimeUpLimit;
		double myMdcHitChi2Cut;

		// --- event T0
		double getEventStartTime();// in ns
		
		// --- get MC truth
		//bool  getMCTruth();
		KalmanFit::Helix getMCHelix();

		// --- get MC charged final states (e, mu, pi, K, p)
		void getMcFinalChargedStates();
		vector<int>        myVecMCTrkId;
		vector<int>        myVecPDG;
		vector<double>     myVecTrkLenFirstHalf;
		vector< vector<double> >  myVecHelix;
		//vector<Event::CgemMcHit*> myVecCgemMcHit;
		//vector<const RecCgemCluster*> myVecCgemXcluster;
		//vector<const RecCgemCluster*> myVecCgemVcluster;
		//vector<const RecCgemCluster*> myVecCgem1DCluster;
		//vector<int> myVecCgemXCluIdx[3][2];
		//vector<int> myVecCgemVCluIdx[3][2];
		vector<const MdcDigi*>        myVecMdcDigi;
		void resetFCSVec();
		void associateDigisToMcParticles();
		
		// --- get digi
		RawDataProviderSvc*           myRawDataProviderSvc;
		const MdcDigi*                myMdcDigiPointer[43][288];
		void  clearMdcDigiPointer();
		vector<const MdcDigi*>        getMdcDigiVec();
		//vector<const RecCgemCluster*> getCgemClusterVec(int view=0);// view = 0(X-view), 1(V-view), 2(X or V-view)

		// --- calibration svc
		MdcCalibFunSvc*           myMdcCalibFunSvc;
		//ICgemCalibFunSvc*         myCgemCalibSvc;

		// --- geometry of MDC
		MdcGeomSvc* myMdcGeomSvc;
		int    myNWire[43];
		int    myWireFlag[43];
		int    myOuterWire[43][288][2]; 
		int    myInnerWire[43][288][2]; 
		double myWirePhi[43][288]; 
		double myRLayer[43];
		
		// --- useful functions
		//IMdcUtilitySvc*  myMdcUtilitySvc;
		double dPhi(double phi1, double phi2);// --- delta phi (-pi, pi)

		// --- for CGEM clusters
		//bool sortCluster(const RecCgemCluster* clusterA , const RecCgemCluster* clusterB);

		// --- global circle/Helix fitter
		DotsHelixFitter myDotsHelixFitter;
		double myChi2CutDiverge;
		void testDotsHelixFitterAllHits();
		void testDotsHelixFitterPartHits();
		double myChiCut_circle;
		int    myNmaxDeact_circle;
		double myChiCut_helix;
		int    myNmaxDeact_helix;

		// --- RecMdcTrackCol
		RecMdcTrackCol* myRecMdcTrackCol;
		RecMdcHitCol*   myRecMdcHitCol;
		bool registerRecMdcTrack();

		// --- store RecMdcTrack from myDotsHelixFitter
		bool saveARecMdcTrack();

		// --- event information
		int myRun;
		int myEvt;
		
		// --- root items for myDotsHelixFitter
		NTuple::Tuple*        myNtHelixFitter;
		NTuple::Item<int>     myRUN;
		NTuple::Item<int>     myEVT;
		NTuple::Item<int>     myPID;
		NTuple::Item<int>     myNPar;
		NTuple::Array<double> myArrayHelixMC;
		NTuple::Array<double> myArrayHelixFitted;
		NTuple::Item<int>     myNHitsCircle;
		NTuple::Array<double> myLayerHitsCircle;
		NTuple::Array<double> myChiHitsCircle;
		NTuple::Item<int>     myNHits;
		NTuple::Array<double> myLayerHits;
		NTuple::Array<double> myChiHits;

};

#endif
