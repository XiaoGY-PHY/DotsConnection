#include "DotsConnection/DotsHelixFitter.h"
#include "DotsConnection/TrkFitFun.h"
//#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "RawEvent/RawDataUtil.h"
#include "Identifier/MdcID.h"
#include "MdcCalibFunSvc/IMdcCalibFunSvc.h"

using namespace TrkFitFun;


DotsHelixFitter::DotsHelixFitter():
	myHelix(NULL),myIsInitialized(false),myMdcGeomSvc(NULL),myCgemGeomSvc(NULL),myMdcCalibFunSvc(NULL),myMdcUtilitySvc(NULL)
{
	myHelix_Ea=HepSymMatrix(5,0);
	myMaxIteration=10;	
	myMinXhitsInCircleFit=3;
	myMinVhitsInHelixFit=2;
	myMinHitsInHelixFit=5;
	myDchi2Converge=5.0;
	myDchi2Diverge=50.0;
	myMaxNChi2Increase=2;
	myChi2Diverge=1000000.0;
	myFitCircleOnly=false;
	myUseAxialHitsOnly=false;
}


DotsHelixFitter::~DotsHelixFitter()
{
	if(myHelix) delete myHelix;
}

void DotsHelixFitter::initialize()
{
	if(myIsInitialized) return;

	// --- get CgemGeomSvc ---
	ISvcLocator* svcLocator = Gaudi::svcLocator();
	ICgemGeomSvc* ISvc;
	StatusCode Cgem_sc=svcLocator->service("CgemGeomSvc", ISvc);
	if(Cgem_sc.isSuccess()) myCgemGeomSvc=dynamic_cast<CgemGeomSvc *>(ISvc);
	else cout<<"DotsHelixFitter::initialize(): can not get CgemGeomSvc"<<endl;
	int nLayerCgem = myCgemGeomSvc->getNumberOfCgemLayer();
	if(nLayerCgem!=3) {
		cout<<"DotsHelixFitter::initialize(): nLayerCgem!=3 !!!"<<endl;
		return;
	}
	for(int i=0; i<nLayerCgem; i++)
	{
		//myCgemLayer[i]=myCgemGeomSvc->getCgemLayer(i);
		CgemGeoLayer* CgemLayer = myCgemGeomSvc->getCgemLayer(i);
		myRmidDGapCgem[i]=(CgemLayer->getMiddleROfGapD())/10.0;//cm
		myR2midDGapCgem[i]=pow(myRmidDGapCgem[i],2);
		myRVCgem[i]=(CgemLayer->getInnerROfAnodeCu1())/10.0;//cm
		myRXCgem[i]=(CgemLayer->getInnerROfAnodeCu2())/10.0;//cm
		myAngStereoCgem[i]=(CgemLayer->getAngleOfStereo())/180.0*CLHEP::pi;// degree -> radians
		myNSheets[i]=CgemLayer->getNumberOfSheet();
		//myRXstrips[i] = myCgemLayer[i]->getInnerROfAnodeCu2();
		//myRVstrips[i] = myCgemLayer[i]->getInnerROfAnodeCu1();
		//if(myPrintFlag) 
		cout<<"CGEM layer "<<i<<": "
			<<" Rmid="<<myRmidDGapCgem[i]<<" cm "
			<<" RV ="<<myRVCgem[i]<<" cm "
			<<" RX ="<<myRXCgem[i]<<" cm "
			<<" , stereo angle = "<<myAngStereoCgem[i]<<" radian "
			<<myNSheets[i]<<" sheets"
			//<<", is reversed "<<IsReverse
			//<<", RX="<<myRXstrips[i]<<", RY="<<myRVstrips[i]
			<<endl;
	}

	// --- get MDCGeomSvc ---
	IMdcGeomSvc* ISvc2;
	StatusCode mdc_sc=svcLocator->service("MdcGeomSvc", ISvc2);
	if(mdc_sc.isSuccess()) myMdcGeomSvc=dynamic_cast<MdcGeomSvc *>(ISvc2);
	else cout<<"DotsHelixFitter::initialize(): can not get MdcGeomSvc"<<endl;
	int nLayer= myMdcGeomSvc->getLayerSize();
	if(nLayer!=43) cout<<"DotsHelixFitter::initialize(): MDC nLayer = "<<nLayer<<" !=43  !!!"<<endl;
	for(int i=0; i<nLayer; i++)
	{
		const MdcGeoLayer* aLayer = myMdcGeomSvc->Layer(i);
		myRWires[i]=0.1*(aLayer->Radius());//  cm
		cout<<"MDC layer "<<i<<" R = "<<myRWires[i]<<endl;
		double nomShift = myMdcGeomSvc->Layer(i)->nomShift();
		if(nomShift>0)      myWireFlag[i]=1;
		else if(nomShift<0) myWireFlag[i]=-1;
		else                myWireFlag[i]=0;
		int nWires=aLayer->NCell();
		myNWire[i]=nWires;
		for(int j=0; j<nWires; j++)
		{
			const MdcGeoWire* aWire = myMdcGeomSvc->Wire(i,j);
			HepPoint3D aPointEast = 0.1 * (aWire->Forward());// in cm
			double* aPointArray = new double[3]{aPointEast.x(),aPointEast.y(),aPointEast.z()};// in cm
			myWestPosWires[i].push_back(aPointArray);
			HepPoint3D aPointWest = 0.1*(aWire->Backward());// in cm
			aPointArray = new double[3]{aPointWest.x(),aPointWest.y(),aPointWest.z()};// in cm
			myEastPosWires[i].push_back(aPointArray);
			HepPoint3D aPointMiddle = 0.5*(aPointEast+aPointWest);
			aPointArray = new double[3]{aPointMiddle.x(),aPointMiddle.y(),aPointMiddle.z()};// in cm
			myPosWires[i].push_back(aPointArray);
			//myPhiWires[i].push_back(aWire->Forward().phi());
			myPhiWires[i].push_back(aPointMiddle.phi());
			myTensionWires[i].push_back(aWire->Tension());;
		}
	}

	// --- MdcCalibFunSvc ---
	IMdcCalibFunSvc* imdcCalibSvc;                                                                                                      
	StatusCode sc = svcLocator->service("MdcCalibFunSvc", imdcCalibSvc);
	if ( sc.isSuccess() ){ 
		myMdcCalibFunSvc = dynamic_cast<MdcCalibFunSvc*>(imdcCalibSvc);
	}
	else {
		cout<<"DotsHelixFitter::initialize(): can not get MdcCalibFunSvc"<<endl;
	} 

	// --- get CgemCalibFunSvc ---
	sc = svcLocator->service("CgemCalibFunSvc", myCgemCalibSvc);
	if ( !(sc.isSuccess()) ) 
	{
		cout<<"DotsHelixFitter::initialize(): can not get CgemCalibFunSvc"<<endl;
	} 
	
	myIsInitialized=true;

	// --- get MdcUtilitySvc
	sc = svcLocator->service("MdcUtilitySvc", myMdcUtilitySvc);
	if ( !(sc.isSuccess()) ) 
	{
		cout<<"DotsHelixFitter::initialize(): can not get MdcUtilitySvc"<<endl;
	}

	return;
}


void DotsHelixFitter::setDChits(vector<const MdcDigi*> aVecMdcDigi, double T0)
{
	myVecMdcDigi=aVecMdcDigi;
	//cout<<"set a vector<const MdcDigi*> with "<<myVecMdcDigi.size()<<" hits"<<endl;
	int nDigi = myVecMdcDigi.size();
	myChiVecMdcDigi.resize(nDigi);
	myAmbiguityMdcDigi.resize(nDigi);
	myMdcDigiIsActive.resize(nDigi);
	myEventT0 = T0;
	
	for(int i=0; i<43; i++) myNumMdcDigiPerLayer[i]=0;
	
	for(int i=0; i<nDigi; i++)
	{
		myChiVecMdcDigi[i]=9999;
		myAmbiguityMdcDigi[i]=0;
		myMdcDigiIsActive[i]=1;
	}
}

void DotsHelixFitter::setCgemClusters(vector<const RecCgemCluster*> aVecCgemCluster)
{
	myVecCgemCluster=aVecCgemCluster;
	int nCluster1D = myVecCgemCluster.size();
	myChiVecCgemCluster.resize(nCluster1D);
	myCgemClusterIsActive.resize(nCluster1D);
	//cout<<"set a vector<RecCgemCluster*> with "<<myVecCgemCluster.size()<<" clusters"<<endl;
	myVecCgemClusterX.clear();
	myVecCgemClusterV.clear();
	vector<const RecCgemCluster*>::iterator iter_cluster = myVecCgemCluster.begin();
	for(int i=0; iter_cluster!=myVecCgemCluster.end(); iter_cluster++, i++)
	{
		int flag = (*iter_cluster)->getflag();
		myChiVecCgemCluster[i]=9999;
		myCgemClusterIsActive[i]=1;
	}
}

int DotsHelixFitter::calculateNewHelix()
{
	bool debug = false; //debug=true;
	
	int state = 0;

	double tension = 9999.;

	int nPar=5;
	if(myFitCircleOnly) nPar=3;
	
	// --- check hits
	//if(debug) cout<<"DotsHelixFitter::calculateNewHelix() starts checking hits ... "<<endl;
	for(int i=0; i<43; i++) myNumMdcDigiPerLayer[i]=0;
	myNActiveHits=0;
	int nXHits[3]={0,0,0}, nVHits[3]={0,0,0};
	int nDigi=myVecMdcDigi.size();
	double* vec_zini = new double[nDigi];
	int maxLayer = 0;
	int i_digi=0;
	vector<const MdcDigi*>::iterator iter_mdcDigi = myVecMdcDigi.begin();
	for(; iter_mdcDigi!=myVecMdcDigi.end(); iter_mdcDigi++, i_digi++)
	{
		// --- initialize vec_zini ---
		vec_zini[i_digi]=0;
		// --- get id, layer, wire ---
		Identifier id = (*iter_mdcDigi)->identify();
		int layer     = MdcID::layer(id);
		int wire      = MdcID::wire(id);
		// --- counting active hits
		if(myMdcDigiIsActive[i_digi])
		{
			myNActiveHits++;

			int hitFlag = myWireFlag[layer];
			if(hitFlag==0) 
			{
				nXHits[0]++;// total
				nXHits[1]++;// MDC
			}
			else 
			{
				nVHits[0]++;
				nVHits[1]++;
			}
		}
		if(maxLayer<layer) maxLayer=layer;
		int nHits = myNumMdcDigiPerLayer[layer];
		if(nHits<2) {
			myIdxMdcDigiNeighbour[layer][nHits]=i_digi;
			myWireIdMdcDigiNeighbour[layer][nHits]=wire;
		}
		myNumMdcDigiPerLayer[layer]++;
	}
	//if(debug) cout<<"DotsHelixFitter::calculateNewHelix() ends checking hits ... "<<endl;
	HepPoint3D farPoint = myHelix->x(M_PI);
	double maxRTrk = farPoint.perp();
	for(int i=0; i<43; i++) 
	{
		if(maxRTrk<myRWires[i]+2) break;// where soft track turning back
		if(myNumMdcDigiPerLayer[i]==2) 
		{
			int wire_1 = myWireIdMdcDigiNeighbour[i][0];
			int wire_2 = myWireIdMdcDigiNeighbour[i][1];
			int delta_n = abs(wire_1-wire_2);
			int ambi_1 = 0;
			int ambi_2 = 0;
			if(delta_n==1) 
			{
				ambi_1 = wire_1>wire_2? 1:-1;
				ambi_2 = wire_1>wire_2? -1:1;
			}
			else if(delta_n==myNWire[i]-1)
			{
				ambi_1 = wire_1<=1? 1:-1;
				ambi_2 = wire_2<=1? 1:-1;
			}
			if(debug)
				cout<<"layer, wire1, wire2, ambi1, ambi2 = "
					<<setw(5)<<i
					<<setw(5)<<wire_1
					<<setw(5)<<wire_2
					<<setw(5)<<ambi_1
					<<setw(5)<<ambi_2
					<<endl;
			// --- fix Ambiguity for neighboured mdc hits
			//myAmbiguityMdcDigi[myIdxMdcDigiNeighbour[i][0]]=ambi_1;
			//myAmbiguityMdcDigi[myIdxMdcDigiNeighbour[i][1]]=ambi_2;
		}
	}
	i_digi=0;
	vector<const RecCgemCluster*>::iterator iter_cluster = myVecCgemCluster.begin();
	for(; iter_cluster!=myVecCgemCluster.end(); iter_cluster++, i_digi++)
	{
		int flag  = (*iter_cluster)->getflag();
		if(flag==0) 
		{
			nXHits[0]++;// total
			nXHits[2]++;// CGEM
		}
		else if(flag==1)
		{
			nVHits[0]++;
			nVHits[2]++;
		}
		myNActiveHits++;
	}
	if(myFitCircleOnly) 
	{
		if(nXHits[0]<myMinXhitsInCircleFit) return 1;
	}
	else 
	{
		if((nXHits[0]+nVHits[0])<myMinHitsInHelixFit) return 2;
		if(nVHits[0]<myMinVhitsInHelixFit)            return 3;
	}




	// ---> start fit
	double lastTotChi2 = 999999;
	int nIterations = 0;
	int n_chi2_increase = 0;

	while(1) // iterations
	{
		myMapFlylenIdx.clear();
		// --- matrix U  P J ---
		HepSymMatrix U(nPar, 0);
		//HepMatrix    P(5,1,0);
		//HepMatrix    J(5,1,0), JT(1,5,0);
		//HepMatrix    J_dmdx(1,3,0), J_dxda(3,5,0);
		HepMatrix    P(nPar,1,0);
		HepMatrix    J(nPar,1,0), JT(1,nPar,0);
		HepMatrix    J_dmdx(1,3,0), J_dxda(3,nPar,0);

		// --- loop MDC hits ---
		double totChi2=0;
		i_digi=0;
		vector<const MdcDigi*>::iterator iter_mdcDigi = myVecMdcDigi.begin();
		for(; iter_mdcDigi!=myVecMdcDigi.end(); iter_mdcDigi++, i_digi++)
		{
			// --- get id, layer, wire ---
			Identifier id = (*iter_mdcDigi)->identify();
			int layer     = MdcID::layer(id);
			int wire      = MdcID::wire(id);
			double charge = (*iter_mdcDigi)->getChargeChannel();

			// --- get doca from track parameters ---
			//tension=myTensionWires[layer][wire];
			double wpos[7]={myEastPosWires[layer][wire][0], myEastPosWires[layer][wire][1], myEastPosWires[layer][wire][2]
				, myWestPosWires[layer][wire][0], myWestPosWires[layer][wire][1], myWestPosWires[layer][wire][2], tension};
			double doca_trk;
			double whitPos[3];// approching position on wire (in cm)
			if(debug) 
			{
				//cout<<"a = "<<myHelix_aVec<<endl;
				cout<<"* layer "<<layer<<", wire "<<wire<<", is active "<<myMdcDigiIsActive[i_digi] <<endl;
				//cout<<"wire positions : ("<<wpos[0]<<", "<<wpos[1]<<", "<<wpos[2]<<")  ("<<wpos[3]<<", "<<wpos[4]<<", "<<wpos[5]<<")"<<endl;
				//cout<<"zini = "<<vec_zini[i_digi]<<endl;
			}
			getDoca(myHelix_a, wpos, doca_trk, whitPos, vec_zini[i_digi]);
			//double doca_trk2 = myMdcUtilitySvc->doca(layer, wire, myHelix_aVec, myHelix_Ea, false, false);
			if(debug) {
				cout<<"doca = "<<doca_trk<<endl;
				//cout<<"doca2 = "<<doca_trk2<<endl;
			}
			int leftRight = 2;
			if(doca_trk!=0) {
				leftRight=int(doca_trk/fabs(doca_trk));
			}
			/*
			if(myAmbiguityMdcDigi[i_digi]!=0) {
				leftRight=myAmbiguityMdcDigi[i_digi];
				//if(debug) cout<<"fix leftRight = "<<leftRight<<endl;
			}*/
			int signDoca=1;
			signDoca=leftRight/fabs(leftRight);
			// --- conversion of left-right into the BESIII convention
			// if(leftRight==-1) leftRight=0;
			if(leftRight==1) leftRight=0;// fixed 2020-11-26
			else leftRight=1;
			//if(debug) cout<<"leftRight = "<<leftRight<<endl;
			
			vec_zini[i_digi] = whitPos[2];// update vec_zini

			// --- get measured doca --- tof in ns, driftTime in ns, T0Walk in ns
			double rawTime  = RawDataUtil::MdcTime((*iter_mdcDigi)->getTimeChannel());
			double tprop  = myMdcCalibFunSvc->getTprop(layer, vec_zini[i_digi]);
			double T0Walk = myMdcCalibFunSvc->getT0(layer,wire) +  myMdcCalibFunSvc->getTimeWalk(layer, charge);
			// --- time of flight (TOF) ---
			KalmanFit::Helix aHelix = *myHelix;
			HepPoint3D aNewPivot(whitPos[0],whitPos[1],whitPos[2]);
			aHelix.pivot(aNewPivot);
			double newPhi0 = aHelix.phi0();
			double dphi = newPhi0-myHelix->phi0();
			while(dphi<-M_PI) dphi+=2*M_PI;
			while(dphi> M_PI) dphi-=2*M_PI;
			double flightLength = fabs(dphi*myHelix->radius())*sqrt(1+myHelix->tanl()*myHelix->tanl());// in cm
			myMapFlylenIdx[flightLength]=i_digi;
			HepLorentzVector p4_pi = myHelix->momentum(dphi, 0.13957);
			double speed = p4_pi.beta()*CC;// cm/second
			double TOF = flightLength/speed*1.e9;// in ns
			// --- drift time ---
			double driftT = rawTime - myEventT0 - TOF - T0Walk - tprop;
			//if(debug) cout<<"myEventT0, driftT = "<<myEventT0<<", "<<driftT<<endl;
			// --- entrance Angle ---
			double phiWire = atan2(whitPos[1],whitPos[0]);
			double phiP    = p4_pi.phi();
			double entranceAngle = phiP-phiWire;
			while(entranceAngle<-M_PI) entranceAngle+=2*M_PI;
			while(entranceAngle> M_PI) entranceAngle-=2*M_PI;
			// --- measured drift distance ---
			double doca_hit = 0.1 * myMdcCalibFunSvc->driftTimeToDist(driftT,layer,wire,leftRight, entranceAngle);// in cm
			// --- get measurement error ---
			double docaErr_hit = 0.1 * myMdcCalibFunSvc->getSigma(layer, leftRight, doca_hit, entranceAngle, myHelix_a[4], vec_zini[i_digi], charge);// in cm
			//if(debug) cout<<"error_doca_hit = "<<docaErr_hit<<endl;

			// --- get derivatives, calculate J, update P, U ---
			doca_hit*=signDoca;
			//if(debug) cout<<"doca_hit = "<<doca_hit<<endl;
			double delD = doca_hit-doca_trk;
			double chi  = delD/docaErr_hit;
			if(debug) 
				cout<<"delta_m, err_m, chi = "<<delD<<", "<<docaErr_hit<<", "<<chi<<endl;
			myChiVecMdcDigi[i_digi]=chi;
			
			if(!myMdcDigiIsActive[i_digi]) continue;
			if(myUseAxialHitsOnly && myWireFlag[layer]!=0) continue;
			
			totChi2+=chi*chi;
			for(int ipar=0; ipar<nPar; ipar++)
			{
				double deri;
				//if(debug) cout<<"ipar = "<<ipar<<endl;
				getDeriLoc(ipar,myHelix_a, deri, wpos, vec_zini[i_digi]);// in cm
				J(ipar+1,1) = deri; 
				//if(debug) cout<<"d(doca)/d(a"<<"_"<<ipar<<") = "<<J(ipar+1,1)<<endl;
				//P(1,ipar+1) += J(1,ipar+1)*(doca_hit-doca_trk)/(docaErr_hit*docaErr_hit); 
				//P(ipar+1,1) += J(ipar+1,1)*chi/docaErr_hit; 
				double scale=1; //if(fabs(chi)>5) scale=fabs(chi*chi);
				P(ipar+1,1) += J(ipar+1,1)*(doca_hit-doca_trk)/(docaErr_hit*docaErr_hit*scale);// large error for large chi
				for(int ipar2=0; ipar2<=ipar; ipar2++)
				{
					//U(ipar+1,ipar2+1)+=J(ipar+1,1)*J(ipar2+1,1)/(docaErr_hit*docaErr_hit);
					U(ipar+1,ipar2+1)+=J(ipar+1,1)*J(ipar2+1,1)/(docaErr_hit*docaErr_hit*scale);// large error for large chi2
				}
			}
		}// loop MDC hits
		if(debug) cout<<"end of MDC hits loop"<<endl;
		// <--- end of MDC hits loop ---

		// --- loop CGEM clusters ---
		i_digi=0;
		vector<const RecCgemCluster*>::iterator iter_cluster = myVecCgemCluster.begin();
		for(; iter_cluster!=myVecCgemCluster.end(); iter_cluster++, i_digi++)
		{
			// --- get layer, cluster flag ---
			int layer = (*iter_cluster)->getlayerid();
			int flag  = (*iter_cluster)->getflag();
			if(myUseAxialHitsOnly && flag!=0) continue;
			int sheet = (*iter_cluster)->getsheetid();
			CgemGeoReadoutPlane* readoutPlane = myCgemGeomSvc->getReadoutPlane(layer,sheet);
			if(debug) 
			{
				cout<<endl<<"layer, sheet, flag = "<<setw(5)<<layer<<setw(5)<<sheet<<setw(5)<<flag<<endl;
			}

			// --- get position & entrance Angle from track parameters ---
			double dphi = IntersectCylinder(myRmidDGapCgem[layer]);
			HepPoint3D pos = myHelix->x(dphi);
			//cout<<"dphi="<<dphi<<", pos="<<pos<<endl;
			double phi_trk = pos.phi();// in radian
			double z_trk = pos.z();// in cm
			J_dxda = dxda_cgem(*myHelix, dphi);
			Hep3Vector p3_trk = myHelix->momentum(dphi);
			double phi_p3trk = p3_trk.phi();
			double incidentPhi = phi_p3trk-phi_trk;
			while(incidentPhi>CLHEP::pi)  incidentPhi-=CLHEP::twopi;
			while(incidentPhi<-CLHEP::pi) incidentPhi+=CLHEP::twopi;

			// --- get X-cluster charge, measured X, derivatives, calculate J, update P, U ---
			double phi_CC(0.0);// in radian
			double X_CC(0.), V_CC(0.);// in cm
			double delta_m, err_m;
			double Q(0.);// in fC
			double T(100);// not valid at present
			int mode(2);
			Q=(*iter_cluster)->getenergydeposit();
			if(flag==0) {
				phi_CC = (*iter_cluster)->getrecphi();
				//phi_CC = (*iter_cluster)->getrecphi_CC();
				double del_phi = phi_CC-phi_trk;
				while(del_phi<-M_PI) del_phi+=2*M_PI;
				while(del_phi> M_PI) del_phi-=2*M_PI;
				X_CC=phi_CC*myRmidDGapCgem[layer];
				if(debug) cout<<"phi_trk, phi_rec = "<<phi_trk<<", "<<phi_CC<<endl;
				//double dX = del_phi*myRmidDGapCgem[layer];
				//double dX = del_phi;
				delta_m=del_phi;
				err_m=myCgemCalibSvc->getSigma(layer,flag,mode,incidentPhi,Q,T)/10.;// in cm
				err_m/=myRmidDGapCgem[layer];// in radian
				J_dmdx(1,1)=-1*pos.y()/myR2midDGapCgem[layer];
				J_dmdx(1,2)=pos.x()/myR2midDGapCgem[layer];
				J_dmdx(1,3)=0.;
			}
			else if(flag==1) {
				double V_trk = readoutPlane->getVFromPhiZ(phi_trk,z_trk*10,false)/10.0;// in cm
				V_CC=(*iter_cluster)->getrecv()/10.;// in cm
				//V_CC=(*iter_cluster)->getrecv_CC()/10.;// in cm
				if(debug) cout<<"V_trk, V_rec = "<<V_trk<<", "<<V_CC<<endl;
				delta_m=V_CC-V_trk;// in cm
				if(fabs(delta_m)>5) {
					/*int nextSheet = myNSheets[layer]-1-sheet;
					  CgemGeoReadoutPlane* nextReadoutPlane = myCgemGeomSvc->getReadoutPlane(layer,nextSheet);
					  double phiminNext = nextReadoutPlane->getPhimin();
					  double V_trk_next = nextReadoutPlane->getVInNextSheetFromV(V_trk*10, phiminNext)/10.0;// cm*10 -> mm
					  double delta_m_next = V_CC-V_trk_next;
					  if(fabs(delta_m_next)<delta_m) delta_m=delta_m_next;
					  if(debug) cout<<"V_trk_next = "<<V_trk_next<<endl;*/
					double V_trk_nearPhiMin = readoutPlane->getVFromPhiZ_nearPhiMin(phi_trk,z_trk*10,false)/10.0;// in cm
					double delta_m_2=V_CC-V_trk_nearPhiMin;// in cm
					if(fabs(delta_m_2)<fabs(delta_m)) delta_m=delta_m_2;
					if(debug) cout<<"V_trk_nearPhiMin= "<<V_trk_nearPhiMin<<endl;
				}
				err_m=myCgemCalibSvc->getSigma(layer,flag,mode,incidentPhi,Q,T)/10.;// in cm
				J_dmdx(1,1)=-myRVCgem[layer]*cos(myAngStereoCgem[layer])*pos.y()/myR2midDGapCgem[layer];
				J_dmdx(1,2)= myRVCgem[layer]*cos(myAngStereoCgem[layer])*pos.x()/myR2midDGapCgem[layer];
				J_dmdx(1,3)= -sin(myAngStereoCgem[layer]);
			}
			else {
				cout<<"flag ="<<flag<<", DotsHelixFitter::calculateNewHelix() is not ready for it!"<<endl;
				continue;// next cluster
			}
			JT=J_dmdx*J_dxda;
			J=JT.T();
			double chi  = delta_m/err_m;
			myChiVecCgemCluster[i_digi]=chi;
			if(debug) 
				cout<<"delta_m, err_m, chi = "<<delta_m<<", "<<err_m<<", "<<chi<<endl;
			totChi2+=chi*chi;
			for(int ipar=0; ipar<nPar; ipar++)
			{
				//if(debug) 
				//	cout<<"d(doca)/d(a"<<"_"<<ipar<<") = "<<J(ipar+1,1)<<endl;
				P(ipar+1,1) += J(ipar+1,1)*chi/err_m; 
				for(int ipar2=0; ipar2<=ipar; ipar2++)
				{
					U(ipar+1,ipar2+1)+=J(ipar+1,1)*J(ipar2+1,1)/(err_m*err_m);
				}
			}
		}// loop myVecCgemCluster
		if(debug) 
			cout<<"end of CGEM cluster loop"<<endl;
		// <--- end of CGEM clusters loop ---

		// --- delta a ---
		//if(debug) 
		//	cout<<"U="<<U<<endl;
		int ifail=99; 
		U.invert(ifail);// also the error matrix
		//if(debug) 
		//{
		//	cout<<"ifail = "<<ifail<<endl;
		//	cout<<"U^-1="<<U<<endl;
		//	cout<<"P="<<P<<endl;
		//}
		HepVector da = U*P;// in cm
		//if(debug) 
		//	cout<<"da = "<<da<<endl;

		// --- update helix ---
		for(int ipar=0; ipar<nPar; ipar++) 
		{
			myHelix_a[ipar]+=da[ipar];
			myHelix_aVec[ipar]=myHelix_a[ipar];
		}
		myHelix->a(myHelix_aVec);
		if(debug) 
		{
			cout<<"aNew = "<<myHelix_aVec
				<<" chi2 = "<<totChi2
				<<endl;
		}
		myChi2=totChi2;
		//if(nIterations==0) cout<<"before fit chi2 = "<<totChi2<<endl;

		// --- converge or too many iterations
		nIterations++;
		if(fabs(lastTotChi2-totChi2)<myDchi2Converge) {
			if(debug) 
				cout<<"DotsHelixFitter::calculateNewHelix(): converge after "<<nIterations<<" iterations"
					<<" with lastTotChi2="<<lastTotChi2<<", totChi2="<<totChi2
					<<endl;
			if(nPar==5) myHelix_Ea = U;
			//cout<<"myHelix_Ea updated "<<endl;
			break;
		}
		else if(totChi2-lastTotChi2>myDchi2Diverge) {
			n_chi2_increase++;
			if(n_chi2_increase>myMaxNChi2Increase||totChi2>myChi2Diverge)
				return 5;// divergence
		}
		if(nIterations>myMaxIteration) {
			if(debug) cout<<"DotsHelixFitter::calculateNewHelix(): "<<nIterations
				<<" iterations, break with lastTotChi2, totChi2 = "<<lastTotChi2<<", "<<totChi2<<endl;
			break;
		}
		lastTotChi2=totChi2;
	}// <--- iterations

	// --- delete vec_zini ---
	delete []vec_zini;

	if(nIterations>myMaxIteration) return 4;
	else return 0;
}


int DotsHelixFitter::deactiveHits(double chi_cut, int nMax)
{
	map<double, int> map_abschi_idx;
	//double maxChi2=0;
	//vector<const MdcDigi*>::iterator iter_mdcDigi_maxChi2=NULL;
	//int i_digi_maxChi2=-1;

	int i_digi=0;
	vector<const MdcDigi*>::iterator iter_mdcDigi = myVecMdcDigi.begin();
	for(; iter_mdcDigi!=myVecMdcDigi.end(); iter_mdcDigi++, i_digi++)
	{
		if(myMdcDigiIsActive[i_digi])
		{
			double abs_chi=fabs(myChiVecMdcDigi[i_digi]);

			if( abs_chi>chi_cut )   // && maxChi2<(abs_chi*abs_chi)
			{
				//maxChi2=abs_chi*abs_chi;
				//iter_mdcDigi_maxChi2=iter_mdcDigi;
				//i_digi_maxChi2=i_digi;
				map_abschi_idx[abs_chi]=i_digi;
			}
		}
	}
	int nDeactive=0;
	if(nMax>0&&map_abschi_idx.size()>0)
	{ 
		map<double, int>::reverse_iterator it_map;
		for(it_map=map_abschi_idx.rbegin(); it_map!=map_abschi_idx.rend(); it_map++)
		{
			myMdcDigiIsActive[it_map->second]=0;
			//cout<<"a hit with |chi|="<<it_map->first<<" deactived"<<endl;
			nDeactive++;
			if(nDeactive>=nMax) break;
		}
	}
	//if(i_digi_maxChi2!=-1) 
	//{
	//	myMdcDigiIsActive[i_digi_maxChi2] = 0;
	//	//cout<<"the worst hit "<<i_digi_maxChi2<<" (chi2="<<maxChi2<<") is deactive"<<endl;
	//	return 1;
	//}
	//else {
	//	//cout<<"no bad hits found!"<<endl;
	//	return 0;
	//}
	return nDeactive;
}


int DotsHelixFitter::activeHits(double chi_cut)
{
	int nActived=0;
	int i_digi=0;
	vector<const MdcDigi*>::iterator iter_mdcDigi = myVecMdcDigi.begin();
	for(; iter_mdcDigi!=myVecMdcDigi.end(); iter_mdcDigi++, i_digi++)
	{
		if(!myMdcDigiIsActive[i_digi])
		{
			double abs_chi=fabs(myChiVecMdcDigi[i_digi]);

			if( abs_chi<chi_cut )
			{
				myMdcDigiIsActive[i_digi] = 1;
				nActived++;
			}
		}
	}
	return nActived;
}


int DotsHelixFitter::deactiveHits(int layer_max, int nHit_max)
{
	int nActive=0;
	int nDeactive=0;

	vector<const MdcDigi*>::iterator iter_mdcDigi_begin = myVecMdcDigi.begin();
	map<double, int>::iterator it_map;
	//cout<<"digi size, flyLen size = "<<myVecMdcDigi.size()<<", "<<myMapFlylenIdx.size()<<endl;
	for(it_map=myMapFlylenIdx.begin(); it_map!=myMapFlylenIdx.end(); it_map++)
	{
		//cout<<"nActive = "<<nActive<<endl;
		vector<const MdcDigi*>::iterator iter_mdcDigi = iter_mdcDigi_begin+(it_map->second);
		if(myMdcDigiIsActive[it_map->second])
		{
			if(nActive>=nHit_max) myMdcDigiIsActive[it_map->second]=0;
			Identifier id = (*iter_mdcDigi)->identify();
			int layer     = MdcID::layer(id);
			if(layer>=layer_max)  myMdcDigiIsActive[it_map->second]=0;
			if(myMdcDigiIsActive[it_map->second]) nActive++;
			else {
				nDeactive++;
				//cout<<"a hit on layer "<<layer<<" is deactived"<<endl;
			}
		}
	}
	//cout<<nDeactive<<" hits deactived "<<endl;

	return nDeactive;
}

void DotsHelixFitter::setInitialHelix(KalmanFit::Helix aHelix)
{
	HepPoint3D pivot(0,0,0);
	myHelix_aVec = aHelix.a();
	if(myHelix!=NULL) delete myHelix;
	myHelix = new KalmanFit::Helix(pivot, myHelix_aVec);
	//cout<<"alpha="<<myHelix->alpha()<<endl;
	for(int i=0; i<5; i++) myHelix_a[i]=myHelix_aVec[i];

}

void DotsHelixFitter::loadOneDcDigi(const MdcDigi* aDcDigi)
{
	// --- get id, layer, wire ---
	Identifier id = aDcDigi->identify();
	myLayer = MdcID::layer(id);
	myWire  = MdcID::wire(id);
	myCharge= aDcDigi->getChargeChannel();

	// --- get myWirePos ---
	double tension = 9999.;
	//tension=myTensionWires[myLayer][myWire];
	myWirePos[0]=myEastPosWires[myLayer][myWire][0];
	myWirePos[1]=myEastPosWires[myLayer][myWire][1];
	myWirePos[2]=myEastPosWires[myLayer][myWire][2]; 
	myWirePos[3]=myWestPosWires[myLayer][myWire][0];
	myWirePos[4]=myWestPosWires[myLayer][myWire][1];
	myWirePos[5]=myWestPosWires[myLayer][myWire][2];
	myWirePos[6]=tension;

}

void DotsHelixFitter::calculateDocaFromTrk(const MdcDigi* aDcDigi)
{
	bool debug = false;

	loadOneDcDigi(aDcDigi);

	// --- get doca from track parameters ---
	getDoca(myHelix_a, myWirePos, myDocaFromTrk, myPosOnWire, myPosOnWire[2]);
	if(debug) 
	{
		cout<<"DotsHelixFitter::UpdateDcDigiInfo(): "<<endl;
		cout<<"doca = "<<myDocaFromTrk<<endl;
		cout<<"point on wire: "<<myPosOnWire[0]<<", "<<myPosOnWire[1]<<", "<<myPosOnWire[2]<<endl;
	}

	// --- flight length ---
	KalmanFit::Helix aHelix = *myHelix;
	HepPoint3D aNewPivot(myPosOnWire[0],myPosOnWire[1],myPosOnWire[2]);
	aHelix.pivot(aNewPivot);
	double newPhi0 = aHelix.phi0();
	double dphi = newPhi0-myHelix->phi0();
	while(dphi<-M_PI) dphi+=2*M_PI;
	while(dphi> M_PI) dphi-=2*M_PI;
	myFlightLength = fabs(dphi*myHelix->radius())*sqrt(1+myHelix->tanl()*myHelix->tanl());// in cm
}

void DotsHelixFitter::updateDcDigiInfo(const MdcDigi* aDcDigi)
{
	bool debug = false;

	loadOneDcDigi(aDcDigi);

	// --- get doca from track parameters ---
	getDoca(myHelix_a, myWirePos, myDocaFromTrk, myPosOnWire, myPosOnWire[2]);
	double phiWire = atan2(myPosOnWire[1],myPosOnWire[0]);
	if(debug) 
	{
		cout<<"DotsHelixFitter::UpdateDcDigiInfo(): "<<endl;
		cout<<"doca = "<<myDocaFromTrk<<endl;
		cout<<"point on wire: "<<myPosOnWire[0]<<", "<<myPosOnWire[1]<<", "<<myPosOnWire[2]<<endl;
	}

	// --- flight length ---
	KalmanFit::Helix aHelix = *myHelix;
	HepPoint3D aNewPivot(myPosOnWire[0],myPosOnWire[1],myPosOnWire[2]);
	aHelix.pivot(aNewPivot);
	double newPhi0 = aHelix.phi0();
	double dphi = newPhi0-myHelix->phi0();
	while(dphi<-M_PI) dphi+=2*M_PI;
	while(dphi> M_PI) dphi-=2*M_PI;
	myFlightLength = fabs(dphi*myHelix->radius())*sqrt(1+myHelix->tanl()*myHelix->tanl());// in cm

	// --- approaching point on Helix
	HepPoint3D posOnTrk = aHelix.x();
	myPosOnTrk[0]=posOnTrk.x();
	myPosOnTrk[1]=posOnTrk.y();
	myPosOnTrk[2]=posOnTrk.z();
	double phiPosOnTrk=atan2(myPosOnTrk[1],myPosOnTrk[0]);

	// --- left or right
	myLeftRight = 2;
	int signDoca=1;
	if(myDocaFromTrk!=0) 
	{
		myLeftRight=int(myDocaFromTrk/fabs(myDocaFromTrk));
		signDoca=myLeftRight;

		// ---> conversion of left-right into the BESIII convention
		//if(myLeftRight==-1) myLeftRight=0;
		if(myLeftRight==1) myLeftRight=0;// fixed 2020-11-26
		else myLeftRight=1;
		double dphiLR=phiWire-phiPosOnTrk;
		while(dphiLR<-M_PI) dphiLR+=2*M_PI;
		while(dphiLR> M_PI) dphiLR-=2*M_PI;
		//if(dphiLR>0)      myLeftRight=0;
		//else if(dphiLR<0) myLeftRight=1;
		//cout<<"myLeftRight, dphiLR = "<<myLeftRight<<", "<<dphiLR<<endl;
	}
	if(debug) cout<<"myLeftRight = "<<myLeftRight<<endl;

	// --- time of flight (TOF) ---
	HepLorentzVector p4_pi = myHelix->momentum(dphi, 0.13957);
	double speed = p4_pi.beta()*CC;// cm/second
	double TOF = myFlightLength/speed*1.e9;// in ns

	// --- get measured doca --- tof in ns, driftTime in ns, T0Walk in ns
	double rawTime  = RawDataUtil::MdcTime(aDcDigi->getTimeChannel());
	double tprop  = myMdcCalibFunSvc->getTprop(myLayer, myPosOnWire[2]);
	double T0Walk = myMdcCalibFunSvc->getT0(myLayer,myWire) +  myMdcCalibFunSvc->getTimeWalk(myLayer, myCharge);

	// --- drift time ---
	myDriftTime = rawTime - myEventT0 - TOF - T0Walk - tprop;
	if(debug) cout<<"driftT = "<<myDriftTime<<endl;
	// --- entrance Angle ---
	double phiP    = p4_pi.phi();
	double entranceAngle = phiP-phiWire;
	while(entranceAngle<-M_PI) entranceAngle+=2*M_PI;
	while(entranceAngle> M_PI) entranceAngle-=2*M_PI;
	myEntranceAngle = entranceAngle;
	// --- measured drift distance ---
	myDocaFromDigi = 0.1 * myMdcCalibFunSvc->driftTimeToDist(myDriftTime,myLayer,myWire,myLeftRight,entranceAngle);// in cm
	myDriftDist[myLeftRight]=myDocaFromDigi;
	myDriftDist[1-myLeftRight]=0.1 * myMdcCalibFunSvc->driftTimeToDist(myDriftTime,myLayer,myWire,1-myLeftRight,entranceAngle);// in cm
	// --- get measurement error ---
	double docaErr_hit = 0.1 * myMdcCalibFunSvc->getSigma(myLayer, myLeftRight, myDocaFromDigi, myEntranceAngle, myHelix_a[4], myPosOnWire[2], myCharge);// in cm
	//if(debug) cout<<"error_doca_hit = "<<docaErr_hit<<endl;
	myDriftDistErr[myLeftRight]=docaErr_hit;
	myDriftDistErr[1-myLeftRight] = 0.1 * myMdcCalibFunSvc->getSigma(myLayer, 1-myLeftRight, myDocaFromDigi, myEntranceAngle, myHelix_a[4], myPosOnWire[2], myCharge);// in cm

	// --- get chi ---
	myDocaFromDigi*=signDoca;
	if(debug) cout<<"doca_hit = "<<myDocaFromDigi<<endl;
	double delD = myDocaFromDigi-myDocaFromTrk;
	myDcChi  = delD/docaErr_hit;
}

RecMdcHit DotsHelixFitter::makeRecMdcHit(const MdcDigi* aDcDigi)
{
	updateDcDigiInfo(aDcDigi);
	RecMdcHit aRecMdcHit; // = new RecMdcHit;
	aRecMdcHit.setDriftDistLeft(-1*myDriftDist[0]); 
	aRecMdcHit.setDriftDistRight(myDriftDist[1]);
	aRecMdcHit.setErrDriftDistLeft(myDriftDistErr[0]);
	aRecMdcHit.setErrDriftDistRight(myDriftDistErr[1]);
	aRecMdcHit.setChisqAdd(myDcChi*myDcChi);
	aRecMdcHit.setFlagLR(myLeftRight);
	//aRecMdcHit.setStat(stat);
	Identifier mdcId = aDcDigi->identify();
	aRecMdcHit.setMdcId(mdcId);
	double tdc = aDcDigi->getTimeChannel();
	double adc = aDcDigi->getChargeChannel();
	aRecMdcHit.setTdc(tdc);
	aRecMdcHit.setAdc(adc);
	aRecMdcHit.setDriftT(myDriftTime);
	// --- doca
	double doca=fabs(myDocaFromTrk);
	if(myLeftRight==0) doca*=-1;
	aRecMdcHit.setDoca(doca);
	aRecMdcHit.setEntra(myEntranceAngle);
	aRecMdcHit.setZhit(myPosOnWire[2]);
	aRecMdcHit.setFltLen(myFlightLength);
	return aRecMdcHit;
}

RecMdcHit DotsHelixFitter::makeRecMdcHit(const MdcDigi* aDcDigi, KalmanFit::Helix aHelix)
{
	setInitialHelix(aHelix);
	return makeRecMdcHit(aDcDigi);
}


vector<RecMdcHit> DotsHelixFitter::makeRecMdcHitVec(int sel)
{
	vector<RecMdcHit> aRecMdcHitVec;
	int i_digi=0;
	vector<const MdcDigi*>::iterator iter_mdcDigi = myVecMdcDigi.begin();
	for(; iter_mdcDigi!=myVecMdcDigi.end(); iter_mdcDigi++, i_digi++)
	{
		if(sel==1&&!myMdcDigiIsActive[i_digi]) continue;// skip inactive hits
		aRecMdcHitVec.push_back(makeRecMdcHit(*iter_mdcDigi));
	}
	return aRecMdcHitVec;
}


double DotsHelixFitter::IntersectCylinder(double r)
{
	double m_rad = myHelix->radius();
	double l = myHelix->center().perp();

	double cosPhi = (m_rad * m_rad + l * l  - r * r) / (2 * m_rad * l);

	if(cosPhi < -1 || cosPhi > 1) return 0;

	double dPhi = myHelix->center().phi() - acos(cosPhi) - myHelix->phi0();

	if(dPhi < -M_PI) dPhi += 2 * M_PI;

	return dPhi;
}

HepMatrix DotsHelixFitter::dxda_cgem(KalmanFit::Helix a, double phi)
{
	HepMatrix dXDA(3,5,0);

	double alpha = a.alpha();
	double dr = a.dr();
	double phi0 = a.phi0();
	double kappa = a.kappa();
	double dz = a.dz();
	double tanl = a.tanl();

	HepPoint3D pos = a.x(phi);
	double x = pos.x();
	double y = pos.y();
	double z = pos.z();
	double r2 = x*x+y*y;

	double cosPhi=cos(phi);
	double sinPhi=sin(phi);

	double dPhiDdr = -(kappa/alpha-cosPhi/(dr+alpha/kappa))/sinPhi;
	//double dPhiDdr2= -kappa*kappa*(2.0*alpha*dr+kappa*(dr*dr+r2))/(2*alpha*(alpha+dr*kappa)*(alpha+dr*kappa)*sinPhi);
	//cout<<"dPhiDdr = "<<dPhiDdr<<endl;
	//cout<<"dPhiDdr2= "<<dPhiDdr2<<endl;
	double dPhiDkappa = -kappa*(dr*dr-r2)*(2*alpha+dr*kappa)/(2*alpha*(alpha+dr*kappa)*(alpha+dr*kappa)*sinPhi);

	double dxDdr = cos(phi0)+alpha/kappa*sin(phi0+phi)*dPhiDdr;
	double dyDdr = sin(phi0)-alpha/kappa*cos(phi0+phi)*dPhiDdr;
	double dxDphi0 = -dr*sin(phi0)+alpha/kappa*(-sin(phi0)+sin(phi0+phi));
	double dyDphi0 = -dr*cos(phi0)+alpha/kappa*( cos(phi0)-cos(phi0+phi));
	double dxDkappa = -alpha/(kappa*kappa)*(cos(phi0)-cos(phi0+phi))+alpha/kappa*sin(phi0+phi)*dPhiDkappa;
	double dyDkappa = -alpha/(kappa*kappa)*(sin(phi0)-sin(phi0+phi))-alpha/kappa*cos(phi0+phi)*dPhiDkappa;
	double dzDdr = -alpha/kappa*tanl*dPhiDdr;
	double dzDkappa = alpha/(kappa*kappa)*tanl*phi-alpha/kappa*tanl*dPhiDkappa;
	double dzDtanl = -alpha/kappa*phi;

	dXDA(1,1) = dxDdr;
	dXDA(1,2) = dxDphi0;
	dXDA(1,3) = dxDkappa;
	dXDA(2,1) = dyDdr;
	dXDA(2,2) = dyDphi0;
	dXDA(2,3) = dyDkappa;
	dXDA(3,1) = dzDdr;
	dXDA(3,3) = dzDkappa;
	dXDA(3,4) = 1.0;
	dXDA(3,5) = dzDtanl;

	return dXDA;
}
