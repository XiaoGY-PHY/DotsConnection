
#include "DotsConnection/DotsConnection.h"
#include "DotsConnection/DotsHelixFitter.h"
#include "DotsConnection/TrkFitFun.h"

#include "GaudiKernel/MsgStream.h"
//#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IPartPropSvc.h"
#include "GaudiKernel/INTupleSvc.h"
//#include "GaudiKernel/IDataManagerSvc.h"
//#include "GaudiKernel/IDataProviderSvc.h"
#include "HepPDT/ParticleDataTable.hh"
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "EventModel/EventHeader.h"
#include "EvTimeEvent/RecEsTime.h"
#include "Identifier/MdcID.h"
#include "RawEvent/RawDataUtil.h"
// #include "MdcCalibFunSvc/IMdcCalibFunSvc.h"

using namespace Event; // for McParticleCol
//using namespace EventModel; // for EventHeader

DotsConnection::DotsConnection(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator)
{
	myMdcCalibFunSvc=NULL;

	// Part 1: Declare the properties
	//declareProperty("MyInt", m_myInt);
	declareProperty("Ntuple", myNtProd=0);
	declareProperty("driftTimeUpLimit",  myDriftTimeUpLimit = 400);
	declareProperty("MdcHitChi2Cut", myMdcHitChi2Cut = 100);
	declareProperty("ChiCut_circle", myChiCut_circle=5);
	declareProperty("NmaxDeact_circle", myNmaxDeact_circle=1);
	declareProperty("ChiCut_helix", myChiCut_helix=5);
	declareProperty("NmaxDeact_helix", myNmaxDeact_helix=1);
	declareProperty("Debug", myDebug=0);
	declareProperty("Chi2CutDiverge", myChi2CutDiverge=99999999);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

StatusCode DotsConnection::initialize(){

	// Part 1: Get the messaging service, print where you are
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << " DotsConnection initialize()" << endreq;

	// Part 2: Print out the property values
	// log << MSG::INFO << "  MyInt =    " << m_myInt << endreq;
	
	// --- Get MdcGeomSvc ---
	IMdcGeomSvc* imdcGeomSvc;
	//StatusCode sc1 =Gaudi::svcLocator()->service("MdcGeomSvc", imdcGeomSvc);
	StatusCode sc = service ("MdcGeomSvc", imdcGeomSvc);
	myMdcGeomSvc =  dynamic_cast<MdcGeomSvc*> (imdcGeomSvc);
	if (sc.isFailure()) {
		return( StatusCode::FAILURE);
	}
	int nWire = myMdcGeomSvc->getWireSize();
	int nLayer= myMdcGeomSvc->getLayerSize();
	int nSuperLayer= myMdcGeomSvc->getSuperLayerSize();
	int nSeg= myMdcGeomSvc->getSegmentNo();
	/*cout<<"nWire = "<<nWire
		<<"  nLayer="<<nLayer
		<<"  nSuperLayer="<<nSuperLayer
		<<"  nSeg="<<nSeg
		<<endl;*/
	int nCellTot=0;
	for(int i=0; i<nLayer; i++)
	{
		const MdcGeoLayer* aLayer = myMdcGeomSvc->Layer(i);
		myNWire[i]=aLayer->NCell();
		myRLayer[i]=aLayer->Radius();
		//cout<<"layer "<<i<<", "<<aLayer->NCell()<<" cells, slant "<<aLayer->Slant()<<", R="<<aLayer->Radius()<<endl;

		double nomShift = aLayer->nomShift();
		if(nomShift>0)      myWireFlag[i]=1;
		else if(nomShift<0) myWireFlag[i]=-1;
		else                myWireFlag[i]=0;

		//nCellTot+=aLayer->NCell();
		nCellTot+=myNWire[i];
		for(int j=0; j<myNWire[i]; j++)
		{
			const MdcGeoWire* aWire = myMdcGeomSvc->Wire(i,j);
			//cout<<"        wire "<<j<<", BWpos ="<<aWire->Backward()
			//	<<", "<<aWire->BWirePos()
			//	<<", FWpos ="<<aWire->Forward()
			//	<<", "<<aWire->FWirePos()
			//	<<endl;
			myWirePhi[i][j]=aWire->Forward().phi();
			int iInnerLayer = i-1;
			if(iInnerLayer>=0) {
				for(int k=0; k<myNWire[iInnerLayer]; k++)
				{
					int k_last = k-1;
					if(k_last<0) k_last=myNWire[iInnerLayer]-1;
					double dphi_last = dPhi(myWirePhi[iInnerLayer][k_last], myWirePhi[i][j]);
					double dphi = dPhi(myWirePhi[iInnerLayer][k], myWirePhi[i][j]);
					if(dphi_last<0&&dphi>0) {
						myInnerWire[i][j][0]=k_last;
						myInnerWire[i][j][1]=k;
						//cout<<"k_last, k ="<<k_last<<", "<<k<<endl;
						break;
					}
				}
			}
		}
	}
	for(int i=0; i<nLayer; i++)
	{
		for(int j=0; j<myNWire[i]; j++)
		{
			int iOuterLayer = i+1;
			if(iOuterLayer<nLayer) {
				for(int k=0; k<myNWire[iOuterLayer]; k++)
				{
					int k_last = k-1;
					if(k_last<0) k_last=myNWire[iOuterLayer]-1;
					double dphi_last = dPhi(myWirePhi[iOuterLayer][k_last], myWirePhi[i][j]);
					double dphi = dPhi(myWirePhi[iOuterLayer][k], myWirePhi[i][j]);
					if(dphi_last<0&&dphi>0) {
						myOuterWire[i][j][0]=k_last;
						myOuterWire[i][j][1]=k;
						//cout<<"k_last, k ="<<k_last<<", "<<k<<endl;
						//cout<<"phi="<<myWirePhi[i][j]<<", outer phi = "<<myWirePhi[iOuterLayer][k_last]<<", "<<myWirePhi[iOuterLayer][k]<<endl;
						break;
					}
				}
			}
		}
	}
	//cout<<"Total "<<nCellTot<<" cells"<<endl;
	// ----------------------
	

	// --- Initialize RawDataProviderSvc ---
	IRawDataProviderSvc* irawDataProviderSvc;
	sc = service ("RawDataProviderSvc", irawDataProviderSvc);
	myRawDataProviderSvc = dynamic_cast<RawDataProviderSvc*> (irawDataProviderSvc);
	if ( sc.isFailure() ){
		log << MSG::FATAL << name()<<" Could not load RawDataProviderSvc!" << endreq;
		return StatusCode::FAILURE;
	} 

	// --- MdcCalibFunSvc ---
	IMdcCalibFunSvc* imdcCalibSvc;                                                                                                      
	sc = service("MdcCalibFunSvc", imdcCalibSvc);
	if ( sc.isSuccess() ){ 
		myMdcCalibFunSvc = dynamic_cast<MdcCalibFunSvc*>(imdcCalibSvc);
	}
	else {
		cout<<"DotsConnection::initialize(): can not get MdcCalibFunSvc"<<endl;
	}
	// --- initialize DotsHelixFitter ---
	myDotsHelixFitter.initialize();
	myDotsHelixFitter.setChi2Diverge(myChi2CutDiverge);
	if(myNtProd&1) 
	{
		NTuplePtr nt_ptr(ntupleSvc(),"TestDotsHelixFitter/trkPar");
		if( nt_ptr ) myNtHelixFitter = nt_ptr;
		else
		{
			myNtHelixFitter = ntupleSvc()->book("TestDotsHelixFitter/trkPar",CLID_ColumnWiseTuple,"trkPar");
			if( myNtHelixFitter ) {
				myNtHelixFitter->addItem       ("run",         myRUN);
				myNtHelixFitter->addItem       ("evt",         myEVT);
				myNtHelixFitter->addItem       ("pid",         myPID);
				myNtHelixFitter->addItem       ("Npar",        myNPar,  0,5);
				myNtHelixFitter->addIndexedItem("HelixMC",     myNPar,  myArrayHelixMC);
				myNtHelixFitter->addIndexedItem("HelixFitted", myNPar,  myArrayHelixFitted);
				myNtHelixFitter->addItem       ("NhitCircle",        myNHitsCircle, 0,100);
				myNtHelixFitter->addIndexedItem("LayerHitsCircle",   myNHitsCircle, myLayerHitsCircle);
				myNtHelixFitter->addIndexedItem("ChiHitsCircle",     myNHitsCircle, myChiHitsCircle);
				myNtHelixFitter->addItem       ("Nhit",        myNHits, 0,100);
				myNtHelixFitter->addIndexedItem("LayerHits",   myNHits, myLayerHits);
				myNtHelixFitter->addIndexedItem("ChiHits",     myNHits, myChiHits);
				//cout<<"myNtHelixFitter added !"<<endl;
			}
		}
	}

	// --- get MdcUtilitySvc
	/*
	sc = service("MdcUtilitySvc", myMdcUtilitySvc);
	if( sc != StatusCode::SUCCESS ){
		log << MSG::FATAL << "can not use MdcUtilitySvc" << endreq;
	}*/


	return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode DotsConnection::execute() {

	// Part 1: Get the messaging service, print where you are
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "DotsConnection execute()" << endreq;

	// Part 2: Print out the different levels of messages
	/*log << MSG::DEBUG << "A DEBUG message" << endreq;
	log << MSG::INFO << "An INFO message" << endreq;
	log << MSG::WARNING << "A WARNING message" << endreq;
	log << MSG::ERROR << "An ERROR message" << endreq;
	log << MSG::FATAL << "A FATAL error message" << endreq;
	*/

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader) {
		log << MSG::WARNING << "Could not find Event Header" << endreq;
		return StatusCode::FAILURE;
	}
	myRun = eventHeader->runNumber();
	myEvt = eventHeader->eventNumber();
	if(myDebug) {
		cout<<endl<<"------------------------------ "<<endl;
		//cout<<"run, evt = "<<myEvt<<", "<<myRun<<endl;
		cout<<"run:"<<myRun<<"  , event: "<<myEvt<<endl;
	}
	
	//if(myRun==763) 
	//testDotsHelixFitterAllHits();
	
	//testDotsHelixFitterPartHits();


	// --- ideal tracking
	getMcFinalChargedStates();
	bool bookTrkCol = registerRecMdcTrack();
	if(!bookTrkCol) 
	{
		cout<<"DotsConnection::execute(): failed to register RecMdcTrackCol!"<<endl;
		return StatusCode::FAILURE;
	}
	associateDigisToMcParticles();



	return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

StatusCode DotsConnection::finalize() {

	// Part 1: Get the messaging service, print where you are
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "DotsConnection finalize()" << endreq;

	return StatusCode::SUCCESS;
}


double DotsConnection::dPhi(double phi1, double phi2)
{
	double dphi=phi1-phi2;
	while(dphi>M_PI)  dphi-=M_PI*2.0;
	while(dphi<-M_PI) dphi+=M_PI*2.0;
	return dphi;
}

KalmanFit::Helix DotsConnection::getMCHelix()
{

	// particle information
	std::string name;
	double charge = 0;
	double pos_x, pos_y, pos_z;
	double px, py, pz;


	// get PartPropSvc
	IPartPropSvc* p_PartPropSvc;
	static const bool CREATEIFNOTTHERE(true);
	StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
	if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
		cout<< " Could not initialize Particle Properties Service" << endl;
	}
	HepPDT::ParticleDataTable* p_particleTable = p_PartPropSvc->PDT();

	// get MC particle collection
	SmartDataPtr<McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
	if (mcParticleCol) {
		McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (;iter_mc != mcParticleCol->end(); iter_mc++ ) {
			if(!(*iter_mc)->primaryParticle()) continue;
			int pid = (*iter_mc)->particleProperty();
			int pid_abs = abs(pid);
			if( p_particleTable->particle(pid) ){
				name = p_particleTable->particle(pid)->name();
				charge= p_particleTable->particle(pid)->charge();
			}else if( p_particleTable->particle(-pid) ){
				name = "anti " + p_particleTable->particle(-pid)->name();
				charge = (-1.)*p_particleTable->particle(-pid)->charge();
			} 
			pos_x = (*iter_mc)->initialPosition().x();
			pos_y = (*iter_mc)->initialPosition().y();
			pos_z = (*iter_mc)->initialPosition().z();
			px    = (*iter_mc)->initialFourMomentum().px();
			py    = (*iter_mc)->initialFourMomentum().py();
			pz    = (*iter_mc)->initialFourMomentum().pz();
			if(pid_abs==11||pid_abs==13||pid_abs==211||pid_abs==321||pid_abs==2212) break;
		}
	}

	// --- get Helix ---
	HepPoint3D posTruth(pos_x, pos_y, pos_z);
	Hep3Vector p3Truth(px, py, pz);
	KalmanFit::Helix helix_truth(posTruth,p3Truth,charge);
	HepPoint3D origin(0, 0, 0);
	helix_truth.pivot(origin);
	if(myDebug) {
		cout<<"DotsConnection::getMCHelix() finds a valid MC particle "<<name<<" at "<<posTruth<<" with p3="<<p3Truth<<endl;
		cout<<"DotsConnection::getMCHelix() Helix = "<<helix_truth.a()<<endl;
	}

	return helix_truth;
}

void DotsConnection::getMcFinalChargedStates()
{
	resetFCSVec();

	// particle information
	std::string name;
	double charge = 0;
	double pos_x, pos_y, pos_z;
	double px, py, pz;
	HepPoint3D origin(0, 0, 0);


	// get PartPropSvc
	IPartPropSvc* p_PartPropSvc;
	static const bool CREATEIFNOTTHERE(true);
	StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
	if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
		cout<< " Could not initialize Particle Properties Service" << endl;
	}
	HepPDT::ParticleDataTable* p_particleTable = p_PartPropSvc->PDT();

	// get MC particle collection
	SmartDataPtr<McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
	if (mcParticleCol) {
		McParticleCol::iterator iter_mc = mcParticleCol->begin();
		if(myDebug) cout<<"MC charged particles: "<<endl;
		for (;iter_mc != mcParticleCol->end(); iter_mc++ ) 
		{
			if( !( (*iter_mc)->primaryParticle() || (*iter_mc)->decayFromGenerator() || (*iter_mc)->decayInFlight() ) ) continue;
			int pid = (*iter_mc)->particleProperty();
			int pid_abs = abs(pid);
			// --- e, mu, pi, K, p
			if(pid_abs==11||pid_abs==13||pid_abs==211||pid_abs==321||pid_abs==2212) 
			{	
				if( p_particleTable->particle(pid) )
				{
					name = p_particleTable->particle(pid)->name();
					charge= p_particleTable->particle(pid)->charge();
				} 
				else if( p_particleTable->particle(-pid) )
				{
					name = "anti " + p_particleTable->particle(-pid)->name();
					charge = (-1.)*p_particleTable->particle(-pid)->charge();
				} 
				pos_x = (*iter_mc)->initialPosition().x();
				pos_y = (*iter_mc)->initialPosition().y();
				pos_z = (*iter_mc)->initialPosition().z();
				px    = (*iter_mc)->initialFourMomentum().px();
				py    = (*iter_mc)->initialFourMomentum().py();
				pz    = (*iter_mc)->initialFourMomentum().pz();
				// --- get Helix ---
				HepPoint3D posTruth(pos_x, pos_y, pos_z);
				Hep3Vector p3Truth(px, py, pz);
				KalmanFit::Helix helix_truth(posTruth,p3Truth,charge);
				helix_truth.pivot(origin);
				// --- get the flight length of the first half in MDC --
				double dphi=helix_truth.IntersectCylinder(myRLayer[42]/10.);// mm -> cm
				if(dphi==0) dphi=M_PI*charge;
				HepPoint3D posOuter = helix_truth.x(dphi);
				double lenOuter = helix_truth.flightLength(posOuter);
				double lenInner = helix_truth.flightLength(posTruth);
				double lenInMdc = lenOuter-lenInner;
				myVecTrkLenFirstHalf.push_back(lenInMdc);
				// --- fill charged final states
				int trkIdx = (*iter_mc)->trackIndex();
				myVecMCTrkId.push_back(trkIdx);
				myVecPDG.push_back(pid);
				vector<double> helix;
				helix.push_back(helix_truth.dr());
				helix.push_back(helix_truth.phi0());
				helix.push_back(helix_truth.kappa());
				helix.push_back(helix_truth.dz());
				helix.push_back(helix_truth.tanl());
				myVecHelix.push_back(helix);
				if(myDebug)
					cout<<" trk idx "<<trkIdx<<", pdg code "<<pid
						<<", p3=("<<px<<", "<<py<<", "<<pz<<")"
						<<", pos=("<<pos_x<<","<<pos_y<<","<<pos_z<<")"
						<<", dphi="<<dphi<<", lenOuter="<<lenOuter<<", lenInner="<<lenInner
						<<", lenInMdc = "<<lenInMdc
						//<<", lenInMdc = "<<lenOuter
						<<endl;
			}
		}
	}
}


void DotsConnection::resetFCSVec()
{
	myVecMCTrkId.clear();
	myVecPDG.clear();
	myVecTrkLenFirstHalf.clear();
	myVecHelix.clear();
	//myVecCgemMcHit.clear();
	//myVecCgemXcluster.clear();
	//myVecCgemVcluster.clear();
}

void DotsConnection::associateDigisToMcParticles()
{
	// --- get event T0 ---
	double T0 = getEventStartTime();// in ns

	// --- get CGEM MC hits
	//SmartDataPtr<Event::CgemMcHitCol> cgemMcHitCol(eventSvc(),"/Event/MC/CgemMcHitCol");
	//if(!cgemMcHitCol) cout<<"DotsConnection::associateDigisToMcParticles() does not find CgemMcHitCol!"<<endl;
	//cout<<"CGEM MC hits obtained"<<endl;

	// --- get CGEM cluster
	//SmartDataPtr<RecCgemClusterCol> aCgemClusterCol(eventSvc(),"/Event/Recon/RecCgemClusterCol");
	//if(!aCgemClusterCol) cout<<"DotsConnection::associateDigisToMcParticles() does not find RecCgemClusterCol!"<<endl;
	//RecCgemClusterCol::iterator iter_cluster_begin=aCgemClusterCol->begin();
	//RecCgemClusterCol::iterator iter_cluster;
	//cout<<"CGEM clusters obtained"<<endl;
	
	// --- get MDC MC hits
	SmartDataPtr<Event::MdcMcHitCol> mdcMcHitCol(eventSvc(),"/Event/MC/MdcMcHitCol");
	if(!mdcMcHitCol) cout<<"DotsConnection::associateDigisToMcParticles() does not find MdcMcHitCol!"<<endl;
	if(myDebug)
		cout<<"MDC MC hits obtained, n="<<mdcMcHitCol->size()<<endl;
	Event::MdcMcHitCol::iterator iter_mdcMcHit = mdcMcHitCol->begin();
	for(; iter_mdcMcHit!=mdcMcHitCol->end(); iter_mdcMcHit++ )
	{
		string creatorProcess = (*iter_mdcMcHit)->getCreatorProcess();
		int trkId = (*iter_mdcMcHit)->getTrackIndex();
		int isSec = (*iter_mdcMcHit)->getIsSecondary();
		double trkLen = (*iter_mdcMcHit)->getFlightLength();
		double px=(*iter_mdcMcHit)->getMomentumX();
		double py=(*iter_mdcMcHit)->getMomentumY();
		double pz=(*iter_mdcMcHit)->getMomentumZ();
		double posx=(*iter_mdcMcHit)->getPositionX();
		double posy=(*iter_mdcMcHit)->getPositionY();
		double posz=(*iter_mdcMcHit)->getPositionZ();
		int pdgCode = (*iter_mdcMcHit)->getCurrentTrackPID();
		int digiIdx = (*iter_mdcMcHit)->getDigiIdx();
		Identifier id = (*iter_mdcMcHit)->identify();
		int layer     = MdcID::layer(id);
		int wire = MdcID::wire(id);
		Hep3Vector pos(posx,posy,0);
		Hep3Vector p3(px,py,0);
		double dotProd = pos*p3;

		if(myDebug)
			cout<<"  MDC MC hit from "<<creatorProcess
				<<", trk "<<trkId
				<<", isSec="<<isSec
				<<", len="<<trkLen/10 // mm -> cm
				<<", p3=("<<px<<", "<<py<<", "<<pz<<")"
				<<", pos=("<<posx<<", "<<posy<<", "<<posz<<")"
				<<", layer "<<layer<<" wire "<<wire
				<<", dotProd="<<dotProd
				<<", pdg code = "<<pdgCode
				<<", digi index = "<<digiIdx
				<<endl;
	}

	// --- get MDC digi vector
	uint32_t getDigiFlag = 0;
	MdcDigiVec mdcDigiVec = myRawDataProviderSvc->getMdcDigiVec(getDigiFlag);
	if(myDebug)
		cout<<"DotsConnection::associateDigisToMcParticles() get "<<mdcDigiVec.size()<<" MDC digis from RawDataProviderSvc"<<endl;
	//vector<const MdcDigi*> aMdcDigiVec;
	clearMdcDigiPointer();
	vector<MdcDigi*>::iterator iter_mdcDigi = mdcDigiVec.begin();
	for(; iter_mdcDigi!=mdcDigiVec.end(); iter_mdcDigi++) 
	{
		// --- get id, layer, wire ---
		Identifier id = (*iter_mdcDigi)->identify();
		int layer     = MdcID::layer(id);
		//if( (layer>=8&&layer<=19) || (layer>=36&&layer<=42) ) // --- only axial hits kept
		int wire = MdcID::wire(id);
		
		double rawTime  = RawDataUtil::MdcTime((*iter_mdcDigi)->getTimeChannel());
		double tprop  = myMdcCalibFunSvc->getTprop(layer, 0);
		double charge = (*iter_mdcDigi)->getChargeChannel();
		double T0Walk = myMdcCalibFunSvc->getT0(layer,wire) +  myMdcCalibFunSvc->getTimeWalk(layer, charge);
		double TOF = 0.;
		double driftT = rawTime - T0 - TOF - T0Walk - tprop;
		if(driftT>myDriftTimeUpLimit) continue;

		myMdcDigiPointer[layer][wire]=*iter_mdcDigi;
	}

	// --- loop MC tracks to get the proper digits (CGEM clusters and MDC hits)
	HepPoint3D origin(0, 0, 0);
	int nMcTrk = myVecMCTrkId.size();
	for(int i=0; i<nMcTrk; i++)
	{
		// --- get track index in mc particle collection
		int trkIdx = myVecMCTrkId[i];
		double trkLenInMdc = myVecTrkLenFirstHalf[i];
		if(myDebug)
			cout<<"MC particle "<<trkIdx<<", pdg code "<<myVecPDG[i]<<endl;

		// --- associate CGEM clusters through CGEM MC hits
		
		// int CgemCluster[3][2]={0,0,0,0,0,0};// [layer][x/v]
		// myVecCgemXcluster.clear();
		// myVecCgemVcluster.clear();
		// myVecCgem1DCluster.clear();
		// for(int l=0; l<3; l++)
		// {
		// 	for(int m=0; m<2; m++)
		// 	{
		// 		myVecCgemXCluIdx[l][m].clear();
		// 		myVecCgemVCluIdx[l][m].clear();
		// 	}
		// }
		// Event::CgemMcHitCol::iterator iter_CgemMcHit = cgemMcHitCol->begin();
		// for(; iter_CgemMcHit!=cgemMcHitCol->end(); iter_CgemMcHit++ )
		// {
		// 	string creatorProcess = (*iter_CgemMcHit)->GetCreatorProcess();
		// 	if(myDebug)
		// 		cout<<"      CGEM MC hit from process "<<creatorProcess;
		// 	if(creatorProcess=="Generator"||creatorProcess=="Decay")// only from generator or decay
		// 	{
		// 		//cout<<", trk id "<<(*iter_CgemMcHit)->GetTrackID()<<endl;
		// 		if((*iter_CgemMcHit)->GetTrackID()!=trkIdx) continue; // keep only matched MC hits
		// 		//cout<<"      is secondary "<<(*iter_CgemMcHit)->GetIsSecondary()<<endl;
		// 		const vector<int> & vec_xCluster = (*iter_CgemMcHit)->GetXclusterIdxVec();// find the X-clusters
		// 		int nXCluster=vec_xCluster.size();
		// 		for(int j=0; j<nXCluster; j++)
		// 		{
		// 			iter_cluster=iter_cluster_begin+vec_xCluster[j];
		// 			int clusterId = (*iter_cluster)->getclusterid();
		// 			int layer=(*iter_cluster)->getlayerid();
		// 			int sheet=(*iter_cluster)->getsheetid();
		// 			if(myDebug)
		// 				cout<<"          find one X-cluster on layer "<<layer;
		// 			if(CgemCluster[layer][0]==0) 
		// 			{
		// 				myVecCgemXcluster.push_back(*iter_cluster);
		// 				myVecCgem1DCluster.push_back(*iter_cluster);
		// 				myVecCgemXCluIdx[layer][sheet].push_back(clusterId);
		// 				CgemCluster[layer][0]++;
		// 				if(myDebug)
		// 					cout<<", associated.   "<<endl;
		// 			}
		// 		}
		// 		//cout<<" X-clusters done "<<endl;
		// 		const vector<int> & vec_vCluster = (*iter_CgemMcHit)->GetVclusterIdxVec();// find the V-clusters
		// 		int nVCluster=vec_vCluster.size();
		// 		for(int j=0; j<nVCluster; j++)
		// 		{
		// 			iter_cluster=iter_cluster_begin+vec_vCluster[j];
		// 			int clusterId = (*iter_cluster)->getclusterid();
		// 			int layer=(*iter_cluster)->getlayerid();
		// 			int sheet=(*iter_cluster)->getsheetid();
		// 			if(myDebug)
		// 				cout<<"          find one V-cluster on layer "<<layer;
		// 			if(CgemCluster[layer][1]==0) 
		// 			{
		// 				myVecCgemVcluster.push_back(*iter_cluster);
		// 				myVecCgem1DCluster.push_back(*iter_cluster);
		// 				myVecCgemVCluIdx[layer][sheet].push_back(clusterId);
		// 				CgemCluster[layer][1]++;
		// 				if(myDebug)
		// 					cout<<", associated.    "<<endl;
		// 			}
		// 		}
		// 	}
		// 	if(myDebug) cout<<endl;
		// }// end of looping CGEM MC hits
		
		
		myVecMdcDigi.clear();
		int nMdcXHits(0), nMdcVHits(0);
		
		// --- associate MDC hits through MC track index
		
		// vector<MdcDigi*>::iterator iter_mdcDigi = mdcDigiVec.begin();
		// for(; iter_mdcDigi!=mdcDigiVec.end(); iter_mdcDigi++) 
		// {
		// 	int mcTrkIdx = (*iter_mdcDigi)->getTrackIndex();
		// 	if(mcTrkIdx==trkIdx) 
		// 	{
		// 		//myVecMdcDigi.push_back((*iter_mdcDigi));
		// 		Identifier id = (*iter_mdcDigi)->identify();
		// 		int layer = MdcID::layer(id);
		// 		int wire  = MdcID::wire(id);
		// 		//cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" associated"<<endl;
		// 		//cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" matched"<<endl;
		// 	}
		// }
		
		// --- associate MDC hits through MdcMcHits 
		Event::MdcMcHitCol::iterator 
		iter_mdcMcHit = mdcMcHitCol->begin();
		for(; iter_mdcMcHit!=mdcMcHitCol->end(); iter_mdcMcHit++ )
		{
			int trkId = (*iter_mdcMcHit)->getTrackIndex();
			if(trkId!=trkIdx) continue;
			if((*iter_mdcMcHit)->getDigiIdx()==-9999) continue;
			Identifier id = (*iter_mdcMcHit)->identify();
			int layer = MdcID::layer(id);
			int wire  = MdcID::wire(id);
			int isSec = (*iter_mdcMcHit)->getIsSecondary();
			if(isSec!=0) {
				/*cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" isSec"
					<<", pid = "<<(*iter_mdcMcHit)->getCurrentTrackPID()
					<<endl;*/
				continue;
			}
			double trkLenMcHit = ((*iter_mdcMcHit)->getFlightLength())/10.;
			if(trkLenMcHit>1.5*trkLenInMdc) {
				//cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" trkLenMcHit>1.5*trkLenInMdc"<<endl;
				continue;
			}
			double px=(*iter_mdcMcHit)->getMomentumX();
			double py=(*iter_mdcMcHit)->getMomentumY();
			double pz=(*iter_mdcMcHit)->getMomentumZ();
			double posx=(*iter_mdcMcHit)->getPositionX();
			double posy=(*iter_mdcMcHit)->getPositionY();
			double posz=(*iter_mdcMcHit)->getPositionZ();
			Hep3Vector pos(posx,posy,0);
			Hep3Vector p3(px,py,0);
			double dotProd = pos*p3;
			if(dotProd<0) {
				//cout<<"      pos="<<pos<<", p3="<<p3<<endl;
				//cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" coming back"<<endl;
				continue;// reject hits from track coming back
			}
			//Identifier id = (*iter_mdcMcHit)->identify();
			//int layer = MdcID::layer(id);
			//int wire  = MdcID::wire(id);
			const MdcDigi* aMdcDigiPt = myMdcDigiPointer[layer][wire];
			if(aMdcDigiPt!=NULL) 
			{
				myVecMdcDigi.push_back(aMdcDigiPt);
				myMdcDigiPointer[layer][wire]=NULL;
				if(myDebug)
					cout<<"      MDC digi on layer "<<layer<<" wire "<<wire<<" associated"<<endl;
				if(layer<8||(layer>19&&layer<36)) nMdcVHits++;
				else nMdcXHits++;
			}
		}

		// --- fit to Helix
		//cout<<" start to fit "<<endl;
		int fitFlag=99;
		//if((myVecCgemVcluster.size()+nMdcVHits)>=2&&(myVecCgemXcluster.size()+nMdcXHits)>=3)
		if((nMdcVHits)>=2&&(nMdcXHits)>=3)
		{
			HepVector a_helix(5,0);
			a_helix(1)=myVecHelix[i][0];
			a_helix(2)=myVecHelix[i][1];
			a_helix(3)=myVecHelix[i][2];
			a_helix(4)=myVecHelix[i][3];
			a_helix(5)=myVecHelix[i][4];
			KalmanFit::Helix  ini_helix(origin,a_helix);
			myDotsHelixFitter.setInitialHelix(ini_helix);
			// myDotsHelixFitter.setCgemClusters(myVecCgem1DCluster);
			myDotsHelixFitter.setDChits(myVecMdcDigi,T0);

			myDotsHelixFitter.fitCircleOnly();
			myDotsHelixFitter.useAxialHitsOnly();
			int nIter=0;
			while(1)
			{
				fitFlag = myDotsHelixFitter.calculateNewHelix();
				nIter++;
				//if(fitFlag!=0) break;
				if(fitFlag==0)
				{
					if(myNtProd&1 && nIter==1) {
						myNHitsCircle=0;
						// --- fill CGEM clusters
						// vector<double> vecChiCgem = myDotsHelixFitter.getVecChiCgemCluster();
						// vector<const RecCgemCluster*>::iterator it_cluster=myVecCgem1DCluster.begin();
						// int i_cluster=0;
						// for(; it_cluster!=myVecCgem1DCluster.end(); it_cluster++, i_cluster++)
						// {
						// 	int flag=(*it_cluster)->getflag();
						// 	if(flag!=0) continue;
						// 	int layer=(*it_cluster)->getlayerid();
						// 	//if(flag==1) layer=-1*(layer+1);// for V-cluster
						// 	if(myNHitsCircle<100) 
						// 	{
						// 		myLayerHitsCircle[myNHitsCircle]=layer;
						// 		myChiHitsCircle[myNHitsCircle]  =vecChiCgem[i_cluster];
						// 		myNHitsCircle++;
						// 	}
						// }

						// --- fill MDC hits
						vector<const MdcDigi*>::iterator it_digi=myVecMdcDigi.begin();
						vector<double> vecChiMdc = myDotsHelixFitter.getVecChiMdcDigi();
						vector<int>    vecIsActMdc = myDotsHelixFitter.getVecMdcDigiIsAct();
						int i_digi=0;
						for(; it_digi!=myVecMdcDigi.end(); it_digi++, i_digi++)
						{
							//if(vecIsActMdc[i_digi]==0) continue;
							Identifier id = (*it_digi)->identify();
							int layer     = MdcID::layer(id);
							if(myWireFlag[layer]!=0) continue;
							if(myNHitsCircle<100) 
							{
								myLayerHitsCircle[myNHitsCircle]=layer;
								myChiHitsCircle[myNHitsCircle]  =vecChiMdc[i_digi];
								myNHitsCircle++;
							}
						}
					}
				}
				else break;
				if(myDotsHelixFitter.deactiveHits(myChiCut_circle,myNmaxDeact_circle)==0) break;
				if(nIter>100) break;
			}
			if(myDebug) {
				cout<<"initial helix "<<ini_helix.a()<<endl;
				cout<<"fitFlag="<<fitFlag<<endl;
				cout<<"nIter="<<nIter<<endl;
				cout<<"fitted circle "<<myDotsHelixFitter.getHelix()<<endl;
				cout<<"circle chi2="<<myDotsHelixFitter.getChi2()<<endl<<endl;
			}
			if(fitFlag==0) {
				myDotsHelixFitter.fitModeHelix();
				nIter=0;
				while(1)
				{
					fitFlag = myDotsHelixFitter.calculateNewHelix();
					if(fitFlag!=0) break;
					if(myDotsHelixFitter.deactiveHits(myChiCut_helix,myNmaxDeact_helix)==0) break;
				}
				while(1)
				{
					if(fitFlag!=0) break;
					if(myDotsHelixFitter.activeHits(myChiCut_helix)==0) break;
					fitFlag = myDotsHelixFitter.calculateNewHelix();
				}
				if(myDebug) {
					cout<<"fitFlag="<<fitFlag<<endl;
					cout<<"fitted helix "<<myDotsHelixFitter.getHelix()<<endl;
					cout<<"helix chi2="<<myDotsHelixFitter.getChi2()<<endl<<endl;
				}
			}
		}

		// --- save good track to RecMdcTrackCol and fill ntuple
		//cout<<"before fill myArrayHelixMC"<<endl;
		if(myNtProd&1) {
			myNHits=0;
			myRUN=myRun;
			myEVT=myEvt;
			//cout<<"fill 1"<<endl;
			myPID=myVecPDG[i];
			//cout<<"fill 2"<<endl;
			myNPar=5;
			//cout<<"fill 3"<<endl;
			myArrayHelixMC[0]=myVecHelix[i][0];
			//cout<<"fill 4"<<endl;
			myArrayHelixMC[1]=myVecHelix[i][1];
			myArrayHelixMC[2]=myVecHelix[i][2];
			myArrayHelixMC[3]=myVecHelix[i][3];
			myArrayHelixMC[4]=myVecHelix[i][4];
		}
		//cout<<"after fill myArrayHelixMC"<<endl;
		if(fitFlag==0)
		{
			saveARecMdcTrack();

			if(myNtProd&1) {
				// --- fill CGEM clusters
				// vector<double> vecChiCgem = myDotsHelixFitter.getVecChiCgemCluster();
				// vector<const RecCgemCluster*>::iterator it_cluster=myVecCgem1DCluster.begin();
				// int i_cluster=0;
				// for(; it_cluster!=myVecCgem1DCluster.end(); it_cluster++, i_cluster++)
				// {
				// 	int flag=(*it_cluster)->getflag();
				// 	int layer=(*it_cluster)->getlayerid();
				// 	if(flag==1) layer=-1*(layer+1);// for V-cluster
				// 	if(myNHits<100) 
				// 	{
				// 		myLayerHits[myNHits]=layer;
				// 		myChiHits[myNHits]  =vecChiCgem[i_cluster];
				// 		myNHits++;
				// 	}
				// }

				// --- fill MDC hits
				vector<const MdcDigi*>::iterator it_digi=myVecMdcDigi.begin();
				vector<double> vecChiMdc = myDotsHelixFitter.getVecChiMdcDigi();
				vector<int>    vecIsActMdc = myDotsHelixFitter.getVecMdcDigiIsAct();
				int i_digi=0;
				for(; it_digi!=myVecMdcDigi.end(); it_digi++, i_digi++)
				{
					if(vecIsActMdc[i_digi]==0) continue;
					Identifier id = (*it_digi)->identify();
					int layer     = MdcID::layer(id);
					if(myNHits<100) 
					{
						myLayerHits[myNHits]=layer;
						myChiHits[myNHits]  =vecChiMdc[i_digi];
						myNHits++;
					}

				}
			}// if myNtProd&1
		}

		// --- write out ntuple
		if(myNtProd&1) myNtHelixFitter->write();

	}// loop MC tracks

	return;
}

vector<const MdcDigi*> DotsConnection::getMdcDigiVec()
{
	uint32_t getDigiFlag = 0;
	/*
	   getDigiFlag += m_maxMdcDigi;
	   if(m_dropHot)     getDigiFlag |= MdcRawDataProvider::b_dropHot;
	   if(m_keepBadTdc)  getDigiFlag |= MdcRawDataProvider::b_keepBadTdc;
	   if(m_keepUnmatch) getDigiFlag |= MdcRawDataProvider::b_keepUnmatch;
	   */
	MdcDigiVec mdcDigiVec = myRawDataProviderSvc->getMdcDigiVec(getDigiFlag);
	if(myDebug) 
		cout<<"DotsConnection::getMdcDigiVec() get "<<mdcDigiVec.size()<<" MDC digis from RawDataProviderSvc"<<endl;
	///*
	// MdcDigiVec aMdcDigiVec;
	vector<const MdcDigi*> aMdcDigiVec;
	vector<MdcDigi*>::iterator iter_mdcDigi = mdcDigiVec.begin();
	for(; iter_mdcDigi!=mdcDigiVec.end(); iter_mdcDigi++) {
		// --- get id, layer, wire ---
		Identifier id = (*iter_mdcDigi)->identify();
		int layer     = MdcID::layer(id);

		// --- only axial hits kept
		//if( (layer>=8&&layer<=19) || (layer>=36&&layer<=42) ) 
		aMdcDigiVec.push_back(*iter_mdcDigi);

	}
	return aMdcDigiVec;
	//*/
	//return mdcDigiVec;
}

double DotsConnection::getEventStartTime()
{
	double T0=0;

	// --- get event start time ---
	SmartDataPtr<RecEsTimeCol> evTimeCol(eventSvc(),"/Event/Recon/RecEsTimeCol");
	if(evTimeCol){
		RecEsTimeCol::iterator iter_evt = evTimeCol->begin();
		if(iter_evt != evTimeCol->end()){
			//T0 = (*iter_evt)->getTest()*1.e-9;//m_bunchT0-s, getTest-ns
			T0 = (*iter_evt)->getTest();// getTest-ns
		}
	}
	else cout<<"error: DotsConnection::getEventStartTime() failed to access event start time"<<endl;

	return T0;
}

void DotsConnection::testDotsHelixFitterAllHits()
{
	// --- set MC helix
	KalmanFit::Helix  ini_helix = getMCHelix();
	myDotsHelixFitter.setInitialHelix(ini_helix);

	// --- get all MDC digi
	vector<const MdcDigi*> vecDigi = getMdcDigiVec();

	// --- get event start time
	double T0 = getEventStartTime();// in ns
	//cout<<"DotsConnection::testDotsHelixFitterAllHits() T0 = "<<T0<<endl;
	//myDotsHelixFitter.setT0(T0);

	// --- set DC digi to DotsHelixFitter
	myDotsHelixFitter.setDChits(vecDigi,T0);

	// --- get DC digi ordered with the fligth length
	vector<const MdcDigi*>::iterator iter_mdcDigi = vecDigi.begin();
	map<double, const MdcDigi*> aMapDcDigi;
	for(; iter_mdcDigi!=vecDigi.end(); iter_mdcDigi++)
	{
		aMapDcDigi[myDotsHelixFitter.getFlightLength(*iter_mdcDigi)] = *iter_mdcDigi;
	}

	// vector<const RecCgemCluster*> vecCgemCluster = getCgemClusterVec();
	// myDotsHelixFitter.setCgemClusters(vecCgemCluster);

	myDotsHelixFitter.fitCircleOnly();
	myDotsHelixFitter.calculateNewHelix();

	// --- fill a NTuple ---
	myNPar=5;
	HepVector a = ini_helix.a();
	HepVector a_new = myDotsHelixFitter.getHelix(); 
	for(int i=0; i<5; i++) 
	{
		myArrayHelixMC[i]=a[i];
		myArrayHelixFitted[i]=a_new[i];
	}
	myNtHelixFitter->write();
}

// vector<const RecCgemCluster*> DotsConnection::getCgemClusterVec(int view)
// {
// 	vector<const RecCgemCluster*> aVecCgemCluster;
// 	SmartDataPtr<RecCgemClusterCol> aCgemClusterCol(eventSvc(),"/Event/Recon/RecCgemClusterCol");
// 	if(aCgemClusterCol)
// 	{
// 		RecCgemClusterCol::iterator iter_cluster=aCgemClusterCol->begin();
// 		int nCluster = aCgemClusterCol->size();
// 		//cout<<"DotsConnection::getCgemClusterVec() finds a RecCgemClusterCol with "<<nCluster<<" clusters!"<<endl;
// 		//		cout<<"~~~~~~~~~~~~~~~~~~~~~~~~ check RecCgemClusterCol:"<<endl;
// 		//		cout    <<setw(10)<<"idx"
// 		//			<<setw(10)<<"layer"
// 		//			<<setw(10)<<"sheet"
// 		//			<<setw(10)<<"XVFlag"
// 		//			<<setw(10)<<"id1 ~"
// 		//			<<setw(10)<<"id2"
// 		//			<<setw(15)<<"phi"
// 		//			<<setw(15)<<"V"
// 		//			<<setw(15)<<"z"
// 		//			<<endl;
// 		for(; iter_cluster!=aCgemClusterCol->end(); iter_cluster++)
// 		{
// 			//			cout    <<setw(10)<<(*iter_cluster)->getclusterid()
// 			//				<<setw(10)<<(*iter_cluster)->getlayerid()
// 			//				<<setw(10)<<(*iter_cluster)->getsheetid()
// 			//				<<setw(10)<<(*iter_cluster)->getflag()
// 			//				<<setw(10)<<(*iter_cluster)->getclusterflagb()
// 			//				<<setw(10)<<(*iter_cluster)->getclusterflage()
// 			//				<<setw(15)<<setprecision(10)<<(*iter_cluster)->getrecphi()
// 			//				<<setw(15)<<setprecision(10)<<(*iter_cluster)->getrecv()
// 			//				<<setw(15)<<setprecision(10)<<(*iter_cluster)->getRecZ()
// 			//				<<endl;
// 			int flag = (*iter_cluster)->getflag();
// 			if(view==0||view==1) 
// 			{
// 				if(flag==view) aVecCgemCluster.push_back(*iter_cluster);
// 			}else if(view=2)
// 			{
// 				if(flag==0||flag==1) aVecCgemCluster.push_back(*iter_cluster);
// 			}
// 		}
// 		//		cout<<"~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
// 	}
// 	else cout<<"DotsConnection::getCgemClusterVec() does not find RecCgemClusterCol!"<<endl;
// 	return aVecCgemCluster;
// }

void DotsConnection::testDotsHelixFitterPartHits()
{
	// --- set MC helix
	KalmanFit::Helix  ini_helix = getMCHelix();
	myDotsHelixFitter.setInitialHelix(ini_helix);

	// --- get all MDC digi
	vector<const MdcDigi*> vecDigi = getMdcDigiVec();

	// --- get event start time
	double T0 = getEventStartTime();// in ns
	//cout<<"DotsConnection::testDotsHelixFitterPartHits() T0 = "<<T0<<endl;
	myDotsHelixFitter.setT0(T0);

	// --- set DC digi to DotsHelixFitter
	//myDotsHelixFitter.setDChits(vecDigi,T0);

	// --- get DC digi ordered with the fligth length
	vector<const MdcDigi*>::iterator iter_mdcDigi = vecDigi.begin();
	map<double, const MdcDigi*> aMapDcDigi;
	for(; iter_mdcDigi!=vecDigi.end(); iter_mdcDigi++)
	{
		aMapDcDigi[myDotsHelixFitter.getFlightLength(*iter_mdcDigi)] = *iter_mdcDigi;
	}

	int layerUdStudy=8;
	int layerUsed=6;
	vector<const MdcDigi*> smallVecDigi;
	map<double, const MdcDigi*>::iterator iter_digi = aMapDcDigi.begin();
	int nLayerCount=-1;
	const MdcDigi* digiUdStudy=NULL;
	for(; iter_digi!=aMapDcDigi.end(); iter_digi++)
	{
		const MdcDigi* aDigi = iter_digi->second;
		int layer = myDotsHelixFitter.getLayer(aDigi);
		if(layer==layerUdStudy) 
		{
			digiUdStudy=aDigi;
			nLayerCount=0;
		}
		if(layer>layerUdStudy&&nLayerCount>=0&&nLayerCount<layerUsed) 
		{
			smallVecDigi.push_back(aDigi);
			nLayerCount++;
		}
		if(nLayerCount==layerUsed) break;
	}
	if( digiUdStudy!=NULL && nLayerCount==layerUsed )
	{
		myDotsHelixFitter.setDChits(smallVecDigi,T0);
		myDotsHelixFitter.fitCircleOnly();
		myDotsHelixFitter.calculateNewHelix();
		double doca = myDotsHelixFitter.getDocaFromTrk(digiUdStudy);
		//cout<<"layer "<<layerUdStudy<<" doca prediction from outer "<<layerUsed<<" digis : "<<doca<<endl;
	}


	// --- fill a NTuple ---
	/*
	   myNPar=5;
	   HepVector a = ini_helix.a();
	   HepVector a_new = myDotsHelixFitter.getHelix(); 
	   for(int i=0; i<5; i++) 
	   {
	   myArrayHelixMC[i]=a[i];
	   myArrayHelixFitted[i]=a_new[i];
	   }
	   myNtHelixFitter->write();
	   */
}


void DotsConnection::clearMdcDigiPointer()
{
	for(int i=0; i<43; i++)
		for(int j=0; j<288; j++)
			myMdcDigiPointer[i][j]=NULL;
}

bool DotsConnection::registerRecMdcTrack()
{
	MsgStream log(msgSvc(), name());
	StatusCode sc;
	IDataProviderSvc* eventSvc = NULL;
	service("EventDataSvc", eventSvc);
	if (eventSvc) {
		log << MSG::INFO << "makeTds:event Svc has been found" << endreq;
	} else {
		log << MSG::FATAL << "makeTds:Could not find eventSvc" << endreq;
		return false;
	}
	IDataManagerSvc *dataManSvc = dynamic_cast<IDataManagerSvc*>(eventSvc);;

	// --- register RecMdcTrack
	myRecMdcTrackCol = new RecMdcTrackCol;
	DataObject *aRecMdcTrackCol;
	eventSvc->findObject("/Event/Recon/RecMdcTrackCol",aRecMdcTrackCol);
	if(aRecMdcTrackCol != NULL) {
		dataManSvc->clearSubTree("/Event/Recon/RecMdcTrackCol");
		eventSvc->unregisterObject("/Event/Recon/RecMdcTrackCol");
	}
	sc =  eventSvc->registerObject("/Event/Recon/RecMdcTrackCol", myRecMdcTrackCol);
	if(sc.isFailure()) {
		log << MSG::FATAL << " Could not register RecMdcTrack collection" <<endreq;
		return false;
	}
	log << MSG::INFO << "RecMdcTrackCol registered successfully!" <<endreq;


	myRecMdcHitCol= new RecMdcHitCol;
	DataObject *aRecMdcHitCol;
	eventSvc->findObject("/Event/Recon/RecMdcHitCol",aRecMdcHitCol);
	if(aRecMdcHitCol != NULL) {
		dataManSvc->clearSubTree("/Event/Recon/RecMdcHitCol");
		eventSvc->unregisterObject("/Event/Recon/RecMdcHitCol");
	}

	sc =  eventSvc->registerObject("/Event/Recon/RecMdcHitCol", myRecMdcHitCol);
	if(sc.isFailure()) {
		log << MSG::FATAL << " Could not register RecMdcHit collection" <<endreq;
		return false;
	}
	log << MSG::INFO << "RecMdcHitCol registered successfully!" <<endreq;


	return true;
}

// bool sortCluster(const RecCgemCluster* clusterA , const RecCgemCluster* clusterB){
// 	return clusterA->getlayerid()<clusterB->getlayerid();
// }

bool DotsConnection::saveARecMdcTrack()
{
	int tkStat =4;

	RecMdcTrack* recMdcTrack = new RecMdcTrack();

	int trackId = myRecMdcTrackCol->size();
	recMdcTrack->setTrackId(trackId);

	double helixPar[5];
	HepVector aHelixVec = myDotsHelixFitter.getHelix();
	helixPar[0]=aHelixVec[0];
	helixPar[1]=aHelixVec[1];
	helixPar[2]=aHelixVec[2];
	helixPar[3]=aHelixVec[3];
	helixPar[4]=aHelixVec[4];//helixPar[4]+=0.1;
	recMdcTrack->setHelix(helixPar);

	// --- for ntuple
	if(myNtProd&1) {
		myArrayHelixFitted[0]=aHelixVec[0];
		myArrayHelixFitted[1]=aHelixVec[1];
		myArrayHelixFitted[2]=aHelixVec[2];
		myArrayHelixFitted[3]=aHelixVec[3];
		myArrayHelixFitted[4]=aHelixVec[4];
	}

	int    q     = helixPar[2]>0? 1:-1;
	KalmanFit::Helix aHelix = myDotsHelixFitter.getClassHelix();
	double pxy   = aHelix.pt();
	double px    = aHelix.momentum(0).x();
	double py    = aHelix.momentum(0).y();
	double pz    = aHelix.momentum(0).z();
	double p     = aHelix.momentum(0).mag();
	double theta = aHelix.direction(0).theta();
	double phi   = aHelix.direction(0).phi();
	HepPoint3D poca  = aHelix.x(0);
	HepPoint3D pivot = aHelix.pivot();
	double r = poca.perp();
	HepSymMatrix Ea = aHelix.Ea();
	//cout<<"Ea="<<Ea<<endl;
	double errorMat[15];
	int k = 0;
	for (int ie = 0 ; ie < 5 ; ie ++){
		for (int je = ie ; je < 5 ; je ++){
			errorMat[k] = Ea[ie][je];
			k++;
		}
	}
	double chisq    = myDotsHelixFitter.getChi2();
	recMdcTrack->setCharge(q);
	recMdcTrack->setPxy(pxy);
	recMdcTrack->setPx(px);
	recMdcTrack->setPy(py);
	recMdcTrack->setPz(pz);
	recMdcTrack->setP(p);
	recMdcTrack->setTheta(theta);
	recMdcTrack->setPhi(phi);
	recMdcTrack->setPoca(poca);
	recMdcTrack->setX(poca.x());//poca
	recMdcTrack->setY(poca.y());
	recMdcTrack->setZ(poca.z());
	recMdcTrack->setR(sqrt(poca.x()*poca.x() + poca.y()*poca.y()));
	HepPoint3D apivot(0.,0.,0.);
	recMdcTrack->setPivot(apivot);
	recMdcTrack->setVX0(0.);//pivot 
	recMdcTrack->setVY0(0.);
	recMdcTrack->setVZ0(0.);
	recMdcTrack->setError(errorMat);
	recMdcTrack->setError(Ea);
	recMdcTrack->setChi2(chisq);
	recMdcTrack->setStat(tkStat);

	int maxLayerId = -1;
	int minLayerId = 43;
	double fiTerm = 0.;
	double fltLen = -0.00001;
	int    layerMaxFltLen=-1;

	// --- CGEM clusters
	// ClusterRefVec clusterRefVec;
	// map<int,int> clusterFitStat;
	// SmartDataPtr<RecCgemClusterCol> aCgemClusterCol(eventSvc(),"/Event/Recon/RecCgemClusterCol");
	// if(!aCgemClusterCol) cout<<"DotsConnection::saveARecMdcTrack() does not find RecCgemClusterCol!"<<endl;
	// RecCgemClusterCol::iterator iter_cluster_begin=aCgemClusterCol->begin();
	// RecCgemClusterCol::iterator iter_cluster=iter_cluster_begin;
	// for(; iter_cluster!=aCgemClusterCol->end(); iter_cluster++)
	// {
	// 	int flag = (*iter_cluster)->getflag();
	// 	if(flag!=2) continue;// skip 1D clusters
	// 	int layer=(*iter_cluster)->getlayerid();
	// 	int sheet=(*iter_cluster)->getsheetid();
	// 	int idXClu = (*iter_cluster)->getclusterflagb();
	// 	bool matchX=false;
	// 	vector<int>::iterator iter=myVecCgemXCluIdx[layer][sheet].begin();
	// 	for(; iter!=myVecCgemXCluIdx[layer][sheet].end(); iter++)
	// 		if((*iter)==idXClu) matchX=true;
	// 	if(!matchX) continue;
	// 	int idVClu = (*iter_cluster)->getclusterflage();
	// 	bool matchV=false;
	// 	iter=myVecCgemVCluIdx[layer][sheet].begin();
	// 	for(; iter!=myVecCgemVCluIdx[layer][sheet].end(); iter++)
	// 		if((*iter)==idVClu) matchV=true;
	// 	if(matchV) 
	// 	{
	// 		const RecCgemCluster* recCgemCluster = (*iter_cluster);
	// 		int clusterid = recCgemCluster->getclusterid();
	// 		clusterRefVec.push_back(recCgemCluster);
	// 		clusterFitStat[clusterid] = 1;
	// 		if(maxLayerId<layer)
	// 		{
	// 			maxLayerId=layer;
	// 		}
	// 	}
	// }


	// --- MDC hits
	int hitId = 0;
	HitRefVec  hitRefVec;
	vector<RecMdcHit> aRecMdcHitVec=myDotsHelixFitter.makeRecMdcHitVec(1);
	int nMdcHits=aRecMdcHitVec.size();
	int nMdcHitsKept=0;
	vector<RecMdcHit>::iterator iter_recMdcHit = aRecMdcHitVec.begin();
	for(; iter_recMdcHit!=aRecMdcHitVec.end(); iter_recMdcHit++)
	{
		if(iter_recMdcHit->getChisqAdd()>myMdcHitChi2Cut) // skip hit with chi2>myMdcHitChi2Cut
			continue;

		RecMdcHit* recMdcHit = new RecMdcHit(*iter_recMdcHit);
		recMdcHit->setId(hitId);
		recMdcHit->setTrkId(trackId);
		recMdcHit->setStat(1);
		myRecMdcHitCol->push_back(recMdcHit);
		SmartRef<RecMdcHit> refHit(recMdcHit);
		hitRefVec.push_back(refHit);
		nMdcHitsKept++;

		Identifier mdcid = recMdcHit->getMdcId();
		int layer = MdcID::layer(mdcid);
		if(layer>maxLayerId)
		{
			maxLayerId = layer;
		}
		if(layer<minLayerId)
		{
			minLayerId = layer;
		}
		if(fltLen<recMdcHit->getFltLen()) {
			fltLen=recMdcHit->getFltLen();
			layerMaxFltLen=layer;
		}
		hitId++;
	}
	if(myDebug)
		cout<<"track "<<trackId<<", "<<nMdcHitsKept<<"/"<<nMdcHits<<" hits kept"<<endl;

	// --- phi term (phi for the outmost hit/cluster)
	// if(maxLayerId>=0&&maxLayerId<3) {
	// 	double rmax=myDotsHelixFitter.getRmidGapCgem(maxLayerId);
	// 	fltLen=aHelix.flightLength(rmax);
	// }
	if(fltLen>0) fiTerm=-fltLen*sin(theta)/aHelix.radius();
	else cout<<"fltLen<0!"<<endl;
	recMdcTrack->setFiTerm(fiTerm);

	// --- check fiTerm
	HepPoint3D posOut = aHelix.x(fiTerm);
	if(myDebug)
		cout<<"fiTerm, layer, fltLen, pos= "<<fiTerm<<", "<<layerMaxFltLen<<", "<<fltLen<<", "<<posOut<<endl;

	// --- setN* functions called in setVec* functions
	recMdcTrack->setVecHits(hitRefVec);
	std::sort(clusterRefVec.begin(),clusterRefVec.end(),sortCluster);
	recMdcTrack->setVecClusters(clusterRefVec,clusterFitStat);
	myRecMdcTrackCol->push_back(recMdcTrack);


	return true;
}

