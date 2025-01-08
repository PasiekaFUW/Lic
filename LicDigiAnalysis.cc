// system include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <bitset>
#include <map>
#include <string>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h" // A lot of content 

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

//#include "DataFormats/MuonDetId/interface/DtDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h" //zad 17
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h" //zad 7
//#PDG
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
//#PDG ^

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TH2D.h"
#include "TH1D.h" //zad 1
#include "TFile.h"
#include <cmath> //zad 5
#include "TH3D.h" //zad 8



const TrackingParticle & ancestor(const TrackingParticle & particle) {

  const TrackingVertexRef&  tpv = particle.parentVertex();
  if (tpv->nSourceTracks() == 0) return particle;
  const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
  return ancestor(parent);
}

std::string print(const TrackingParticle & tp) {
  std::stringstream ss;
  ss <<" pid: "<<tp.pdgId()
     <<" pt: "<<tp.pt()
     <<" eta: "<<tp.eta()
     <<" phi: "<<tp.phi()
     <<" vtx[r,z]:  ["<<tp.parentVertex()->position().Rho() <<", "<<tp.parentVertex()->position().z()<<"]"
     <<" time: "<<tp.parentVertex()->position().T()
     ;
  return ss.str();
}
 





class LicDigiAnalysis : public edm::one::EDAnalyzer<> {
public:
  LicDigiAnalysis (const edm::ParameterSet & cfg);
  virtual ~LicDigiAnalysis(){ 
    TFile f("licDigiHistos.root","RECREATE");
    f.cd();
    hLicExample->Write(); // Example
    hPt->Write(); //zad 1
    hVz->Write(); //zad 2
    hEta->Write(); //zad 3
    hVxy->Write(); //zad 4
    hT->Write(); //zad 5
    hPSimHit->Write(); //zad 7.1
    hPSimTrackerHit->Write(); //zad 7.2
    hPSimHitVector->Write(); //zad 7.3
    hPSimHitXYZ->Write(); //zad 8
    hPPGvS->Write(); //zad 11
    hPSimHitRZ->Write(); //zad 13
    hPSimHitXY->Write(); //zad 16
    hVPPGPT1->Write(); //zad 17.1
    hVPPGPT2->Write(); //zad 17.2
    hVSPT1->Write(); //zad 17.3
    hVSPT2->Write(); //zad 17.4
    h1Dtest->Write(); //zad 17.5
    h2Dtest->Write(); //zad 18.0
    hLandau->Write(); //zad 18
    hPhiComp->Write(); //zad 19
    hPhiB_st1->Write(); //zad 27
    hPhiB_st2->Write(); //zad 20
    hPhiBCompSt1->Write(); //zad 21
    hPhiBCompSt2->Write(); //zad 26
    hHowMany1->Write(); //zad 22
    hCheck0->Write(); //zad 22.5
    hPhiCompSt1->Write(); //zad 25.1
    hPhiCompSt2->Write(); //zad 25.2
    hDeltaPhiB1->Write(); //zad 28.1
    hDeltaPhiB2->Write(); //zad 28.2
    hDeltaPhi1->Write(); //zad 29.1
    hDeltaPhi2->Write(); //zad 29.2
    hDeltaBCodeSt1->Write(); //zad 30.1
    hDeltaBCodeSt2->Write(); //zad 30.2


    f.Write();


    std::cout << "AT THE END: "; printStat();
  }

  virtual void analyze(const edm::Event &ev, const edm::EventSetup& es) {
    theEventCnt++;
    debug = 0; //Debug Flag 0 - off, 1 - on
    analyzeDT(ev,es);

    if (theEventCnt%100==1) printStat(); 
  }
  void analyzeDT(const edm::Event&, const edm::EventSetup& es);
  void printStat();

private:
  edm::EDGetTokenT<L1MuDTChambPhContainer> inputDTPh_leg;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputDTTh_leg;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
//  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<std::vector<PSimHit> > inputPSimHit; //zad 7

  //#PPG
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeomteryToken;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> theFieldToken;
  edm::ESGetToken<Propagator, TrackingComponentsRecord> thePropagatorToken;

  bool debug;
  unsigned int theEventCnt;
  unsigned int theAllDtPDigisCnt, theAllDtTDigisCnt;
  TH2D *hLicExample; //Example
  TH1D *hPt; //zad 1
  TH1D *hVz; //zad 2
  TH2D *hEta; //zad 3
  TH2D *hVxy; //zad 4
  TH2D *hT; //zad 5
  TH1D *hPSimHit; //zad 7
  TH1D *hPSimTrackerHit; //zad 7
  TH1D *hPSimHitVector; //zad 7
  TH3D *hPSimHitXYZ; //zad 8
  TH1D *hPPGvS; //zad 11
  TH2D *hPSimHitRZ; //zad 13
  TH2D *hPSimHitXY; //zad 16
  TH2D *hVPPGPT1; //zad 17.1
  TH2D *hVPPGPT2; //zad 17.2
  TH2D *hVSPT1; //zad 17.3
  TH2D *hVSPT2; //zad 17.4
  TH1D *h1Dtest; //zad 17.5
  TH1D *h2Dtest; //zad 18.0
  TH1D *hLandau; //zad 18
  TH2D *hPhiComp; //zad 19
  TH2D *hPhiB_st1; //zad 27
  TH2D *hPhiB_st2; //zad 20
  TH2D *hPhiBCompSt1; //zad 21
  TH2D *hPhiBCompSt2; //zad 26
  TH1D *hHowMany1; //zad 22
  TH2D *hCheck0; //zad 22.5
  TH2D *hPhiCompSt1; //zad 25.1
  TH2D *hPhiCompSt2; //zad 25.2
  TH1D *hDeltaPhiB1; //zad 28.1
  TH1D *hDeltaPhiB2; //zad 28.2
  TH1D *hDeltaPhi1; //zad 29.1
  TH1D *hDeltaPhi2; //zad 29.2
  TH2D *hDeltaBCodeSt1; //zad 30.1
  TH2D *hDeltaBCodeSt2; //zad 30.2


};

void LicDigiAnalysis::printStat()
{
   std::cout<<std::dec <<"=========> Analyzed #"<<theEventCnt
            <<" DtPh: "<<theAllDtPDigisCnt <<"  DtTh: "<<theAllDtTDigisCnt <<std::endl; 
}

LicDigiAnalysis::LicDigiAnalysis(const edm::ParameterSet & cfg) 
  : debug (false), theEventCnt(0),
    theAllDtPDigisCnt(0),
    theAllDtTDigisCnt(0)
{
  inputDTPh_leg = consumes<L1MuDTChambPhContainer>(cfg.getParameter<edm::InputTag>("srcDTPh_leg"));
  inputDTTh_leg = consumes<L1MuDTChambThContainer>(cfg.getParameter<edm::InputTag>("srcDTTh_leg"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputPSimHit = consumes<std::vector<PSimHit> >(edm::InputTag("g4SimHits", "MuonDTHits")); //zad 7 

  //#PPG
  theGeomteryToken=esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>();
  theFieldToken=esConsumes<MagneticField, IdealMagneticFieldRecord>();   
  thePropagatorToken=esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("","SteppingHelixPropagatorAny")); 
 

  hLicExample = new TH2D("hLicExample","hLicExample", 12,0.5,12.5, 8,-0.5,7.5); //Example
  hPt = new TH1D("hPt", "Muon Transverse Momentum", 105, 0, 105); //zad 1
  hVz = new TH1D("hVz", "Muon Vertex Z", 30, -15, 15); //zad 2
  hEta = new TH2D("hEta", "Muon Pseudorapidity", 100, -2, 2, 100, -2, 2); //zad 3
  hVxy = new TH2D("hVxy", "Muon Vertex X and Y", 1000, -0.005, 0.005, 1000, -0.005, 0.005); //zad 4
  hT = new TH2D("hT", "Determining the time's unit - ns", 1000, 0, 1e-9, 1000, 0, 1e-9); //zad 5 
  hPSimHit = new TH1D("hPSimHit", "Size of numberOfHits", 100, 0, 100); //zad 7
  hPSimTrackerHit = new TH1D("hPSimTrackerHit", "Size of numberOfTrackerHits", 100, 0, 100); //zad 7
  hPSimHitVector = new TH1D("hPSimHitVector", "Size of PSimHit DT-Collection from Vector", 200, 0, 200); //zad 7
  hPSimHitXYZ = new TH3D("hPSimHitXYZ", "Position of PSimHit", 240, -120, 120, 200, -100, 100, 200, -1, 1); //zad 8
  hPPGvS = new TH1D("hPPGvS", "The difference between positions obtained from propagation and simulation", 40, 0, 400); //zad 11
  hPSimHitRZ = new TH2D("hPSimHitRZ", "PSimHit global position at entry", 400, -1000, 1000, 400, -1000, -1000); //zad 13
  hPSimHitXY = new TH2D("hPSimHitXY", "PSimHit global position at entry", 400, -1000, 1000, 400, -1000, 1000); //zad 16
  hVPPGPT1 = new TH2D("hVPPGPT1", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.1
  hVPPGPT2 = new TH2D("hVPPGPT2", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.2
  hVSPT1 = new TH2D("hVSPT1", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.3
  hVSPT2 = new TH2D("hVSPT2", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.4
  h1Dtest = new TH1D("h1Dtest", "Transverse momentum for specific pt station 1", 100, 0.75, 1); //zad 17.5
  h2Dtest = new TH1D("h2Dtest", "Transverse momentum for specific pt station 2", 100, 0.75, 1); //zad zad 18.0
  hLandau = new TH1D("hLandau", "Landau Curve of tp loss between Stations 1-2", 60, 0, 3); //zad 18
  hPhiComp = new TH2D("hPhiComp", "Comparison of Phi and PhiB", 100, -3.15, 3.15, 50, -1, 1); //zad 19
  hPhiB_st1 = new TH2D("hPhiB_st1", "PhiB in a function of pt st 1", 4096, 0, 100, 2000, -0.5, 0.5); //zad 27
  hPhiB_st2 = new TH2D("hPhiB_st2", "PhiB in a function of pt st 2", 4096, 0, 100, 2000, -0.5, 0.5); //zad 20
  hPhiBCompSt1 = new TH2D("hPhiBCompSt1", "PhiB Comparison at station 1", 1000, -1, 1, 1000, -1, 1); //zad 21
  hHowMany1 = new TH1D("hHowMany1", "Checking", 4096, -2500, 1596); //zad 22
  hCheck0 = new TH2D("hCheck0", "Checking", 10, 0, 10, 10, 0, 10); //zad 22.5
  hPhiCompSt1 = new TH2D("hPhiCompSt1", "Comparing phi at station 1", 4096, -4, 4, 4096, -4, 4); //zad 25.1
  hPhiCompSt2 = new TH2D("hPhiCompSt2", "Comparing phi at station 2", 4096, -4, 4, 4096, -4, 4); //zad 25.2
  hPhiBCompSt2 = new TH2D("hPhiBCompSt2", "PhiB Comparison at station 2", 1000, -1, 1, 1000, -1, 1); //zad 26
  hDeltaPhiB1 = new TH1D("hDeltaPhiB1", "Delta PhiB at station 1 entry", 2000, -0.1, 0.1); //zad 28.1
  hDeltaPhiB2 = new TH1D("hDeltaPhiB2", "Delta PhiB at station 2 entry", 2000, -0.1, 0.1); //zad 28.2
  hDeltaPhi1 = new TH1D("hDeltaPhi1", "Delta Phi at station 1 entry", 2000, -0.1, 0.1); //zad 29.1
  hDeltaPhi2 = new TH1D("hDeltaPhi2", "Delta Phi at station 2 entry", 2000, -0.1, 0.1); //zad 29.2
  hDeltaBCodeSt1 = new TH2D("hDeltaBCodeSt1", "Delta PhiB in the chDigi.code() variable function, st1", 10, 0, 9, 2000, -0.1, 0.1); //zad 30.1
  hDeltaBCodeSt2 = new TH2D("hDeltaBCodeSt2", "Delta PhiB in the chDigi.code() variable function, st2", 10, 0, 9, 2000, -0.1, 0.1); //zad 30.2

 
}
int FailedPPG = 0; //Debuging variable
void LicDigiAnalysis::analyzeDT( const edm::Event &ev, const edm::EventSetup& es) {
  if (debug) std::cout << "-------- Tracking Particles -----------" << std::endl;
  //#PPG
  GlobalTrajectoryParameters gtp;
  const MagneticField * field = &es.getData(theFieldToken);
  const GlobalTrackingGeometry & globalGeometry = es.getData(theGeomteryToken);


  edm::Handle<TrackingParticleCollection> tpColl;
  ev.getByToken(inputTP, tpColl);
  const TrackingParticleCollection & myTP = *(tpColl.product());
  const std::vector<PSimHit> & myPSimHits = ev.get(inputPSimHit); 

  if (debug){ 
    std::cout<<" TRACKING PARTICLES: " << myTP.size() << std::endl;
    std::cout << "myPSimHits Size: " << myPSimHits.size() << std::endl;
    for (const auto & ah: myPSimHits) {
      hPSimHitXYZ->Fill(ah.localPosition().x(), ah.localPosition().y(), ah.localPosition().z()); //zad 8
      std::cout << ah.trackId() << std::endl;
    }
  }

  const TrackingParticle *myMuon = 0;
  for (const auto & tp : myTP) {
    if ( abs( tp.pdgId())==13  && tp.pt() > 1.) {
      myMuon = & tp;
      break;
    }
  } 
  const TrackingParticle & tp = *myMuon;
  const auto & vtx = tp.parentVertex()->position();
  const auto & mom = tp.momentum();
  gtp=GlobalTrajectoryParameters(GlobalPoint(vtx.x(), vtx.y(), vtx.z()), GlobalVector(mom.x(), mom.y(), mom.z()), tp.charge(), field);
  FreeTrajectoryState fts(gtp);

  int i_hits = 0; 
  int firsttime = 0; //zad 17
  double firstmom = 0; //zad 18
  int count = 0; //zad 19
  int ft = 0; //zad 20
  int ftst2 = 0; //zad 27
  int ftt = 0; //zad 21
  //int hits = 0; //zad 22
  double PhiB_Sim_St1 = 0; //zad 21
  double Phi_Sim_St1 = 0; //zad 25.1
  double Phi_Sim_St2 = 0; //zad 25.2
  int ftt2 = 0; //zad 25.2
  double PhiB_Sim_St2 = 0; // zad 26


  hEta->Fill(tp.eta() , tp.eta() * tp.charge()); //zad 3 

  for(const auto & ah: myPSimHits) {
    int station_P = 0; //zad 17
    int station_S = 0; //zad 17
    i_hits++;
     //std::cout << "ile hitow " << i_hits << std::endl;
    if (ah.trackId() != 1) continue;
    //GlobalPoint trackPosition = globalGeometry.idToDet(ah.detUnitId())->position(); //not track but rather a detector
    //if(debug) std::cout << "trackPosition: " << trackPosition << std::endl; 
    const GeomDet * geomDet = globalGeometry.idToDet(ah.detUnitId()); //geographicalId() -> z DTChamberID.h DTChamberId
    DTChamberId chamber(geomDet->geographicalId()); //zad 17
    if(debug) std::cout << "Station_P: " << chamber.station() << std::endl;
    const Propagator & propagator = es.getData(thePropagatorToken);
    TrajectoryStateOnSurface stateAtDet =  propagator.propagate(fts, geomDet->surface());
    if(stateAtDet.isValid() == 0) {
        ++FailedPPG;
        if(debug) std::cout<< " Number of failed propagations: " << FailedPPG << std::endl;            
        continue;
    } 
    if(debug) {
      std::cout <<" Is PROPAGATION VALID: "<<stateAtDet.isValid() <<std::endl;
      if(stateAtDet.isValid() == 1){
        std::cout<<" Local position from PSimHit: " << ah.localPosition() << std::endl;
        std::cout<<" Local position after propagation: "<< stateAtDet.localPosition() <<std::endl; 
        std::cout<<" The difference: "<< ah.localPosition() - stateAtDet.localPosition() <<std::endl;
      }
    }
    
    hPPGvS->Fill(sqrt( pow(ah.localPosition().x() - stateAtDet.localPosition().x(), 2) + pow(ah.localPosition().y() - stateAtDet.localPosition().y(), 2)));//, ah.localPosition().z() - stateAtDet.localPosition().z()); //zad 11
    
    DTChamberId dtChamberId(ah.detUnitId()); //zad 17
    
    if(station_P!=chamber.station() && station_S!=dtChamberId.station()){
      station_P=chamber.station();
      station_S=dtChamberId.station();
      if(station_P!=1){
        //hPHPT->Fill(pt_P - sqrt(pow(stateAtDet.localMomentum().x(), 2) + pow(stateAtDet.localMomentum().y(), 2) + pow(stateAtDet.localMomentum().z(), 2)), pt_S - ah.pabs());
      }
      //pt_P=sqrt(pow(stateAtDet.localMomentum().x(), 2) + pow(stateAtDet.localMomentum().y(), 2) + pow(stateAtDet.localMomentum().z(), 2));
      //pt_S=ah.pabs();
    }

    //zad 17.1 tp.numberOfTrackerHits() == 1
    if(i_hits==1 && station_P == 1) {
    hVPPGPT1->Fill(tp.pt(), sqrt(pow(stateAtDet.globalMomentum().x(), 2) + pow(stateAtDet.globalMomentum().y(), 2)  ));
    
    //+ pow(stateAtDet.localMomentum().z(), 2))); 
    }

    //zad 17.2
    if(i_hits==1 && station_P == 2) {
    hVPPGPT2->Fill(tp.pt(), sqrt(pow(stateAtDet.globalMomentum().x(), 2) + pow(stateAtDet.globalMomentum().y(), 2))); 
    }
    //zad 17.3 
    GlobalVector globalMomentumPSimHit = geomDet->toGlobal(ah.momentumAtEntry());
    if(i_hits==1 && station_S == 1) {
    //hVSPT1->Fill(tp.pt(), ah.pabs()); 
    hVSPT1->Fill(tp.pt(), sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) );
    }
    //zad 17.4 
    if(i_hits==1 && station_S == 2) {
    //hVSPT2->Fill(tp.pt(), ah.pabs()); 
    hVSPT2->Fill(tp.pt(), sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) );
    }
    //zad 17.5
    if(firsttime==0 && station_S == 1 && 25 <  tp.pt() && tp.pt() < 45 ) {
      h1Dtest->Fill(sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) / tp.pt() );
      firstmom = sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) );
      firsttime++;
    }
   
   
    //zad 18
    if(firsttime==1 && station_S == 2) {
      h2Dtest->Fill(sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) /tp.pt() );
      hLandau->Fill( firstmom - sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) );
      firsttime++;
    }
    
    GlobalPoint entry = geomDet->toGlobal(ah.entryPoint());
    GlobalPoint position = geomDet->toGlobal(ah.localPosition());

    //GlobalPoint exit = geomDet->toGlobal(ah.exitPoint());
    //GlobalPoint globalHit = geomDet->toGlobal(ah.localPosition());
    //zad 19 
    if(tp.pt() < 50 && station_S == 1 && count == 0) {
      count++;
      //hPhiComp -> Fill(ah.phiAtEntry(), globalMomentumPSimHit.phi());
      //hPhiComp -> Fill(globalHit.phi(), abs(exit.phi() - entry.phi()));
      //std::cout << "zad 19: phi= " << entry.phi() << "phiB= " << globalMomentumPSimHit.phi() - entry.phi() << std::endl;
      hPhiComp -> Fill(entry.phi(), globalMomentumPSimHit.phi() - entry.phi());

    }

       //zad 27
    if(station_S==1 && ftst2 == 0){
      ftst2++; 
      hPhiB_st1 -> Fill(tp.pt(), globalMomentumPSimHit.phi() - entry.phi() );
    }

    //zad 20
    if(station_S==2 && ft == 0){
      ft++; 
      //hPhiB -> Fill( sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ), abs(exit.phi() - entry.phi()) );
      hPhiB_st2 -> Fill(tp.pt(), globalMomentumPSimHit.phi() - entry.phi() );

    }

    if(ftt == 0 && station_S == 1){
      PhiB_Sim_St1 = globalMomentumPSimHit.phi() - entry.phi(); //zad 21
      Phi_Sim_St1 = position.phi(); //zad 25.1
      ftt++;
    }

    if(ftt2 == 0 && station_S ==2){
      ftt2++;
      PhiB_Sim_St2 = globalMomentumPSimHit.phi() - entry.phi(); //zad 26
      Phi_Sim_St2 = position.phi(); //zad 25.2
    }




    

    if(debug) std::cout << "Station: " << dtChamberId.station() << std::endl;
    if(debug) {
      std::cout<<" Simulated momentum at entry: " << ah.pabs() <<  std::endl; //zad 12
      std::cout<<" Propagated momentum at entry: " << sqrt( pow(stateAtDet.localMomentum().x(), 2) + pow(stateAtDet.localMomentum().y(), 2) + pow(stateAtDet.localMomentum().z(), 2)) << std::endl; 
    }
    GlobalPoint globalEntryHisto = geomDet->toGlobal(ah.localPosition());
    hPSimHitRZ->Fill(sqrt(pow(globalEntryHisto.x(), 2) + pow(globalEntryHisto.y(), 2)), globalEntryHisto.z()); //zad 13 
    hPSimHitXY->Fill(globalEntryHisto.x(), globalEntryHisto.y()); //zad 16 
        
  }


  double PhiB_Rec_St1 = 0; //zad 21
  double Phi_Rec_St1 = 0; //zad 25.1
  double Phi_Rec_St2 = 0; //zad 25.2
  int ftt_rec_1 = 0; //zad 21
  int ftt_rec_2 = 0; //zad 25.2
  double PhiB_Rec_St2 = 0; //zad 26
  int codeSt1 = 0; //zad 30.1
  int codeSt2 = 0; //zad 30.2
  if (debug) std::cout << "-------- HERE DIGI COMPARE DT ---------" << std::endl;
  //std::cout << "gtp w digi = " << gtp << std::endl;
  edm::Handle<L1MuDTChambPhContainer> digiCollectionDTPh_leg;
  ev.getByToken(inputDTPh_leg, digiCollectionDTPh_leg);
  const L1MuDTChambPhContainer& dtphDigisLeg= *digiCollectionDTPh_leg.product();
  if (debug) std::cout <<" DTPh digis from BMTF " << dtphDigisLeg.getContainer()->size()<< std::endl; //Container Size
  for (const auto &  chDigi : *dtphDigisLeg.getContainer() ) {
    //zad 22
    if(chDigi.stNum()==1 && chDigi.scNum()==2) {
      hHowMany1->Fill(chDigi.phi());
    }


    if(chDigi.stNum() == 1 && ftt_rec_1 ==0){
      ftt_rec_1++;
      PhiB_Rec_St1 = chDigi.phiB(); //zad 21
      PhiB_Rec_St1 /= 512.;
      Phi_Rec_St1 = chDigi.phi();
      Phi_Rec_St1 /= 4096.;
      //ScNum nie powinien być modyfikowany, a całość nie jest przesunięta o Pi/12
      Phi_Rec_St1 += M_PI/6. * (chDigi.scNum() + 0) ; //zad 25.1
      if(Phi_Rec_St1 > M_PI){
        Phi_Rec_St1 -= 2.*M_PI;
      }
      //std::cout << chDigi.scNum() << std::endl;
      //scNum() from 0 to 11
      codeSt1 = chDigi.code(); //zad 30.1
    }
    if(chDigi.stNum() == 2 && ftt_rec_2 == 0){
      ftt_rec_2++;
      PhiB_Rec_St2 = chDigi.phiB(); // / 512; //zad 26
      PhiB_Rec_St2 /= 512.;
      Phi_Rec_St2 = chDigi.phi();
      Phi_Rec_St2 /= 4096.;
      Phi_Rec_St2 += M_PI/6. * (chDigi.scNum() + 0) ;//zad 25.2
      if(Phi_Rec_St2 > M_PI){
        Phi_Rec_St2 -= 2.*M_PI;
      }
      codeSt2 = chDigi.code(); //zad 30.2
    }
    //std::cout << chDigi.phiB() << "phiB" << std::endl;
    if (abs(chDigi.whNum()) != 2) continue;
    if (chDigi.stNum() ==4) continue;
    if (chDigi.bxNum() != 0) continue;
    if (chDigi.code()==7) continue;
    // DTChamberId chId(chDigi.whNum(),chDigi.stNum(),chDigi.scNum()+1);
    theAllDtPDigisCnt++;
    if (debug) std::cout <<"DtDataWord64 BMTF    " 
        <<" bxNum: "<<chDigi.bxNum()
        <<" whNum: "<<chDigi.whNum()
        <<" station: "<< chDigi.stNum()
        <<" sector: "<< chDigi.scNum()
        <<" phi:   "<<chDigi.phi() //pozycja
        <<" phiB:   "<<chDigi.phiB() // kierunek
        <<" code(q): "<< chDigi.code()
        << std::endl;
    hLicExample->Fill(chDigi.scNum()+1,chDigi.code()); //Example


  }
  //zad 22.5
  //hCheck0->Fill(PhiB_Rec.size(), PhiB_Sim.size());
  //zad 21 dla stacji 1 i 2 osobno, phiB zdefiniowac jako phi z pedu - phi z polozenia
  //for(size_t k = 0; k < std::min(PhiB_Sim.size(), PhiB_Rec.size()); ++k){
    //hPhiBComp -> Fill(PhiB_Rec[k], PhiB_Sim[k]); //zad 21
  //}
  if(Phi_Rec_St1 != 0 && Phi_Sim_St1 != 0) {
    hPhiCompSt1 -> Fill(Phi_Sim_St1, Phi_Rec_St1); //zad 25.1
    //hDeltaPhi1 -> Fill(Phi_Sim_St1 - Phi_Rec_St1); //zad 29.1
    hDeltaPhi1 -> Fill(reco::deltaPhi(Phi_Sim_St1, Phi_Rec_St1)); //zad 29.1

  }
  if(Phi_Rec_St2 != 0 && Phi_Sim_St2 != 0) {
    hPhiCompSt2 -> Fill(Phi_Sim_St2, Phi_Rec_St2); //zad 25.2
    hDeltaPhi2 -> Fill(reco::deltaPhi(Phi_Sim_St2, Phi_Rec_St2)); //zad 29.2
  }
 if(PhiB_Rec_St1 != 0 && PhiB_Sim_St1 != 0) {
    hPhiBCompSt1 -> Fill(PhiB_Sim_St1, PhiB_Rec_St1); //zad 21
   // hDeltaPhiB1 -> Fill(PhiB_Sim_St1 - PhiB_Rec_St1); //zad 28.1
    hDeltaPhiB1 -> Fill(reco::deltaPhi(PhiB_Sim_St1, PhiB_Rec_St1)); //zad 29.1
    hDeltaBCodeSt1 -> Fill(codeSt1, reco::deltaPhi(PhiB_Sim_St1, PhiB_Rec_St1) ); //zad 30.1

  }
  if(PhiB_Rec_St2 != 0 && PhiB_Sim_St2 != 0) {
   hPhiBCompSt2 -> Fill(PhiB_Sim_St2, PhiB_Rec_St2); //zad 26
   hDeltaPhiB2 -> Fill(reco::deltaPhi(PhiB_Sim_St2, PhiB_Rec_St2)); //zad 28.2
   hDeltaBCodeSt2 -> Fill(codeSt2, reco::deltaPhi(PhiB_Sim_St2, PhiB_Rec_St2) ); //zad 30.2
  }
 
 

}

DEFINE_FWK_MODULE(LicDigiAnalysis);
