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
    hPvSX1->Write(); //zad 18.1
    hPvSX2->Write(); //zad 18.2
    hPvSY1->Write(); //zad 18.3
    hPvSY2->Write(); //zad 18.4
    hHelp->Write(); //zad 19.0
    hPtX->Write(); //zad 19
    hSPt->Write(); //zad 20
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
  TH1D *hEta; //zad 3
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
  TH1D *hPvSX1; //zad 18.1
  TH1D *hPvSX2; //zad 18.2
  TH1D *hPvSY1; //zad 18.3
  TH1D *hPvSY2; //zad 18.4
  TH1D *hHelp; //zad 19.0
  TH2D *hPtX; //zad 19
  TH2D *hSPt; //zad 20

 
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
  hEta = new TH1D("hEta", "Muon Pseudorapidity", 100, -2, 2); //zad 3
  hVxy = new TH2D("hVxy", "Muon Vertex X and Y", 1000, -0.005, 0.005, 1000, -0.005, 0.005); //zad 4
  hT = new TH2D("hT", "Determining the time's unit - ns", 1000, 0, 1e-9, 1000, 0, 1e-9); //zad 5 
  hPSimHit = new TH1D("hPSimHit", "Size of numberOfHits", 100, 0, 100); //zad 7
  hPSimTrackerHit = new TH1D("hPSimTrackerHit", "Size of numberOfTrackerHits", 100, 0, 100); //zad 7
  hPSimHitVector = new TH1D("hPSimHitVector", "Size of PSimHit DT-Collection from Vector", 200, 0, 200); //zad 7
  hPSimHitXYZ = new TH3D("hPSimHitXYZ", "Position of PSimHit", 240, -120, 120, 200, -100, 100, 200, -1, 1); //zad 8
  hPPGvS = new TH1D("hPPGvS", "The difference between positions obtained from propagation and simulation", 40, 0, 400); //zad 11
  hPSimHitRZ = new TH2D("hPSimHitRZ", "PSimHit global position at entry", 300, 400, 700, 800, -200, -1000); //zad 13
  hPSimHitXY = new TH2D("hPSimHitXY", "PSimHit global position at entry", 400, -1000, 1000, 400, -1000, 1000); //zad 16
  hVPPGPT1 = new TH2D("hVPPGPT1", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.1
  hVPPGPT2 = new TH2D("hVPPGPT2", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.2
  hVSPT1 = new TH2D("hVSPT1", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.3
  hVSPT2 = new TH2D("hVSPT2", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.4
  h1Dtest = new TH1D("h1Dtest", "Transverse momentum for specific pt", 100, 0.75, 1); //zad 17.5
  hPvSX1 = new TH1D("hPvSX1", "The difference in X plane position obtained from propagation and simulation for 5-6 GeV Pt", 30, -60, 60); //zad 18.1
  hPvSX2 = new TH1D("hPvSX2", "The difference in X plane position obtained from propagation and simulation for 20-25 GeV Pt", 60, -60, 60); //zad 18.2
  hPvSY1 = new TH1D("hPvSY1", "The difference in Y plane position obtained from propagation and simulation for 5-6 GeV Pt", 30, -60, 60); //zad 18.3
  hPvSY2 = new TH1D("hPvSY2", "The difference in Y plane position obtained from propagation and simulation for 20-25 GeV Pt", 60, -60, 60); //zad 18.4
  hHelp = new TH1D("hHelp", "", 60, -60, 60); //zad 19.0
  hPtX = new TH2D("hPtX", "Full width at half maximum of X histogram in relation to Pt value", 100, 0, 100, 10, -5, 5); //zad 19
  hSPt = new TH2D("hSPt", "Spread in relation to Pt", 100, 0, 100, 100, 0, 100); //zad 20
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
  for(const auto & ah: myPSimHits) {
    int station_P = 0; //zad 17
    int station_S = 0; //zad 17
    i_hits++;
    int i_iterations = 0;
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
    if(i_hits==1 && station_S == 1 && 25 <  tp.pt() && tp.pt() < 45 && i_iterations == 0) {
      //i_mom = tp.pt();
      h1Dtest->Fill(sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ) / tp.pt() );
      i_iterations++;
    }
    //std::cout << "i_mom = " << i_mom << std::endl;
    //zad 18.1
    if(tp.pt()>5 && tp.pt()<6){
      if(chamber.station()==2 && dtChamberId.station()==2){
        hPvSX1->Fill(ah.localPosition().x()-stateAtDet.localPosition().x());
      }
    }
    //zad 18.2
    if(tp.pt()>20 && tp.pt()<25){
      if(chamber.station()==2 && dtChamberId.station()==2){
        hPvSX2->Fill(ah.localPosition().x()-stateAtDet.localPosition().x());
      }
    }
    //zad 18.3
    if(tp.pt()>5 && tp.pt()<6){
      if(chamber.station()==2 && dtChamberId.station()==2){
        hPvSY1->Fill(ah.localPosition().y()-stateAtDet.localPosition().y());
      }
    }
    //zad 18.4
    if(tp.pt()>20 && tp.pt()<25){
      if(chamber.station()==2 && dtChamberId.station()==2){
        hPvSY2->Fill(ah.localPosition().y()-stateAtDet.localPosition().y());
      }
    }
    
    //zad 19
    hHelp->Fill(ah.localPosition().x()-stateAtDet.localPosition().x());
    
    

    if(debug) std::cout << "Station: " << dtChamberId.station() << std::endl;
    if(debug) {
      std::cout<<" Simulated momentum at entry: " << ah.pabs() <<  std::endl; //zad 12
      std::cout<<" Propagated momentum at entry: " << sqrt( pow(stateAtDet.localMomentum().x(), 2) + pow(stateAtDet.localMomentum().y(), 2) + pow(stateAtDet.localMomentum().z(), 2)) << std::endl; 
    }
    GlobalPoint globalEntryHisto = geomDet->toGlobal(ah.localPosition());
    hPSimHitRZ->Fill(sqrt(pow(globalEntryHisto.x(), 2) + pow(globalEntryHisto.y(), 2)), globalEntryHisto.z()); //zad 13 
    hPSimHitXY->Fill(globalEntryHisto.x(), globalEntryHisto.y()); //zad 16 
        
  }


  /*
  for (const auto & tp : myTP) {
    int station_P = 0; //zad 17
    int station_S = 0; //zad 17
    //int pt_P = 0; //zad 17
    //int pt_S = 0; //zad 17

    //#PPG
    if ( abs( tp.pdgId())==13  && tp.pt() > 1. && gtp.charge()==0) {
      int i_hits = 0;
      int i_iterations = 0;
      double i_mom = 0;
      for(const auto & ah: myPSimHits) {
        if (ah.trackId() != 1) continue;
        i_hits++;

        const auto & vtx = tp.parentVertex()->position();
        const auto & mom = tp.momentum();
        gtp=GlobalTrajectoryParameters(GlobalPoint(vtx.x(), vtx.y(), vtx.z()), GlobalVector(mom.x(), mom.y(), mom.z()), tp.charge(), field);
        //std::cout << "gtp = " << gtp << std::endl;
        FreeTrajectoryState fts(gtp);
        
      }
    } 
    //zad 19
    int nBins = hHelp->GetNbinsX();
    double halfMax = hHelp->GetMaximum() / 2;
    int leftBin = -1;
    int rightBin = -1;
    for (int i = 1; i <= nBins; ++i) {
        if (hHelp->GetBinContent(i) >= halfMax) {
            leftBin = i;
        }
        if (hHelp->GetBinContent(i) >= halfMax) {
            rightBin = i;
        }
    }
    double FWHM = hHelp->GetBinCenter(rightBin) - hHelp->GetBinCenter(leftBin);
    hPtX->Fill(tp.pt(), FWHM);
    hHelp->Reset(); 

    if(abs(tp.pdgId())==13){
      const auto & simTracks = tp.g4Tracks();
      if (simTracks.size() > 0) {
        if(debug) std::cout << simTracks.size() << ", muon Track Id: " << simTracks[0].trackId() << std::endl;
      }
      hPt->Fill(tp.pt()); //zad 1
      hVz->Fill(tp.parentVertex()->position().z()); //zad 2
      hEta->Fill(tp.eta()); //zad 3
      hVxy->Fill(tp.parentVertex()->position().x(), tp.parentVertex()->position().y()); //zad 4 
    }
    if(myTP.size()==2){ //zad 5 
      double spdl = 29979245800; //speed of light cm/s
      double z_muon = 0;
      double rho_e = 0;
      double z_e = 0;
      if(abs(tp.pdgId())==13){
        z_muon = tp.parentVertex()->position().z();
      }
      if(abs(tp.pdgId())==11){
        rho_e = tp.parentVertex()->position().Rho();
        z_e = tp.parentVertex()->position().z();
      }
      
      hT->Fill(tp.parentVertex()->position().T(), sqrt(pow(rho_e, 2)+pow(z_muon - z_e, 2))/ spdl); 
    }
    hPSimHit->Fill(tp.numberOfHits()); //zad 7 Equivalent to trackPSimHit().size()
    hPSimTrackerHit->Fill(tp.numberOfTrackerHits()); //zad 7 Equivalent to trackPSimHit(DetId::Tracker).size()
    hPSimHitVector->Fill(myPSimHits.size()); //zad 7

  
//    if ( abs( tp.pdgId())!=13  || tp.pt() < 1. || tp.parentVertex()->position().Rho()>200. ||  fabs(tp.parentVertex()->position().T()) > 24.) continue;
    
    if (debug) std::cout << print(tp)<<std::endl; 
  }
  */
  if (debug) std::cout << "-------- HERE DIGI COMPARE DT ---------" << std::endl;
  //std::cout << "gtp w digi = " << gtp << std::endl;
  edm::Handle<L1MuDTChambPhContainer> digiCollectionDTPh_leg;
  ev.getByToken(inputDTPh_leg, digiCollectionDTPh_leg);
  const L1MuDTChambPhContainer& dtphDigisLeg= *digiCollectionDTPh_leg.product();
  if (debug) std::cout <<" DTPh digis from BMTF " << dtphDigisLeg.getContainer()->size()<< std::endl; //Container Size
  for (const auto &  chDigi : *dtphDigisLeg.getContainer() ) {
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

    //Implementation of the Propagation Code #PPG
  /*
    if (gtp.charge() !=0) {  //petla po simhit powinna byc
      std::cout<<" Generated starting Global position: "<< gtp.position() << std::endl;
      FreeTrajectoryState fts(gtp);
      DTChamberId dtChamberId(chDigi.whNum(),chDigi.stNum(),chDigi.scNum()+1);
      GlobalPoint detPosition = globalGeometry.idToDet(dtChamberId)->position(); //zamiast dtChamberId->detUnitId 
      std::cout<<" Position of the Detector: "<< detPosition << std::endl;
      //const DTChamber * chamber = geomDt.chamber(dtChamberId);
      const GeomDet * geomDet = globalGeometry.idToDet(dtChamberId);
      const Propagator & propagator = es.getData(thePropagatorToken);
      TrajectoryStateOnSurface stateAtDet =  propagator.propagate(fts, geomDet->surface());
      std::cout <<" Is PROPAGATION VALID: "<<stateAtDet.isValid() <<std::endl;
      LocalPoint localPositionBefore;
      GlobalPoint globalPositionBefore;
      if (stateAtDet.isValid() == 1) {
        localPositionBefore = stateAtDet.localPosition();
        globalPositionBefore = stateAtDet.globalPosition();
        std::cout<<" Local position after propagation: "<< localPositionBefore <<std::endl;
        std::cout<<" Global position after propagation: "<< globalPositionBefore <<std::endl; 
        // std::cout<<" The difference between Global position from stateAtDet and the one calculated from toGlobal Function: "<< geomDet->toGlobal(stateAtDet.localPosition()) - stateAtDet.globalPosition() << std::endl; //zad 9
      //   for(int i=0; i<200; ++i) { //zad 10 mozna usunac, porownujemy w ramach okreslonego detektora, trzeba sprawdzac czy trackId=1
      //     TrajectoryStateOnSurface stateAtDet =  propagator.propagate(fts, geomDet->surface());
      //     //std::cout<< "Debug global: " << stateAtDet.globalPosition() << std::endl;
      //     if(stateAtDet.isValid() == 1) continue;
      //   }
      //  std::cout<<" Local position change after 200 propagations: "<< stateAtDet.localPosition() - localPositionBefore <<std::endl;
      //  std::cout<<" Global position change after 200 propagations: "<< stateAtDet.globalPosition() - globalPositionBefore <<std::endl;  
      }
      int i = 0;
      int j=0;
      for (const auto & ah: myPSimHits) {
        if(ah.trackId()==1){
        j++;
        }
      }
      //std::cout<< j << std::endl;
      GlobalPoint simulatedPosition;
      for (const auto & ah: myPSimHits) {
        if(ah.trackId() == 1){
          auto globalEntryHisto = geomDet->toGlobal(ah.entryPoint());
          hPSimHitRZ->Fill(sqrt(pow(globalEntryHisto.x(), 2) + pow(globalEntryHisto.y(), 2)), globalEntryHisto.z()); //zad 13 
        }
        if(ah.trackId()==1 && i==0) {
          std::cout<< "PSimHit Local Entry Point: " << ah.entryPoint() << std::endl;
          std::cout<< "PSimHit Global Entry Point: " << geomDet->toGlobal(ah.entryPoint()) << std::endl;
        }
        if(ah.trackId()==1 && i==j-1) {
          std::cout<< "PSimHit Local Exit Point: " << ah.exitPoint() << std::endl;
          std::cout<< "PSimHit Global Exit Point: " << geomDet->toGlobal(ah.exitPoint()) << std::endl;
          simulatedPosition = geomDet->toGlobal(ah.exitPoint());
        }
        i++;
      }
      //std::cout<< "test= " << simulatedPosition << std::endl;
      hPPGvS->Fill(globalPositionBefore.x() - simulatedPosition.x(), globalPositionBefore.y() - simulatedPosition.y(), globalPositionBefore.z() - simulatedPosition.z()); //zad 11
      double simulatedEnergyLoss = 0;
      for (const auto & ah: myPSimHits) {
      simulatedEnergyLoss = simulatedEnergyLoss + ah.energyLoss(); //to nie sa wszystkie kroki, straty energie brac z pedu pbabs
      //std::cout<<"The energy deposit in the PSimHit: "<< ah.energyLoss() << std::endl;
      }
      for(const auto & tp: myTP) {
      std::cout << "Brief energy from the first SimTrack: " << tp.energy() << std::endl; 
      }
      std::cout << "Simulated total energy loss: " << simulatedEnergyLoss << std::endl; //zad 12
      std::cout << "Propagateted total energy loss: " << "For now lack of such data." << std::endl;
      //end of function
    }*/
  }


/*
  edm::Handle<L1MuDTChambThContainer> digiCollectionDTTh_BMTF;
  ev.getByToken(inputDTTh_BMTF, digiCollectionDTTh_BMTF);
  const L1MuDTChambThContainer& dtthDigisBMTF = *digiCollectionDTTh_BMTF.product();
  if (debug) std::cout <<" DTTh digis from BMTF " << dtthDigisBMTF.getContainer()->size()<< std::endl;
  for (const auto &  chDigi : *dtthDigisBMTF.getContainer() ) {
    unsigned int eta = 0;
    unsigned int etaQ = 0;
    for (unsigned int ipos=0; ipos <7; ipos++) {
     if (chDigi.position(ipos) >1 ) if (debug) std::cout <<" HERE PROBLEM !!!!" << std::endl;
     if (chDigi.position(ipos)==1) eta |= (1 <<ipos);
     if (chDigi.quality(ipos)==1) etaQ |= (1 <<ipos);
    }
    if (chDigi.bxNum() != 0) continue;
//    if (abs(chDigi.bxNum()) >2) continue;
    if (eta==0 || etaQ==0) continue;
    if (abs(chDigi.whNum()) != 2) continue;
    if (debug) std::cout <<"DtDataWord64 BMTF TH " 
        <<" bxNum: "<<chDigi.bxNum()
        <<" whNum: "<<chDigi.whNum()
        <<" station: "<< chDigi.stNum()
        <<" sector: "<< chDigi.scNum()
        <<" eta: " << eta
        <<" etaQ: " << etaQ
        << std::endl;
  }
*/
}

DEFINE_FWK_MODULE(LicDigiAnalysis);
