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
    hPhiB->Write(); //zad 20
    hPhiBComp->Write(); //zad 21
    hHowMany1->Write(); //zad 22
    hCheck0->Write(); //zad 22.5
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
  TH2D *hPhiB; //zad 20
  TH2D *hPhiBComp; //zad 21
  TH1D *hHowMany1; //zad 22
  TH2D *hCheck0; //zad 22.5
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
  hPSimHitRZ = new TH2D("hPSimHitRZ", "PSimHit global position at entry", 300, 400, 700, 800, -200, -1000); //zad 13
  hPSimHitXY = new TH2D("hPSimHitXY", "PSimHit global position at entry", 400, -1000, 1000, 400, -1000, 1000); //zad 16
  hVPPGPT1 = new TH2D("hVPPGPT1", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.1
  hVPPGPT2 = new TH2D("hVPPGPT2", "Comparison of Tranverse Momentum from Propagation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.2
  hVSPT1 = new TH2D("hVSPT1", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 1 entry", 80, 0, 80, 80, 0, 80); //zad 17.3
  hVSPT2 = new TH2D("hVSPT2", "Comparison of Tranverse Momentum from Simulation and Vertex at Station 2 entry", 80, 0, 80, 80, 0, 80); //zad 17.4
  h1Dtest = new TH1D("h1Dtest", "Transverse momentum for specific pt station 1", 100, 0.75, 1); //zad 17.5
  h2Dtest = new TH1D("h2Dtest", "Transverse momentum for specific pt station 2", 100, 0.75, 1); //zad zad 18.0
  hLandau = new TH1D("hLandau", "Landau Curve of tp loss between Stations 1-2", 60, 0, 3); //zad 18
  hPhiComp = new TH2D("hPhiComp", "Comparison of Phi and PhiB", 100, -3.15, 3.15, 50, -1, 1); //zad 19
  hPhiB = new TH2D("hPhiB", "PhiB in a function of pt", 100, 0, 100, 200, -1, 1); //zad 20
  hPhiBComp = new TH2D("hPhiBComp", "PhiB Comparison", 1000, -150, 400, 200, -0.2, 0.8); //zad 21
  hHowMany1 = new TH1D("hHowMany1", "Number of hits at station 1", 100, 0, 100); //zad 22
  hCheck0 = new TH2D("hCheck0", "Checking vectors size", 10, 0, 10, 10, 0, 10); //zad 22.5

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
  int ftt = 0; //zad 21
  int hits = 0; //zad 22
  std::vector<double> PhiB_Sim; //zad 21

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

    //zad 20
    if(station_S==2 && ft == 0){
      ft++; 
      //hPhiB -> Fill( sqrt( pow(globalMomentumPSimHit.x(), 2) + pow(globalMomentumPSimHit.y(), 2) ), abs(exit.phi() - entry.phi()) );
      hPhiB -> Fill(tp.pt(), globalMomentumPSimHit.phi() - entry.phi() );

    }

    if(ftt == 0 && station_S == 1){
      PhiB_Sim.push_back(globalMomentumPSimHit.phi() - entry.phi()); //zad 21
      ftt++;
    }

    //zad 22
    if(station_S == 1 && tp.pt() > 50){
      hits++;
      hHowMany1 -> Fill(hits);
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


  std::vector<double> PhiB_Rec; //zad 21
  int ftt_rec_1 = 0; //zad 21
  if (debug) std::cout << "-------- HERE DIGI COMPARE DT ---------" << std::endl;
  //std::cout << "gtp w digi = " << gtp << std::endl;
  edm::Handle<L1MuDTChambPhContainer> digiCollectionDTPh_leg;
  ev.getByToken(inputDTPh_leg, digiCollectionDTPh_leg);
  const L1MuDTChambPhContainer& dtphDigisLeg= *digiCollectionDTPh_leg.product();
  if (debug) std::cout <<" DTPh digis from BMTF " << dtphDigisLeg.getContainer()->size()<< std::endl; //Container Size
  for (const auto &  chDigi : *dtphDigisLeg.getContainer() ) {
    if(chDigi.stNum() == 1 && ftt_rec_1 ==0){
      ftt_rec_1++;
      PhiB_Rec.push_back(chDigi.phiB()); //zad 21
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
  hCheck0->Fill(PhiB_Rec.size(), PhiB_Sim.size());
  //zad 21 dla stacji 1 i 2 osobno, phiB zdefiniowac jako phi z pedu - phi z polozenia
  for(size_t k = 0; k < std::min(PhiB_Sim.size(), PhiB_Rec.size()); ++k){
    hPhiBComp -> Fill(PhiB_Rec[k], PhiB_Sim[k]); //zad 21
  }
  

}

DEFINE_FWK_MODULE(LicDigiAnalysis);
