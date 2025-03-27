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

#include "DataFormats/MuonDetId/interface/DTChamberId.h" //zad 17
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h" //Ph2
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTThContainer.h" //Ph2

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h" //zad 7

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"



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
    //hLicExample->Write(); // Example
    hPhiB_st1->Write(); //zad 27
    hPhiB_st2->Write(); //zad 20
    hPhiBCompSt1->Write(); //zad 21
    hPhiBCompSt2->Write(); //zad 26
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
  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputDTPh_upg; //Ph2
  edm::EDGetTokenT<L1Phase2MuDTThContainer> inputDTTh_upg; //Ph2

  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
//  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<std::vector<PSimHit> > inputPSimHit; //zad 7
  
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeomteryToken;

  bool debug;
  unsigned int theEventCnt;
  unsigned int theAllDtPDigisCnt, theAllDtTDigisCnt;
  //TH2D *hLicExample; //Example

  TH2D *hPhiB_st1; //zad 27
  TH2D *hPhiB_st2; //zad 20
  TH2D *hPhiBCompSt1; //zad 21
  TH2D *hPhiBCompSt2; //zad 26
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
  inputDTPh_upg = consumes<L1Phase2MuDTPhContainer>(cfg.getParameter<edm::InputTag>("srcDTPh_upg")); //Ph2
  inputDTTh_upg = consumes<L1Phase2MuDTThContainer>(cfg.getParameter<edm::InputTag>("srcDTPh_upg")); //Ph2
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputPSimHit = consumes<std::vector<PSimHit> >(edm::InputTag("g4SimHits", "MuonDTHits")); //zad 7 

  theGeomteryToken=esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>();

 

  //hLicExample = new TH2D("hLicExample","hLicExample", 12,0.5,12.5, 8,-0.5,7.5); //Example
  hPhiB_st1 = new TH2D("hPhiB_st1", "PhiB in a function of pt st 1", 4096, 0, 100, 2000, -0.5, 0.5); //zad 27
  hPhiB_st2 = new TH2D("hPhiB_st2", "PhiB in a function of pt st 2", 4096, 0, 100, 2000, -0.5, 0.5); //zad 20
  hPhiBCompSt1 = new TH2D("hPhiBCompSt1", "PhiB Comparison at station 1", 1000, -1, 1, 1000, -1, 1); //zad 21
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
void LicDigiAnalysis::analyzeDT( const edm::Event &ev, const edm::EventSetup& es) {
  if (debug) std::cout << "-------- Tracking Particles -----------" << std::endl;
  edm::Handle<TrackingParticleCollection> tpColl;
  const GlobalTrackingGeometry & globalGeometry = es.getData(theGeomteryToken);
  ev.getByToken(inputTP, tpColl);
  const TrackingParticleCollection & myTP = *(tpColl.product());
  const std::vector<PSimHit> & myPSimHits = ev.get(inputPSimHit); 

  if (debug){ 
    std::cout<<" TRACKING PARTICLES: " << myTP.size() << std::endl;
    std::cout << "myPSimHits Size: " << myPSimHits.size() << std::endl;
    for (const auto & ah: myPSimHits) {
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
  int i_hits = 0; 
  int ft = 0; //zad 20
  int ftst2 = 0; //zad 27
  int ftt = 0; //zad 21
  double PhiB_Sim_St1 = 0; //zad 21
  double Phi_Sim_St1 = 0; //zad 25.1
  double Phi_Sim_St2 = 0; //zad 25.2
  int ftt2 = 0; //zad 25.2
  double PhiB_Sim_St2 = 0; // zad 26

  for(const auto & ah: myPSimHits) {
    int station_P = 0; //zad 17
    int station_S = 0; //zad 17
    i_hits++;
    if (ah.trackId() != 1) continue;
    const GeomDet * geomDet = globalGeometry.idToDet(ah.detUnitId()); //geographicalId() -> z DTChamberID.h DTChamberId
    DTChamberId chamber(geomDet->geographicalId()); //zad 17
    if(debug) std::cout << "Station_P: " << chamber.station() << std::endl;

    DTChamberId dtChamberId(ah.detUnitId()); //zad 17
    
    if(station_P!=chamber.station() && station_S!=dtChamberId.station()){
      station_P=chamber.station();
      station_S=dtChamberId.station();
    }

    GlobalVector globalMomentumPSimHit = geomDet->toGlobal(ah.momentumAtEntry());   
    GlobalPoint entry = geomDet->toGlobal(ah.entryPoint());
    GlobalPoint position = geomDet->toGlobal(ah.localPosition());

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

  edm::Handle<L1MuDTChambPhContainer> digiCollectionDTPh_leg;
  ev.getByToken(inputDTPh_leg, digiCollectionDTPh_leg);
  const L1MuDTChambPhContainer& dtphDigisLeg= *digiCollectionDTPh_leg.product();
  if (debug) std::cout <<" DTPh digis from BMTF " << dtphDigisLeg.getContainer()->size()<< std::endl; //Container Size
  for (const auto &  chDigi : *dtphDigisLeg.getContainer() ) {

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
    //hLicExample->Fill(chDigi.scNum()+1,chDigi.code()); //Example

  }

  edm::Handle<L1Phase2MuDTPhContainer> digiCollectionDTPh_upg;
  ev.getByToken(inputDTPh_upg, digiCollectionDTPh_upg);
  const L1Phase2MuDTPhContainer& dtphDigisUpg= *digiCollectionDTPh_upg.product();
  if (debug) std::cout <<" Upgrade DTPh digis from BMTF " << dtphDigisUpg.getContainer()->size()<< std::endl;
  for (const auto &  chDigi : *dtphDigisUpg.getContainer() ) {
//    if (abs(chDigi.whNum()) != 2) continue;
//    if (chDigi.stNum() ==4) continue;
//    if (chDigi.bxNum() != 0) continue;
    if (debug) std::cout <<"DtDataWord64 BMTF    " 
        <<" bxNum: "<<chDigi.bxNum()
        <<" whNum: "<<chDigi.whNum()
        <<" station: "<< chDigi.stNum()
        <<" sector: "<< chDigi.scNum()
        <<" phi:   "<<chDigi.phi()
        << std::endl;
  }

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
