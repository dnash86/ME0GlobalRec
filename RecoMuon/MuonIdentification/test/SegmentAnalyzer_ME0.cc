#include "FWCore/Framework/interface/Event.h"

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/MuonReco/interface/EmulatedME0Segment.h>
#include <DataFormats/MuonReco/interface/EmulatedME0SegmentCollection.h>

#include <DataFormats/MuonReco/interface/RealME0Muon.h>
#include <DataFormats/MuonReco/interface/RealME0MuonCollection.h>

// #include "CLHEP/Matrix/SymMatrix.h"
// #include "CLHEP/Matrix/Matrix.h"
// #include "CLHEP/Vector/ThreeVector.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "TRandom3.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"
//#include <deltaR.h>
#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "TH1.h" 
#include <TH2.h>
#include "TFile.h"
#include <TProfile.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include <Geometry/GEMGeometry/interface/ME0EtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/ME0DetId.h>



class SegmentAnalyzer_ME0 : public edm::EDAnalyzer {
public:
  explicit SegmentAnalyzer_ME0(const edm::ParameterSet&);
  ~SegmentAnalyzer_ME0();


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginJob();

  //protected:
  
  private:
  TFile* histoFile;
  TH1F *Segment_Eta;    TH1F *Segment_Phi;    TH1F *Segment_R;  TH2F *Segment_Pos;  
  TH1F *Rechit_Eta;    TH1F *Rechit_Phi;    TH1F *Rechit_R;  TH2F *Rechit_Pos;  
  TH1F  *GenMuon_Eta;   TH1F *GenMuon_Phi;    TH1F *GenMuon_R;  TH2F *GenMuon_Pos;    TH1F *GenMuon_Pt;
 

//Removing this
};

SegmentAnalyzer_ME0::SegmentAnalyzer_ME0(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
}



void SegmentAnalyzer_ME0::beginJob()
{
  

  Segment_Eta = new TH1F("Segment_Eta"      , "Segment #eta"   , 40, 2.4, 4.0 );
  Segment_Phi = new TH1F("Segment_Phi"      , "Segment #phi"   , 60, -3, 3. );
  Segment_R = new TH1F("Segment_R"      , "Segment r"   , 30, 0, 150 );
  Segment_Pos = new TH2F("Segment_Pos"      , "Segment x,y"   ,100,-100.,100., 100,-100.,100. );

  Rechit_Eta = new TH1F("Rechit_Eta"      , "Rechit #eta"   , 40, 2.4, 4.0 );
  Rechit_Phi = new TH1F("Rechit_Phi"      , "Rechit #phi"   , 60, -3, 3. );
  Rechit_R = new TH1F("Rechit_R"      , "Rechit r"   , 30, 0, 150 );
  Rechit_Pos = new TH2F("Rechit_Pos"      , "Rechit x,y"   ,100,-100.,100., 100,-100.,100. );

  GenMuon_Eta = new TH1F("GenMuon_Eta"      , "GenMuon #eta"   , 40, 2.4, 4.0 );
  GenMuon_Phi = new TH1F("GenMuon_Phi"      , "GenMuon #phi"   , 60, -3, 3. );
  GenMuon_R = new TH1F("GenMuon_R"      , "GenMuon r"   , 30, 0, 150 );
  GenMuon_Pos = new TH2F("GenMuon_Pos"      , "GenMuon x,y"   ,100,-100.,100., 100,-100.,100. );
  GenMuon_Pt = new TH1F("GenMuon_Pt"      , "Muon p_{T}"   , 80,0 , 8. );


}


SegmentAnalyzer_ME0::~SegmentAnalyzer_ME0(){}

void
SegmentAnalyzer_ME0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;

  //run_ = (int)iEvent.id().run();
  //event_ = (int)iEvent.id().event();


    //David's functionality
    

  using namespace reco;

  

  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel<GenParticleCollection>("genParticles", genParticles);

  Handle<ME0SegmentCollection> OurSegments;
  iEvent.getByLabel("me0Segments","",OurSegments);


  edm::ESHandle<ME0Geometry> me0Geom;
  iSetup.get<MuonGeometryRecord>().get(me0Geom);

 
  unsigned int gensize=genParticles->size();
  for(unsigned int i=0; i<gensize; ++i) {
    const reco::GenParticle& CurrentParticle=(*genParticles)[i];
    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  

      GenMuon_Eta->Fill(CurrentParticle.eta());
      GenMuon_Phi->Fill(CurrentParticle.phi());
      GenMuon_Pt->Fill(CurrentParticle.pt());

    }
  }


  //================  For Segment Plotting
  //std::cout <<"Number of Segments "<<OurSegments.size()<<std::endl;
  for (auto thisSegment = OurSegments->begin(); thisSegment != OurSegments->end(); 
       ++thisSegment){
    ME0DetId id = thisSegment->me0DetId();
    //std::cout<<"ME0DetId =  "<<id<<std::endl;

    auto roll = me0Geom->etaPartition(id); 
    auto GlobVect(roll->toGlobal(thisSegment->localPosition()));
    Segment_Eta->Fill(GlobVect.eta());
    Segment_Phi->Fill(GlobVect.phi());
    Segment_R->Fill(GlobVect.perp());
    Segment_Pos->Fill(GlobVect.x(),GlobVect.y());

    auto theseRecHits = thisSegment->specificRecHits();
    std::cout <<"ME0 Ensemble Det Id "<<id<<"  Number of RecHits "<<theseRecHits.size()<<std::endl;
    
    for (auto thisRecHit = theseRecHits.begin(); thisRecHit!= theseRecHits.end(); thisRecHit++){
      auto me0id = thisRecHit->me0Id();
      auto rollForRechit = me0Geom->etaPartition(me0id);
      
      auto thisRecHitGlobalPoint = rollForRechit->toGlobal(thisRecHit->localPosition()); 
      
      Rechit_Eta->Fill(thisRecHitGlobalPoint.eta());
      Rechit_Phi->Fill(thisRecHitGlobalPoint.phi());
      Rechit_R->Fill(thisRecHitGlobalPoint.perp());
      Rechit_Pos->Fill(thisRecHitGlobalPoint.x(),thisRecHitGlobalPoint.y());
      
    }
  }
  //==================
  

  
}

void SegmentAnalyzer_ME0::endJob() 
{
  histoFile->cd();
  TCanvas *c1 = new TCanvas("c1", "canvas" );

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  


  Segment_Phi->Write();   Segment_Phi->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Segment_Phi.png");
  Segment_R->Write();   Segment_R->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Segment_R.png");
  Segment_Pos->Write();   Segment_Pos->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Segment_Pos.png");


  Rechit_Phi->Write();   Rechit_Phi->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Rechit_Phi.png");
  Rechit_R->Write();   Rechit_R->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Rechit_R.png");
  Rechit_Pos->Write();   Rechit_Pos->Draw();  c1->Print("OutputSegmentAnalyzerPlots/Rechit_Pos.png");


  GenMuon_Eta->Write();   GenMuon_Eta->Draw();  c1->Print("OutputSegmentAnalyzerPlots/GenMuon_Eta.png");

  GenMuon_Pt->Write();   GenMuon_Pt->Draw();  c1->Print("OutputSegmentAnalyzerPlots/GenMuon_Pt.png");

  c1->SetLogy();
  Segment_Eta->Write();   Segment_Eta->Draw();  GenMuon_Eta->SetLineColor(2);GenMuon_Eta->Draw("SAME"); c1->Print("OutputSegmentAnalyzerPlots/Segment_Eta.png");
  Rechit_Eta->Write();   Rechit_Eta->Draw();   c1->Print("OutputSegmentAnalyzerPlots/Rechit_Eta.png");
  delete histoFile; histoFile = 0;
}

DEFINE_FWK_MODULE(SegmentAnalyzer_ME0);
