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



class TestAnalyzer_ME0 : public edm::EDAnalyzer {
public:
  explicit TestAnalyzer_ME0(const edm::ParameterSet&);
  ~TestAnalyzer_ME0();


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginJob();

  //protected:
  
  private:
  TFile* histoFile;
  TH1F *Candidate_Eta;  TH1F *Mass_h; 
  TH1F *Segment_Eta;    TH1F *Segment_Phi;    TH1F *Segment_R;  TH2F *Segment_Pos;  
  TH1F *Rechit_Eta;    TH1F *Rechit_Phi;    TH1F *Rechit_R;  TH2F *Rechit_Pos;  
     TH1F *GenMuon_Phi;    TH1F *GenMuon_R;  TH2F *GenMuon_Pos;  
  TH1F *Track_Eta; TH1F *Track_Pt;  TH1F *ME0Muon_Eta; TH1F *ME0Muon_Pt; 
  TH1F *UnmatchedME0Muon_Eta; TH1F *UnmatchedME0Muon_Pt; 
  TH1F *TracksPerSegment_h;  TH2F *TracksPerSegment_s;  TProfile *TracksPerSegment_p;
  TH2F *ClosestDelR_s; TProfile *ClosestDelR_p;
  TH2F *PtDiff_s; TProfile *PtDiff_p; TH1F *PtDiff_h; TH1F *QOverPtDiff_h; TH1F *PtDiff_rms;
  TH1F *VertexDiff_h;
  TH2F *PDiff_s; TProfile *PDiff_p; TH1F *PDiff_h;
  TH1F *FakeTracksPerSegment_h;  TH2F *FakeTracksPerSegment_s;  TProfile *FakeTracksPerSegment_p;
  TH1F *FakeTracksPerAssociatedSegment_h;  TH2F *FakeTracksPerAssociatedSegment_s;  TProfile *FakeTracksPerAssociatedSegment_p;
  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt;   TH1F *MatchedME0Muon_Eta; TH1F *MatchedME0Muon_Pt; 
  TH1F *MuonRecoEff_Eta;  TH1F *MuonRecoEff_Pt;
  TH1F *MuonAllTracksEff_Eta;  TH1F *MuonAllTracksEff_Pt;
  TH1F *MuonUnmatchedTracksEff_Eta;  TH1F *MuonUnmatchedTracksEff_Pt; TH1F *FractionMatched_Eta;


//Removing this
};

TestAnalyzer_ME0::TestAnalyzer_ME0(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
}



void TestAnalyzer_ME0::beginJob()
{
  Candidate_Eta = new TH1F("Candidate_Eta"      , "Candidate #eta"   , 40, 2.4, 4.0 );

  Track_Eta = new TH1F("Track_Eta"      , "Track #eta"   , 40, 2.4, 4.0 );
  Track_Pt = new TH1F("Track_Pt"      , "Muon p_{T}"   , 80,0 , 8. );

  Segment_Eta = new TH1F("Segment_Eta"      , "Segment #eta"   , 40, 2.4, 4.0 );
  Segment_Phi = new TH1F("Segment_Phi"      , "Segment #phi"   , 60, -3, 3. );
  Segment_R = new TH1F("Segment_R"      , "Segment r"   , 30, 0, 150 );
  Segment_Pos = new TH2F("Segment_Pos"      , "Segment x,y"   ,100,-100.,100., 100,-100.,100. );

  Rechit_Eta = new TH1F("Rechit_Eta"      , "Rechit #eta"   , 40, 2.4, 4.0 );
  Rechit_Phi = new TH1F("Rechit_Phi"      , "Rechit #phi"   , 60, -3, 3. );
  Rechit_R = new TH1F("Rechit_R"      , "Rechit r"   , 30, 0, 150 );
  Rechit_Pos = new TH2F("Rechit_Pos"      , "Rechit x,y"   ,100,-100.,100., 100,-100.,100. );

  //  GenMuon_Eta = new TH1F("GenMuon_Eta"      , "GenMuon #eta"   , 40, 2.4, 4.0 );
  GenMuon_Phi = new TH1F("GenMuon_Phi"      , "GenMuon #phi"   , 60, -3, 3. );
  GenMuon_R = new TH1F("GenMuon_R"      , "GenMuon r"   , 30, 0, 150 );
  GenMuon_Pos = new TH2F("GenMuon_Pos"      , "GenMuon x,y"   ,100,-100.,100., 100,-100.,100. );

  ME0Muon_Eta = new TH1F("ME0Muon_Eta"      , "Muon #eta"   , 40, 2.4, 4.0 );
  ME0Muon_Pt = new TH1F("ME0Muon_Pt"      , "Muon p_{T}"   , 80,0 , 8. );

  GenMuon_Eta = new TH1F("GenMuon_Eta"      , "Muon #eta"   , 40, 2.4, 4.0 );
  GenMuon_Pt = new TH1F("GenMuon_Pt"      , "Muon p_{T}"   , 80,0 , 8. );

  MatchedME0Muon_Eta = new TH1F("MatchedME0Muon_Eta"      , "Muon #eta"   , 40, 2.4, 4.0 );
  MatchedME0Muon_Pt = new TH1F("MatchedME0Muon_Pt"      , "Muon p_{T}"   , 8,0 , 40 );

  UnmatchedME0Muon_Eta = new TH1F("UnmatchedME0Muon_Eta"      , "Muon #eta"   , 40, 2.4, 4.0 );
  UnmatchedME0Muon_Pt = new TH1F("UnmatchedME0Muon_Pt"      , "Muon p_{T}"   , 8,0 , 40 );

  Mass_h = new TH1F("Mass_h"      , "Mass"   , 100, 0., 200 );

  MuonRecoEff_Eta = new TH1F("MuonRecoEff_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,40, 2.4, 4.0  );
  MuonRecoEff_Pt = new TH1F("MuonRecoEff_Pt"      , "Fraction of ME0Muons matched to gen muons"   ,8, 0,40  );

  MuonAllTracksEff_Eta = new TH1F("MuonAllTracksEff_Eta"      , "All ME0Muons over all tracks"   ,40, 2.4, 4.0  );
  MuonAllTracksEff_Pt = new TH1F("MuonAllTracksEff_Pt"      , "All ME0Muons over all tracks"   ,8, 0,40  );

  MuonUnmatchedTracksEff_Eta = new TH1F("MuonUnmatchedTracksEff_Eta"      , "Unmatched ME0Muons over all ME0Muons"   ,40, 2.4, 4.0  );
  MuonUnmatchedTracksEff_Pt = new TH1F("MuonUnmatchedTracksEff_Pt"      , "Unmatched ME0Muons over all ME0Muons"   ,8, 0,40  );

  TracksPerSegment_h = new TH1F("TracksPerSegment_h", "Number of tracks", 60,0.,60.);
  TracksPerSegment_s = new TH2F("TracksPerSegment_s" , "Tracks per segment vs |#eta|", 40, 2.4, 4.0, 60,0.,60.);
  TracksPerSegment_p = new TProfile("TracksPerSegment_p" , "Tracks per segment vs |#eta|", 40, 2.4, 4.0, 0.,60.);

  FakeTracksPerSegment_h = new TH1F("FakeTracksPerSegment_h", "Number of fake tracks", 60,0.,60.);
  FakeTracksPerSegment_s = new TH2F("FakeTracksPerSegment_s" , "Fake tracks per segment", 10, 2.4, 4.0, 100,0.,60.);
  FakeTracksPerSegment_p = new TProfile("FakeTracksPerSegment_p" , "Average N_{tracks}/segment not matched to genmuons", 10, 2.4, 4.0, 0.,60.);

  FakeTracksPerAssociatedSegment_h = new TH1F("FakeTracksPerAssociatedSegment_h", "Number of fake tracks", 60,0.,60.);
  FakeTracksPerAssociatedSegment_s = new TH2F("FakeTracksPerAssociatedSegment_s" , "Fake tracks per segment", 10, 2.4, 4.0, 100,0.,60.);
  FakeTracksPerAssociatedSegment_p = new TProfile("FakeTracksPerAssociatedSegment_p" , "Average N_{tracks}/segment not matched to genmuons", 10, 2.4, 4.0, 0.,60.);

  ClosestDelR_s = new TH2F("ClosestDelR_s" , "#Delta R", 40, 2.4, 4.0, 15,0.,0.15);
  ClosestDelR_p = new TProfile("ClosestDelR_p" , "#Delta R", 40, 2.4, 4.0, 0.,0.15);
  
  FractionMatched_Eta = new TH1F("FractionMatched_Eta"      , "Fraction of ME0Muons that end up successfully matched (matched/all)"   ,40, 2.4, 4.0  );

  PtDiff_s = new TH2F("PtDiff_s" , "Relative pt difference", 40, 2.4, 4.0, 50,0.,0.5);
  PtDiff_h = new TH1F("PtDiff_s" , "pt resolution", 100,-0.5,0.5);
  QOverPtDiff_h = new TH1F("PtDiff_s" , "q/pt resolution", 100,-0.5,0.5);
  PtDiff_p = new TProfile("PtDiff_p" , "pt resolution vs. #eta", 40, 2.4, 4.0, -1.0,1.0,"s");

  PtDiff_rms    = new TH1F( "PtDiff_rms",    "RMS", 40, 2.4, 4.0 ); 

  PDiff_s = new TH2F("PDiff_s" , "Relative p difference", 40, 2.4, 4.0, 50,0.,0.5);
  PDiff_h = new TH1F("PDiff_s" , "Relative p difference", 50,0.,0.5);
  PDiff_p = new TProfile("PDiff_p" , "Relative p difference", 40, 2.4, 4.0, 0.,1.0,"s");

  VertexDiff_h = new TH1F("VertexDiff_h", "Difference in vertex Z", 50, 0, 0.2);

}


TestAnalyzer_ME0::~TestAnalyzer_ME0(){}

void
TestAnalyzer_ME0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;

  //run_ = (int)iEvent.id().run();
  //event_ = (int)iEvent.id().event();


    //David's functionality
    

  using namespace reco;

  // Handle <ME0MuonCollection > OurMuons;
  // iEvent.getByLabel <ME0MuonCollection> ("me0SegmentMatcher", OurMuons);

  Handle <std::vector<RecoChargedCandidate> > OurCandidates;
  iEvent.getByLabel <std::vector<RecoChargedCandidate> > ("me0MuonConverter", OurCandidates);

  //Handle<std::vector<EmulatedME0Segment> > OurSegments;
  //iEvent.getByLabel<std::vector<EmulatedME0Segment> >("me0SegmentProducer", OurSegments);

  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel<GenParticleCollection>("genParticles", genParticles);

  Handle <TrackCollection > generalTracks;
  iEvent.getByLabel <TrackCollection> ("generalTracks", generalTracks);

  Handle <std::vector<RealME0Muon> > OurMuons;
  iEvent.getByLabel <std::vector<RealME0Muon> > ("me0SegmentMatcher", OurMuons);

  Handle<ME0SegmentCollection> OurSegments;
  iEvent.getByLabel("me0Segments","",OurSegments);


  edm::ESHandle<ME0Geometry> me0Geom;
  iSetup.get<MuonGeometryRecord>().get(me0Geom);

  //=====Finding RealME0Muons that match gen muons, plotting the closest of those
  //    -----First, make a vector of bools for each RealME0Muon

  std::vector<bool> IsMatched;
  std::vector<int> SegIdForMatch;
  for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
       thisMuon != OurMuons->end(); ++thisMuon){
    IsMatched.push_back(false);
    SegIdForMatch.push_back(-1);
  }
  std::cout<<IsMatched.size()<<" total me0muons"<<std::endl;
  //   -----Now, loop over each gen muon to compare it to each RealME0Muon
  //   -----For each gen muon, take the closest RealME0Muon that is a match within delR 0.15
  //   -----Each time a match on an RealME0Muon is made, change the IsMatched bool corresponding to it to true
  //   -----Also, each time a match on an RealME0Muon is made, we plot the pt and eta of the gen muon it was matched to

  // unsigned int gensize=genParticles->size();
  // for(unsigned int i=0; i<gensize; ++i) {
  //   const reco::GenParticle& CurrentParticle=(*genParticles)[i];
  //   if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  
  //     for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
  // 	   thisTrack != generalTracks->end();++thisTrack){
  // 	VertexDiff_h->Fill(fabs(thisTrack->vz()-CurrentParticle.vz()));
  // 	PtDiff_h->Fill(fabs(thisTrack->pt() - CurrentParticle.pt())/CurrentParticle.pt());
  // 	PtDiff_s->Fill(CurrentParticle.eta(),fabs(thisTrack->pt() - CurrentParticle.pt())/CurrentParticle.pt());
  // 	PtDiff_p->Fill(CurrentParticle.eta(),fabs(thisTrack->pt() - CurrentParticle.pt())/CurrentParticle.pt());
  // 	PDiff_h->Fill(fabs(thisTrack->p() - CurrentParticle.p())/CurrentParticle.p());
  // 	PDiff_s->Fill(CurrentParticle.eta(),fabs(thisTrack->p() - CurrentParticle.p())/CurrentParticle.p());
  // 	PDiff_p->Fill(CurrentParticle.eta(),fabs(thisTrack->p() - CurrentParticle.p())/CurrentParticle.p());
  //     }
  //   }
  // }
  //--------------------------------------------------------------

  std::vector<int> MatchedSegIds;
  unsigned int gensize=genParticles->size();
  for(unsigned int i=0; i<gensize; ++i) {
    const reco::GenParticle& CurrentParticle=(*genParticles)[i];
    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  

      GenMuon_Eta->Fill(CurrentParticle.eta());
      GenMuon_Phi->Fill(CurrentParticle.phi());

      //auto GlobVect(roll->toGlobal(Seg.localPosition()));
      //GenMuon_R->Fill(CurrentParticle.());
      
      if ( ((CurrentParticle.eta()) > 2.4) && ((CurrentParticle.eta()) < 3.2) ) GenMuon_Pt->Fill(CurrentParticle.pt());

      double LowestDelR = 9999;
      double thisDelR = 9999;
      int MatchedID = -1;
      int RealME0MuonID = 0;

      std::cout<<"Size = "<<OurMuons->size()<<std::endl;
      for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
	   thisMuon != OurMuons->end(); ++thisMuon){
	TrackRef tkRef = thisMuon->innerTrack();
	SegIdForMatch.push_back(thisMuon->me0segid());
	thisDelR = reco::deltaR(CurrentParticle,*tkRef);
	if (thisDelR < 0.15 ){
	  if (thisDelR < LowestDelR){
	    LowestDelR = thisDelR;
	    //if (fabs(tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt() < 0.50) MatchedID = RealME0MuonID;
	    MatchedID = RealME0MuonID;
	  }
	}
	VertexDiff_h->Fill(fabs(tkRef->vz()-CurrentParticle.vz()));
	PtDiff_h->Fill((tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt());	
	QOverPtDiff_h->Fill(( (tkRef->charge() /tkRef->pt()) - (CurrentParticle.charge()/CurrentParticle.pt() ) )/  (CurrentParticle.charge()/CurrentParticle.pt() ) );
	PtDiff_s->Fill(CurrentParticle.eta(),fabs(tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt());
	PtDiff_p->Fill(CurrentParticle.eta(),(tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt());
	for(int i=1; i<=PtDiff_p->GetNbinsX(); ++i) {
	  PtDiff_rms->SetBinContent(i, PtDiff_p->GetBinError(i)); 
	  }
	PDiff_h->Fill(fabs(tkRef->p() - CurrentParticle.p())/CurrentParticle.p());
	PDiff_s->Fill(CurrentParticle.eta(),fabs(tkRef->p() - CurrentParticle.p())/CurrentParticle.p());
	PDiff_p->Fill(CurrentParticle.eta(),fabs(tkRef->p() - CurrentParticle.p())/CurrentParticle.p());
	RealME0MuonID++;
      }
      if (MatchedID != -1){
	IsMatched[MatchedID] = true;
	MatchedME0Muon_Eta->Fill(CurrentParticle.eta());
	MatchedSegIds.push_back(SegIdForMatch[MatchedID]);	

	if ( ((CurrentParticle.eta()) > 2.4) && ((CurrentParticle.eta()) < 3.8) ) {
	  MatchedME0Muon_Pt->Fill(CurrentParticle.pt());
	  
	}
      }
    }
  }

  //Del R study ===========================
  for(unsigned int i=0; i<gensize; ++i) {
    const reco::GenParticle& CurrentParticle=(*genParticles)[i];
    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  

      double LowestDelR = 9999;
      double thisDelR = 9999;

      for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
	   thisMuon != OurMuons->end(); ++thisMuon){
	TrackRef tkRef = thisMuon->innerTrack();
	thisDelR = reco::deltaR(CurrentParticle,*tkRef);
	if (thisDelR < LowestDelR) LowestDelR = thisDelR;
      }
    
    ClosestDelR_s->Fill(CurrentParticle.eta(), LowestDelR);
    ClosestDelR_p->Fill(CurrentParticle.eta(), LowestDelR);
    }
  }

  //====================================

  //   -----Finally, we loop over all the RealME0Muons in the event
  //   -----Before, we plotted the gen muon pt and eta for the efficiency plot of matches
  //   -----Now, each time a match failed, we plot the RealME0Muon pt and eta
  int RealME0MuonID = 0;
  for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
       thisMuon != OurMuons->end(); ++thisMuon){
    if (!IsMatched[RealME0MuonID]){
      TrackRef tkRef = thisMuon->innerTrack();
      UnmatchedME0Muon_Eta->Fill(tkRef->eta());
      if ( (TMath::Abs(tkRef->eta()) > 2.4) && (TMath::Abs(tkRef->eta()) < 4.0) ) UnmatchedME0Muon_Pt->Fill(tkRef->pt());
    }
    RealME0MuonID++;
  }
  


  // for (std::vector<ME0Segment>::const_iterator thisSegment = OurSegments->begin();
  //      thisSegment != OurSegments->end();++thisSegment){
  //   LocalVector TempVect(thisSegment->localDirection().x(),thisSegment->localDirection().y(),thisSegment->localDirection().z());
  //   Segment_Eta->Fill(TempVect.eta());
  // }

  
  for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
       thisTrack != generalTracks->end();++thisTrack){
    Track_Eta->Fill(thisTrack->eta());
    if ( (TMath::Abs(thisTrack->eta()) > 2.4) && (TMath::Abs(thisTrack->eta()) < 4.0) ) Track_Pt->Fill(thisTrack->pt());
  }

  // for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
  //      thisMuon != OurMuons->end(); ++thisMuon){
  //   TrackRef tkRef = thisMuon->innerTrack();
  //   ME0Muon_Eta->Fill(tkRef->eta());
  //   if ( (TMath::Abs(tkRef->eta()) > 2.4) && (TMath::Abs(tkRef->eta()) < 4.0) ) ME0Muon_Pt->Fill(tkRef->pt());
  // }
  
  std::vector<double> SegmentEta, SegmentPhi, SegmentR, SegmentX, SegmentY;
  // std::vector<const ME0Segment*> Ids;
  // std::vector<const ME0Segment*> Ids_NonGenMuons;
  // std::vector<const ME0Segment*> UniqueIdList;
   std::vector<int> Ids;
   std::vector<int> Ids_NonGenMuons;
   std::vector<int> UniqueIdList;
   int TrackID=0;

   for (std::vector<RealME0Muon>::const_iterator thisMuon = OurMuons->begin();
	thisMuon != OurMuons->end(); ++thisMuon){
    TrackRef tkRef = thisMuon->innerTrack();
    //ME0Segment segRef = thisMuon->me0segment();
    //const ME0Segment* SegId = segRef->get();

    ME0Segment Seg = thisMuon->me0segment();
    ME0DetId id =Seg.me0DetId();
    auto roll = me0Geom->etaPartition(id); 
    auto GlobVect(roll->toGlobal(Seg.localPosition()));

    int SegId=thisMuon->me0segid();

    //std::cout<<SegId<<std::endl;

  
    bool IsNew = true;
    for (unsigned int i =0; i < Ids.size(); i++){
      if (SegId == Ids[i]) IsNew=false;
    }

    if (IsNew) {
      UniqueIdList.push_back(SegId);
      //std::cout<<"New SegId = "<<SegId<<std::endl;
      //std::cout<<GlobVect<<std::endl;
      SegmentEta.push_back(GlobVect.eta());
      SegmentPhi.push_back(GlobVect.phi());
      SegmentR.push_back(GlobVect.perp());
      SegmentX.push_back(GlobVect.x());
      SegmentY.push_back(GlobVect.y());
    }
    Ids.push_back(SegId);
    if (!IsMatched[TrackID]) Ids_NonGenMuons.push_back(SegId);

    ME0Muon_Eta->Fill(tkRef->eta());
    if ( (TMath::Abs(tkRef->eta()) > 2.4) && (TMath::Abs(tkRef->eta()) < 4.0) ) ME0Muon_Pt->Fill(tkRef->pt());

    TrackID++;
  }
  
   std::cout<<UniqueIdList.size()<<" unique segments per event"<<std::endl;
  for (unsigned int i = 0; i < UniqueIdList.size(); i++){
    int Num_Total=0, Num_Fake = 0, Num_Fake_Associated = 0;
    for (unsigned int j = 0; j < Ids.size(); j++){
      if (Ids[j] == UniqueIdList[i]) Num_Total++;
    }

    for (unsigned int j = 0; j < Ids_NonGenMuons.size(); j++){
      if (Ids_NonGenMuons[j] == UniqueIdList[i]) Num_Fake++;
      bool AssociatedWithMatchedSegment = false;
      for (unsigned int isegid=0;isegid < MatchedSegIds.size();isegid++){
	if (MatchedSegIds[isegid]==Ids_NonGenMuons[j]) AssociatedWithMatchedSegment=true;
      }
      if (AssociatedWithMatchedSegment) Num_Fake_Associated++;
    }

    TracksPerSegment_h->Fill((double)Num_Total);
    TracksPerSegment_s->Fill(SegmentEta[i], (double)Num_Total);
    TracksPerSegment_p->Fill(SegmentEta[i], (double)Num_Total);

    FakeTracksPerSegment_h->Fill((double)Num_Fake);
    FakeTracksPerSegment_s->Fill(SegmentEta[i], (double)Num_Fake);
    FakeTracksPerSegment_p->Fill(SegmentEta[i], (double)Num_Fake);

    FakeTracksPerAssociatedSegment_h->Fill((double)Num_Fake_Associated);
    FakeTracksPerAssociatedSegment_s->Fill(SegmentEta[i], (double)Num_Fake_Associated);
    FakeTracksPerAssociatedSegment_p->Fill(SegmentEta[i], (double)Num_Fake_Associated);

    // if (SegmentEta[i] > 2.4){
    //   Segment_Eta->Fill(SegmentEta[i]);
    //   Segment_Phi->Fill(SegmentPhi[i]);
    //   Segment_R->Fill(SegmentR[i]);
    //   Segment_Pos->Fill(SegmentX[i],SegmentY[i]);
    // }
  }

  //================  For Segment Plotting
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
    //std::cout <<"ME0 Ensemble Det Id "<<id<<"  Number of RecHits "<<theseRecHits.size()<<std::endl;
    
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
  

  //std::cout<<recosize<<std::endl;
  for (std::vector<RecoChargedCandidate>::const_iterator thisCandidate = OurCandidates->begin();
       thisCandidate != OurCandidates->end(); ++thisCandidate){
    TLorentzVector CandidateVector;
    CandidateVector.SetPtEtaPhiM(thisCandidate->pt(),thisCandidate->eta(),thisCandidate->phi(),0);
    //std::cout<<"On a muon"<<std::endl;
    //std::cout<<thisCandidate->eta()<<std::endl;
    Candidate_Eta->Fill(thisCandidate->eta());
  }

  if (OurCandidates->size() == 2){
    TLorentzVector CandidateVector1,CandidateVector2;
    CandidateVector1.SetPtEtaPhiM((*OurCandidates)[0].pt(),(*OurCandidates)[0].eta(),(*OurCandidates)[0].phi(),0);
    CandidateVector2.SetPtEtaPhiM((*OurCandidates)[1].pt(),(*OurCandidates)[1].eta(),(*OurCandidates)[1].phi(),0);
    Double_t Mass = (CandidateVector1+CandidateVector2).M();
    Mass_h->Fill(Mass);
  }
  
  
}

void TestAnalyzer_ME0::endJob() 
{
  histoFile->cd();
  TCanvas *c1 = new TCanvas("c1", "canvas" );

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
     
  Candidate_Eta->Write();   Candidate_Eta->Draw();  c1->Print("OutputAnalyzerPngs/Candidate_Eta.png");
  Track_Eta->Write();   Track_Eta->Draw();  c1->Print("OutputAnalyzerPngs/Track_Eta.png");
  Track_Pt->Write();   Track_Pt->Draw();  c1->Print("OutputAnalyzerPngs/Track_Pt.png");

  Segment_Eta->Write();   Segment_Eta->Draw();  GenMuon_Eta->SetLineColor(2);GenMuon_Eta->Draw("SAME"); c1->Print("OutputAnalyzerPngs/Segment_Eta.png");
  Segment_Phi->Write();   Segment_Phi->Draw();  c1->Print("OutputAnalyzerPngs/Segment_Phi.png");
  Segment_R->Write();   Segment_R->Draw();  c1->Print("OutputAnalyzerPngs/Segment_R.png");
  Segment_Pos->Write();   Segment_Pos->Draw();  c1->Print("OutputAnalyzerPngs/Segment_Pos.png");

  Rechit_Eta->Write();   Rechit_Eta->Draw();   c1->Print("OutputAnalyzerPngs/Rechit_Eta.png");
  Rechit_Phi->Write();   Rechit_Phi->Draw();  c1->Print("OutputAnalyzerPngs/Rechit_Phi.png");
  Rechit_R->Write();   Rechit_R->Draw();  c1->Print("OutputAnalyzerPngs/Rechit_R.png");
  Rechit_Pos->Write();   Rechit_Pos->Draw();  c1->Print("OutputAnalyzerPngs/Rechit_Pos.png");

  ME0Muon_Eta->Write();   ME0Muon_Eta->Draw();  c1->Print("OutputAnalyzerPngs/ME0Muon_Eta.png");
  //c1->SetLogy();
  ME0Muon_Pt->Write();   ME0Muon_Pt->Draw();  c1->Print("OutputAnalyzerPngs/ME0Muon_Pt.png");

  GenMuon_Eta->Write();   GenMuon_Eta->Draw();  c1->Print("OutputAnalyzerPngs/GenMuon_Eta.png");

  GenMuon_Pt->Write();   GenMuon_Pt->Draw();  c1->Print("OutputAnalyzerPngs/GenMuon_Pt.png");

  MatchedME0Muon_Eta->Write();   MatchedME0Muon_Eta->Draw();  c1->Print("OutputAnalyzerPngs/MatchedME0Muon_Eta.png");
  MatchedME0Muon_Pt->Write();   MatchedME0Muon_Pt->Draw();  c1->Print("OutputAnalyzerPngs/MatchedME0Muon_Pt.png");

  UnmatchedME0Muon_Eta->Write();   UnmatchedME0Muon_Eta->Draw();  c1->Print("OutputAnalyzerPngs/UnmatchedME0Muon_Eta.png");
  UnmatchedME0Muon_Pt->Write();   UnmatchedME0Muon_Pt->Draw();  c1->Print("OutputAnalyzerPngs/UnmatchedME0Muon_Pt.png");

  Mass_h->Write();   Mass_h->Draw();  c1->Print("OutputAnalyzerPngs/Mass_h.png");
  TracksPerSegment_s->SetMarkerStyle(1);
  TracksPerSegment_s->SetMarkerSize(3.0);
  TracksPerSegment_s->Write();     TracksPerSegment_s->Draw();  c1->Print("OutputAnalyzerPngs/TracksPerSegment_s.png");

  TracksPerSegment_h->Write();     TracksPerSegment_h->Draw();  c1->Print("OutputAnalyzerPngs/TracksPerSegment_h.png");

  TracksPerSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  TracksPerSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  TracksPerSegment_p->Write();     TracksPerSegment_p->Draw();  c1->Print("OutputAnalyzerPngs/TracksPerSegment_p.png");

  ClosestDelR_s->SetMarkerStyle(1);
  ClosestDelR_s->SetMarkerSize(3.0);
  ClosestDelR_s->Write();     ClosestDelR_s->Draw();  c1->Print("OutputAnalyzerPngs/ClosestDelR_s.png");

  ClosestDelR_p->GetXaxis()->SetTitle("Gen Muon #eta");
  ClosestDelR_p->GetYaxis()->SetTitle("Average closest #Delta R track");
  std::cout<<"  ClosestDelR_p values:"<<std::endl;
  for (int i=1; i<=ClosestDelR_p->GetNbinsX(); ++i){
    std::cout<<2.4+(double)i*((4.0-2.4)/40.)<<","<<ClosestDelR_p->GetBinContent(i)<<std::endl;
  }
  ClosestDelR_p->Write();     ClosestDelR_p->Draw();  c1->Print("OutputAnalyzerPngs/ClosestDelR_p.png");

  FakeTracksPerSegment_s->SetMarkerStyle(1);
  FakeTracksPerSegment_s->SetMarkerSize(3.0);
  FakeTracksPerSegment_s->Write();     FakeTracksPerSegment_s->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerSegment_s.png");

  FakeTracksPerSegment_h->Write();     FakeTracksPerSegment_h->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerSegment_h.png");

  FakeTracksPerSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  FakeTracksPerSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  FakeTracksPerSegment_p->Write();     FakeTracksPerSegment_p->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerSegment_p.png");

  FakeTracksPerAssociatedSegment_s->SetMarkerStyle(1);
  FakeTracksPerAssociatedSegment_s->SetMarkerSize(3.0);
  FakeTracksPerAssociatedSegment_s->Write();     FakeTracksPerAssociatedSegment_s->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerAssociatedSegment_s.png");

  FakeTracksPerAssociatedSegment_h->Write();     FakeTracksPerAssociatedSegment_h->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerAssociatedSegment_h.png");

  FakeTracksPerAssociatedSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  FakeTracksPerAssociatedSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  FakeTracksPerAssociatedSegment_p->Write();     FakeTracksPerAssociatedSegment_p->Draw();  c1->Print("OutputAnalyzerPngs/FakeTracksPerAssociatedSegment_p.png");

  GenMuon_Eta->Sumw2();  MatchedME0Muon_Eta->Sumw2();
  GenMuon_Pt->Sumw2();  MatchedME0Muon_Pt->Sumw2();

  Track_Eta->Sumw2();  ME0Muon_Eta->Sumw2();
  Track_Pt->Sumw2();  ME0Muon_Pt->Sumw2();

  UnmatchedME0Muon_Eta->Sumw2();
  UnmatchedME0Muon_Pt->Sumw2();
  
  MuonRecoEff_Eta->Divide(MatchedME0Muon_Eta, GenMuon_Eta, 1, 1, "B");
  std::cout<<"GenMuon_Eta =  "<<GenMuon_Eta->Integral()<<std::endl;
  std::cout<<"MatchedME0Muon_Eta =  "<<MatchedME0Muon_Eta->Integral()<<std::endl;
  MuonRecoEff_Eta->GetXaxis()->SetTitle("Gen Muon #eta");
  MuonRecoEff_Eta->GetYaxis()->SetTitle("Matching Efficiency");
  MuonRecoEff_Eta->Write();   MuonRecoEff_Eta->Draw();  c1->Print("OutputAnalyzerPngs/MuonRecoEff_Eta.png");

  std::cout<<"  MuonRecoEff_Eta values:"<<std::endl;
  for (int i=1; i<=MuonRecoEff_Eta->GetNbinsX(); ++i){
    std::cout<<2.4+(double)i*((4.0-2.4)/40.)<<","<<MuonRecoEff_Eta->GetBinContent(i)<<std::endl;
  }
  

  // MuonRecoEff_Pt->Divide(MatchedME0Muon_Pt, GenMuon_Pt, 1, 1, "B");
  // MuonRecoEff_Pt->GetXaxis()->SetTitle("Gen Muon p_{T}");
  // MuonRecoEff_Pt->GetYaxis()->SetTitle("Matching Efficiency");
  // MuonRecoEff_Pt->SetMinimum(.85);
  // MuonRecoEff_Pt->Write();   MuonRecoEff_Pt->Draw();  c1->Print("OutputAnalyzerPngs/MuonRecoEff_Pt.png");

  MuonAllTracksEff_Eta->Divide(ME0Muon_Eta, Track_Eta, 1, 1, "B");
  MuonAllTracksEff_Eta->Write();   MuonAllTracksEff_Eta->Draw();  c1->Print("OutputAnalyzerPngs/MuonAllTracksEff_Eta.png");

  // MuonAllTracksEff_Pt->Divide(ME0Muon_Pt, Track_Pt, 1, 1, "B");
  // MuonAllTracksEff_Pt->Write();   MuonAllTracksEff_Pt->Draw();  c1->Print("OutputAnalyzerPngs/MuonAllTracksEff_Pt.png");

  // MuonUnmatchedTracksEff_Eta->Divide(UnmatchedME0Muon_Eta, ME0Muon_Eta, 1, 1, "B");
  // MuonUnmatchedTracksEff_Eta->Write();   Candidate_Eta->Draw();  c1->Print("OutputAnalyzerPngs/Candidate_Eta.png");

  // MuonUnmatchedTracksEff_Pt->Divide(UnmatchedME0Muon_Pt, ME0Muon_Pt, 1, 1, "B");
  // MuonUnmatchedTracksEff_Pt->Write();   Candidate_Eta->Draw();  c1->Print("OutputAnalyzerPngs/Candidate_Eta.png");
  FractionMatched_Eta->Divide(MatchedME0Muon_Eta, ME0Muon_Eta, 1, 1, "B");
  FractionMatched_Eta->GetXaxis()->SetTitle("Gen Muon #eta");
  FractionMatched_Eta->GetYaxis()->SetTitle("Matched/All ME0Muons");
  FractionMatched_Eta->Write();   FractionMatched_Eta->Draw();  c1->Print("OutputAnalyzerPngs/FractionMatched_Eta.png");

  gStyle->SetOptStat(1);
  PtDiff_h->GetXaxis()->SetTitle("(pt track-ptgen)/ptgen");
  PtDiff_h->Write();     PtDiff_h->Draw();  c1->Print("OutputAnalyzerPngs/PtDiff_h.png");

  QOverPtDiff_h->GetXaxis()->SetTitle("(q/pt track-q/ptgen)/(q/ptgen)");
  QOverPtDiff_h->Write();     QOverPtDiff_h->Draw();  c1->Print("OutputAnalyzerPngs/QOverPtDiff_h.png");

  gStyle->SetOptStat(0);
  PtDiff_s->SetMarkerStyle(1);
  PtDiff_s->SetMarkerSize(3.0);
  PtDiff_s->GetXaxis()->SetTitle("Gen Muon #eta");
  PtDiff_s->GetYaxis()->SetTitle("|(pt track-ptgen)/ptgen|");
  PtDiff_s->Write();     PtDiff_s->Draw();  c1->Print("OutputAnalyzerPngs/PtDiff_s.png");
  
  PtDiff_p->SetMarkerStyle(1);
  PtDiff_p->SetMarkerSize(3.0);
  PtDiff_p->GetXaxis()->SetTitle("Gen Muon #eta");
  PtDiff_p->GetYaxis()->SetTitle("Average (pt track-ptgen)/ptgen");
  PtDiff_p->Write();     PtDiff_p->Draw();  c1->Print("OutputAnalyzerPngs/PtDiff_p.png");

  //PtDiff_rms->SetMarkerStyle(1);
  //PtDiff_rms->SetMarkerSize(3.0);

  PtDiff_rms->SetMarkerStyle(22); 
  PtDiff_rms->SetMarkerSize(1.2); 
  PtDiff_rms->SetMarkerColor(kBlue); 
  //PtDiff_rms->SetLineColor(kRed); 
  
  //PtDiff_rms->Draw("PL"); 

  PtDiff_rms->GetXaxis()->SetTitle("Gen Muon #eta");
  PtDiff_rms->GetYaxis()->SetTitle("RMS of (pt track-ptgen)/ptgen");
  PtDiff_rms->Write();     PtDiff_rms->Draw("P");  c1->Print("OutputAnalyzerPngs/PtDiff_rms.png");

  PDiff_h->GetXaxis()->SetTitle("|(p track-pgen)/pgen|");
  PDiff_h->Write();     PDiff_h->Draw();  c1->Print("OutputAnalyzerPngs/PDiff_h.png");

  PDiff_s->SetMarkerStyle(1);
  PDiff_s->SetMarkerSize(3.0);
  PDiff_s->GetXaxis()->SetTitle("Gen Muon #eta");
  PDiff_s->GetYaxis()->SetTitle("|(p track-pgen)/pgen|");
  PDiff_s->Write();     PDiff_s->Draw();  c1->Print("OutputAnalyzerPngs/PDiff_s.png");
  
  PDiff_p->SetMarkerStyle(1);
  PDiff_p->SetMarkerSize(3.0);
  PDiff_p->GetXaxis()->SetTitle("Gen Muon #eta");
  PDiff_p->GetYaxis()->SetTitle("Average |(p track-pgen)/pgen|");
  PDiff_p->Write();     PDiff_p->Draw();  c1->Print("OutputAnalyzerPngs/PDiff_p.png");

  VertexDiff_h->Write();     VertexDiff_h->Draw();  c1->Print("OutputAnalyzerPngs/VertexDiff_h.png");

  delete histoFile; histoFile = 0;
}

DEFINE_FWK_MODULE(TestAnalyzer_ME0);
