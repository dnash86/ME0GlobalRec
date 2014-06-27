/** \file ME0SegmentMatcher.cc
 *
 * \author David Nash
 */

#include <RecoMuon/MuonIdentification/src/ME0SegmentMatcher.h>

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <DataFormats/MuonReco/interface/RealME0Muon.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "DataFormats/GeometrySurface/interface/LocalError.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"


//For Plots ---------------Remove later....

#include <TLegend.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
//---------------------------------------

ME0SegmentMatcher::ME0SegmentMatcher(const edm::ParameterSet& pas) : iev(0){
	
  produces<std::vector<reco::RealME0Muon> >();  //May have to later change this to something that makes more sense, OwnVector, RefVector, etc


}

ME0SegmentMatcher::~ME0SegmentMatcher() {}

void ME0SegmentMatcher::produce(edm::Event& ev, const edm::EventSetup& setup) {

    LogDebug("ME0SegmentMatcher") << "start producing segments for " << ++iev << "th event ";

   

    //Getting the objects we'll need

    
    using namespace edm;
    ESHandle<MagneticField> bField;
    setup.get<IdealMagneticFieldRecord>().get(bField);
    ESHandle<Propagator> shProp;
    setup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);


    
    using namespace reco;

    // Handle<std::vector<ME0Segment> > OurSegments;
    // ev.getByLabel<std::vector<ME0Segment> >("me0SegmentProducer", OurSegments);

    Handle<ME0SegmentCollection> OurSegments;
    ev.getByLabel("me0Segments","",OurSegments);


    // Handle<CSCSegmentCollection> TestSegments;
    // ev.getByLabel("cscSegments","",TestSegments);


    std::auto_ptr<std::vector<RealME0Muon> > oc( new std::vector<RealME0Muon> ); 
    std::vector<RealME0Muon> TempStore; 

    Handle <TrackCollection > generalTracks;
    ev.getByLabel <TrackCollection> ("generalTracks", generalTracks);


    
    int TrackNumber = 0;
    std::vector<int> TkMuonNumbers, TkIndex, TkToKeep;
    std::vector<GlobalVector>FinalTrackPosition;


    for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
	 thisTrack != generalTracks->end(); ++thisTrack,++TrackNumber){
      //Initializing our plane
     
      //FIXME: Remove later
      if (fabs(thisTrack->eta()) < 1.8) continue;
      //if (fabs(thisTrack->pt()) < 0.6) continue;


      float zSign  = thisTrack->pz()/fabs(thisTrack->pz());

      //float zValue = 560. * zSign;
      float zValue = 526.75 * zSign;
      Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
      //Getting the initial variables for propagation
      int chargeReco = thisTrack->charge(); 
      GlobalVector p3reco, r3reco;

      p3reco = GlobalVector(thisTrack->outerPx(), thisTrack->outerPy(), thisTrack->outerPz());
      r3reco = GlobalVector(thisTrack->outerX(), thisTrack->outerY(), thisTrack->outerZ());

      AlgebraicSymMatrix66 covReco;
      //This is to fill the cov matrix correctly
      AlgebraicSymMatrix55 covReco_curv;
      covReco_curv = thisTrack->outerStateCovariance();
      FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
      getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);

      //Now we propagate and get the propagated variables from the propagated state
      SteppingHelixStateInfo startrecostate(initrecostate);
      SteppingHelixStateInfo lastrecostate;

      const SteppingHelixPropagator* ThisshProp = 
	dynamic_cast<const SteppingHelixPropagator*>(&*shProp);
	
      lastrecostate = ThisshProp->propagate(startrecostate, *plane);
	
      FreeTrajectoryState finalrecostate;
      lastrecostate.getFreeState(finalrecostate);

      AlgebraicSymMatrix66 covFinalReco;
      GlobalVector p3FinalReco_glob, r3FinalReco_globv;
      getFromFTS(finalrecostate, p3FinalReco_glob, r3FinalReco_globv, chargeReco, covFinalReco);
      FinalTrackPosition.push_back(r3FinalReco_globv);

      
      //To transform the global propagated track to local coordinates
      int SegmentNumber = 0;
     

      for (auto thisSegment = OurSegments->begin(); thisSegment != OurSegments->end(); 
	   ++thisSegment,++SegmentNumber){


	ME0DetId id = thisSegment->me0DetId();
	auto roll = me0Geom->etaPartition(id); 


	if ( zSign * roll->toGlobal(thisSegment->localPosition()).z() < 0 ) continue;

	GlobalPoint r3FinalReco_glob(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z());
	LocalPoint r3FinalReco = roll->toLocal(r3FinalReco_glob);
	LocalVector p3FinalReco=roll->toLocal(p3FinalReco_glob);

	//GlobalPoint thisPosition(roll->toGlobal(thisSegment->localPosition()));
	LocalPoint thisPosition(thisSegment->localPosition());
	LocalVector thisDirection(thisSegment->localDirection().x(),thisSegment->localDirection().y(),thisSegment->localDirection().z());  //FIXME
	//The same goes for the error
	AlgebraicMatrix thisCov(4,4,0);   

	for (int i = 1; i <=4; i++){
	  for (int j = 1; j <=4; j++){
	    thisCov(i,j) = thisSegment->parametersError()(i,j);
	  }
	}



	//Only taking the diagonal terms

	LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,chargeReco);
	JacobianCartesianToLocal jctl(roll->surface(),ltp);
	AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian(); 

	AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc); 
	AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
	for(int i=0; i<5; ++i) {
	  for(int j=0; j<5; ++j) {
	    C[i][j] = Ctmp[i][j]; 

	  }
	}  



	Double_t sigmax = sqrt(C[3][3]+thisSegment->localPositionError().xx() );      
	Double_t sigmay = sqrt(C[4][4]+thisSegment->localPositionError().yy() );

	bool X_MatchFound = false, Y_MatchFound = false;
	

	if ( (fabs(thisPosition.x()-r3FinalReco.x()) < (3.0 * sigmax))  ) X_MatchFound = true;
	if ( (fabs(thisPosition.y()-r3FinalReco.y()) < (3.0 * sigmay))  ) Y_MatchFound = true;

	
	if (X_MatchFound && Y_MatchFound) {
	  TrackRef thisTrackRef(generalTracks,TrackNumber);
	 
	  TempStore.push_back(reco::RealME0Muon(thisTrackRef,(*thisSegment),SegmentNumber));
	  TkIndex.push_back(TrackNumber);
	}
      }
    }

    for (unsigned int i = 0; i < TkIndex.size(); ++i){     //Now we construct a vector of unique TrackNumbers of tracks that have been stored
      bool AlreadyStoredInTkMuonNumbers = false;
      for (unsigned int j = 0; j < TkMuonNumbers.size(); ++j){
	if (TkMuonNumbers[j]==TkIndex[i]) AlreadyStoredInTkMuonNumbers = true;
      }
      if (!AlreadyStoredInTkMuonNumbers) TkMuonNumbers.push_back(TkIndex[i]);
    }

    for (unsigned int i = 0; i < TkMuonNumbers.size(); ++i){            //Now we loop over each TrackNumber that has been stored
      int ReferenceMuonNumber = TkMuonNumbers[i];          // The muon number of the track, starts at 0 and increments
      double RefDelR = 99999.9, ComparisonIndex = 0;
      int WhichTrackToKeep=-1;
      for (std::vector<RealME0Muon>::const_iterator thisMuon = TempStore.begin();    //Now we have the second nested loop, over the RealME0Muons
	   thisMuon != TempStore.end(); ++thisMuon, ++ComparisonIndex){
	
	int thisMuonNumber = TkIndex[ComparisonIndex];    //The track number of the muon we are currently looking at
	if (thisMuonNumber == ReferenceMuonNumber){        //This means we're looking at the same track in the ME0Muon and TkMuons
	  ME0Segment Seg = thisMuon->me0segment();
	  ME0DetId id = Seg.me0DetId();
	  auto roll = me0Geom->etaPartition(id); 

	  TrackRef TkRef = thisMuon->innerTrack();
	  //Here LocalPoint is used, although the local frame and global frame coincide, hence all calculations are made in global coordinates
	  //  NOTE: Correct this when making the change to "real" ME0Segments, since these will be in real local coordinates

	  //LocalPoint SegPos(Seg.localPosition().x(),Seg.localPosition().y(),Seg.localPosition().z());
	  GlobalPoint SegPos(roll->toGlobal(Seg.localPosition()));

	  //LocalPoint TkPos(TkRef->vx(),TkRef->vy(),TkRef->vz());
	  //LocalPoint TkPos(FinalTrackPosition[thisMuonNumber].x(),FinalTrackPosition[thisMuonNumber].y(),FinalTrackPosition[thisMuonNumber].z());	  
	  GlobalPoint TkPos(FinalTrackPosition[thisMuonNumber].x(),FinalTrackPosition[thisMuonNumber].y(),FinalTrackPosition[thisMuonNumber].z());
	  double delR = reco::deltaR(SegPos,TkPos);


	  if (delR < RefDelR) {
	    WhichTrackToKeep = ComparisonIndex;  //Storing a list of the vector indices of tracks to keep
	                                         //Note: These are not the same as the "Track Numbers"
	    RefDelR=delR;
	  }
	}
      }
      if (WhichTrackToKeep != -1) TkToKeep.push_back(WhichTrackToKeep);
    }

    for (unsigned int i = 0; i < TkToKeep.size(); ++i){
      int thisKeepIndex = TkToKeep[i];

      oc->push_back(TempStore[thisKeepIndex]);    //Filling the collection


    }
  	
    // put collection in event
    ev.put(oc);
}

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix55& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix66& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void ME0SegmentMatcher::getFromFTS(const FreeTrajectoryState& fts,
				    GlobalVector& p3, GlobalVector& r3, 
				    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  //Yikes, was setting this to p3T instead of r3T!?!
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}


void ME0SegmentMatcher::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(me0Geom);

}


void ME0SegmentMatcher::endRun()
{
  
}

 DEFINE_FWK_MODULE(ME0SegmentMatcher);
