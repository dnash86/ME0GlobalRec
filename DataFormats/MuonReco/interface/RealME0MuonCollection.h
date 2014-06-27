#ifndef DataFormats_RealME0MuonCollection_H
#define DataFormats_RealME0MuonCollection_H

/** \class RealME0MuonCollection
 *
 * The collection of RealME0Muon's. See \ref RealME0MuonCollection.h for details.
 *
 *  $Date: 2010/03/12 13:08:15 $
 *  \author David Nash
 */

#include "DataFormats/MuonReco/interface/RealME0Muon.h"
#include "DataFormats/Common/interface/Ref.h"

/// collection of RealME0Muons
typedef std::vector<reco::RealME0Muon> RealME0MuonCollection;

/// persistent reference to a RealME0Muon
typedef edm::Ref<RealME0MuonCollection> RealME0MuonRef;

//typedef std::vector<ME0Muon> ME0MuonCollection; 
	
#endif
