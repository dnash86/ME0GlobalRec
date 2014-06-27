import FWCore.ParameterSet.Config as cms

process = cms.Process("ME0SegmentMatching")

#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#process.load("Configuration.StandardSequences.MagneticField_38T_cff")

#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")

#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')


#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
#process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')



#process.load('Geometry.GEMGeometry.cmsExtendedGeometryPostLS1plusGEMr08v01XML_cfi')


## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCal_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon4EtaReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon4Eta_cff')
#process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

## global tag for 2019 upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" }

)

#process.me0SegmentProducer = cms.EDProducer("EmulatedME0SegmentProducer")

process.me0SegmentMatcher = cms.EDProducer("ME0SegmentMatcher"
)

process.me0MuonConverter = cms.EDProducer("ME0MuonConverter")

process.p = cms.Path(process.me0SegmentMatcher+
                     process.me0MuonConverter
                     )

process.PoolSource.fileNames = [
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/step3_SingleMuPt20Extended.root'

]


process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
#                              process.AODSIMEventContent,
                              fileName = cms.untracked.string(
    '/tmp/dnash/test20.root'

                                                              )
                              )

process.outpath = cms.EndPath(process.o1)
