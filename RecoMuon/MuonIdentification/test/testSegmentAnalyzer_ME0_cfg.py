import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")


#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.Geometry.GeometryExtended2023Muon4EtaReco_cff')

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" }

)

process.Test = cms.EDAnalyzer("SegmentAnalyzer_ME0",
                              HistoFile = cms.string('OutputTestHistos_SegmentsOnly.root'),
                              
)

process.p = cms.Path(process.Test)
process.PoolSource.fileNames = [
    #'file:FirstTest.root'
    #'file:/afs/cern.ch/work/d/dnash/ME0Segments/CommitToCMSSW/CMSSW_6_1_2_SLHC8/src/RecoMuon/MuonIdentification/test/FirstTest.root'
    #'file:/tmp/dnash/Zmumu_FlatMuonPt_SLHC8.root'
    #'file:/afs/cern.ch/work/d/dnash/ME0Segments/FullSimPixel/CMSSW_6_2_0_SLHC8/src/'
    #'file:/afs/cern.ch/work/d/dnash/ME0Segments/FullSimPixel/CMSSW_6_2_0_SLHC8/src/ZMMTest.root'
    #'file:/tmp/dnash/ZMMTest_Again.root'
    #'file:/tmp/dnash/step3.root'
    #'file:/tmp/dnash/FourMu_ME0Muons.root'
    #'file:/tmp/dnash/step3.root',
    #'file:/tmp/dnash/step3_copied.root',
    #'file:/tmp/dnash/FourMu_5kevts_output.root'
    #'file:/tmp/dnash/step3_TryWithoutMET.root'
    #'file:/tmp/dnash/step3_trackingfix.root'
    #'file:/tmp/dnash/RunOnstep3_TryWithoutMET_WithRealSegments.root'
    #'file:/tmp/dnash/TestingRealSegments.root'
    #'file:/tmp/dnash/test.root'
    #'file:/tmp/dnash/RealME0Muons_3Sigma_SingleMuPt10Extended.root'
    #'file:/tmp/dnash/RealME0Muons_3Sigma_testjustme0Seg.root'
    #'root://eoscms//eos/cms//store/group/upgrade/muon/ME0GlobalReco/RealME0Muons_3Sigma_SingleMuPt5Extended.root'
    #'root://eoscms//eos/cms//store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Matched_Muplus_Pt5-gun_9.root'
    #'file:/tmp/dnash/RunOnstep3_TryWithoutMET_WithRealSegments.root'
    #'file:/tmp/dnash/step3_Post130.root',
    #'file:/tmp/dnash/step3_Post262.root',
    #'file:/tmp/dnash/step3_Post600.root'
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_6.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_11.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_12.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_13.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_14.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_15.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_16.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_17.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_18.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_19.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muminus_Pt5-gun_20.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_1.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_10.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_2.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_3.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_4.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_5.root',

    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_7.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_8.root',
    'root://eoscms//eos/cms/store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/Muplus_Pt5-gun_9.root'

    
]
