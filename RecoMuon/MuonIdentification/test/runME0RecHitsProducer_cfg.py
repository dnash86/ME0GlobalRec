import FWCore.ParameterSet.Config as cms

process = cms.Process("LocalReco")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:/tmp/dnash/Muplus_Pt5_gun.root'
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/0202C8CE-B2DB-E311-B903-002354EF3BD2.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/040E3CF1-9ADB-E311-9965-0025905A6092.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/04742CE5-69DB-E311-8763-002618943882.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/0698C5AF-CBDB-E311-B2CB-00248C0BE018.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/08E2FD7A-BDDB-E311-963C-00304867BED8.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/08ECDF7A-63DB-E311-9BAF-0026189438BF.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/0A38A807-98DB-E311-844F-0025905A60D6.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/0C0FE2E2-6DDB-E311-B7AE-0026189438D6.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/0EE9D039-65DB-E311-9FAC-002354EF3BDB.root',
        #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-RECO/PU140bx25_PH2_1K_FB_V2-v1/00000/101D8E73-B5DB-E311-9D7E-00261894396B.root'
        #'TEMPLATEIN'
    'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt5-gun/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V2-v1/00000/001C83EB-36DB-E311-B3A6-002618943870.root'
    #'root://xrootd.unl.edu//store/mc/Muon2023Upg14DR/Muplus_Pt50-gun/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V2-v1/00000/02479FFF-2ADA-E311-809E-0025905A611C.root'
    )
)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")

process.load('RecoLocalMuon.GEMRecHit.me0RecHits_cfi')
process.load('RecoLocalMuon.GEMRecHit.me0Segments_cfi')

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(
    '/tmp/dnash/test.root'
    #'TEMPLATEOUT'
    ),
    outputCommands = cms.untracked.vstring(
        'keep  *_*_*_*',
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('rechit_step')
    )
)


process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")
process.rechit_step  = cms.Path(process.me0RecHits*process.me0Segments)
#process.rechit_step  = cms.Path(process.me0Segments)
process.endjob_step  = cms.Path(process.endOfProcess)
process.out_step     = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.rechit_step,
    process.endjob_step,
    process.out_step
)
