import sys 

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
        #filename
        #'file:/eos/cms/store/relval/CMSSW_12_4_0_pre4/RelValDisplacedSUSY_14TeV/MINIAODSIM/124X_mcRun3_2021_realistic_v1-v1/2580000/4454d87c-cf0e-41ce-bf05-746bbd02631f.root'
        #'file:/eos/cms/store/relval/CMSSW_12_4_0_pre4/RelValDisplacedSUSY_14TeV/MINIAODSIM/PU_124X_mcRun3_2021_realistic_v1-v1/2580000/55b84045-4ec8-45b6-899e-31f151e24f62.root'
        'file:/eos/cms/store/relval/CMSSW_12_4_0_pre4/RelValDisplacedSUSY_14TeV/MINIAODSIM/PU_124X_mcRun3_2021_realistic_v1-v1/2580000/e4802125-09fa-42a6-b79f-c3379e768a9d.root'

    )
)

process.ntuple = cms.EDAnalyzer('Analyzer',

   
           #globaltracks =cms.InputTag("globalMuons"),

           displacedglobaltracks =cms.InputTag("displacedGlobalMuons"),

           #standalonetracks =cms.InputTag("standAloneMuons"),

           displacedstandalonetracks =cms.InputTag("displacedStandAloneMuons"),

           displacedtracks =cms.InputTag("displacedTracks"),

           #displacedoutintracks = cms.InputTag("muonSeededTracksOutInDisplaced"),

           #earlydisplacedmuons = cms.InputTag("earlyDisplacedMuons"),

           #genparticle =cms.InputTag("genParticles"),

           #recovertex =cms.InputTag("offlinePrimaryVertices"),

           #muons = cms.InputTag("muons"),

           #displacedmuons = cms.InputTag("displacedMuons"),

           #tracks = cms.InputTag("tracks")

           genparticle =cms.InputTag("prunedGenParticles"),

           TruthMatchMuonMaxR = cms.untracked.double(0.3), # [eta-phi]


         )

outputfilename = 'file:truthmatch_added_with_pu_2.root'

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(outputfilename)
                                   )
process.p = cms.Path(process.ntuple)

