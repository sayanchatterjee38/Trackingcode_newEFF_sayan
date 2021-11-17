import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.EventContent.EventContent_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')

# Add PbPb collision event selection                                                                                                                                                                       
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('TrackingCode.HIRun2015Ana.HITrackCorrectionAnalyzer_cfi')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('alltrk_crosscheck_Nov17_2021.root')
)

process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpRecoAssocGeneralTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("generalTracks")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

#process.load("SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi")

# Input source
process.source = cms.Source("PoolSource",
                        duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames =  cms.untracked.vstring(
                                
                                #'/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_668.root'
                                '/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_1.root'
                                #'/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_10.root'
                                #'/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_101.root'
                                #'/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_102.root'
                            )
)

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSelect = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSelect.HLTPaths = [
 "HLT_HIMinimumBias_*",
 ]
process.hltSelect.andOr = cms.bool(True)
process.hltSelect.throw = cms.bool(False)



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic_hi', '') 
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

# __________________ Event selection _________________                                                                                                                                                     
# Load centrality producer for centrality calculation                                                                                                                                                      
# Add PbPb centrality                                                                                                                                                                                      
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.newCentralityBin = process.centralityBin.clone()
process.centralityBin.nonDefaultGlauberModel = cms.string("")


process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
            tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x02_mc"),
            #tag = cms.string("CentralityTable_HFtowers200_HydjetCymbal5Ev8_v1020x01_mc"),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
      label = cms.untracked.string("HFtowers")
   ),
])

process.p = cms.Path(process.hiCentrality * process.centralityBin)


### Track cuts ###
# input collections
process.HITrackCorrections.generator = cms.InputTag("generator")
process.HITrackCorrections.centralitySrc = cms.InputTag("centralityBin","HFtowers")
process.HITrackCorrections.trackSrc = cms.InputTag("generalTracks")
process.HITrackCorrections.qualityString = cms.string("highPurity")
process.HITrackCorrections.vertexSrc = cms.InputTag("offlinePrimaryVertices")
process.HITrackCorrections.pfCandSrc = cms.InputTag("particleFlow")
process.HITrackCorrections.jetSrc = cms.InputTag("ak4CaloJets")
# options
process.HITrackCorrections.useCentrality = True
process.HITrackCorrections.applyTrackCuts = True
process.HITrackCorrections.fillNTuples = False
process.HITrackCorrections.applyVertexZCut = True
process.HITrackCorrections.doVtxReweighting = False
process.HITrackCorrections.doCaloMatched = True
# cut values
process.HITrackCorrections.dxyErrMax = 3.0
process.HITrackCorrections.dzErrMax = 3.0
process.HITrackCorrections.ptErrMax = 0.1
process.HITrackCorrections.nhitsMin = 11
process.HITrackCorrections.chi2nMax = 0.18
process.HITrackCorrections.reso = 0.5
#process.HITrackCorrections.crossSection = 1.0 #1.0 is no reweigh

#algo 
#process.HITrackCorrections.algoParameters = cms.vint32(4,5,6,7)
# vertex reweight parameters
process.HITrackCorrections.vtxWeightParameters = cms.vdouble(0.0306789, 0.427748, 5.16555, 0.0228019, -0.02049, 7.01258 )
###


process.ana_step = cms.Path(
process.offlinePrimaryVerticesRecovery )



process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.newCentralityBin = process.centralityBin.clone()
###
process.p = cms.Path(
    process.tpClusterProducer *
    process.quickTrackAssociatorByHits *
    process.tpRecoAssocGeneralTracks *
    process.centralityBin *
    process.offlinePrimaryVerticesRecovery *
    process.hltSelect *
    process.hfCoincFilter2Th4 *           
    process.primaryVertexFilter *         
    process.clusterCompatibilityFilter *  
    process.HITrackCorrections
    )


# peripheral pv recovery                                                                                                                                                                                   
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
