import FWCore.ParameterSet.Config as cms


HITrackCorrections = cms.EDAnalyzer('HITrackCorrectionAnalyzer',
    trackSrc = cms.InputTag("generalTracks"),
    mvaSrc = cms.InputTag('generalTracks','MVAValues'),
    vertexSrc = cms.InputTag("offlinePrimaryVertices"),
    #pfCandSrc = cms.InputTag("particleFlowTmp"),
    pfCandSrc = cms.InputTag("particleFlow"),
    #jetSrc = cms.InputTag("akPu4CaloJets"),
    jetSrc = cms.InputTag("ak4CaloJets"),
    tpEffSrc = cms.InputTag('mix','MergedTrackTruth'),
    tpFakSrc = cms.InputTag('mix','MergedTrackTruth'),
    associatorMap = cms.InputTag('tpRecoAssocGeneralTracks'),
    tracksgen = cms.InputTag("genParticles"),
    
    ptBins = cms.vdouble(
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.05, 1.1, 1.15, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0,
        20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 
        180.0, 300.0, 500.0
    ),
    etaBins = cms.vdouble( 
#       -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
 #       0.4, 0.8, 1.2, 1.6, 2.0, 2.4
       -2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0,
        0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4
    ),
    occBins = cms.vdouble(
#	0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,  70.0, 80.0, 90.0, 100.0, 
#	110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0
        0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 
        140.0, 200.0
    ),

    phiBins = cms.vdouble(
#        -3.3, -3.0, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3
#         -3.3, -3.0, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3 
         -3.15, -3.05, -2.95, -2.85, -2.75, -2.65, -2.55, -2.45, -2.35, -2.25, -2.15, -2.05, -1.95, -1.85, -1.75, -1.65, -1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15

   ),


    vtxWeightParameters = cms.vdouble( 4.49636e-02, 1.36629e-01, 5.30010e+00,
                                       2.50170e-02, 4.59123e-01, 9.64888e+00 ),
    algoParameters = cms.vint32(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    doCaloMatched = cms.bool(True),
    reso = cms.double(0.5),
    crossSection = cms.double(1.0),
    doVtxReweighting = cms.bool(False),
    applyVertexZCut = cms.bool(True),
    vertexZMax = cms.double(15.),
    applyTrackCuts = cms.bool(True),
    qualityString = cms.string("highPurity"),
    dxyErrMax = cms.double(3.0),
    dzErrMax = cms.double(3.0),
    ptErrMax = cms.double(0.1),
    nhitsMin = cms.int32(11),
    chi2nMax = cms.double(0.18), #sayan ~ 26.09.2021
    doMomRes = cms.bool(False),
    fillNTuples = cms.bool(False),
    useCentrality = cms.bool(True), #sayan ~ 26.09.2021
    centralitySrc = cms.InputTag("centralityBin","HFTowers")
)
