import FWCore.ParameterSet.Config as cms

diTauCandidateSVfitHistManager = cms.PSet(
    pluginName = cms.string('diTauCandidateSVfitHistManager'),
    pluginType = cms.string('DiCandidatePairSVfitHistManager'),

    diTauCandidateSource = cms.InputTag(''),

    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    genParticleSource = cms.InputTag('genParticles'),

    dqmDirectory_store = cms.string('DiTauCandidateSVfitQuantities'),

    SVfitAlgorithms = cms.VPSet(
        cms.PSet(
            name = cms.string("psKine")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt_ptBalance")
        ),
        #cms.PSet(
            #name = cms.string("polKine"),
            #polarizationHypotheses = cms.vstring("LL", "LR", "RL", "RR", "best{LL,RR}", "best{LR,RL}", "best{LL,LR,RL,RR}"),
            #massHypotheses = cms.vdouble(90.)
        #),
        #cms.PSet(
            #name = cms.string("polKine_MEt"),
            #polarizationHypotheses = cms.vstring("LL", "LR", "RL", "RR", "best{LL,RR}", "best{LR,RL}", "best{LL,LR,RL,RR}"),
            #massHypotheses = cms.vdouble(90.)
        #),
        #cms.PSet(
            #name = cms.string("polKine_MEt_ptBalance"),
            #polarizationHypotheses = cms.vstring("LL", "LR", "RL", "RR", "best{LL,RR}", "best{LR,RL}", "best{LL,LR,RL,RR}"),
            #massHypotheses = cms.vdouble(90.)
        #),
        ##cms.PSet(
        ##    name = cms.string("polKine_MEt_ptBalance_Zprod"),
        ##    polarizationHypotheses = cms.vstring("LL", "LR", "RL", "RR", "best{LL,RR}", "best{LR,RL}", "best{LL,LR,RL,RR}"),
        ##    massHypotheses = cms.vdouble(90.)
        ##),
    ),

    vertexPtThresholds = cms.vdouble(5., 10., 15., 20.),

    #requireGenMatch = cms.bool(True),
    requireGenMatch = cms.bool(False),

    #normalization = cms.string("diTauCandidates"),
    normalization = cms.string("events")
)
