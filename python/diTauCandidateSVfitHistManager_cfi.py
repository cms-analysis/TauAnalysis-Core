import FWCore.ParameterSet.Config as cms

diTauCandidateSVfitHistManager = cms.PSet(    
    pluginName = cms.string('diTauCandidateSVfitHistManager'),
    pluginType = cms.string('DiCandidatePairSVfitHistManager'),
      
    diTauCandidateSource = cms.InputTag(''),

    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    genParticleSource = cms.InputTag('genParticles'),

    dqmDirectory_store = cms.string('DiTauCandidateSVfitQuantities'),

    #requireGenMatch = cms.bool(True),
    requireGenMatch = cms.bool(False),

    #normalization = cms.string("diTauCandidates"),
    normalization = cms.string("events")
)
