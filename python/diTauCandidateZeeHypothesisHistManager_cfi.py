import FWCore.ParameterSet.Config as cms

diTauCandidateZeeHypothesisHistManager = cms.PSet(    
    pluginName = cms.string('diTauCandidateZeeHypothesisHistManager'),
    pluginType = cms.string('PATElecTauPairZeeHypothesisHistManager'),
      
    ZeeHypothesisSource = cms.InputTag('elecTauPairZeeHypotheses'),

    dqmDirectory_store = cms.string('DiTauCandidateZeeHypothesisQuantities'),

    #normalization = cms.string("diTauCandidates")
    normalization = cms.string("events")
)
