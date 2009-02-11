import FWCore.ParameterSet.Config as cms

metHistManager = cms.PSet(
  name = cms.string('metHistManager'),
  type = cms.string('MEtHistManager'),
      
  metSource = cms.InputTag('allLayer1METs'),

  dqmDirectory_store = cms.string('MEtQuantities')
)
