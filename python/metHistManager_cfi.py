import FWCore.ParameterSet.Config as cms

metHistManager = cms.PSet(
  pluginName = cms.string('metHistManager'),
  pluginType = cms.string('MEtHistManager'),
      
  metSource = cms.InputTag('layer1METs'),
  metSignificanceSource = cms.InputTag('met'),

  dqmDirectory_store = cms.string('MEtQuantities')
)
