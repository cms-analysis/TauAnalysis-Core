import FWCore.ParameterSet.Config as cms

metHistManager = cms.PSet(
  pluginName = cms.string('metHistManager'),
  pluginType = cms.string('MEtHistManager'),
      
  metSource = cms.InputTag('layer1METs'),
  #metSignificanceSource = cms.InputTag('met'),
  metSignificanceSource = cms.InputTag('metsignificance'),

  dqmDirectory_store = cms.string('MEtQuantities')
)
