import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodAnalyzer = cms.PSet(    
    pluginName = cms.string('svFitLikelihoodAnalyzer'),
    pluginType = cms.string('SVfitLikelihoodMuTauPairAnalyzer'),
      
    diTauCandidateSource = cms.InputTag('selectedMuTauPairsForAHtoMuTauCollinearApproxZmassVetoCumulative'),

    svFitAlgorithms = cms.VPSet(
        cms.PSet(
            name = cms.string("psKine")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt_ptBalance")
        )
    ),

    svFitLikelihoodFunctions = cms.VPSet(
        svFitLikelihoodDiTauKinematicsPhaseSpace,
        svFitLikelihoodMEt,
        svFitLikelihoodDiTauPtBalance
    )
)
