import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodAnalyzer = cms.PSet(
    pluginName = cms.string('svFitLikelihoodAnalyzer'),
    pluginType = cms.string('SVfitLikelihoodDiTauAnalyzer'),

    diTauCandidateSource = cms.InputTag(''),

    svFitAlgorithms = cms.VPSet(
        cms.PSet(
            name = cms.string("psKine")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt_ptBalance")
        ),
        cms.PSet(
            name = cms.string("psKine_MEt_ptBalance2")
        ),
        #cms.PSet(
            #name = cms.string("psKine_MEt_Track_ptBalance")
        #),
        #cms.PSet(
            #name = cms.string("polKine")
        #),
        #cms.PSet(
            #name = cms.string("polKine_MEt")
        #),
        #cms.PSet(
            #name = cms.string("polKine_MEt_ptBalance")
        #)
    ),

    svFitLikelihoodFunctions = cms.VPSet(
        svFitLikelihoodDiTauKinematicsPhaseSpace,
        svFitLikelihoodDiTauKinematicsPolarized,
        svFitLikelihoodDiTauMEt,
        svFitLikelihoodDiTauPtBalance
    )
)
