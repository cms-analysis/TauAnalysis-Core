#include "TauAnalysis/Core/plugins/SVfitLikelihoodAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <iostream>

template<typename T1, typename T2>
SVfitLikelihoodAnalyzer<T1,T2>::SVfitLikelihoodAnalyzer(const edm::ParameterSet& cfg)
{
  diTauCandidateSrc_ = cfg.getParameter<edm::InputTag>("diTauCandidateSource");

  typedef std::vector<edm::ParameterSet> vParameterSet;

  vParameterSet cfgSVfitAlgorithms = cfg.getParameter<vParameterSet>("svFitAlgorithms");
  for ( vParameterSet::const_iterator cfgSVfitAlgorithm = cfgSVfitAlgorithms.begin();
	cfgSVfitAlgorithm != cfgSVfitAlgorithms.end(); ++cfgSVfitAlgorithm ) {
    svFitAlgorithmType svfitAlgorithm;
    svfitAlgorithm.algorithmName_ = cfgSVfitAlgorithm->getParameter<std::string>("name");    
    if ( cfgSVfitAlgorithm->exists("polarizationHypotheses") ) {
      svfitAlgorithm.polarizationHypotheses_ = cfgSVfitAlgorithm->getParameter<vstring>("polarizationHypotheses");
    } else {
      svfitAlgorithm.polarizationHypotheses_.push_back(std::string("Unknown"));
    }
    svFitAlgorithms_.push_back(svfitAlgorithm);
  }

  vParameterSet cfgLikelihoodFunctions = cfg.getParameter<vParameterSet>("svFitLikelihoodFunctions");
  for ( vParameterSet::const_iterator cfgLikelihoodFunction = cfgLikelihoodFunctions.begin();
	cfgLikelihoodFunction != cfgLikelihoodFunctions.end(); ++cfgLikelihoodFunction ) {
    std::string pluginType = cfgLikelihoodFunction->getParameter<std::string>("pluginType");
    typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<T1,T2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
    SVfitDiTauLikelihoodBase<T1,T2>* svFitLikelihoodFunction 
	= SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgLikelihoodFunction);
    svFitLikelihoodFunctions_.push_back(svFitLikelihoodFunction);
  }
}

template<typename T1, typename T2>
SVfitLikelihoodAnalyzer<T1,T2>::~SVfitLikelihoodAnalyzer()
{
  for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::iterator it = svFitLikelihoodFunctions_.begin();
	it != svFitLikelihoodFunctions_.end(); ++it ) {
    delete (*it);
  }
}

template<typename T1, typename T2>
void SVfitLikelihoodAnalyzer<T1,T2>::analyze(const edm::Event& evt, const edm::EventSetup& es, double weight)
{
  std::cout << "<SVfitLikelihoodAnalyzer::analyze>:" << std::endl;

  for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator svFitLikelihoodFunction = svFitLikelihoodFunctions_.begin();
	svFitLikelihoodFunction != svFitLikelihoodFunctions_.end(); ++svFitLikelihoodFunction ) {
    (*svFitLikelihoodFunction)->beginEvent(evt, es);
  }

  typedef std::vector<CompositePtrCandidateT1T2MEt<T1,T2> > CompositePtrCandidateCollection;
  edm::Handle<CompositePtrCandidateCollection> diTauCandidates;
  pf::fetchCollection(diTauCandidates, diTauCandidateSrc_, evt);

  unsigned iDiTauCandidate = 0;
  for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {
    std::cout << "DiTauCandidate(" << iDiTauCandidate << "):" << std::endl;    
    std::cout << " Pt = " << diTauCandidate->pt() << std::endl;
    std::cout << " genPt = " << diTauCandidate->p4gen().pt() << std::endl;
    std::cout << " theta = " << diTauCandidate->theta()*180./TMath::Pi() 
	      << " (eta = " << diTauCandidate->eta() << ")" << std::endl;
    std::cout << " phi = " << diTauCandidate->phi()*180./TMath::Pi() << std::endl;
    
    for ( typename std::vector<svFitAlgorithmType>::const_iterator svFitAlgorithm = svFitAlgorithms_.begin();
	  svFitAlgorithm != svFitAlgorithms_.end(); ++svFitAlgorithm ) {
      for ( vstring::const_iterator polarizationHypothesis = svFitAlgorithm->polarizationHypotheses_.begin();
	    polarizationHypothesis != svFitAlgorithm->polarizationHypotheses_.end(); ++polarizationHypothesis ) {
	std::cout << "SVfit algorithm = " << svFitAlgorithm->algorithmName_ << "," 
		  << " polarization = " << (*polarizationHypothesis) << ":" << std::endl;
	const SVfitDiTauSolution* solution = diTauCandidate->svFitSolution(svFitAlgorithm->algorithmName_, *polarizationHypothesis);
	if ( !solution ) {
	  edm::LogWarning ("SVfitLikelihoodAnalyzer::analyze") 
	    << " Failed to find SVfit solution; algorithm name = " << svFitAlgorithm->algorithmName_ << ","
	    << " polarization hypothesis = " << (*polarizationHypothesis) << " !!";
	  continue;
	}
	
	for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator svFitLikelihoodFunction = svFitLikelihoodFunctions_.begin();
	      svFitLikelihoodFunction != svFitLikelihoodFunctions_.end(); ++svFitLikelihoodFunction ) {
	  std::cout << " nLL(" << (*svFitLikelihoodFunction)->name() << ") = " 
		    << (**svFitLikelihoodFunction)(*diTauCandidate, *solution) << std::endl;
	}
      }
    }

    ++iDiTauCandidate;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SVfitLikelihoodAnalyzer<reco::Candidate, reco::Candidate> SVfitLikelihoodDiCandidatePairAnalyzer;
typedef SVfitLikelihoodAnalyzer<pat::Electron, pat::Tau> SVfitLikelihoodElecTauPairAnalyzer;
typedef SVfitLikelihoodAnalyzer<pat::Muon, pat::Tau> SVfitLikelihoodMuTauPairAnalyzer;
typedef SVfitLikelihoodAnalyzer<pat::Tau, pat::Tau> SVfitLikelihoodDiTauPairAnalyzer;
typedef SVfitLikelihoodAnalyzer<pat::Electron, pat::Muon> SVfitLikelihoodElecMuPairAnalyzer;

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, SVfitLikelihoodDiCandidatePairAnalyzer, "SVfitLikelihoodDiCandidatePairAnalyzer");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, SVfitLikelihoodElecTauPairAnalyzer, "SVfitLikelihoodElecTauPairAnalyzer");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, SVfitLikelihoodMuTauPairAnalyzer, "SVfitLikelihoodMuTauPairAnalyzer");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, SVfitLikelihoodDiTauPairAnalyzer, "SVfitLikelihoodDiTauPairAnalyzer");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, SVfitLikelihoodElecMuPairAnalyzer, "SVfitLikelihoodElecMuPairAnalyzer");
  
