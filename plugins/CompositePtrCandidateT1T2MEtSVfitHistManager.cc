#include "TauAnalysis/Core/plugins/CompositePtrCandidateT1T2MEtSVfitHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"

#include <TMath.h>
#include <TFile.h>

const double epsilon = 0.01;

using namespace svFitHistogramManager_namespace;

template<typename T1, typename T2>
bool matchesGenCandidatePair(const CompositePtrCandidateT1T2MEt<T1,T2>& compositePtrCandidate)
{
  bool isGenMatched = false;
// not implemented yet...
  return isGenMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

svFitHistManagerType::svFitHistManagerType(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  algorithmName_ = cfg.getParameter<std::string>("algorithmName");    
  polarizationHypothesis_ = cfg.exists("polarizationHypothesis") ? 
    cfg.getParameter<std::string>("polarizationHypothesis") : "Unknown";
  
  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  dqmDirectory_store_ = dqmDirectoryName(dqmDirectory_store_).append(algorithmName_);
  if ( polarizationHypothesis_ != "Unknown" ) 
    dqmDirectory_store_ = dqmDirectoryName(dqmDirectory_store_).append(polarizationHypothesis_);
}

void svFitHistManagerType::bookHistogramsImp()
{
  hX1_ = book1D("X1", "X_{1}", 51, -0.01, 1.01);
  hX2_ = book1D("X2", "X_{2}", 51, -0.01, 1.01);
/*
  hX1vsGenX1_ = book2D("X1vsGenX1", "X_{1} vs. gen. X_{1}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX1vsGenX1Profile_ = bookProfile1D("X1vsGenX1Profile", "X_{1} vs. gen. X_{1}", 51, -0.01, 1.01);
  hX2vsGenX2_ = book2D("X2vsGenX2", "X_{2} vs. gen. X_{2}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX2vsGenX2Profile_ = bookProfile1D("X2vsGenX2Profile", "X_{2} vs. gen. X_{2}", 51, -0.01, 1.01);
 */
  hMass_ = book1D("Mass", "Mass", 50, 0., 250.);
  hMassActErr_ = book1D("MassActErr", "rec. - gen. Mass", 100, -125., +125.);
  hMassEstErr_ = book1D("MassEstErr", "estimated Uncertainty on rec. Mass", 100, -125., +125.);
  hMassRes_ = book1D("MassRes", "Mass Resolution", 100, -2.5, +2.5);
  hMassPull_ = book1D("MassPull", "(rec. - gen. Mass)/estimated Uncertainty", 100, -5.0, +5.0);
  hGenLeg1RecLeg2Mass_ = book1D("GenLeg1RecLeg2Mass", "gen. leg_{1} + rec. leg_{2} Invariant Mass", 50, 0., 250.);
  hRecLeg1GenLeg2Mass_ = book1D("RecLeg1GenLeg2Mass", "rec. leg_{1} + gen. leg_{2} Invariant Mass", 50, 0., 250.);

  hMassVsLogLikelihood_ = book2D("MassVsLogLikelihood", "SVfit Mass vs. log-Likelihood", 50, 0., 250., 20, -35., 25.);

  hLogLikelihood_ = book1D("LogLikelihood", "SVfit log-Likelihood", 100, -50., 50.);

  hDecayTimeLeg1_ = book1D("DecayTimeLeg1", "SVfit leg_{1} Decay eigentime", 100, 0., 1000.);
  hDecayTimeLeg2_ = book1D("DecayTimeLeg2", "SVfit leg_{2} Decay eigentime", 100, 0., 1000.);

  hSVfitStatus_ = book1D("SVfitStatus", "SVfit Status", 10, -2.5, 7.5);
}

template<typename T1, typename T2>
void svFitHistManagerType::customFillHistograms(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate, double weight)
{
  const SVfitDiTauSolution* svFitSolution = diTauCandidate.svFitSolution(algorithmName_, polarizationHypothesis_);
  
  if ( !svFitSolution ) {
    edm::LogError("<svFitHistManagerType::customFillHistograms>")
      << " No SVfitDiTauSolution object reconstructed for algorithm = " << algorithmName_ << "," 
      << " polarizaton hypothesis = " << polarizationHypothesis_ << " --> histograms will NOT be filled !!";
    return;
  }
  
  if ( svFitSolution->isValidSolution() ) {
    hX1_->Fill(svFitSolution->leg1().x(), weight);
    hX2_->Fill(svFitSolution->leg2().x(), weight);
/*
    if ( diTauCandidate->p4Leg1gen().energy() > epsilon && 
      diTauCandidate->p4Leg2gen().energy() > epsilon ) {
      hX1vsGenX1_->Fill(diTauCandidate->x1gen(), svFitSolution->leg1().x(), weight);
      hX1vsGenX1Profile_->getTProfile()->Fill(diTauCandidate->x1gen(), svFitSolution->leg1().x(), weight);
      hX2vsGenX2_->Fill(diTauCandidate->x2gen(), svFitSolution->leg2().x(), weight);
      hX2vsGenX2Profile_->getTProfile()->Fill(diTauCandidate->x2gen(), svFitSolution->leg2().x(), weight);
    }
 */	
    double genMass = diTauCandidate.p4gen().mass();
    double recMass = svFitSolution->mass();
    hMass_->Fill(recMass, weight);
    if ( genMass > 0. ) {
      hMassActErr_->Fill(recMass - genMass, weight);
      hMassRes_->Fill((recMass - genMass)/genMass, weight);
      if ( svFitSolution->hasErrorEstimates() ) {
	double recMassErr = svFitSolution->massErr();
	hMassEstErr_->Fill(recMassErr, weight);
	if ( recMassErr > 0. ) hMassPull_->Fill((recMass - genMass)/recMassErr, weight);
      }
    }
    hGenLeg1RecLeg2Mass_->Fill((diTauCandidate.p4Leg1gen() + svFitSolution->leg2().p4()).mass(), weight);
    hRecLeg1GenLeg2Mass_->Fill((svFitSolution->leg1().p4() + diTauCandidate.p4Leg2gen()).mass(), weight);

    hMassVsLogLikelihood_->Fill(svFitSolution->p4().mass(), svFitSolution->logLikelihood(), weight);

    hLogLikelihood_->Fill(svFitSolution->logLikelihood(), weight);

    hDecayTimeLeg1_->Fill(compDecayEigenTime(svFitSolution->leg1DecayVertex(),
					     svFitSolution->eventVertexPosSVrefitted(), svFitSolution->leg1().p4().energy()), weight);
    hDecayTimeLeg2_->Fill(compDecayEigenTime(svFitSolution->leg2DecayVertex(),
					     svFitSolution->eventVertexPosSVrefitted(), svFitSolution->leg2().p4().energy()), weight);
  }
    
  hSVfitStatus_->Fill(svFitSolution->minuitStatus(), weight);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

template<typename T1, typename T2>
CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::CompositePtrCandidateT1T2MEtSVfitHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::CompositePtrCandidateT1T2MEtSVfitHistManager>:" << std::endl;

  diTauCandidateSrc_ = cfg.getParameter<edm::InputTag>("diTauCandidateSource");
  //std::cout << " diTauCandidateSrc = " << diTauCandidateSrc_ << std::endl;

  diTauLeg1WeightExtractors_ = getTauJetWeightExtractors<T1>(cfg, "diTauLeg1WeightSource");
  diTauLeg2WeightExtractors_ = getTauJetWeightExtractors<T2>(cfg, "diTauLeg2WeightSource");

  requireGenMatch_ = cfg.getParameter<bool>("requireGenMatch");
  //std::cout << " requireGenMatch = " << requireGenMatch_ << std::endl;

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "diTauCandidates");

  std::string dqmDirectory_store = cfg.getParameter<std::string>("dqmDirectory_store");

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgSVfitAlgorithms = cfg.getParameter<vParameterSet>("SVfitAlgorithms");
  for ( vParameterSet::const_iterator cfgSVfitAlgorithm = cfgSVfitAlgorithms.begin();
	cfgSVfitAlgorithm != cfgSVfitAlgorithms.end(); ++cfgSVfitAlgorithm ) {
    std::string name = cfgSVfitAlgorithm->getParameter<std::string>("name");
    if ( cfgSVfitAlgorithm->exists("polarizationHypotheses") ) {
      typedef std::vector<std::string> vstring;
      vstring polarizationHypotheses = cfgSVfitAlgorithm->getParameter<vstring>("polarizationHypotheses");
      for ( vstring::const_iterator polarizationHypothesis = polarizationHypotheses.begin();
	    polarizationHypothesis != polarizationHypotheses.end(); ++polarizationHypothesis ) {
	edm::ParameterSet cfgSVfitAlgorithm_customized = (*cfgSVfitAlgorithm);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("algorithmName", name);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("polarizationHypothesis", *polarizationHypothesis);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("dqmDirectory_store", dqmDirectory_store);
	svFitAlgorithmHistManagers_.push_back(new svFitHistManagerType(cfgSVfitAlgorithm_customized));
      }
    } else {
      edm::ParameterSet cfgSVfitAlgorithm_customized = (*cfgSVfitAlgorithm);
      cfgSVfitAlgorithm_customized.addParameter<std::string>("algorithmName", name);
      cfgSVfitAlgorithm_customized.addParameter<std::string>("dqmDirectory_store", dqmDirectory_store);
      svFitAlgorithmHistManagers_.push_back(new svFitHistManagerType(cfgSVfitAlgorithm_customized));
    }
  }
}

template<typename T1, typename T2>
CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::~CompositePtrCandidateT1T2MEtSVfitHistManager()
{
  for ( typename std::vector<FakeRateJetWeightExtractor<T1>*>::iterator it = diTauLeg1WeightExtractors_.begin();
	it != diTauLeg1WeightExtractors_.end(); ++it ) {
    delete (*it);
  }

  for ( typename std::vector<FakeRateJetWeightExtractor<T2>*>::iterator it = diTauLeg2WeightExtractors_.begin();
	it != diTauLeg2WeightExtractors_.end(); ++it ) {
    delete (*it);
  }

  for ( std::vector<svFitHistManagerType*>::iterator it = svFitAlgorithmHistManagers_.begin();
	it != svFitAlgorithmHistManagers_.end(); ++it ) {
    delete (*it);
  }
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::bookHistogramsImp()
{
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::bookHistogramsImp>:" << std::endl;

  for ( std::vector<svFitHistManagerType*>::iterator svFitAlgorithmHistManager = svFitAlgorithmHistManagers_.begin();
	svFitAlgorithmHistManager != svFitAlgorithmHistManagers_.end(); ++svFitAlgorithmHistManager ) {
    (*svFitAlgorithmHistManager)->beginJob();
  }

  bookWeightHistograms(*dqmStore_, "DiTauCandidateWeight", "Composite Weight", 
		       hDiTauCandidateWeightPosLog_, hDiTauCandidateWeightNegLog_, hDiTauCandidateWeightZero_, 
		       hDiTauCandidateWeightLinear_);
}

template<typename T1, typename T2>
double CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::getDiTauCandidateWeight(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate)
{
  double diTauLeg1Weight = getTauJetWeight<T1>(*diTauCandidate.leg1(), diTauLeg1WeightExtractors_);
  double diTauLeg2Weight = getTauJetWeight<T2>(*diTauCandidate.leg2(), diTauLeg2WeightExtractors_);
  return (diTauLeg1Weight*diTauLeg2Weight);
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::fillHistogramsImp(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::fillHistogramsImp>:" << std::endl; 
  
  typedef std::vector<CompositePtrCandidateT1T2MEt<T1,T2> > CompositePtrCandidateCollection;
  edm::Handle<CompositePtrCandidateCollection> diTauCandidates;
  getCollection(evt, diTauCandidateSrc_, diTauCandidates);
  
  //std::cout << " diTauCandidateSrc = " << diTauCandidateSrc_.label() << ":" 
  //	      << " diTauCandidates.size = " << diTauCandidates->size() << std::endl;
  
  double diTauCandidateWeightSum = 0.;
  for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {
    if ( requireGenMatch_ && !matchesGenCandidatePair(*diTauCandidate) ) continue;
    
    diTauCandidateWeightSum += getDiTauCandidateWeight(*diTauCandidate);
  }
  
  for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {
    
    //bool isGenMatched = matchesGenCandidatePair(*diTauCandidate);
    //std::cout << " Pt = " << diTauCandidate->pt() << ", phi = " << diTauCandidate->phi() << "," 
    //          << " visMass = " << diTauCandidate->p4Vis().mass() << std::endl;
    //std::cout << " isGenMatched = " << isGenMatched << std::endl;
    
    if ( requireGenMatch_ && !matchesGenCandidatePair(*diTauCandidate) ) continue;
    
    double diTauCandidateWeight = getDiTauCandidateWeight(*diTauCandidate);
    double weight = getWeight(evtWeight, diTauCandidateWeight, diTauCandidateWeightSum);
      
    for ( std::vector<svFitHistManagerType*>::iterator svFitAlgorithmHistManager = svFitAlgorithmHistManagers_.begin();
	  svFitAlgorithmHistManager != svFitAlgorithmHistManagers_.end(); ++svFitAlgorithmHistManager ) {
      (*svFitAlgorithmHistManager)->customFillHistograms<T1,T2>(*diTauCandidate, weight);
    }

    fillWeightHistograms(hDiTauCandidateWeightPosLog_, hDiTauCandidateWeightNegLog_, hDiTauCandidateWeightZero_, 
			 hDiTauCandidateWeightLinear_, diTauCandidateWeight);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef CompositePtrCandidateT1T2MEtSVfitHistManager<reco::Candidate, reco::Candidate> DiCandidatePairSVfitHistManager;
typedef CompositePtrCandidateT1T2MEtSVfitHistManager<pat::Electron, pat::Tau> PATElecTauPairSVfitHistManager;
typedef CompositePtrCandidateT1T2MEtSVfitHistManager<pat::Muon, pat::Tau> PATMuTauPairSVfitHistManager;
typedef CompositePtrCandidateT1T2MEtSVfitHistManager<pat::Tau, pat::Tau> PATDiTauPairSVfitHistManager;
typedef CompositePtrCandidateT1T2MEtSVfitHistManager<pat::Electron, pat::Muon> PATElecMuPairSVfitHistManager;

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, DiCandidatePairSVfitHistManager, "DiCandidatePairSVfitHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, DiCandidatePairSVfitHistManager, "DiCandidatePairSVfitHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATElecTauPairSVfitHistManager, "PATElecTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATElecTauPairSVfitHistManager, "PATElecTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATMuTauPairSVfitHistManager, "PATMuTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATMuTauPairSVfitHistManager, "PATMuTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATDiTauPairSVfitHistManager, "PATDiTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATDiTauPairSVfitHistManager, "PATDiTauPairSVfitHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATElecMuPairSVfitHistManager, "PATElecMuPairSVfitHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATElecMuPairSVfitHistManager, "PATElecMuPairSVfitHistManager");
  
#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<DiCandidatePairSVfitHistManager> DiCandidatePairSVfitAnalyzer;
typedef HistManagerAdapter<PATElecTauPairSVfitHistManager> PATElecTauPairSVfitAnalyzer;
typedef HistManagerAdapter<PATMuTauPairSVfitHistManager> PATMuTauPairSVfitAnalyzer;
typedef HistManagerAdapter<PATDiTauPairSVfitHistManager> PATDiTauPairSVfitAnalyzer;
typedef HistManagerAdapter<PATElecMuPairSVfitHistManager> PATElecMuPairSVfitAnalyzer;

DEFINE_FWK_MODULE(DiCandidatePairSVfitAnalyzer);
DEFINE_FWK_MODULE(PATElecTauPairSVfitAnalyzer);
DEFINE_FWK_MODULE(PATMuTauPairSVfitAnalyzer);
DEFINE_FWK_MODULE(PATDiTauPairSVfitAnalyzer);
DEFINE_FWK_MODULE(PATElecMuPairSVfitAnalyzer);

