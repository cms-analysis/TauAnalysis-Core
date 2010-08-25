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

#include <TMath.h>
#include <TFile.h>

const double epsilon = 0.01;

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
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::bookHistogramsImp()
{
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::bookHistogramsImp>:" << std::endl;

  hX1_ = book1D("X1", "X_{1}", 51, -0.01, 1.01);
  hX2_ = book1D("X2", "X_{2}", 51, -0.01, 1.01);
/*
  hX1vsGenX1_ = book2D("X1vsGenX1", "X_{1} vs. gen. X_{1}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX1vsGenX1Profile_ = bookProfile1D("X1vsGenX1Profile", "X_{1} vs. gen. X_{1}", 51, -0.01, 1.01);
  hX2vsGenX2_ = book2D("X2vsGenX2", "X_{2} vs. gen. X_{2}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX2vsGenX2Profile_ = bookProfile1D("X2vsGenX2Profile", "X_{2} vs. gen. X_{2}", 51, -0.01, 1.01);
 */
  hMass_ = book1D("Mass", " Mass", 50, 0., 250.);
  hMass1stSolution_ = book1D("Mass1stSolution", "SVfit Mass (1st Solution)", 50, 0., 250.);
  hGenLeg1RecLeg2Mass1stSolution_ = book1D("GenLeg1RecLeg2Mass1stSolution", 
					   "gen. leg_{1} + rec. leg_{2} Invariant Mass (1st Solution)", 50, 0., 250.);
  hRecLeg1GenLeg2Mass1stSolution_ = book1D("RecLeg1GenLeg2Mass1stSolution", 
					   "rec. leg_{1} + gen. leg_{2} Invariant Mass (1st Solution)", 50, 0., 250.);
  hMass2ndSolution_ = book1D("Mass2ndSolution", "SVfit  Mass (2nd Solution)", 50, 0., 250.);
  hGenLeg1RecLeg2Mass2ndSolution_ = book1D("GenLeg1RecLeg2Mass2ndSolution", 
					   "gen. leg_{1} + rec. leg_{2} Invariant Mass (2nd Solution)", 50, 0., 250.);
  hRecLeg1GenLeg2Mass2ndSolution_ = book1D("RecLeg1GenLeg2Mass2ndSolution", 
					   "rec. leg_{1} + gen. leg_{2} Invariant Mass (2nd Solution)", 50, 0., 250.);
  hMass3rdSolution_ = book1D("Mass3rdSolution", "SVfit  Mass (3rd Solution)", 50, 0., 250.);
  hGenLeg1RecLeg2Mass3rdSolution_ = book1D("GenLeg1RecLeg2Mass3rdSolution", 
					   "gen. leg_{1} + rec. leg_{2} Invariant Mass (3rd Solution)", 50, 0., 250.);
  hRecLeg1GenLeg2Mass3rdSolution_ = book1D("RecLeg1GenLeg2Mass3rdSolution", 
					   "rec. leg_{1} + gen. leg_{2} Invariant Mass (3rd Solution)", 50, 0., 250.);
  hMass4thSolution_ = book1D("Mass4thSolution", "SVfit  Mass (4th Solution)", 50, 0., 250.);  
  hGenLeg1RecLeg2Mass4thSolution_ = book1D("GenLeg1RecLeg2Mass4thSolution", 
					   "gen. leg_{1} + rec. leg_{2} Invariant Mass (4th Solution)", 50, 0., 250.);
  hRecLeg1GenLeg2Mass4thSolution_ = book1D("RecLeg1GenLeg2Mass4thSolution", 
					   "rec. leg_{1} + gen. leg_{2} Invariant Mass (4th Solution)", 50, 0., 250.);
  hMassBestMatch_ = book1D("MassBestMatch", "SVfit  Mass best matching gen. Mass", 50, 0., 250.);
  hMassAverage_ = book1D("MassAverage", "SVfit  Mass (Average of valid Solutions)", 50, 0., 250.);
  hMassNumSolutionsAveraged_ = book1D("MassNumSolutionsAveraged", "SVfit Num. Mass Solutions incl. in Average", 5, -0.5, 4.5);
  hMassVsLogLikelihood_ = book2D("MassVsLogLikelihood", "SVfit Mass vs. log-Likelihood", 50, 0., 250., 20, -35., 25.);
  hLogLikelihood_ = book1D("LogLikelihood", "SVfit log-Likelihood", 100, -50., 50.);
  hDecayTimeLeg1_ = book1D("DecayTimeLeg1", "SVfit leg_{1} Decay eigentime", 100, 0., 1000.);
  hDecayTimeLeg2_ = book1D("DecayTimeLeg2", "SVfit leg_{2} Decay eigentime", 100, 0., 1000.);
  hSVfitStatus_ = book1D("SVfitStatus", "SVfit Status", 10, -2.5, 7.5);

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
void fillSVmassRecoSolutionHistogram(unsigned iSolution, 
				     MonitorElement* hRecLeg1RecLeg2, MonitorElement* hGenLeg1RecLeg2, MonitorElement* hRecLeg1GenLeg2,  
				     const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate, double weight)
{
  const std::vector<SVmassRecoSolution>& svFitSolutions = diTauCandidate.svFitSolutions();
  if ( iSolution >= 0 && iSolution < svFitSolutions.size() ) {
    if ( svFitSolutions[iSolution].isValidSolution() ) {
      hRecLeg1RecLeg2->Fill(svFitSolutions[iSolution].p4().mass(), weight);
      hGenLeg1RecLeg2->Fill((diTauCandidate.p4Leg1gen() + svFitSolutions[iSolution].p4Leg2()).mass(), weight);
      hRecLeg1GenLeg2->Fill((svFitSolutions[iSolution].p4Leg1() + diTauCandidate.p4Leg2gen()).mass(), weight);
    }
  }
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
    //std::cout << " Pt = " << diTauCandidate->pt() << ", phi = " << diTauCandidate->phi() << ", visMass = " << diTauCandidate->p4Vis().mass() << std::endl;
    //std::cout << " isGenMatched = " << isGenMatched << std::endl;
    
    if ( requireGenMatch_ && !matchesGenCandidatePair(*diTauCandidate) ) continue;
    
    double diTauCandidateWeight = getDiTauCandidateWeight(*diTauCandidate);
    double weight = getWeight(evtWeight, diTauCandidateWeight, diTauCandidateWeightSum);
          
    double minLogLikelihood = 1.e+6;
      
    const std::vector<SVmassRecoSolution>& svFitSolutions = diTauCandidate->svFitSolutions();
    //std::cout << "svFitSolutions.size = " << svFitSolutions.size() << std::endl;
    for ( std::vector<SVmassRecoSolution>::const_iterator svFitSolution = svFitSolutions.begin();
	  svFitSolution != svFitSolutions.end(); ++svFitSolution ) {
      //std::cout << " svFitSolution->isValidSolution[" << svFitSolution - svFitSolutions.begin() << "]" 
      //	  << " = " << svFitSolution->isValidSolution() << std::endl;
      if ( svFitSolution->isValidSolution() ) {
	
	hX1_->Fill(svFitSolution->x1(), weight);
	hX2_->Fill(svFitSolution->x2(), weight);
/*
	if ( diTauCandidate->p4Leg1gen().energy() > epsilon && 
	     diTauCandidate->p4Leg2gen().energy() > epsilon ) {
	  hX1vsGenX1_->Fill(diTauCandidate->x1gen(), svFitSolution->x1(), weight);
	  hX1vsGenX1Profile_->getTProfile()->Fill(diTauCandidate->x1gen(), svFitSolution->x1(), weight);
	  hX2vsGenX2_->Fill(diTauCandidate->x2gen(), svFitSolution->x2(), weight);
	  hX2vsGenX2Profile_->getTProfile()->Fill(diTauCandidate->x2gen(), svFitSolution->x2(), weight);
	}
 */	
	hMass_->Fill(svFitSolution->p4().mass(), weight);
	hMassVsLogLikelihood_->Fill(svFitSolution->p4().mass(), svFitSolution->logLikelihood(), weight);
	hLogLikelihood_->Fill(svFitSolution->logLikelihood(), weight);
	double leg1TotEnergy = ( svFitSolution->x1() > 0 && svFitSolution->x1() <= 1 ) ?
	  svFitSolution->p4VisLeg1().energy()/svFitSolution->x1() : svFitSolution->p4VisLeg1().energy();
	hDecayTimeLeg1_->Fill(compDecayEigenTime(svFitSolution->decayVertexPosLeg1(), 
						 svFitSolution->primaryVertexPosSVrefitted(), leg1TotEnergy), weight);
	double leg2TotEnergy = ( svFitSolution->x2() > 0 && svFitSolution->x2() <= 1 ) ?
	  svFitSolution->p4VisLeg2().energy()/svFitSolution->x2() : svFitSolution->p4VisLeg2().energy();
	hDecayTimeLeg2_->Fill(compDecayEigenTime(svFitSolution->decayVertexPosLeg2(), 
						 svFitSolution->primaryVertexPosSVrefitted(), leg2TotEnergy), weight);
	
	if ( svFitSolution->logLikelihood() < minLogLikelihood ) minLogLikelihood = svFitSolution->logLikelihood();
      }
      
      hSVfitStatus_->Fill(svFitSolution->svFitStatus(), weight);
    }
    
    fillSVmassRecoSolutionHistogram(0, hMass1stSolution_, hGenLeg1RecLeg2Mass1stSolution_, hRecLeg1GenLeg2Mass1stSolution_, 
				    *diTauCandidate, weight);
    fillSVmassRecoSolutionHistogram(1, hMass2ndSolution_, hGenLeg1RecLeg2Mass2ndSolution_, hRecLeg1GenLeg2Mass2ndSolution_, 
				    *diTauCandidate, weight);
    fillSVmassRecoSolutionHistogram(2, hMass3rdSolution_, hGenLeg1RecLeg2Mass3rdSolution_, hRecLeg1GenLeg2Mass3rdSolution_, 
				    *diTauCandidate, weight);
    fillSVmassRecoSolutionHistogram(3, hMass4thSolution_, hGenLeg1RecLeg2Mass4thSolution_, hRecLeg1GenLeg2Mass4thSolution_, 
				    *diTauCandidate, weight);
    
    double genDiTauMass = diTauCandidate->p4gen().mass();
    double svFitMassBestMatch = 1.e+6;
    double svFitMassAverage = 0.;
    unsigned svFitMassNumSolutionsAveraged = 0;
    double svFitMassAverageNumSigmaCut = 2.;
    
    for ( std::vector<SVmassRecoSolution>::const_iterator svFitSolution = svFitSolutions.begin();
	  svFitSolution != svFitSolutions.end(); ++svFitSolution ) {
      if ( svFitSolution->isValidSolution() && svFitSolution->svFitStatus() == 0 ) {
	if ( TMath::Abs(svFitSolution->p4().mass() - genDiTauMass) < TMath::Abs(svFitMassBestMatch - genDiTauMass) ) {
	  svFitMassBestMatch = svFitSolution->p4().mass();
	}
	
	if ( svFitSolution->logLikelihood() < (0.5*svFitMassAverageNumSigmaCut*svFitMassAverageNumSigmaCut) ) {
	  svFitMassAverage += svFitSolution->p4().mass();
	  ++svFitMassNumSolutionsAveraged;
	}
      }
    }
    
    if ( svFitMassBestMatch != 1.e+6 ) hMassBestMatch_->Fill(svFitMassBestMatch, weight);
    
    if ( svFitMassNumSolutionsAveraged > 0 ) hMassAverage_->Fill(svFitMassAverage/svFitMassNumSolutionsAveraged, weight);
    hMassNumSolutionsAveraged_->Fill(svFitMassNumSolutionsAveraged, weight);
    
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

