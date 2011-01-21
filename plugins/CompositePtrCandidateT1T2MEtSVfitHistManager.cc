#include "TauAnalysis/Core/plugins/CompositePtrCandidateT1T2MEtSVfitHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>

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

class SVfitHistManagerEntryBase : public HistManagerBase
{
 public:
  SVfitHistManagerEntryBase(const edm::ParameterSet&);
  ~SVfitHistManagerEntryBase();

  template<typename T1, typename T2>
  void customFillHistograms(const CompositePtrCandidateT1T2MEt<T1,T2>&, const reco::VertexCollection&, double);

 protected:
  std::string algorithmName_;
  std::string polarizationHypothesis_;

  typedef std::vector<std::string> vstring;
  vstring bestPolarizationHypothesis_;

  virtual void bookHistogramsImp();

  template<typename T1, typename T2>
  const SVfitDiTauSolution* getSVfitSolution(const CompositePtrCandidateT1T2MEt<T1,T2>&);

//-- dummy implementation of fillHistogramsImp function
//   declared as purely virtual function in HistManagerBase class
  void fillHistogramsImp(const edm::Event&, const edm::EventSetup&, double) {}

 private:
//--- histograms
  MonitorElement* hX1_;
  MonitorElement* hX2_;
/*
  MonitorElement* hX1vsGenX1_;
  MonitorElement* hX1vsGenX1Profile_;
  MonitorElement* hX2vsGenX2_;
  MonitorElement* hX2vsGenX2Profile_;
 */
  MonitorElement* hX1res_;
  MonitorElement* hX2res_;

  MonitorElement* hMass_;
  MonitorElement* hMassL_;
  MonitorElement* hMassXL_;
  MonitorElement* hMassActErr_;
  MonitorElement* hMassEstErr_;
  MonitorElement* hMassRes_;
  MonitorElement* hMassPull_;

  typedef std::vector<MonitorElement*> vMonitorElement;
  vMonitorElement hMassVsNumVertices_;
  typedef std::vector<double> vdouble;
  vdouble vertexPtThresholds_;

  MonitorElement* hPolarizationHypothesis_;

  MonitorElement* hGenLeg1RecLeg2Mass_;
  MonitorElement* hRecLeg1GenLeg2Mass_;

  MonitorElement* hMassVsNegLogLikelihood_;

  MonitorElement* hNegLogLikelihood_;

  MonitorElement* hDecayTimeLeg1_;
  MonitorElement* hDecayTimeLeg2_;

  MonitorElement* hSVfitStatus_;
};

SVfitHistManagerEntryBase::SVfitHistManagerEntryBase(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<SVfitHistManagerEntryBase::SVfitHistManagerEntryBase>:" << std::endl;

  algorithmName_ = cfg.getParameter<std::string>("algorithmName");
  //std::cout << " algorithmName = " << algorithmName_ << std::endl;
  polarizationHypothesis_ = cfg.exists("polarizationHypothesis") ?
    cfg.getParameter<std::string>("polarizationHypothesis") : "Unknown";
  //std::cout << " polarizationHypothesis = " << polarizationHypothesis_ << std::endl;

//--- check if polarizationHypothesis represents just a single hypothesis
//    or the "best" (i.e. the one with the highest likelihood) out of a set of different solutions
  TPRegexp regexpParser_bestPolarizationHypothesis("best\\{([LR]{2})(,[LR]{2})*\\}");

  TString polarizationHypothesis_tstring = polarizationHypothesis_.data();
  if ( regexpParser_bestPolarizationHypothesis.Match(polarizationHypothesis_tstring) >= 1 ) {
    TObjArray* subStrings = regexpParser_bestPolarizationHypothesis.MatchS(polarizationHypothesis_tstring);

    std::string subString = ((TObjString*)subStrings->At(0))->GetString().Data();
    //std::cout << " subString = " << subString << std::endl;

    size_t pos_start = subString.find("{", 0);
    size_t pos_end = subString.find("}", pos_start);
    subString = std::string(subString, pos_start + 1, pos_end - pos_start - 1);

    pos_start = 0;
    do {
      pos_end = subString.find(",", pos_start);
      bestPolarizationHypothesis_.push_back(std::string(subString, pos_start, pos_end - pos_start));
      pos_start = pos_end + 1;
    } while ( pos_end != std::string::npos );

    //std::cout << " bestPolarizationHypothesis = " << format_vstring(bestPolarizationHypothesis_)
    //	        << " (" << bestPolarizationHypothesis_.size() << " entries)" << std::endl;

//--- replace "special" characters { '{', ',', '}' } by underscores,
//    in order to compose a valid dqmDirectory name
    int errorFlag = 0;
    polarizationHypothesis_ = replace_string(polarizationHypothesis_, "{", "_", 0,  1, errorFlag);
    polarizationHypothesis_ = replace_string(polarizationHypothesis_, ",", "_", 0, 10, errorFlag);
    polarizationHypothesis_ = replace_string(polarizationHypothesis_, "}", "",  0,  1, errorFlag);
    //std::cout << " polarizationHypothesis(modified) = " << polarizationHypothesis_ << std::endl;
  }

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  dqmDirectory_store_ = dqmDirectoryName(dqmDirectory_store_).append(algorithmName_);
  if ( polarizationHypothesis_ != "Unknown" ) {
    dqmDirectory_store_ = dqmDirectoryName(dqmDirectory_store_).append(polarizationHypothesis_);
  }
}

SVfitHistManagerEntryBase::~SVfitHistManagerEntryBase()
{
// nothing to be done yet...
}

void setAxisLabelPolarizationHypothesis(TAxis* axis)
{
  axis->SetBinLabel(1, "Unknown");
  axis->SetBinLabel(3, "LL");
  axis->SetBinLabel(4, "LR");
  axis->SetBinLabel(5, "RL");
  axis->SetBinLabel(6, "RR");
}

void SVfitHistManagerEntryBase::bookHistogramsImp()
{
  hX1_ = book1D("X1", "X_{1}", 51, -0.01, 1.01);
  hX2_ = book1D("X2", "X_{2}", 51, -0.01, 1.01);
/*
  hX1vsGenX1_ = book2D("X1vsGenX1", "X_{1} vs. gen. X_{1}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX1vsGenX1Profile_ = bookProfile1D("X1vsGenX1Profile", "X_{1} vs. gen. X_{1}", 51, -0.01, 1.01);
  hX2vsGenX2_ = book2D("X2vsGenX2", "X_{2} vs. gen. X_{2}", 21, -0.01, 1.01, 100, -2.5, +2.5);
  hX2vsGenX2Profile_ = bookProfile1D("X2vsGenX2Profile", "X_{2} vs. gen. X_{2}", 51, -0.01, 1.01);
 */
  hX1res_ = book1D("X1res", "rec. X_{1} - gen. X_{1}", 201, -1.005, + 1.005);
  hX2res_ = book1D("X2res", "rec. X_{2} - gen. X_{2}", 201, -1.005, + 1.005);

  hMass_ = book1D("Mass", "Mass", 50, 0., 250.);
  hMassL_ = book1D("MassL", "Mass", 100, 0., 500.);
  hMassXL_ = book1D("MassXL", "Mass", 150, 0., 750.);
  hMassActErr_ = book1D("MassActErr", "rec. - gen. Mass", 100, -125., +125.);
  hMassEstErr_ = book1D("MassEstErr", "estimated Uncertainty on rec. Mass", 100, -125., +125.);
  hMassRes_ = book1D("MassRes", "Mass Resolution", 100, -2.5, +2.5);
  hMassPull_ = book1D("MassPull", "(rec. - gen. Mass)/estimated Uncertainty", 100, -5.0, +5.0);

  for ( vdouble::const_iterator vertexPtThreshold = vertexPtThresholds_.begin();
	vertexPtThreshold != vertexPtThresholds_.end(); ++vertexPtThreshold ) {
    std::ostringstream meName_ostringstream;
    meName_ostringstream << "MassVsNumVerticesPtGt" << std::fixed << std::setprecision(1) << (*vertexPtThreshold);
    int errorFlag = 0;
    std::string meName_string = replace_string(meName_ostringstream.str(), ".", "_", 0, 1, errorFlag);
    MonitorElement* me = book2D(meName_string.data(), meName_string.data(), 50, 0., 250., 10, -0.5, 9.5);
    hMassVsNumVertices_.push_back(me);
  }

  hPolarizationHypothesis_ = book1D("PolarizationHypothesis", "(best) Polarization hypothesis", 6, -0.5, 5.5);
  setAxisLabelPolarizationHypothesis(hPolarizationHypothesis_->getTH1()->GetXaxis());

  hGenLeg1RecLeg2Mass_ = book1D("GenLeg1RecLeg2Mass", "gen. leg_{1} + rec. leg_{2} Invariant Mass", 50, 0., 250.);
  hRecLeg1GenLeg2Mass_ = book1D("RecLeg1GenLeg2Mass", "rec. leg_{1} + gen. leg_{2} Invariant Mass", 50, 0., 250.);

  hMassVsNegLogLikelihood_ = book2D("MassVsNegLogLikelihood", "SVfit Mass vs. -log(Likelihood)", 50, 0., 250., 20, -35., 25.);

  hNegLogLikelihood_ = book1D("NegLogLikelihood", "SVfit -log(Likelihood)", 100, -50., 50.);

  hDecayTimeLeg1_ = book1D("DecayTimeLeg1", "SVfit leg_{1} Decay eigentime", 100, 0., 1000.);
  hDecayTimeLeg2_ = book1D("DecayTimeLeg2", "SVfit leg_{2} Decay eigentime", 100, 0., 1000.);

  hSVfitStatus_ = book1D("SVfitStatus", "SVfit Status", 10, -2.5, 7.5);
}

template<typename T1, typename T2>
const SVfitDiTauSolution* SVfitHistManagerEntryBase::getSVfitSolution(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate)
{
  const SVfitDiTauSolution* svFitSolution_best = 0;
  double negLogLikelihood_best = -1.;
  if ( bestPolarizationHypothesis_.size() > 0 ) {
    for ( vstring::const_iterator polarizationHypothesis_i = bestPolarizationHypothesis_.begin();
	  polarizationHypothesis_i != bestPolarizationHypothesis_.end(); ++polarizationHypothesis_i ) {
      const SVfitDiTauSolution* svFitSolution_i = diTauCandidate.svFitSolution(algorithmName_, *polarizationHypothesis_i);

      if ( !svFitSolution_i ) continue;

      double negLogLikelihood_i = svFitSolution_i->negLogLikelihood();
      if ( svFitSolution_best == 0 || negLogLikelihood_i < negLogLikelihood_best ) {
	svFitSolution_best = svFitSolution_i;
	negLogLikelihood_best = negLogLikelihood_i;
      }
    }
  } else {
    svFitSolution_best = diTauCandidate.svFitSolution(algorithmName_, polarizationHypothesis_);
  }

  return svFitSolution_best;
}

template<typename T1, typename T2>
void SVfitHistManagerEntryBase::customFillHistograms(
       const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate, const reco::VertexCollection& recoVertices, double weight)
{
  const SVfitDiTauSolution* svFitSolution = getSVfitSolution<T1,T2>(diTauCandidate);

  if ( !svFitSolution ) {
    edm::LogError("<SVfitHistManagerEntryBase::customFillHistograms>")
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
      hX1vsGenX1_->Fill(diTauCandidate.x1gen(), svFitSolution->leg1().x(), weight);
      hX1vsGenX1Profile_->getTProfile()->Fill(diTauCandidate.x1gen(), svFitSolution->leg1().x(), weight);
      hX2vsGenX2_->Fill(diTauCandidate.x2gen(), svFitSolution->leg2().x(), weight);
      hX2vsGenX2Profile_->getTProfile()->Fill(diTauCandidate.x2gen(), svFitSolution->leg2().x(), weight);
    }
 */
    hX1res_->Fill(svFitSolution->leg1().x() - diTauCandidate.x1gen(), weight);
    hX2res_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);

    double genMass = diTauCandidate.p4gen().mass();
    double recMass = svFitSolution->mass();
    hMass_->Fill(recMass, weight);
    hMassL_->Fill(recMass, weight);
    hMassXL_->Fill(recMass, weight);
    if ( genMass > 0. ) {
      hMassActErr_->Fill(recMass - genMass, weight);
      hMassRes_->Fill((recMass - genMass)/genMass, weight);
      if ( svFitSolution->hasErrorEstimates() ) {
	double recMassErr = svFitSolution->massErr();
	hMassEstErr_->Fill(recMassErr, weight);
	if ( recMassErr > 0. ) hMassPull_->Fill((recMass - genMass)/recMassErr, weight);
      }
    }

    std::vector<double> trackPtSums = compTrackPtSums(recoVertices);
    assert(vertexPtThresholds_.size() == hMassVsNumVertices_.size());
    size_t numVertexPtThresholds = vertexPtThresholds_.size();
    for ( size_t iVertexPtThreshold = 0; iVertexPtThreshold < numVertexPtThresholds; ++iVertexPtThreshold ) {
      size_t numVertices = getNumVerticesPtGtThreshold(trackPtSums, vertexPtThresholds_[iVertexPtThreshold]);
      hMassVsNumVertices_[iVertexPtThreshold]->Fill(recMass, numVertices, weight);
    }

    hPolarizationHypothesis_->getTH1()->Fill(svFitSolution->polarizationHypothesisName().data(), weight);

    hGenLeg1RecLeg2Mass_->Fill((diTauCandidate.p4Leg1gen() + svFitSolution->leg2().p4()).mass(), weight);
    hRecLeg1GenLeg2Mass_->Fill((svFitSolution->leg1().p4() + diTauCandidate.p4Leg2gen()).mass(), weight);

    hMassVsNegLogLikelihood_->Fill(svFitSolution->p4().mass(), svFitSolution->negLogLikelihood(), weight);

    hNegLogLikelihood_->Fill(svFitSolution->negLogLikelihood(), weight);

    hDecayTimeLeg1_->Fill(compDecayEigenTime(svFitSolution->leg1DecayDistance(), svFitSolution->leg1().p4().energy()), weight);
    hDecayTimeLeg2_->Fill(compDecayEigenTime(svFitSolution->leg2DecayDistance(), svFitSolution->leg2().p4().energy()), weight);
  }

  hSVfitStatus_->Fill(svFitSolution->minuitStatus(), weight);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

template<typename T1, typename T2>
class SVfitHistManagerEntryTemplateSpecific : public SVfitHistManagerEntryBase
{
 public:
  SVfitHistManagerEntryTemplateSpecific(const edm::ParameterSet& cfg)
    : SVfitHistManagerEntryBase(cfg)
  {}

  void customFillHistograms(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate,
			    const reco::VertexCollection& recoVertices, double weight)
  {
    SVfitHistManagerEntryBase::customFillHistograms<T1,T2>(diTauCandidate, recoVertices, weight);
  }

 private:
  void bookHistogramsImp()
  {
    SVfitHistManagerEntryBase::bookHistogramsImp();
  }
};

template <>
class SVfitHistManagerEntryTemplateSpecific<pat::Muon,pat::Tau> : public SVfitHistManagerEntryBase
{
 public:
  SVfitHistManagerEntryTemplateSpecific(const edm::ParameterSet& cfg)
    : SVfitHistManagerEntryBase(cfg)
  {}

  void customFillHistograms(const CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau>& diTauCandidate,
			    const reco::VertexCollection& recoVertices, double weight)
  {
    SVfitHistManagerEntryBase::customFillHistograms<pat::Muon,pat::Tau>(diTauCandidate, recoVertices, weight);

    const SVfitDiTauSolution* svFitSolution = getSVfitSolution<pat::Muon,pat::Tau>(diTauCandidate);

    if ( !svFitSolution ) {
      edm::LogError("<SVfitHistManagerEntryTemplateSpecific::customFillHistograms>")
	<< " No SVfitDiTauSolution object reconstructed for algorithm = " << algorithmName_ << ","
	<< " polarizaton hypothesis = " << polarizationHypothesis_ << " --> histograms will NOT be filled !!";
      return;
    }

    if ( svFitSolution->isValidSolution() && diTauCandidate.leg2()->genJet() != 0 ) {
      double recMass = svFitSolution->mass();

      int recTauDecayMode = diTauCandidate.leg2()->decayMode();
      std::string recTauDecayMode_string = getTauDecayModeName(recTauDecayMode);
      std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(*diTauCandidate.leg2()->genJet());

      if ( recTauDecayMode_string == genTauDecayMode ) {
	if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) {
	  hMassMuonOneProngNoPi0s_->Fill(recMass, weight);
	  hX2resOneProngNoPi0s_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);
	}
	if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero ) {
	  hMassMuonOneProngOnePi0_->Fill(recMass, weight);
	  hX2resOneProngOnePi0_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);
	}
	if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero ) {
	  hMassMuonOneProngTwoPi0s_->Fill(recMass, weight);
	  hX2resOneProngTwoPi0s_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);
	}
	if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero ) {
	  hMassMuonThreeProngNoPi0s_->Fill(recMass, weight);
	  hX2resThreeProngNoPi0s_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);
	}
	if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion1PiZero ) {
	  hMassMuonThreeProngOnePi0_->Fill(recMass, weight);
	  hX2resThreeProngOnePi0_->Fill(svFitSolution->leg2().x() - diTauCandidate.x2gen(), weight);
	}
      }
    }
  }

 private:
  void bookHistogramsImp()
  {
    SVfitHistManagerEntryBase::bookHistogramsImp();

    hMassMuonOneProngNoPi0s_   = book1D("MassMuonOneProngNoPi0s", "MassMuonOneProngNoPi0s", 50, 0., 250.);
    hMassMuonOneProngOnePi0_   = book1D("MassMuonOneProngOnePi0", "MassMuonOneProngOnePi0", 50, 0., 250.);
    hMassMuonOneProngTwoPi0s_  = book1D("MassMuonOneProngTwoPi0s", "MassMuonOneProngTwoPi0s", 50, 0., 250.);
    hMassMuonThreeProngNoPi0s_ = book1D("MassMuonThreeProngNoPi0s", "MassMuonThreeProngNoPi0s", 50, 0., 250.);
    hMassMuonThreeProngOnePi0_ = book1D("MassMuonThreeProngOnePi0", "MassMuonThreeProngOnePi0", 50, 0., 250.);

    hX2resOneProngNoPi0s_      = book1D("X2resOneProngNoPi0s", "X2resOneProngNoPi0s", 201, -1.005, + 1.005);
    hX2resOneProngOnePi0_      = book1D("X2resOneProngOnePi0", "X2resOneProngOnePi0", 201, -1.005, + 1.005);
    hX2resOneProngTwoPi0s_     = book1D("X2resOneProngTwoPi0s", "X2resOneProngTwoPi0s", 201, -1.005, + 1.005);
    hX2resThreeProngNoPi0s_    = book1D("X2resThreeProngNoPi0s", "X2resThreeProngNoPi0s", 201, -1.005, + 1.005);
    hX2resThreeProngOnePi0_    = book1D("X2resThreeProngOnePi0", "X2resThreeProngOnePi0", 201, -1.005, + 1.005);
  }

//--- histograms
  MonitorElement* hMassMuonOneProngNoPi0s_;
  MonitorElement* hMassMuonOneProngOnePi0_;
  MonitorElement* hMassMuonOneProngTwoPi0s_;
  MonitorElement* hMassMuonThreeProngNoPi0s_;
  MonitorElement* hMassMuonThreeProngOnePi0_;

  MonitorElement* hX2resOneProngNoPi0s_;
  MonitorElement* hX2resOneProngOnePi0_;
  MonitorElement* hX2resOneProngTwoPi0s_;
  MonitorElement* hX2resThreeProngNoPi0s_;
  MonitorElement* hX2resThreeProngOnePi0_;
};

//
//-----------------------------------------------------------------------------------------------------------------------
//

template<typename T1, typename T2>
CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::CompositePtrCandidateT1T2MEtSVfitHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::CompositePtrCandidateT1T2MEtSVfitHistManager>:" << std::endl;

  diTauCandidateSrc_ = cfg.getParameter<edm::InputTag>("diTauCandidateSource");
  //std::cout << " diTauCandidateSrc = " << diTauCandidateSrc_.label() << std::endl;

  genParticleSrc_ = cfg.getParameter<edm::InputTag>("genParticleSource");
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  vertexSrc_ = cfg.getParameter<edm::InputTag>("vertexSource");
  //std::cout << " vertexSrc = " << vertexSrc_ << std::endl;

  diTauLeg1WeightExtractors_ = getTauJetWeightExtractors<T1>(cfg, "diTauLeg1WeightSource");
  diTauLeg2WeightExtractors_ = getTauJetWeightExtractors<T2>(cfg, "diTauLeg2WeightSource");

  requireGenMatch_ = cfg.getParameter<bool>("requireGenMatch");
  //std::cout << " requireGenMatch = " << requireGenMatch_ << std::endl;

  pdgIdsElectron_.push_back(+11);
  pdgIdsElectron_.push_back(-11);
  pdgIdsMuon_.push_back(+13);
  pdgIdsMuon_.push_back(-13);
  pdgIdsPhoton_.push_back(22);
  for ( int iQuarkType = 1; iQuarkType <= 6; ++iQuarkType ) {
    pdgIdsJet_.push_back(+iQuarkType);
    pdgIdsJet_.push_back(-iQuarkType);
  }
  pdgIdsJet_.push_back(21);

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "diTauCandidates");

  std::string dqmDirectory_store = cfg.getParameter<std::string>("dqmDirectory_store");

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgSVfitAlgorithms = cfg.getParameter<vParameterSet>("SVfitAlgorithms");
  for ( vParameterSet::const_iterator cfgSVfitAlgorithm = cfgSVfitAlgorithms.begin();
	cfgSVfitAlgorithm != cfgSVfitAlgorithms.end(); ++cfgSVfitAlgorithm ) {
    std::string name = cfgSVfitAlgorithm->getParameter<std::string>("name");
    algorithmNames_.push_back(name);

    if ( cfgSVfitAlgorithm->exists("polarizationHypotheses") ) {
      vstring polarizationHypotheses = cfgSVfitAlgorithm->getParameter<vstring>("polarizationHypotheses");
      for ( vstring::const_iterator polarizationHypothesis = polarizationHypotheses.begin();
	    polarizationHypothesis != polarizationHypotheses.end(); ++polarizationHypothesis ) {
	edm::ParameterSet cfgSVfitAlgorithm_customized = (*cfgSVfitAlgorithm);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("algorithmName", name);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("polarizationHypothesis", *polarizationHypothesis);
	cfgSVfitAlgorithm_customized.addParameter<std::string>("dqmDirectory_store", dqmDirectory_store);
	svFitAlgorithmHistManagers_.push_back(new SVfitHistManagerEntryTemplateSpecific<T1,T2>(cfgSVfitAlgorithm_customized));

	if ( polarizationHypothesis->find("best") == std::string::npos ) polarizationHypotheses_[name].push_back(*polarizationHypothesis);
      }
    } else {
      edm::ParameterSet cfgSVfitAlgorithm_customized = (*cfgSVfitAlgorithm);
      cfgSVfitAlgorithm_customized.addParameter<std::string>("algorithmName", name);
      cfgSVfitAlgorithm_customized.addParameter<std::string>("dqmDirectory_store", dqmDirectory_store);
      svFitAlgorithmHistManagers_.push_back(new SVfitHistManagerEntryTemplateSpecific<T1,T2>(cfgSVfitAlgorithm_customized));
    }

    if ( cfgSVfitAlgorithm->exists("massHypotheses") ) {
      massHypotheses_[name] = cfgSVfitAlgorithm->getParameter<vdouble>("massHypotheses");
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

  for ( typename SVfitHistManagerEntryCollection::iterator it = svFitAlgorithmHistManagers_.begin();
	it != svFitAlgorithmHistManagers_.end(); ++it ) {
    delete (*it);
  }

  for ( typename std::map<std::string, std::vector<massHypothesisEntryType*> >::iterator it1 = massHypothesisEntries_.begin();
	it1 != massHypothesisEntries_.end(); ++it1 ) {
    for ( typename std::vector<massHypothesisEntryType*>::iterator it2 = it1->second.begin();
	  it2 != it1->second.end(); ++it2 ) {
      delete (*it2);
    }
  }
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtSVfitHistManager<T1,T2>::bookHistogramsImp()
{
  //std::cout << "<CompositePtrCandidateT1T2MEtSVfitHistManager::bookHistogramsImp>:" << std::endl;

  for ( typename SVfitHistManagerEntryCollection::iterator svFitAlgorithmHistManager = svFitAlgorithmHistManagers_.begin();
	svFitAlgorithmHistManager != svFitAlgorithmHistManagers_.end(); ++svFitAlgorithmHistManager ) {
    (*svFitAlgorithmHistManager)->beginJob();
  }

  for ( vstring::const_iterator algorithmName = algorithmNames_.begin();
	algorithmName != algorithmNames_.end(); ++algorithmName ) {
    std::string dqmDirectory_store_algorithm = dqmDirectoryName(dqmDirectory_store_).append(*algorithmName);
    dqmStore_->setCurrentFolder(dqmDirectory_store_algorithm);

    hMassLRvsRLbestLR_[*algorithmName] = book2D("MassLRvsRLbestLR", "Mass for LR vs. RL (-log(Likelihood_{LR} < -log(Likelihood_{RL}",
						50, 0., 250., 50, 0., 250.);
    hMassLRvsRLbestRL_[*algorithmName] = book2D("MassLRvsRLbestRL", "Mass for LR vs. RL (-log(Likelihood_{RL} < -log(Likelihood_{LR}",
						50, 0., 250., 50, 0., 250.);
    hMassLLvsRRbestLL_[*algorithmName] = book2D("MassLLvsRRbestLL", "Mass for LL vs. RR (-log(Likelihood_{LL} < -log(Likelihood_{RR}",
						50, 0., 250., 50, 0., 250.);
    hMassLLvsRRbestRR_[*algorithmName] = book2D("MassLLvsRRbestRR", "Mass for LL vs. RR (-log(Likelihood_{RR} < -log(Likelihood_{LL}",
						50, 0., 250., 50, 0., 250.);


    for ( vdouble::const_iterator massHypothesis = massHypotheses_[*algorithmName].begin();
	  massHypothesis != massHypotheses_[*algorithmName].end(); ++massHypothesis ) {
      std::ostringstream massHypothesisString;
      massHypothesisString << (*massHypothesis);

      massHypothesisEntryType* massHypothesisEntry = new massHypothesisEntryType();

      std::string hMassName = std::string("Mass").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMass_ = book1D(hMassName, hMassName, 50, 0., 250.);
      std::string hMassLname = std::string("MassL").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassL_ = book1D(hMassLname, hMassLname, 100, 0., 500.);
      std::string hMassXLname = std::string("MassXL").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassXL_ = book1D(hMassXLname, hMassXLname, 150, 0., 750.);

      std::string hMassGenLeg2ElectronName = std::string("MassGenLeg2Electron").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassGenLeg2Electron_ = book1D(hMassGenLeg2ElectronName, hMassGenLeg2ElectronName, 50, 0., 250.);
      std::string hMassGenLeg2MuonName = std::string("MassGenLeg2Muon").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassGenLeg2Muon_ = book1D(hMassGenLeg2MuonName, hMassGenLeg2MuonName, 50, 0., 250.);
      std::string hMassGenLeg2PhotonName = std::string("MassGenLeg2Photon").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassGenLeg2Photon_ = book1D(hMassGenLeg2PhotonName, hMassGenLeg2PhotonName, 50, 0., 250.);
      std::string hMassGenLeg2JetName = std::string("MassGenLeg2Jet").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hMassGenLeg2Jet_ = book1D(hMassGenLeg2JetName, hMassGenLeg2JetName, 50, 0., 250.);

      std::string hPolarizationHypothesisName = std::string("PolarizationHypothesis").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hPolarizationHypothesis_ = book1D(hPolarizationHypothesisName, hPolarizationHypothesisName, 6, -0.5, 5.5);
      setAxisLabelPolarizationHypothesis(massHypothesisEntry->hPolarizationHypothesis_->getTH1()->GetXaxis());

      std::string hX1resName = std::string("X1res").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hX1res_ = book1D(hX1resName, hX1resName, 201, -1.005, + 1.005);

      std::string hX2resName = std::string("X2res").append("_").append(massHypothesisString.str());
      massHypothesisEntry->hX2res_ = book1D(hX2resName, hX2resName, 201, -1.005, + 1.005);

      massHypothesisEntries_[*algorithmName].push_back(massHypothesisEntry);
    }
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

  edm::Handle<reco::GenParticleCollection> genParticles;
  if ( genParticleSrc_.label() != "" ) evt.getByLabel(genParticleSrc_, genParticles);

  edm::Handle<reco::VertexCollection> recoVertices;
  evt.getByLabel(vertexSrc_, recoVertices);

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

    for ( typename SVfitHistManagerEntryCollection::iterator svFitAlgorithmHistManager = svFitAlgorithmHistManagers_.begin();
	  svFitAlgorithmHistManager != svFitAlgorithmHistManagers_.end(); ++svFitAlgorithmHistManager ) {
      (*svFitAlgorithmHistManager)->customFillHistograms(*diTauCandidate, *recoVertices, weight);
    }

    int errorFlag;
    for ( vstring::const_iterator algorithmName = algorithmNames_.begin();
	  algorithmName != algorithmNames_.end(); ++algorithmName ) {
      const SVfitDiTauSolution* svFitSolution_LR = diTauCandidate->svFitSolution(*algorithmName, "LR", &errorFlag);
      const SVfitDiTauSolution* svFitSolution_RL = diTauCandidate->svFitSolution(*algorithmName, "RL", &errorFlag);
      if ( svFitSolution_LR && svFitSolution_RL ) {
	double mass_LR = svFitSolution_LR->mass();
	double negLogLikelihood_LR = svFitSolution_LR->negLogLikelihood();
	double mass_RL = svFitSolution_RL->mass();
	double negLogLikelihood_RL = svFitSolution_RL->negLogLikelihood();
	if      ( negLogLikelihood_LR < negLogLikelihood_RL) hMassLRvsRLbestLR_[*algorithmName]->Fill(mass_RL, mass_LR,  weight);
	else if ( negLogLikelihood_RL < negLogLikelihood_LR) hMassLRvsRLbestRL_[*algorithmName]->Fill(mass_RL, mass_LR,  weight);
      }

      const SVfitDiTauSolution* svFitSolution_LL = diTauCandidate->svFitSolution(*algorithmName, "LL", &errorFlag);
      const SVfitDiTauSolution* svFitSolution_RR = diTauCandidate->svFitSolution(*algorithmName, "RR", &errorFlag);
      if ( svFitSolution_LL && svFitSolution_RR ) {
	double mass_LL = svFitSolution_LL->mass();
	double negLogLikelihood_LL = svFitSolution_LL->negLogLikelihood();
	double mass_RR = svFitSolution_RR->mass();
	double negLogLikelihood_RR = svFitSolution_RR->negLogLikelihood();
	if      ( negLogLikelihood_LL < negLogLikelihood_RR) hMassLLvsRRbestLL_[*algorithmName]->Fill(mass_RR, mass_LL,  weight);
	else if ( negLogLikelihood_RR < negLogLikelihood_LL) hMassLLvsRRbestRR_[*algorithmName]->Fill(mass_RR, mass_LL,  weight);
      }

      vstring polarizationHypotheses = polarizationHypotheses_[*algorithmName];
      if ( polarizationHypotheses.size() > 0 ) {
	vdouble massHypotheses = massHypotheses_[*algorithmName];
	std::vector<massHypothesisEntryType*> massHypothesisEntries = massHypothesisEntries_[*algorithmName];
	assert(massHypotheses.size() == massHypothesisEntries.size());

	unsigned numMassHypotheses = massHypotheses.size();
	for ( unsigned iMassHypothesis = 0; iMassHypothesis < numMassHypotheses; ++iMassHypothesis ) {
	  double mass0 = massHypotheses[iMassHypothesis];

	  massHypothesisEntryType* massHypothesisEntry_i = massHypothesisEntries[iMassHypothesis];

	  const SVfitDiTauSolution* svFitSolution_best = 0;
	  double dMass_best = -1.;
	  for ( vstring::const_iterator polarizationHypothesis_i = polarizationHypotheses.begin();
		polarizationHypothesis_i != polarizationHypotheses.end(); ++polarizationHypothesis_i ) {
	    const SVfitDiTauSolution* svFitSolution_i = diTauCandidate->svFitSolution(*algorithmName, *polarizationHypothesis_i);
	    if ( !svFitSolution_i ) continue;

	    double dMass_i = TMath::Abs(mass0 - svFitSolution_i->mass());
	    if ( svFitSolution_best == 0 || dMass_i < dMass_best ) {
	      svFitSolution_best = svFitSolution_i;
	      dMass_best = dMass_i;
	    }
	  }

	  if ( svFitSolution_best ) {
	    double svFitMass = svFitSolution_best->mass();
	    massHypothesisEntry_i->hMass_->Fill(svFitMass, weight);
	    massHypothesisEntry_i->hMassL_->Fill(svFitMass, weight);
	    massHypothesisEntry_i->hMassXL_->Fill(svFitMass, weight);
	    if ( genParticles.isValid() ) {
	      fillHistogramGenMatch(massHypothesisEntry_i->hMassGenLeg2Electron_, svFitMass,
				    diTauCandidate->leg2()->p4(),  *genParticles, pdgIdsElectron_, weight);
	      fillHistogramGenMatch(massHypothesisEntry_i->hMassGenLeg2Muon_, svFitMass,
				    diTauCandidate->leg2()->p4(),  *genParticles, pdgIdsMuon_, weight);
	      fillHistogramGenMatch(massHypothesisEntry_i->hMassGenLeg2Photon_,svFitMass,
				    diTauCandidate->leg2()->p4(),  *genParticles, pdgIdsPhoton_, weight);
	      fillHistogramGenMatch(massHypothesisEntry_i->hMassGenLeg2Jet_, svFitMass,
				    diTauCandidate->leg2()->p4(),  *genParticles, pdgIdsJet_, weight);
	    }
	    std::string polarizationHypothesisName_best = svFitSolution_best->polarizationHypothesisName();
	    massHypothesisEntry_i->hPolarizationHypothesis_->getTH1()->Fill(polarizationHypothesisName_best.data(), weight);
	    massHypothesisEntry_i->hX1res_->Fill(svFitSolution_best->leg1().x() - diTauCandidate->x1gen(), weight);
	    massHypothesisEntry_i->hX2res_ ->Fill(svFitSolution_best->leg2().x() - diTauCandidate->x2gen(), weight);
	  }
	}
      }
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

