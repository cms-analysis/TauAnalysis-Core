#ifndef TauAnalysis_Core_histManagerAuxFunctions_h
#define TauAnalysis_Core_histManagerAuxFunctions_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "TauAnalysis/Core/interface/FakeRateJetWeightExtractor.h"

void fillLeptonIsoDepositHistograms(const pat::IsoDeposit*, MonitorElement*, MonitorElement*, MonitorElement*, double);

void clearIsoParam(reco::isodeposit::AbsVetos&);

int getMatchingGenParticlePdgId(const reco::Particle::LorentzVector&, edm::Handle<reco::GenParticleCollection>&);

double getDeltaRnearestJet(const reco::Particle::LorentzVector&, edm::Handle<pat::JetCollection>&, double = -1.);

template<typename T>
std::vector<FakeRateJetWeightExtractor<T>*> getTauJetWeightExtractors(const edm::ParameterSet& cfg, const char* label)
{
  std::vector<FakeRateJetWeightExtractor<T>*> tauJetWeightExtractors;
  
  if ( cfg.exists(label) ) {
    typedef std::vector<std::string> vstring;
    vstring tauJetWeightSources = cfg.getParameter<vstring>(label);
    for ( vstring::const_iterator tauJetWeightSource = tauJetWeightSources.begin();
	  tauJetWeightSource != tauJetWeightSources.end(); ++tauJetWeightSource ) {
      tauJetWeightExtractors.push_back(new FakeRateJetWeightExtractor<T>(*tauJetWeightSource));
    }
  }

  return tauJetWeightExtractors;
}

template<typename T>
double getTauJetWeight(const T& patTau, const std::vector<FakeRateJetWeightExtractor<T>*>& tauJetWeightExtractors)
{
  double tauJetWeight = 1.;
  for ( typename std::vector<FakeRateJetWeightExtractor<T>*>::const_iterator tauJetWeightExtractor = tauJetWeightExtractors.begin();
	tauJetWeightExtractor != tauJetWeightExtractors.end(); ++tauJetWeightExtractor ) {
    tauJetWeight *= (**tauJetWeightExtractor)(patTau);
  }

  return tauJetWeight;
}

#endif
