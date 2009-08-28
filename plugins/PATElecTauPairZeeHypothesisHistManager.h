#ifndef TauAnalysis_Core_PATElecTauPairZeeHypothesisHistManager_h  
#define TauAnalysis_Core_PATElecTauPairZeeHypothesisHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "TauAnalysis/Core/interface/FakeRateJetWeightExtractor.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesisFwd.h"

class PATElecTauPairZeeHypothesisHistManager : public HistManagerBase 
{
 public:  
  explicit PATElecTauPairZeeHypothesisHistManager(const edm::ParameterSet&);
  ~PATElecTauPairZeeHypothesisHistManager();
  
 private:
//--- histogram booking and filling functions 
//    inherited from HistManagerBase class
  void bookHistograms();

  double getTauWeight(const PATElecTauPairZeeHypothesis&);
  
  void fillHistograms(const edm::Event&, const edm::EventSetup&, double);

//--- configuration parameters
  edm::InputTag ZeeHypothesisSrc_;

  std::string tauJetWeightSrc_;

  std::string dqmDirectory_store_;

//--- "helper" class for accessing weight values
//    associated to second tau decay products
//    (efficiency/fake-rate with which the tau-jet passes the tau id. criteria)
  FakeRateJetWeightExtractor<pat::Tau>* tauJetWeightExtractor_;

//--- histograms
  MonitorElement* hGenElectron1Pt_;
  MonitorElement* hGenElectron1Eta_;
  MonitorElement* hGenElectron1Phi_;

  MonitorElement* hGenElectron2Pt_;
  MonitorElement* hGenElectron2Eta_;
  MonitorElement* hGenElectron2Phi_;

  MonitorElement* hGenVisMass_;

  MonitorElement* hElectron1bestMatchPt_;
  MonitorElement* hElectron1bestMatchEta_;
  MonitorElement* hElectron1bestMatchPhi_;
  MonitorElement* hElectron1bestMatchType_;

  MonitorElement* hElectron2bestMatchPt_;
  MonitorElement* hElectron2bestMatchEta_;
  MonitorElement* hElectron2bestMatchPhi_;
  MonitorElement* hElectron2bestMatchType_;

  MonitorElement* hElectron1bestMatchPtRes_;
  MonitorElement* hElectron1bestMatchEtaRes_;
  MonitorElement* hElectron1bestMatchPhiRes_;

  MonitorElement* hElectron2bestMatchPtRes_;
  MonitorElement* hElectron2bestMatchEtaRes_;
  MonitorElement* hElectron2bestMatchPhiRes_;

  MonitorElement* hVisMassBestMach_;
  
  MonitorElement* hVisMassFromCaloJets_;
  MonitorElement* hVisMassFromPFJets_;
  MonitorElement* hVisMassFromTracks_;
  MonitorElement* hVisMassFromGsfElectrons_;
  MonitorElement* hVisMassFromGsfTracks_;

  int dqmError_;
};

#endif  


