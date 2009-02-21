#ifndef TauAnalysis_Core_JetHistManager_h  
#define TauAnalysis_Core_JetHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>
#include <string>

class JetHistManager : public HistManagerBase 
{
 public:
  
  explicit JetHistManager(const edm::ParameterSet&);
  ~JetHistManager();
  
  void bookHistograms(const edm::EventSetup&);
  void fillHistograms(const edm::Event&, const edm::EventSetup&);

 private:

//--- private functions
  void bookJetHistograms(DQMStore&, MonitorElement*&, MonitorElement*&, MonitorElement*&, const char*);
  
  void fillJetHistograms(const pat::Jet&, MonitorElement*, MonitorElement*, MonitorElement*);
  void fillNumCentralJetsToBeVetoesHistograms(const std::vector<pat::Jet>&, MonitorElement*, double, double, double);

//--- configuration parameters
  edm::InputTag jetSrc_;

  std::string dqmDirectory_store_;

  bool requireGenJetMatch_;

  typedef std::vector<double> vdouble;
  vdouble centralJetsToBeVetoedEtMin_;
  vdouble centralJetsToBeVetoedEtaMax_;
  vdouble centralJetsToBeVetoedAlphaMin_;

//--- histograms
  MonitorElement* hNumJets_;

  MonitorElement* hJetPt_;
  MonitorElement* hJetEta_;
  MonitorElement* hJetPtVsEta_;
  MonitorElement* hJetPhi_;
  
  MonitorElement* hJetAlpha_;
  MonitorElement* hJetNumTracks_;
  MonitorElement* hJetTrkPt_;
  MonitorElement* hJetLeadTrkPt_;

  std::vector<MonitorElement*> hNumCentralJetsToBeVetoed_;
};

#endif  


