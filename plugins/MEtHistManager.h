#ifndef TauAnalysis_Core_MEtHistManager_h  
#define TauAnalysis_Core_MEtHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <vector>
#include <string>

class MEtHistManager : public HistManagerBase 
{
 public:
  
  explicit MEtHistManager(const edm::ParameterSet&);
  ~MEtHistManager();
  
  void bookHistograms(const edm::EventSetup&);
  void fillHistograms(const edm::Event&, const edm::EventSetup&);

 private:

//--- configuration parameters
  edm::InputTag metSrc_;

  std::string dqmDirectory_store_;

//--- histograms
  MonitorElement* hMEtPt_;
  MonitorElement* hMEtPhi_;
  MonitorElement* hMEtPtCompGen_;
  MonitorElement* hMEtPtRecVsGen_;
  MonitorElement* hMEtPhiCompGen_;
  MonitorElement* hMEtPhiRecVsGen_;
  MonitorElement* hMEtGenPt_;
  MonitorElement* hMEtGenPhi_;
};

#endif  


