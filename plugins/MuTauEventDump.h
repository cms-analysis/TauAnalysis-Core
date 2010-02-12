#ifndef TauAnalysis_Core_MuTauEventDump_h  
#define TauAnalysis_Core_MuTauEventDump_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TauAnalysis/Core/interface/GenericEventDump.h"
#include "TauAnalysis/Core/interface/ObjectDumpBase.h"

class MuTauEventDump : public GenericEventDump
{
 public:  
  explicit MuTauEventDump(const edm::ParameterSet&);
  ~MuTauEventDump();

  void printDiTauCandidateInfo(const edm::Event& evt) const {
    printDiTauCandidateInfoImp<pat::Muon, pat::Tau>(evt);
  }
  
 protected:
  void print(const edm::Event&, const edm::EventSetup&, 
	     const filterResults_type&, const filterResults_type&, double) const;

  ObjectDumpBase* muonDump_;
  ObjectDumpBase* tauDump_;
  ObjectDumpBase* muTauDump_;

  void printMuTauZmumuHypothesisInfo(const edm::Event&) const;
  void printDiMuZmumuHypothesisInfo(const edm::Event&) const;

  edm::InputTag muTauZmumuHypothesisSource_;
  edm::InputTag diMuZmumuHypothesisSource_;
};

#endif  


