#include "TauAnalysis/Core/plugins/MEtHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"

#include <TMath.h>

MEtHistManager::MEtHistManager(const edm::ParameterSet& cfg)
{
  //std::cout << "<MEtHistManager::MEtHistManager>:" << std::endl;

  metSrc_ = cfg.getParameter<edm::InputTag>("metSource");
  //std::cout << " metSrc = " << metSrc_ << std::endl;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;
}

MEtHistManager::~MEtHistManager()
{
//--- nothing to be done yet...
}

void MEtHistManager::bookHistograms(const edm::EventSetup& setup)
{
  //std::cout << "<MEtHistManager::bookHistograms>:" << std::endl;

  if ( edm::Service<DQMStore>().isAvailable() ) {
    DQMStore& dqmStore = (*edm::Service<DQMStore>());

    dqmStore.setCurrentFolder(dqmDirectory_store_);

    hMEtPt_ = dqmStore.book1D("MEtPt", "MEtPt", 75, 0., 150.);
    hMEtPhi_ = dqmStore.book1D("MEtPhi", "MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hMEtPtCompGen_ = dqmStore.book1D("MEtPtCompToGen", "MEtPtCompToGen", 100, -5.0, +5.0);
    hMEtPtRecVsGen_ = dqmStore.book2D("MEtPtRecVsGen", "MEtPtRecVsGen", 75, 0., 150., 75, 0., 150.);
    hMEtPhiCompGen_ = dqmStore.book1D("MEtPhiCompToGen", "MEtPhiCompToGen", 72, -TMath::Pi(), +TMath::Pi());
    hMEtPhiRecVsGen_ = dqmStore.book2D("MEtPhiRecVsGen", "MEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
    hMEtGenPt_ = dqmStore.book1D("MEtGenPt", "MEtGenPt", 75, 0., 150.);
    hMEtGenPhi_ = dqmStore.book1D("MEtGenPhi", "MEtGenPhi", 36, -TMath::Pi(), +TMath::Pi());
  }
}

void MEtHistManager::fillHistograms(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{  
  //std::cout << "<MEtHistManager::fillHistograms>:" << std::endl; 

  edm::Handle<std::vector<pat::MET> > patMETs;
  iEvent.getByLabel(metSrc_, patMETs);
  if ( patMETs->size() == 1 ) {
    const pat::MET& theEventMET = (*patMETs->begin());

    hMEtPt_->Fill(theEventMET.pt());
    hMEtPhi_->Fill(theEventMET.phi());

    if ( theEventMET.genMET() ) {
      hMEtPtCompGen_->Fill((theEventMET.pt() - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()));
      hMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), theEventMET.pt());
      hMEtPhiCompGen_->Fill(theEventMET.phi() - theEventMET.genMET()->phi());
      hMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), theEventMET.phi());
      hMEtGenPt_->Fill(theEventMET.genMET()->pt());
      hMEtGenPhi_->Fill(theEventMET.genMET()->phi());
    }
  } else {
    edm::LogError ("MEtHistManager::fillHistograms") << " Exactly one MET object expected per event --> skipping !!";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(HistManagerPluginFactory, MEtHistManager, "MEtHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<MEtHistManager> MEtAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(MEtAnalyzer);
