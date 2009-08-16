#include "TauAnalysis/Core/plugins/PATElecTauPairZeeHypothesisHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesisFwd.h"

#include <TMath.h>

//
//-----------------------------------------------------------------------------------------------------------------------
//

PATElecTauPairZeeHypothesisHistManager::PATElecTauPairZeeHypothesisHistManager(const edm::ParameterSet& cfg)
  : tauJetWeightExtractor_(0),
    dqmError_(0)
{
  //std::cout << "<PATElecTauPairZeeHypothesisHistManager::PATElecTauPairZeeHypothesisHistManager>:" << std::endl;

  ZeeHypothesisSrc_ = cfg.getParameter<edm::InputTag>("ZeeHypothesisSource");
  //std::cout << " ZeeHypothesisSrc = " << ZeeHypothesis_ << std::endl;

  if ( cfg.exists("tauJetWeightSource") ) {
    tauJetWeightSrc_ = cfg.getParameter<std::string>("tauJetWeightSource");
    tauJetWeightExtractor_ = new FakeRateJetWeightExtractor<pat::Tau>(tauJetWeightSrc_);
  }

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;
}

PATElecTauPairZeeHypothesisHistManager::~PATElecTauPairZeeHypothesisHistManager()
{
//--- nothing to be done yet...
}

void PATElecTauPairZeeHypothesisHistManager::bookHistograms()
{
  //std::cout << "<PATElecTauPairZeeHypothesisHistManager::bookHistograms>:" << std::endl;

  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("bookHistograms") << " Failed to access dqmStore --> histograms will NOT be booked !!";
    dqmError_ = 1;
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  
  dqmStore.setCurrentFolder(dqmDirectory_store_);

  hGenElectron1Pt_ = dqmStore.book1D("GenElectron1Pt", "gen. P_{T}^{e1}", 75, 0., 150.);
  hGenElectron1Eta_ = dqmStore.book1D("GenElectron1Eta", "gen. #eta_{e1}", 60, -3., +3.);
  hGenElectron1Phi_ = dqmStore.book1D("GenElectron1Phi", "gen. #phi_{e1}", 36, -TMath::Pi(), +TMath::Pi());

  hGenElectron2Pt_ = dqmStore.book1D("GenElectron2Pt", "gen. P_{T}^{e2}", 75, 0., 150.);
  hGenElectron2Eta_ = dqmStore.book1D("GenElectron2Eta", "gen. #eta_{e2}", 60, -3., +3.);
  hGenElectron2Phi_ = dqmStore.book1D("GenElectron2Phi", "gen. #phi_{e2}", 36, -TMath::Pi(), +TMath::Pi());

  hGenVisMass_ = dqmStore.book1D("GenVisMass", "gen. Z^{0} Mass", 40, 0., 200.);

  hElectron1bestMatchPt_ = dqmStore.book1D("Electron1bestMatchPt", "P_{T} of rec. Object best matching e_{1}", 75, 0., 150.);
  hElectron1bestMatchEta_= dqmStore.book1D("Electron1bestMatchEta", "#eta of rec. Object best matching e_{1}", 60, -3., +3.);
  hElectron1bestMatchPhi_ = dqmStore.book1D("Electron1bestMatchPhi", "#phi of rec. Object best matching e_{1}", 36, -TMath::Pi(), +TMath::Pi());
  hElectron1bestMatchType_ = dqmStore.book1D("Electron1bestMatchType", "Type of rec. Object best matching e_{1}", 10, -0.5, 9.5);

  hElectron2bestMatchPt_ = dqmStore.book1D("Electron2bestMatchPt", "P_{T} of rec. Object best matching e_{2}", 75, 0., 150.);
  hElectron2bestMatchEta_= dqmStore.book1D("Electron2bestMatchEta", "#eta of rec. Object best matching e_{2}", 60, -3., +3.);
  hElectron2bestMatchPhi_ = dqmStore.book1D("Electron2bestMatchPhi", "#phi of rec. Object best matching e_{2}", 36, -TMath::Pi(), +TMath::Pi());
  hElectron2bestMatchType_ = dqmStore.book1D("Electron2bestMatchType", "Type of rec. Object best matching e_{2}", 10, -0.5, 9.5);

  hElectron1bestMatchPtRes_ = dqmStore.book1D("Electron1bestMatchPtRes", "gen. P_{T}^{e1} - P_{T} of best matching rec. Object", 50, -25., +25.);
  hElectron1bestMatchEtaRes_ = dqmStore.book1D("Electron1bestMatchEtaRes", "gen. #eta_{e1} - #eta of best matching rec. Object", 100, -0.25, +0.25);
  hElectron1bestMatchPhiRes_ = dqmStore.book1D("Electron1bestMatchPhiRes", "gen. #phi_{e1} - #phi of best matching rec. Object", 100, -0.25, +0.25);

  hElectron2bestMatchPtRes_ = dqmStore.book1D("Electron2bestMatchPtRes", "gen. P_{T}^{e2} - P_{T} of best matching rec. Object", 50, -25., +25.);
  hElectron2bestMatchEtaRes_ = dqmStore.book1D("Electron2bestMatchEtaRes", "gen. #eta_{e2} - #eta of best matching rec. Object", 100, -0.25, +0.25);
  hElectron2bestMatchPhiRes_ = dqmStore.book1D("Electron2bestMatchPhiRes", "gen. #phi_{e2} - #phi of best matching rec. Object", 100, -0.25, +0.25);

  hVisMassBestMach_ = dqmStore.book1D("VisMassBestMach", "Z #rightarrow e^{+} e^{-} Mass hypothesis", 40, 0., 200.);

  hVisMassFromCaloJets_ = dqmStore.book1D("VisMassFromCaloJets", "hypothetic Z^{0} Mass from calo. Jets", 40, 0., 200.);
  hVisMassFromPFJets_ = dqmStore.book1D("VisMassFromPFJets", "hypothetic Z^{0} Mass from particle-flow Jets", 40, 0., 200.);
  hVisMassFromTracks_ = dqmStore.book1D("VisMassFromTracks", "hypothetic Z^{0} Mass from Tracks", 40, 0., 200.);
  hVisMassFromGsfElectrons_ = dqmStore.book1D("VisMassFromGsfElectrons", "hypothetic Z^{0} Mass from GSF Electrons", 40, 0., 200.);
  hVisMassFromGsfTracks_ = dqmStore.book1D("VisMassFromGsfTracks", "hypothetic Z^{0} Mass from GSF Tracks", 40, 0., 200.);
}

void PATElecTauPairZeeHypothesisHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<PATElecTauPairZeeHypothesisHistManager::fillHistograms>:" << std::endl; 

  if ( dqmError_ ) {
    edm::LogError ("fillHistograms") << " Failed to access dqmStore --> histograms will NOT be filled !!";
    return;
  }

  edm::Handle<PATElecTauPairZeeHypothesisCollection> ZeeHypotheses;
  evt.getByLabel(ZeeHypothesisSrc_, ZeeHypotheses);

  double tauJetWeightSum = 0.;
  if ( tauJetWeightExtractor_ ) {
    for ( PATElecTauPairZeeHypothesisCollection::const_iterator ZeeHypothesis = ZeeHypotheses->begin();
	  ZeeHypothesis != ZeeHypotheses->end(); ++ZeeHypothesis ) {
      tauJetWeightSum += (*tauJetWeightExtractor_)(*ZeeHypothesis->elecTauPair()->leg2());
    }
  }

  for ( PATElecTauPairZeeHypothesisCollection::const_iterator ZeeHypothesis = ZeeHypotheses->begin();
	ZeeHypothesis != ZeeHypotheses->end(); ++ZeeHypothesis ) {
    
    double weight = evtWeight;
    if ( tauJetWeightExtractor_ ) {
      double tauJetWeight = (*tauJetWeightExtractor_)(*ZeeHypothesis->elecTauPair()->leg2());
      weight *= (tauJetWeight/tauJetWeightSum);
    }

    hElectron1bestMatchPt_->Fill(ZeeHypothesis->p4Elec1bestMatch().pt(), weight);
    hElectron1bestMatchEta_->Fill(ZeeHypothesis->p4Elec1bestMatch().eta(), weight);
    hElectron1bestMatchPhi_->Fill(ZeeHypothesis->p4Elec1bestMatch().phi(), weight);
    hElectron1bestMatchType_->Fill(ZeeHypothesis->typeElec1bestMatch(), weight);

    hElectron2bestMatchPt_->Fill(ZeeHypothesis->p4Elec2bestMatch().pt(), weight);
    hElectron2bestMatchEta_->Fill(ZeeHypothesis->p4Elec2bestMatch().eta(), weight);
    hElectron2bestMatchPhi_->Fill(ZeeHypothesis->p4Elec2bestMatch().phi(), weight);
    hElectron2bestMatchType_->Fill(ZeeHypothesis->typeElec2bestMatch(), weight);

    if ( ZeeHypothesis->genElec1().isAvailable() ) {
      hGenElectron1Pt_->Fill(ZeeHypothesis->genElec1()->pt(), weight);
      hGenElectron1Eta_->Fill(ZeeHypothesis->genElec1()->eta(), weight);
      hGenElectron1Phi_->Fill(ZeeHypothesis->genElec1()->phi(), weight);

      hElectron1bestMatchPtRes_->Fill(ZeeHypothesis->genElec1()->pt() - ZeeHypothesis->p4Elec1bestMatch().pt(), weight);
      hElectron1bestMatchEtaRes_->Fill(ZeeHypothesis->genElec1()->eta() - ZeeHypothesis->p4Elec1bestMatch().eta(), weight);
      hElectron1bestMatchPhiRes_->Fill(ZeeHypothesis->genElec1()->phi() - ZeeHypothesis->p4Elec1bestMatch().phi(), weight);
    }

    if ( ZeeHypothesis->genElec2().isAvailable() ) {
      hGenElectron2Pt_->Fill(ZeeHypothesis->genElec2()->pt(), weight);
      hGenElectron2Eta_->Fill(ZeeHypothesis->genElec2()->eta(), weight);
      hGenElectron2Phi_->Fill(ZeeHypothesis->genElec2()->phi(), weight);

      hElectron2bestMatchPtRes_->Fill(ZeeHypothesis->genElec2()->pt() - ZeeHypothesis->p4Elec2bestMatch().pt(), weight);
      hElectron2bestMatchEtaRes_->Fill(ZeeHypothesis->genElec2()->eta() - ZeeHypothesis->p4Elec2bestMatch().eta(), weight);
      hElectron2bestMatchPhiRes_->Fill(ZeeHypothesis->genElec2()->phi() - ZeeHypothesis->p4Elec2bestMatch().phi(), weight);
    }

    if ( ZeeHypothesis->genElec1().isAvailable() &&
	 ZeeHypothesis->genElec2().isAvailable() ) {
      hGenVisMass_->Fill((ZeeHypothesis->genElec2()->p4() + ZeeHypothesis->genElec2()->p4()).mass(), weight);
    }

    hVisMassBestMach_->Fill(ZeeHypothesis->p4Z0bestMatch().mass(), weight);

    if ( ZeeHypothesis->elec1matchedCaloJet().isAvailable() &&
	 ZeeHypothesis->elec2matchedCaloJet().isAvailable() ) {
      hVisMassFromCaloJets_->Fill((ZeeHypothesis->elec1matchedCaloJet()->p4() + ZeeHypothesis->elec2matchedCaloJet()->p4()).mass(), weight);
    }

    if ( ZeeHypothesis->elec1matchedPFJet().isAvailable() &&
	 ZeeHypothesis->elec2matchedPFJet().isAvailable() ) {
      hVisMassFromPFJets_->Fill((ZeeHypothesis->elec1matchedPFJet()->p4() + ZeeHypothesis->elec2matchedPFJet()->p4()).mass(), weight);
    }

    if ( ZeeHypothesis->elec1matchedTrack().isAvailable() &&
	 ZeeHypothesis->elec2matchedTrack().isAvailable() ) {
      const edm::Ptr<reco::Track> track1 = ZeeHypothesis->elec1matchedTrack();
      reco::Particle::LorentzVector track1Momentum(track1->px(), track1->py(), track1->pz(), track1->p());
      const edm::Ptr<reco::Track> track2 = ZeeHypothesis->elec2matchedTrack();
      reco::Particle::LorentzVector track2Momentum(track2->px(), track2->py(), track2->pz(), track2->p());
      hVisMassFromTracks_->Fill((track1Momentum + track2Momentum).mass(), weight);
    }

    if ( ZeeHypothesis->elec1matchedGsfElectron().isAvailable() &&
	 ZeeHypothesis->elec2matchedGsfElectron().isAvailable() ) {
      hVisMassFromCaloJets_->Fill((ZeeHypothesis->elec1matchedGsfElectron()->p4() + ZeeHypothesis->elec2matchedGsfElectron()->p4()).mass(), weight);
    }

    if ( ZeeHypothesis->elec1matchedGsfTrack().isAvailable() &&
	 ZeeHypothesis->elec2matchedGsfTrack().isAvailable() ) {
      const edm::Ptr<reco::Track> gsfTrack1 = ZeeHypothesis->elec1matchedGsfTrack();
      reco::Particle::LorentzVector gsfTrack1Momentum(gsfTrack1->px(), gsfTrack1->py(), gsfTrack1->pz(), gsfTrack1->p());
      const edm::Ptr<reco::Track> gsfTrack2 = ZeeHypothesis->elec2matchedGsfTrack();
      reco::Particle::LorentzVector gsfTrack2Momentum(gsfTrack2->px(), gsfTrack2->py(), gsfTrack2->pz(), gsfTrack2->p());
      hVisMassFromTracks_->Fill((gsfTrack1Momentum + gsfTrack2Momentum).mass(), weight);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATElecTauPairZeeHypothesisHistManager, "PATElecTauPairZeeHypothesisHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATElecTauPairZeeHypothesisHistManager, "PATElecTauPairZeeHypothesisHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<PATElecTauPairZeeHypothesisHistManager> PATElecTauPairZeeHypothesisAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(PATElecTauPairZeeHypothesisAnalyzer);
