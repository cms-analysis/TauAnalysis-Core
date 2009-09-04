#include "TauAnalysis/Core/plugins/MuonHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"

#include <TMath.h>

#include <stdlib.h>

bool matchesGenMuon(const pat::Muon& patMuon)
{
  //std::cout << "<matchesGenMuon>:" << std::endl;

  bool isGenMuonMatched = false;
  std::vector<reco::GenParticleRef> associatedGenParticles = patMuon.genParticleRefs();
  for ( std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin(); 
	it != associatedGenParticles.end(); ++it ) {
    if ( it->isAvailable() ) {
      const reco::GenParticleRef& genParticle = (*it);
      if ( genParticle->pdgId() == -13 || genParticle->pdgId() == +13 ) isGenMuonMatched = true;
    } else {
      edm::LogWarning("matchesGenMuon") << " edm::Ref of genParticle associated to pat::Muon is invalid !!";
    }
  }
  return isGenMuonMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

MuonHistManager::MuonHistManager(const edm::ParameterSet& cfg)
  : dqmError_(0)
{
  //std::cout << "<MuonHistManager::MuonHistManager>:" << std::endl;

  muonSrc_ = cfg.getParameter<edm::InputTag>("muonSource");
  //std::cout << " muonSrc = " << muonSrc_ << std::endl;

  vertexSrc_ = ( cfg.exists("vertexSource") ) ? cfg.getParameter<edm::InputTag>("vertexSource") : edm::InputTag();
  if ( vertexSrc_.label() == "" ) {
    edm::LogWarning("MuonHistManager") << " Configuration parameter 'vertexSource' not specified" 
				       << " --> Impact Parameter histograms will NOT be plotted !!";
  }
  //std::cout << " vertexSrc = " << vertexSrc_ << std::endl;

  genParticleSrc_ = ( cfg.exists("genParticleSource") ) ? cfg.getParameter<edm::InputTag>("genParticleSource") : edm::InputTag();
  if ( genParticleSrc_.label() == "" ) {
    edm::LogWarning("MuonHistManager") << " Configuration parameter 'genParticleSource' not specified" 
				       << " --> matching gen. Particle PdgId histogram will NOT be plotted !!";
  }
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;

  requireGenMuonMatch_ = cfg.getParameter<bool>("requireGenMuonMatch");
  //std::cout << " requireGenMuonMatch = " << requireGenMuonMatch_ << std::endl;

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "muons");

  makeIsoPtCtrlHistograms_ = ( cfg.exists("makeIsoPtCtrlHistograms") ) ? 
    cfg.getParameter<bool>("makeIsoPtCtrlHistograms") : false;

  makeIsoPtConeSizeDepHistograms_ = ( cfg.exists("makeIsoPtConeSizeDepHistograms") ) ? 
    cfg.getParameter<bool>("makeIsoPtConeSizeDepHistograms") : false;

  numMuonIsoConeSizes_ = 15;
  muonIsoConeSizeIncr_ = 0.1;
  numMuonIsoPtThresholds_ = 4;
  muonIsoPtThresholdIncr_ = 0.5;

//--- create "veto" objects for computation of IsoDeposit sums
  muonTrkIsoParam_.push_back(IsoDepositVetoFactory::make("0.02"));
  muonTrkIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(1.0)"));
  
  muonEcalIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonEcalIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.0)"));

  muonHcalIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonHcalIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.0)"));

  muonParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.5)"));
}

MuonHistManager::~MuonHistManager()
{
//--- delete "veto" objects for computation of IsoDeposit sums
  clearIsoParam(muonTrkIsoParam_);
  clearIsoParam(muonEcalIsoParam_);
  clearIsoParam(muonHcalIsoParam_);
  clearIsoParam(muonParticleFlowIsoParam_);
}

void MuonHistManager::bookHistograms()
{
  //std::cout << "<MuonHistManager::bookHistograms>:" << std::endl;

  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("bookHistograms") << " Failed to access dqmStore --> histograms will NOT be booked !!";
    dqmError_ = 1;
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_store_);

//--- book histogram for number of muons in each event
  hNumMuons_ = dqmStore.book1D("NumMuons", "NumMuons", 10, -0.5, 9.5);

//--- book histograms for Pt, eta and phi distributions
//    of muons passing all id. and isolation selections
  bookMuonHistograms(dqmStore, hMuonPt_, hMuonEta_, hMuonPhi_, "Muon");
  hMuonPtVsEta_ = dqmStore.book2D("MuonPtVsEta", "MuonPtVsEta", 24, -3., +3., 30, 0., 150.);
  hMuonCharge_ = dqmStore.book1D("MuonCharge", "Muon Charge", 3, -1.5, +1.5);
  
  hMuonPtCompToGen_ = dqmStore.book1D("MuonPtCompToGen", "MuonPtCompToGen", 200, -0.10, +0.10);
  hMuonThetaCompToGen_ = dqmStore.book1D("MuonThetaCompToGen", "MuonThetaCompToGen", 200, -0.010, +0.010);
  hMuonPhiCompToGen_ = dqmStore.book1D("MuonPhiCompToGen", "MuonPhiCompToGen", 200, -0.010, +0.010);
  
  hMuonMatchingGenParticlePdgId_ = dqmStore.book1D("MuonMatchingGenParticlePdgId", "matching gen. Particle PdgId", 26, -1.5, 24.5);

  hMuonTrackIPxy_ = dqmStore.book1D("MuonTrackIPxy", "MuonTrackIPxy", 100, -0.100, 0.100);
  hMuonTrackIPz_ = dqmStore.book1D("MuonTrackIPz", "MuonTrackIPz", 100, -1.0, 1.0);
  
  hMuonEcalDeposits_ = dqmStore.book1D("MuonEcalDeposits", "MuonEcalDeposits", 100, 0., 20.);
  hMuonHcalDeposits_ = dqmStore.book1D("MuonHcalDeposits", "MuonHcalDeposits", 100, 0., 20.);
  hMuonCaloDeposits_ = dqmStore.book1D("MuonCaloDeposits", "MuonCaloDeposits", 100, 0., 20.);
  hMuonCaloCompatibility_ = dqmStore.book1D("MuonCaloCompatibility", "MuonCaloCompatibility", 102, -0.01, 1.01);
  
  hMuonNumberOfChambers_ = dqmStore.book1D("MuonNumberOfChambers", "MuonNumberOfChambers", 25, -0.5, 24.5);
  hMuonSegmentCompatibility_ = dqmStore.book1D("MuonSegmentCompatibility", "MuonSegmentCompatibility", 102, -0.01, 1.01);
  
  hMuonTrkIsoPt_ = dqmStore.book1D("MuonTrkIsoPt", "MuonTrkIsoPt", 100, 0., 20.);    
  hMuonEcalIsoPt_ = dqmStore.book1D("MuonEcalIsoPt", "MuonEcalIsoPt", 100, 0., 20.);
  hMuonHcalIsoPt_ = dqmStore.book1D("MuonHcalIsoPt", "MuonHcalIsoPt", 100, 0., 20.);
  hMuonIsoSumPt_ = dqmStore.book1D("MuonIsoSumPt", "MuonIsoSumPt", 100, 0., 20.);
  hMuonTrkIsoPtRel_ = dqmStore.book1D("MuonTrkIsoPtRel", "MuonTrkIsoPtRel", 200, 0., 2.);    
  hMuonEcalIsoPtRel_ = dqmStore.book1D("MuonEcalIsoPtRel", "MuonEcalIsoPtRel", 200, 0., 2.);
  hMuonHcalIsoPtRel_ = dqmStore.book1D("MuonHcalIsoPtRel", "MuonHcalIsoPtRel", 200, 0., 2.);
  hMuonIsoSumPtRel_ = dqmStore.book1D("MuonIsoSumPtRel", "MuonIsoSumPtRel", 200, 0., 2.);

  hMuonParticleFlowIsoPt_ = dqmStore.book1D("MuonParticleFlowIsoPt", "MuonParticleFlowIsoPt", 100, 0., 20.);    
  hMuonPFChargedHadronIsoPt_ = dqmStore.book1D("MuonPFChargedHadronIsoPt", "MuonPFChargedHadronIsoPt", 100, 0., 20.);   
  hMuonPFNeutralHadronIsoPt_ = dqmStore.book1D("MuonPFNeutralHadronIsoPt", "MuonPFNeutralHadronIsoPt", 100, 0., 20.);   
  hMuonPFGammaIsoPt_ = dqmStore.book1D("MuonPFGammaIsoPt", "MuonPFGammaIsoPt", 100, 0., 20.);  
  
//--- book "control" histograms to check agreement between muon isolation variables
//    computed by PAT-level IsoDeposits based on particle flow 
//    with values computed on AOD level, based on ECAL recHits/CaloTowers and reco::Tracks
  if ( makeIsoPtCtrlHistograms_ ) {
    hMuonPFChargedHadronIsoPtCtrl_ = dqmStore.book2D("MuonPFChargedHadronIsoPtCtrl", "PFChargedHadron vs. reco::Track Isolation P_{T}", 40, 0., 20., 40, 0., 20.);
    hMuonPFGammaIsoPtCtrl_ = dqmStore.book2D("MuonPFGammaIsoPtCtrl", "PFGamma vs. ECAL recHit/CaloTower Isolation P_{T}", 40, 0., 20., 40, 0., 20.);
  } else {
    hMuonPFChargedHadronIsoPtCtrl_ = 0;
    hMuonPFGammaIsoPtCtrl_ = 0;
  } 

  hMuonTrkIsoValProfile_ = dqmStore.book1D("MuonTrkIsoValProfile", "MuonTrkIsoValProfile", 100, 0., 10.);
  hMuonTrkIsoEtaDistProfile_ = dqmStore.book1D("MuonTrkIsoEtaDistProfile", "MuonTrkIsoEtaDistProfile", 15, 0., 1.5);
  hMuonTrkIsoPhiDistProfile_ = dqmStore.book1D("MuonTrkIsoPhiDistProfile", "MuonTrkIsoPhiDistProfile", 15, 0., 1.5);
  
  hMuonEcalIsoValProfile_ = dqmStore.book1D("MuonEcalIsoValProfile", "MuonEcalIsoValProfile", 100, 0., 10.);
  hMuonEcalIsoEtaDistProfile_ = dqmStore.book1D("MuonEcalIsoEtaDistProfile", "MuonEcalIsoEtaDistProfile", 15, 0., 1.5);
  hMuonEcalIsoPhiDistProfile_ = dqmStore.book1D("MuonEcalIsoPhiDistProfile", "MuonEcalIsoPhiDistProfile", 15, 0., 1.5);
  
  hMuonHcalIsoValProfile_ = dqmStore.book1D("MuonHcalIsoValProfile", "MuonHcalIsoValProfile", 100, 0., 10.);
  hMuonHcalIsoEtaDistProfile_ = dqmStore.book1D("MuonHcalIsoEtaDistProfile", "MuonHcalIsoEtaDistProfile", 15, 0., 1.5);
  hMuonHcalIsoPhiDistProfile_  = dqmStore.book1D("MuonHcalIsoPhiDistProfile", "MuonHcalIsoPhiDistProfile", 15, 0., 1.5);
  
  if ( makeIsoPtConeSizeDepHistograms_ ) bookMuonIsoConeSizeDepHistograms(dqmStore);
}

double MuonHistManager::getMuonWeight(const pat::Muon& patMuon)
{
  return 1.;
}

void MuonHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<MuonHistManager::fillHistograms>:" << std::endl; 

  if ( dqmError_ ) {
    edm::LogError ("fillHistograms") << " Failed to access dqmStore --> histograms will NOT be filled !!";
    return;
  }

  edm::Handle<std::vector<pat::Muon> > patMuons;
  evt.getByLabel(muonSrc_, patMuons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(genParticleSrc_, genParticles);

  //std::cout << " patMuons.size = " << patMuons->size() << std::endl;
  hNumMuons_->Fill(patMuons->size(), evtWeight);

  double muonWeightSum = 0.;
  for ( std::vector<pat::Muon>::const_iterator patMuon = patMuons->begin(); 
	patMuon != patMuons->end(); ++patMuon ) {
    if ( requireGenMuonMatch_ && !matchesGenMuon(*patMuon) ) continue;

    muonWeightSum += getMuonWeight(*patMuon);
  }

  for ( std::vector<pat::Muon>::const_iterator patMuon = patMuons->begin(); 
	patMuon != patMuons->end(); ++patMuon ) {
  
    //bool isGenMuonMatched = matchesGenMuon(*patMuon);
    //std::cout << " Pt = " << patMuon->pt() << ", eta = " << patMuon->eta() << ", phi = " << patMuon->phi() << std::endl;
    //std::cout << " isGenMuonMatched = " << isGenMuonMatched << std::endl;
    
    if ( requireGenMuonMatch_ && !matchesGenMuon(*patMuon) ) continue;
    
    double weight = ( normMethod_ == kNormEvents ) ? evtWeight*(getMuonWeight(*patMuon)/muonWeightSum) : getMuonWeight(*patMuon);

    fillMuonHistograms(*patMuon, hMuonPt_, hMuonEta_, hMuonPhi_, weight);
    hMuonPtVsEta_->Fill(patMuon->eta(), patMuon->pt(), weight);
    hMuonCharge_->Fill(patMuon->charge(), weight);

//--- compare reconstructed muon to generator level one;
//    normalize difference between reconstructed and generated Pt
//    to expected Pt dependence of resolution
    if ( patMuon->genLepton() ) {
      hMuonPtCompToGen_->Fill((patMuon->pt() - patMuon->genLepton()->pt())/patMuon->genLepton()->pt(), weight);
      hMuonThetaCompToGen_->Fill(patMuon->theta() - patMuon->genLepton()->theta(), weight);
      hMuonPhiCompToGen_->Fill(patMuon->phi() - patMuon->genLepton()->phi(), weight);
    }

    int matchingGenParticlePdgId = getMatchingGenParticlePdgId(patMuon->p4(), genParticles);
    if ( matchingGenParticlePdgId == -1 ) {
      hMuonMatchingGenParticlePdgId_->Fill(-1, weight);
    } else if ( abs(matchingGenParticlePdgId) > 22 ) {
      hMuonMatchingGenParticlePdgId_->Fill(24, weight);
    } else {
      hMuonMatchingGenParticlePdgId_->Fill(abs(matchingGenParticlePdgId), weight);
    }

    if ( vertexSrc_.label() != "" && patMuon->track().isAvailable() && patMuon->track().isNonnull() ) {
      edm::Handle<std::vector<reco::Vertex> > recoVertices;
      evt.getByLabel(vertexSrc_, recoVertices);
      if ( recoVertices->size() >= 1 ) {
	const reco::Vertex& thePrimaryEventVertex = (*recoVertices->begin());
	hMuonTrackIPxy_->Fill(patMuon->track()->dxy(thePrimaryEventVertex.position()), weight);
	hMuonTrackIPz_->Fill(patMuon->track()->dz(thePrimaryEventVertex.position()), weight);
      }
    }

    const reco::Muon* recoMuon = dynamic_cast<const reco::Muon*>(patMuon->originalObject());
    if ( recoMuon != NULL ) {
      if ( recoMuon->isEnergyValid() ) {
	hMuonEcalDeposits_->Fill(recoMuon->calEnergy().em, weight);
	hMuonHcalDeposits_->Fill(recoMuon->calEnergy().had + recoMuon->calEnergy().ho, weight);
	hMuonCaloDeposits_->Fill(recoMuon->calEnergy().em + recoMuon->calEnergy().had + recoMuon->calEnergy().ho, weight);
      }
      if ( recoMuon->isCaloCompatibilityValid() ) hMuonCaloCompatibility_->Fill(recoMuon->caloCompatibility(), weight);
      
      hMuonNumberOfChambers_->Fill(recoMuon->numberOfChambers(), weight);
      double segmentCompatibility = muon::segmentCompatibility(*recoMuon);
      hMuonSegmentCompatibility_->Fill(segmentCompatibility, weight);
    } else {
      edm::LogError("analyze") << " Failed to access reco::Muon linked to pat::Muon object --> some histograms will NOT be filled !!";
    }

    fillMuonIsoHistograms(*patMuon, weight);
    if ( makeIsoPtConeSizeDepHistograms_ ) fillMuonIsoConeSizeDepHistograms(*patMuon, weight);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void MuonHistManager::bookMuonHistograms(DQMStore& dqmStore, MonitorElement*& hMuonPt, MonitorElement*& hMuonEta, MonitorElement*& hMuonPhi, const char* histoSetName)
{
  std::string hMuonPtName = std::string(histoSetName).append("Pt");
  hMuonPt = dqmStore.book1D(hMuonPtName, hMuonPtName, 75, 0., 150.);
  
  std::string hMuonEtaName = std::string(histoSetName).append("Eta");
  hMuonEta = dqmStore.book1D(hMuonEtaName, hMuonEtaName, 60, -3., +3.);
  
  std::string hMuonPhiName = std::string(histoSetName).append("Phi");
  hMuonPhi = dqmStore.book1D(hMuonPhiName, hMuonPhiName, 36, -TMath::Pi(), +TMath::Pi());
}

void MuonHistManager::bookMuonIsoConeSizeDepHistograms(DQMStore& dqmStore)
{
  for ( unsigned iConeSize = 1; iConeSize <= numMuonIsoConeSizes_; ++iConeSize ) {
    std::ostringstream iConeSizeString;
    iConeSizeString << std::setfill('0') << std::setw(2) << iConeSize;
    
    std::string hMuonTrkIsoPtConeSizeDepName_i 
      = std::string("MuonTrkIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonTrkIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonTrkIsoPtConeSizeDepName_i, 
							hMuonTrkIsoPtConeSizeDepName_i, 40, 0., 10.));
    std::string hMuonEcalIsoPtConeSizeDepName_i 
      = std::string("MuonEcalIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonEcalIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonEcalIsoPtConeSizeDepName_i, 
							 hMuonEcalIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonHcalIsoPtConeSizeDepName_i 
      = std::string("MuonHcalIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonHcalIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonHcalIsoPtConeSizeDepName_i, 
							 hMuonHcalIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonParticleFlowIsoPtConeSizeDepName_i 
      = std::string("MuonParticleFlowIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonParticleFlowIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonParticleFlowIsoPtConeSizeDepName_i, 
								 hMuonParticleFlowIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFChargedHadronIsoPtConeSizeDepName_i 
      = std::string("MuonChargedHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFChargedHadronIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonPFChargedHadronIsoPtConeSizeDepName_i, 
								    hMuonPFChargedHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFNeutralHadronIsoPtConeSizeDepName_i 
      = std::string("MuonPFNeutralHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFNeutralHadronIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonPFNeutralHadronIsoPtConeSizeDepName_i, 
								    hMuonPFNeutralHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFGammaIsoPtConeSizeDepName_i 
      = std::string("MuonPFGammaIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFGammaIsoPtConeSizeDep_.push_back(dqmStore.book1D(hMuonPFGammaIsoPtConeSizeDepName_i, 
							    hMuonPFGammaIsoPtConeSizeDepName_i, 40, 0., 10.));
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void MuonHistManager::fillMuonHistograms(const pat::Muon& patMuon, 
					 MonitorElement* hMuonPt, MonitorElement* hMuonEta, MonitorElement* hMuonPhi, 
					 double weight)
{
  //std::cout << "<MuonHistManager::fillMuonHistograms>:" << std::endl;

  hMuonPt->Fill(patMuon.pt(), weight);
  hMuonEta->Fill(patMuon.eta(), weight);
  hMuonPhi->Fill(patMuon.phi(), weight);
}

void MuonHistManager::fillMuonIsoHistograms(const pat::Muon& patMuon, double weight)
{
  //std::cout << "<MuonHistManager::fillMuonIsoHistograms>:" << std::endl;

  hMuonTrkIsoPt_->Fill(patMuon.trackIso(), weight);
  hMuonEcalIsoPt_->Fill(patMuon.ecalIso(), weight);
  hMuonHcalIsoPt_->Fill(patMuon.hcalIso(), weight);
  hMuonIsoSumPt_->Fill(patMuon.trackIso() + patMuon.ecalIso() + patMuon.hcalIso(), weight);
  hMuonTrkIsoPtRel_->Fill(patMuon.trackIso()/patMuon.pt(), weight);
  hMuonEcalIsoPtRel_->Fill(patMuon.ecalIso()/patMuon.pt(), weight);
  hMuonHcalIsoPtRel_->Fill(patMuon.hcalIso()/patMuon.pt(), weight);
  hMuonIsoSumPtRel_->Fill((patMuon.trackIso() + patMuon.ecalIso() + patMuon.hcalIso())/patMuon.pt(), weight);

  //std::cout << " particleIso = " << patMuon.particleIso() << std::endl;
  //std::cout << " chargedParticleIso = " << patMuon.chargedParticleIso() << std::endl;
  //std::cout << " neutralParticleIso = " << patMuon.neutralParticleIso() << std::endl;
  //std::cout << " gammaParticleIso = " << patMuon.gammaParticleIso() << std::endl;
  
  hMuonParticleFlowIsoPt_->Fill(patMuon.particleIso(), weight);
  hMuonPFChargedHadronIsoPt_->Fill(patMuon.chargedParticleIso(), weight);
  hMuonPFNeutralHadronIsoPt_->Fill(patMuon.neutralParticleIso(), weight);
  hMuonPFGammaIsoPt_->Fill(patMuon.gammaParticleIso(), weight);
  
  if ( makeIsoPtCtrlHistograms_ ) {
    hMuonPFChargedHadronIsoPtCtrl_->Fill(patMuon.trackIso(), patMuon.chargedParticleIso(), weight);
    hMuonPFGammaIsoPtCtrl_->Fill(patMuon.ecalIso(), patMuon.gammaParticleIso(), weight);
  }

  fillLeptonIsoDepositHistograms(patMuon.trackerIsoDeposit(), 
				 hMuonTrkIsoValProfile_, hMuonTrkIsoEtaDistProfile_, hMuonTrkIsoPhiDistProfile_, weight);
  fillLeptonIsoDepositHistograms(patMuon.ecalIsoDeposit(), 
				 hMuonEcalIsoValProfile_, hMuonEcalIsoEtaDistProfile_, hMuonEcalIsoPhiDistProfile_, weight);
  fillLeptonIsoDepositHistograms(patMuon.hcalIsoDeposit(), 
				 hMuonHcalIsoValProfile_, hMuonHcalIsoEtaDistProfile_, hMuonHcalIsoPhiDistProfile_, weight);
}

void MuonHistManager::fillMuonIsoConeSizeDepHistograms(const pat::Muon& patMuon, double weight)
{
  //std::cout << "<MuonHistManager::fillMuonIsoConeSizeDepHistograms>:" << std::endl;

  for ( unsigned iConeSize = 1; iConeSize <= numMuonIsoConeSizes_; ++iConeSize ) {
    double isoConeSize_i = iConeSize*muonIsoConeSizeIncr_;
    
    double muonTrkIsoDeposit_i = patMuon.trackerIsoDeposit()->countWithin(isoConeSize_i, muonTrkIsoParam_, false);
    hMuonTrkIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonTrkIsoDeposit_i, weight);
    
    double muonEcalIsoDeposit_i = patMuon.ecalIsoDeposit()->countWithin(isoConeSize_i, muonEcalIsoParam_, false);
    hMuonEcalIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonEcalIsoDeposit_i, weight);
    
    double muonHcalIsoDeposit_i = patMuon.hcalIsoDeposit()->countWithin(isoConeSize_i, muonHcalIsoParam_, false);
    hMuonHcalIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonHcalIsoDeposit_i, weight);

    if ( patMuon.isoDeposit(pat::ParticleIso) ) {
      double muonParticleFlowIsoDeposit_i 
	= patMuon.isoDeposit(pat::ParticleIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonParticleFlowIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonParticleFlowIsoDeposit_i, weight);
    }
    
    if ( patMuon.isoDeposit(pat::ChargedParticleIso) ) {
      double muonPFChargedHadronIsoDeposit_i 
	= patMuon.isoDeposit(pat::ChargedParticleIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFChargedHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFChargedHadronIsoDeposit_i, weight);
    }
    
    if ( patMuon.isoDeposit(pat::NeutralParticleIso) ) {
      double muonPFNeutralHadronIsoDeposit_i 
	= patMuon.isoDeposit(pat::NeutralParticleIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFNeutralHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFNeutralHadronIsoDeposit_i, weight);
    }

    if ( patMuon.isoDeposit(pat::GammaParticleIso) ) {
      double muonPFGammaIsoDeposit_i 
	= patMuon.isoDeposit(pat::GammaParticleIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFGammaIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFGammaIsoDeposit_i, weight);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, MuonHistManager, "MuonHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, MuonHistManager, "MuonHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<MuonHistManager> MuonAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(MuonAnalyzer);
