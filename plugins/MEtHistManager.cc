#include "TauAnalysis/Core/plugins/MEtHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"

#include <TMath.h>

MEtHistManager::MEtHistManager(const edm::ParameterSet& cfg)
  : dqmError_(0)
{
  //std::cout << "<MEtHistManager::MEtHistManager>:" << std::endl;

  metSrc_ = cfg.getParameter<edm::InputTag>("metSource");
  //std::cout << " metSrc = " << metSrc_ << std::endl;
  
  metSignificanceSrc_ = cfg.getParameter<edm::InputTag>("metSignificanceSource");
  //std::cout << " metSignificanceSrc = " << metSignificanceSrc_ << std::endl;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;
}

MEtHistManager::~MEtHistManager()
{
//--- nothing to be done yet...
}

void MEtHistManager::bookHistograms()
{
  //std::cout << "<MEtHistManager::bookHistograms>:" << std::endl;

  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("bookHistograms") << " Failed to access dqmStore --> histograms will NOT be booked !!";
    dqmError_ = 1;
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  
  dqmStore.setCurrentFolder(dqmDirectory_store_);
  
  hRAWplusJESplusMUONplusTAU_MEtPt_ = dqmStore.book1D("RAWplusJESplusMUONplusTAU_MEtPt", "RAWplusJESplusMUONplusTAU_MEtPt", 75, 0., 150.);
  hRAWplusJESplusMUONplusTAU_MEtPhi_ = dqmStore.book1D("RAWplusJESplusMUONplusTAU_MEtPhi", "RAWplusJESplusMUONplusTAU_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUONplusTAU_MEtPx_ = dqmStore.book1D("RAWplusJESplusMUONplusTAU_MEtPx", "RAWplusJESplusMUONplusTAU_MEtPx", 150, -150., 150.);
  hRAWplusJESplusMUONplusTAU_MEtPy_ = dqmStore.book1D("RAWplusJESplusMUONplusTAU_MEtPy", "RAWplusJESplusMUONplusTAU_MEtPy", 150, -150., 150.);
  
  hRAWplusJESplusMUON_MEtPt_ = dqmStore.book1D("RAWplusJESplusMUON_MEtPt", "RAWplusJESplusMUON_MEtPt", 75, 0., 150.);
  hRAWplusJESplusMUON_MEtPhi_ = dqmStore.book1D("RAWplusJESplusMUON_MEtPhi", "RAWplusJESplusMUON_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUON_MEtPx_ = dqmStore.book1D("RAWplusJESplusMUON_MEtPx", "RAWplusJESplusMUON_MEtPx", 150, -150., 150.);
  hRAWplusJESplusMUON_MEtPy_ = dqmStore.book1D("RAWplusJESplusMUON_MEtPy", "RAWplusJESplusMUON_MEtPy", 150, -150., 150.);
  
  hRAWplusMUONplusTAU_MEtPt_ = dqmStore.book1D("RAWplusMUONplusTAU_MEtPt", "RAWplusMUONplusTAU_MEtPt", 75, 0., 150.);
  hRAWplusMUONplusTAU_MEtPhi_ = dqmStore.book1D("RAWplusMUONplusTAU_MEtPhi", "RAWplusMUONplusTAU_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
  hRAWplusMUONplusTAU_MEtPx_ = dqmStore.book1D("RAWplusMUONplusTAU_MEtPx", "RAWplusMUONplusTAU_MEtPx", 150, -150., 150.);
  hRAWplusMUONplusTAU_MEtPy_ = dqmStore.book1D("RAWplusMUONplusTAU_MEtPy", "RAWplusMUONplusTAU_MEtPy", 150, -150., 150.);
  
  hRAW_MEtPt_ = dqmStore.book1D("RAW_MEtPt", "RAW_MEtPt", 75, 0., 150.);
  hRAW_MEtPhi_ = dqmStore.book1D("RAW_MEtPhi", "RAW_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
  hRAW_MEtPx_ = dqmStore.book1D("RAW_MEtPx", "RAW_MEtPx", 150, -150., 150.);
  hRAW_MEtPy_ = dqmStore.book1D("RAW_MEtPy", "RAW_MEtPy", 150, -150., 150.);
  
  hRAW_MEtSignificance_ = dqmStore.book1D("RAW_MEtSignificance", "RAW_MEtSignificance", 101, -0.5, 100.05);

  hRAWplusJES_MEtPt_ = dqmStore.book1D("RAWplusJES_MEtPt", "RAWplusJES_MEtPt", 75, 0., 150.);
  hRAWplusJES_MEtPhi_ = dqmStore.book1D("RAWplusJES_MEtPhi", "RAWplusJES_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
  hRAWplusJES_MEtPx_ = dqmStore.book1D("RAWplusJES_MEtPx", "RAWplusJES_MEtPx", 150, -150., 150.);
  hRAWplusJES_MEtPy_ = dqmStore.book1D("RAWplusJES_MEtPy", "RAWplusJES_MEtPy", 150, -150., 150.);
  
  hMUON_MExCorrection_ = dqmStore.book1D("MUON_MExCorrection", "MUON_MExCorrection", 150, -150., 150.);
  hMUON_MEyCorrection_ = dqmStore.book1D("MUON_MEyCorrection", "MUON_MEyCorrection", 150, -150., 150.);
  hTAU_MExCorrection_ = dqmStore.book1D("TAU_MExCorrection", "TAU_MExCorrection", 150, -150., 150.);
  hTAU_MEyCorrection_ = dqmStore.book1D("TAU_MEyCorrection", "TAU_MEyCorrection", 150, -150., 150.);
  hJES_MExCorrection_ = dqmStore.book1D("JES_MExCorrection", "JES_MExCorrection", 150, -150., 150.);
  hJES_MEyCorrection_ = dqmStore.book1D("JES_MEyCorrection", "JES_MEyCorrection", 150, -150., 150.);
  
  hRAWplusJESplusMUONplusTAUMEtPtCompGen_ = dqmStore.book1D("RAWplusJESplusMUONplusTAUMEtPtCompGen", "RAWplusJESplusMUONplusTAUMEtPtCompGen", 100, -5.0, +5.0);
  hRAWplusJESplusMUONplusTAUMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONplusTAUMEtPtRecVsGen", "RAWplusJESplusMUONplusTAUMEtPtRecVsGen", 30, 0., 150., 30, 0., 150.);
  hRAWplusJESplusMUONplusTAUMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESplusMUONplusTAUMEtPhiCompGen", "RAWplusJESplusMUONplusTAUMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUONplusTAUMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONplusTAUMEtPhiRecVsGen", "RAWplusJESplusMUONplusTAUMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWplusJESplusMUONMEtPtCompGen_ = dqmStore.book1D("RAWplusJESplusMUONMEtPtCompGen", "RAWplusJESplusMUONMEtPtCompGen", 100, -5.0, +5.0);
  hRAWplusJESplusMUONMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONMEtPtRecVsGen", "RAWplusJESplusMUONMEtPtRecVsGen", 30, 0., 150., 30, 0., 150.);
  hRAWplusJESplusMUONMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESplusMUONMEtPhiCompGen", "RAWplusJESplusMUONMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUONMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONMEtPhiRecVsGen", "RAWplusJESplusMUONMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWplusJESMEtPtCompGen_ = dqmStore.book1D("RAWplusJESMEtPtCompGen", "RAWplusJESMEtPtCompGen", 100, -5.0, +5.0);
  hRAWplusJESMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESMEtPtRecVsGen", "RAWplusJESMEtPtRecVsGen", 30, 0., 150., 30, 0., 150.);
  hRAWplusJESMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESMEtPhiCompGen", "RAWplusJESMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESMEtPhiRecVsGen", "RAWplusJESMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWMEtPtCompGen_ = dqmStore.book1D("RAWMEtPtCompGen", "RAWMEtPtCompGen", 100, -5.0, +5.0);
  hRAWMEtPtRecVsGen_ = dqmStore.book2D("RAWMEtPtRecVsGen", "RAWMEtPtRecVsGen", 30, 0., 150., 30, 0., 150.);
  hRAWMEtPhiCompGen_ = dqmStore.book1D("RAWMEtPhiCompGen", "RAWMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWMEtPhiRecVsGen_ = dqmStore.book2D("RAWMEtPhiRecVsGen", "RAWMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Phi_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Phi", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWplusJESplusMUONMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Pt", "GenMEtDeltaRAWplusJESplusMUONMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONMEt_Phi_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Phi", "GenMEtDeltaRAWplusJESplusMUONMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
  hGenMEtDeltaRAWplusJESplusMUONMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Px", "GenMEtDeltaRAWplusJESplusMUONMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Py", "GenMEtDeltaRAWplusJESplusMUONMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWplusJESMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Pt", "GenMEtDeltaRAWplusJESMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESMEt_Phi_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Phi", "GenMEtDeltaRAWplusJESMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
  hGenMEtDeltaRAWplusJESMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Px", "GenMEtDeltaRAWplusJESMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Py", "GenMEtDeltaRAWplusJESMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Pt", "GenMEtDeltaRAWMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWMEt_Phi_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Phi", "GenMEtDeltaRAWMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
  hGenMEtDeltaRAWMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Px", "GenMEtDeltaRAWMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Py", "GenMEtDeltaRAWMEt_Py", 150, -150.0, 150.0);
  
  hGenMEt_Pt_ = dqmStore.book1D("GenMEt_Pt", "GenMEt_Pt", 75, 0., 150.);
  hGenMEt_Phi_ = dqmStore.book1D("GenMEt_Phi", "GenMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
}

void MEtHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)

{  
  //std::cout << "<MEtHistManager::fillHistograms>:" << std::endl; 

  if ( dqmError_ ) {
    edm::LogError ("fillHistograms") << " Failed to access dqmStore --> histograms will NOT be filled !!";
    return;
  }

  edm::Handle<std::vector<pat::MET> > patMETs;
  evt.getByLabel(metSrc_, patMETs);

  edm::Handle<std::vector<reco::CaloMET> > recoMETs;
  evt.getByLabel(metSignificanceSrc_, recoMETs);
  double RAW_MEtSignificance = ( recoMETs->size() == 1 ) ? recoMETs->begin()->metSignificance() : -1.;

  if ( patMETs->size() == 1 ) {
    const pat::MET& theEventMET = (*patMETs->begin());

    double RAWplusJESplusMUONplusTAU_MEtPt = theEventMET.pt();
    double RAWplusJESplusMUONplusTAU_MEtPhi = theEventMET.phi();
    double RAWplusJESplusMUONplusTAU_MEtPx = RAWplusJESplusMUONplusTAU_MEtPt * cos(RAWplusJESplusMUONplusTAU_MEtPhi);
    double RAWplusJESplusMUONplusTAU_MEtPy = RAWplusJESplusMUONplusTAU_MEtPt * sin(RAWplusJESplusMUONplusTAU_MEtPhi);
    math::XYZTLorentzVector RAWplusJESplusMUONplusTAU_MetVector(RAWplusJESplusMUONplusTAU_MEtPx , RAWplusJESplusMUONplusTAU_MEtPy , 0 , 0);
    hRAWplusJESplusMUONplusTAU_MEtPt_->Fill(RAWplusJESplusMUONplusTAU_MEtPt, evtWeight);
    hRAWplusJESplusMUONplusTAU_MEtPhi_->Fill(RAWplusJESplusMUONplusTAU_MEtPhi, evtWeight);
    hRAWplusJESplusMUONplusTAU_MEtPx_->Fill(RAWplusJESplusMUONplusTAU_MEtPx, evtWeight);
    hRAWplusJESplusMUONplusTAU_MEtPy_->Fill(RAWplusJESplusMUONplusTAU_MEtPy, evtWeight);

    double RAWplusJESplusMUON_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrTAU);
    double RAWplusJESplusMUON_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrTAU);
    double RAWplusJESplusMUON_MEtPx = RAWplusJESplusMUON_MEtPt * cos(RAWplusJESplusMUON_MEtPhi);
    double RAWplusJESplusMUON_MEtPy = RAWplusJESplusMUON_MEtPt * sin(RAWplusJESplusMUON_MEtPhi);
    math::XYZTLorentzVector RAWplusJESplusMUON_MetVector(RAWplusJESplusMUON_MEtPx , RAWplusJESplusMUON_MEtPy , 0 , 0);
    hRAWplusJESplusMUON_MEtPt_->Fill(RAWplusJESplusMUON_MEtPt, evtWeight);
    hRAWplusJESplusMUON_MEtPhi_->Fill(RAWplusJESplusMUON_MEtPhi, evtWeight);
    hRAWplusJESplusMUON_MEtPx_->Fill(RAWplusJESplusMUON_MEtPx, evtWeight);
    hRAWplusJESplusMUON_MEtPy_->Fill(RAWplusJESplusMUON_MEtPy, evtWeight);

    double RAWplusMUONplusTAU_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrJES);
    double RAWplusMUONplusTAU_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrJES);
    double RAWplusMUONplusTAU_MEtPx = RAWplusMUONplusTAU_MEtPt * cos(RAWplusMUONplusTAU_MEtPhi);
    double RAWplusMUONplusTAU_MEtPy = RAWplusMUONplusTAU_MEtPt * sin(RAWplusMUONplusTAU_MEtPhi);
    math::XYZTLorentzVector RAWplusMUONplusTAU_MetVector(RAWplusMUONplusTAU_MEtPx , RAWplusMUONplusTAU_MEtPy , 0 , 0);
    hRAWplusMUONplusTAU_MEtPt_->Fill(RAWplusMUONplusTAU_MEtPt, evtWeight);
    hRAWplusMUONplusTAU_MEtPhi_->Fill(RAWplusMUONplusTAU_MEtPhi, evtWeight);
    hRAWplusMUONplusTAU_MEtPx_->Fill(RAWplusMUONplusTAU_MEtPx, evtWeight);
    hRAWplusMUONplusTAU_MEtPy_->Fill(RAWplusMUONplusTAU_MEtPy, evtWeight);

    double RAW_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrALL);
    double RAW_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrALL);
    double RAW_MEtPx = RAW_MEtPt * cos(RAW_MEtPhi);
    double RAW_MEtPy = RAW_MEtPt * sin(RAW_MEtPhi);
    math::XYZTLorentzVector RAW_MetVector(RAW_MEtPx , RAW_MEtPy , 0 , 0);
    hRAW_MEtPt_->Fill(RAW_MEtPt, evtWeight);
    hRAW_MEtPhi_->Fill(RAW_MEtPhi, evtWeight);
    hRAW_MEtPx_->Fill(RAW_MEtPx, evtWeight);
    hRAW_MEtPy_->Fill(RAW_MEtPy, evtWeight);

    hRAW_MEtSignificance_->Fill(RAW_MEtSignificance, evtWeight);

    math::XYZTLorentzVector deltaTAU_MetVector = RAWplusJESplusMUONplusTAU_MetVector - RAWplusJESplusMUON_MetVector;
    math::XYZTLorentzVector deltaJES_MetVector = RAWplusJESplusMUONplusTAU_MetVector - RAWplusMUONplusTAU_MetVector;
    math::XYZTLorentzVector deltaMUON_MetVector = RAWplusJESplusMUONplusTAU_MetVector - deltaJES_MetVector - deltaTAU_MetVector - RAW_MetVector;

    math::XYZTLorentzVector RAWplusJES_MetVector = RAW_MetVector + deltaJES_MetVector;
    double RAWplusJES_MEtPt = RAWplusJES_MetVector.Pt();
    double RAWplusJES_MEtPhi = RAWplusJES_MetVector.Phi();
    double RAWplusJES_MEtPx = RAWplusJES_MetVector.Px();
    double RAWplusJES_MEtPy = RAWplusJES_MetVector.Py();
    hRAWplusJES_MEtPt_->Fill(RAWplusJES_MEtPt, evtWeight);
    hRAWplusJES_MEtPhi_->Fill(RAWplusJES_MEtPhi, evtWeight);
    hRAWplusJES_MEtPx_->Fill(RAWplusJES_MEtPx, evtWeight);
    hRAWplusJES_MEtPy_->Fill(RAWplusJES_MEtPy, evtWeight);

    hMUON_MExCorrection_->Fill(deltaMUON_MetVector.Px(), evtWeight);
    hMUON_MEyCorrection_->Fill(deltaMUON_MetVector.Py(), evtWeight);
    hTAU_MExCorrection_->Fill(deltaTAU_MetVector.Px(), evtWeight);
    hTAU_MEyCorrection_->Fill(deltaTAU_MetVector.Py(), evtWeight);
    hJES_MExCorrection_->Fill(deltaJES_MetVector.Px(), evtWeight);
    hJES_MEyCorrection_->Fill(deltaJES_MetVector.Py(), evtWeight);

    //MUONplusTAU_MEtPhi = (MUONplusTAU_MEtPhi>=0) ? MUONplusTAU_MEtPhi : (MUONplusTAU_MEtPhi + (2.0 * TMath::Pi()));

    if ( theEventMET.genMET() ) {

      hRAWplusJESplusMUONplusTAUMEtPtCompGen_->Fill((RAWplusJESplusMUONplusTAU_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()), evtWeight);
      hRAWplusJESplusMUONplusTAUMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJESplusMUONplusTAU_MEtPt, evtWeight);
      hRAWplusJESplusMUONplusTAUMEtPhiCompGen_->Fill(RAWplusJESplusMUONplusTAU_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hRAWplusJESplusMUONplusTAUMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJESplusMUONplusTAU_MEtPhi, evtWeight);
      hRAWplusJESplusMUONMEtPtCompGen_->Fill((RAWplusJESplusMUON_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()), evtWeight);
      hRAWplusJESplusMUONMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJESplusMUON_MEtPt, evtWeight);
      hRAWplusJESplusMUONMEtPhiCompGen_->Fill(RAWplusJESplusMUON_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hRAWplusJESplusMUONMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJESplusMUON_MEtPhi, evtWeight);
      hRAWplusJESMEtPtCompGen_->Fill((RAWplusJES_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()), evtWeight);
      hRAWplusJESMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJES_MEtPt, evtWeight);
      hRAWplusJESMEtPhiCompGen_->Fill(RAWplusJES_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hRAWplusJESMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJES_MEtPhi, evtWeight);
      hRAWMEtPtCompGen_->Fill((RAW_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()), evtWeight);
      hRAWMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAW_MEtPt, evtWeight);
      hRAWMEtPhiCompGen_->Fill(RAW_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hRAWMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAW_MEtPhi, evtWeight);

      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt_->Fill(RAWplusJESplusMUONplusTAU_MEtPt - theEventMET.genMET()->pt(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px_->Fill(RAWplusJESplusMUONplusTAU_MEtPx - theEventMET.genMET()->px(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py_->Fill(RAWplusJESplusMUONplusTAU_MEtPy - theEventMET.genMET()->py(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONMEt_Pt_->Fill(RAWplusJESplusMUON_MEtPt - theEventMET.genMET()->pt(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONMEt_Px_->Fill(RAWplusJESplusMUON_MEtPx - theEventMET.genMET()->px(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONMEt_Py_->Fill(RAWplusJESplusMUON_MEtPy - theEventMET.genMET()->py(), evtWeight);
      hGenMEtDeltaRAWplusJESMEt_Pt_->Fill(RAWplusJES_MEtPt - theEventMET.genMET()->pt(), evtWeight);
      hGenMEtDeltaRAWplusJESMEt_Px_->Fill(RAWplusJES_MEtPx - theEventMET.genMET()->px(), evtWeight);
      hGenMEtDeltaRAWplusJESMEt_Py_->Fill(RAWplusJES_MEtPy - theEventMET.genMET()->py(), evtWeight);
      hGenMEtDeltaRAWMEt_Pt_->Fill(RAW_MEtPt - theEventMET.genMET()->pt(), evtWeight);
      hGenMEtDeltaRAWMEt_Px_->Fill(RAW_MEtPx - theEventMET.genMET()->px(), evtWeight);
      hGenMEtDeltaRAWMEt_Py_->Fill(RAW_MEtPy - theEventMET.genMET()->py(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Phi_->Fill(RAWplusJESplusMUONplusTAU_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hGenMEtDeltaRAWplusJESplusMUONMEt_Phi_->Fill(RAWplusJESplusMUON_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hGenMEtDeltaRAWplusJESMEt_Phi_->Fill(RAWplusJES_MEtPhi - theEventMET.genMET()->phi(), evtWeight);
      hGenMEtDeltaRAWMEt_Phi_->Fill(RAW_MEtPhi - theEventMET.genMET()->phi(), evtWeight);

      hGenMEt_Pt_->Fill(theEventMET.genMET()->pt(), evtWeight);
      hGenMEt_Phi_->Fill(theEventMET.genMET()->phi(), evtWeight);
    }
  } else {
    edm::LogError ("MEtHistManager::fillHistograms") << " Exactly one MET object expected per event --> skipping !!";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, MEtHistManager, "MEtHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, MEtHistManager, "MEtHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<MEtHistManager> MEtAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(MEtAnalyzer);

