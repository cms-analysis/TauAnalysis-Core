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
  hRAWplusJESplusMUONplusTAUMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONplusTAUMEtPtRecVsGen", "RAWplusJESplusMUONplusTAUMEtPtRecVsGen", 75, 0., 150., 75, 0., 150.);
  hRAWplusJESplusMUONplusTAUMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESplusMUONplusTAUMEtPhiCompGen", "RAWplusJESplusMUONplusTAUMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUONplusTAUMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONplusTAUMEtPhiRecVsGen", "RAWplusJESplusMUONplusTAUMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWplusJESplusMUONMEtPtCompGen_ = dqmStore.book1D("RAWplusJESplusMUONMEtPtCompGen", "RAWplusJESplusMUONMEtPtCompGen", 100, -5.0, +5.0);
  hRAWplusJESplusMUONMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONMEtPtRecVsGen", "RAWplusJESplusMUONMEtPtRecVsGen", 75, 0., 150., 75, 0., 150.);
  hRAWplusJESplusMUONMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESplusMUONMEtPhiCompGen", "RAWplusJESplusMUONMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESplusMUONMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESplusMUONMEtPhiRecVsGen", "RAWplusJESplusMUONMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWplusJESMEtPtCompGen_ = dqmStore.book1D("RAWplusJESMEtPtCompGen", "RAWplusJESMEtPtCompGen", 100, -5.0, +5.0);
  hRAWplusJESMEtPtRecVsGen_ = dqmStore.book2D("RAWplusJESMEtPtRecVsGen", "RAWplusJESMEtPtRecVsGen", 75, 0., 150., 75, 0., 150.);
  hRAWplusJESMEtPhiCompGen_ = dqmStore.book1D("RAWplusJESMEtPhiCompGen", "RAWplusJESMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWplusJESMEtPhiRecVsGen_ = dqmStore.book2D("RAWplusJESMEtPhiRecVsGen", "RAWplusJESMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hRAWMEtPtCompGen_ = dqmStore.book1D("RAWMEtPtCompGen", "RAWMEtPtCompGen", 100, -5.0, +5.0);
  hRAWMEtPtRecVsGen_ = dqmStore.book2D("RAWMEtPtRecVsGen", "RAWMEtPtRecVsGen", 75, 0., 150., 75, 0., 150.);
  hRAWMEtPhiCompGen_ = dqmStore.book1D("RAWMEtPhiCompGen", "RAWMEtPhiCompGen", 72, -TMath::Pi(), +TMath::Pi());
  hRAWMEtPhiRecVsGen_ = dqmStore.book2D("RAWMEtPhiRecVsGen", "RAWMEtPhiRecVsGen", 36, -TMath::Pi(), +TMath::Pi(), 36, -TMath::Pi(), +TMath::Pi());
  
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", "GenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWplusJESplusMUONMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Pt", "GenMEtDeltaRAWplusJESplusMUONMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Px", "GenMEtDeltaRAWplusJESplusMUONMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESplusMUONMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESplusMUONMEt_Py", "GenMEtDeltaRAWplusJESplusMUONMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWplusJESMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Pt", "GenMEtDeltaRAWplusJESMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Px", "GenMEtDeltaRAWplusJESMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWplusJESMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWplusJESMEt_Py", "GenMEtDeltaRAWplusJESMEt_Py", 150, -150.0, 150.0);
  
  hGenMEtDeltaRAWMEt_Pt_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Pt", "GenMEtDeltaRAWMEt_Pt", 150, -150.0, 150.0);
  hGenMEtDeltaRAWMEt_Px_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Px", "GenMEtDeltaRAWMEt_Px", 150, -150.0, 150.0);
  hGenMEtDeltaRAWMEt_Py_ = dqmStore.book1D("GenMEtDeltaRAWMEt_Py", "GenMEtDeltaRAWMEt_Py", 150, -150.0, 150.0);
  
  hGenMEt_Pt_ = dqmStore.book1D("GenMEt_Pt", "GenMEt_Pt", 75, 0., 150.);
  hGenMEt_Phi_ = dqmStore.book1D("GenMEt_Phi", "GenMEt_Phi", 36, -TMath::Pi(), +TMath::Pi());
}

void MEtHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es)

{  
  //std::cout << "<MEtHistManager::fillHistograms>:" << std::endl; 

  if ( dqmError_ ) {
    edm::LogError ("fillHistograms") << " Failed to access dqmStore --> histograms will NOT be filled !!";
    return;
  }

  edm::Handle<std::vector<pat::MET> > patMETs;
  evt.getByLabel(metSrc_, patMETs);
  if ( patMETs->size() == 1 ) {
    const pat::MET& theEventMET = (*patMETs->begin());

    double RAWplusJESplusMUONplusTAU_MEtPt = theEventMET.pt();
    double RAWplusJESplusMUONplusTAU_MEtPhi = theEventMET.phi();
    double RAWplusJESplusMUONplusTAU_MEtPx = RAWplusJESplusMUONplusTAU_MEtPt * cos(RAWplusJESplusMUONplusTAU_MEtPhi);
    double RAWplusJESplusMUONplusTAU_MEtPy = RAWplusJESplusMUONplusTAU_MEtPt * sin(RAWplusJESplusMUONplusTAU_MEtPhi);
    math::XYZTLorentzVector RAWplusJESplusMUONplusTAU_MetVector(RAWplusJESplusMUONplusTAU_MEtPx , RAWplusJESplusMUONplusTAU_MEtPy , 0 , 0);
    hRAWplusJESplusMUONplusTAU_MEtPt_->Fill(RAWplusJESplusMUONplusTAU_MEtPt);
    hRAWplusJESplusMUONplusTAU_MEtPhi_->Fill(RAWplusJESplusMUONplusTAU_MEtPhi);
    hRAWplusJESplusMUONplusTAU_MEtPx_->Fill(RAWplusJESplusMUONplusTAU_MEtPx);
    hRAWplusJESplusMUONplusTAU_MEtPy_->Fill(RAWplusJESplusMUONplusTAU_MEtPy);

    double RAWplusJESplusMUON_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrTAU);
    double RAWplusJESplusMUON_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrTAU);
    double RAWplusJESplusMUON_MEtPx = RAWplusJESplusMUON_MEtPt * cos(RAWplusJESplusMUON_MEtPhi);
    double RAWplusJESplusMUON_MEtPy = RAWplusJESplusMUON_MEtPt * sin(RAWplusJESplusMUON_MEtPhi);
    math::XYZTLorentzVector RAWplusJESplusMUON_MetVector(RAWplusJESplusMUON_MEtPx , RAWplusJESplusMUON_MEtPy , 0 , 0);
    hRAWplusJESplusMUON_MEtPt_->Fill(RAWplusJESplusMUON_MEtPt);
    hRAWplusJESplusMUON_MEtPhi_->Fill(RAWplusJESplusMUON_MEtPhi);
    hRAWplusJESplusMUON_MEtPx_->Fill(RAWplusJESplusMUON_MEtPx);
    hRAWplusJESplusMUON_MEtPy_->Fill(RAWplusJESplusMUON_MEtPy);

    double RAWplusMUONplusTAU_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrJES);
    double RAWplusMUONplusTAU_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrJES);
    double RAWplusMUONplusTAU_MEtPx = RAWplusMUONplusTAU_MEtPt * cos(RAWplusMUONplusTAU_MEtPhi);
    double RAWplusMUONplusTAU_MEtPy = RAWplusMUONplusTAU_MEtPt * sin(RAWplusMUONplusTAU_MEtPhi);
    math::XYZTLorentzVector RAWplusMUONplusTAU_MetVector(RAWplusMUONplusTAU_MEtPx , RAWplusMUONplusTAU_MEtPy , 0 , 0);
    hRAWplusMUONplusTAU_MEtPt_->Fill(RAWplusMUONplusTAU_MEtPt);
    hRAWplusMUONplusTAU_MEtPhi_->Fill(RAWplusMUONplusTAU_MEtPhi);
    hRAWplusMUONplusTAU_MEtPx_->Fill(RAWplusMUONplusTAU_MEtPx);
    hRAWplusMUONplusTAU_MEtPy_->Fill(RAWplusMUONplusTAU_MEtPy);

    double RAW_MEtPt = theEventMET.uncorrectedPt(pat::MET::uncorrALL);
    double RAW_MEtPhi = theEventMET.uncorrectedPhi(pat::MET::uncorrALL);
    double RAW_MEtPx = RAW_MEtPt * cos(RAW_MEtPhi);
    double RAW_MEtPy = RAW_MEtPt * sin(RAW_MEtPhi);
    math::XYZTLorentzVector RAW_MetVector(RAW_MEtPx , RAW_MEtPy , 0 , 0);
    hRAW_MEtPt_->Fill(RAW_MEtPt);
    hRAW_MEtPhi_->Fill(RAW_MEtPhi);
    hRAW_MEtPx_->Fill(RAW_MEtPx);
    hRAW_MEtPy_->Fill(RAW_MEtPy);

    math::XYZTLorentzVector deltaTAU_MetVector = RAWplusJESplusMUONplusTAU_MetVector - RAWplusJESplusMUON_MetVector;
    math::XYZTLorentzVector deltaJES_MetVector = RAWplusJESplusMUONplusTAU_MetVector - RAWplusMUONplusTAU_MetVector;
    math::XYZTLorentzVector deltaMUON_MetVector = RAWplusJESplusMUONplusTAU_MetVector - deltaJES_MetVector - deltaTAU_MetVector - RAW_MetVector;

    math::XYZTLorentzVector RAWplusJES_MetVector = RAW_MetVector + deltaJES_MetVector;
    double RAWplusJES_MEtPt = RAWplusJES_MetVector.Pt();
    double RAWplusJES_MEtPhi = RAWplusJES_MetVector.Phi();
    double RAWplusJES_MEtPx = RAWplusJES_MetVector.Px();
    double RAWplusJES_MEtPy = RAWplusJES_MetVector.Py();
    hRAWplusJES_MEtPt_->Fill(RAWplusJES_MEtPt);
    hRAWplusJES_MEtPhi_->Fill(RAWplusJES_MEtPhi);
    hRAWplusJES_MEtPx_->Fill(RAWplusJES_MEtPx);
    hRAWplusJES_MEtPy_->Fill(RAWplusJES_MEtPy);

    hMUON_MExCorrection_->Fill(deltaMUON_MetVector.Px());
    hMUON_MEyCorrection_->Fill(deltaMUON_MetVector.Py());
    hTAU_MExCorrection_->Fill(deltaTAU_MetVector.Px());
    hTAU_MEyCorrection_->Fill(deltaTAU_MetVector.Py());
    hJES_MExCorrection_->Fill(deltaJES_MetVector.Px());
    hJES_MEyCorrection_->Fill(deltaJES_MetVector.Py());

    //MUONplusTAU_MEtPhi = (MUONplusTAU_MEtPhi>=0) ? MUONplusTAU_MEtPhi : (MUONplusTAU_MEtPhi + (2.0 * TMath::Pi()));

    if ( theEventMET.genMET() ) {
      hRAWplusJESplusMUONplusTAUMEtPtCompGen_->Fill((RAWplusJESplusMUONplusTAU_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()));
      hRAWplusJESplusMUONplusTAUMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJESplusMUONplusTAU_MEtPt);
      hRAWplusJESplusMUONplusTAUMEtPhiCompGen_->Fill(RAWplusJESplusMUONplusTAU_MEtPhi - theEventMET.genMET()->phi());
      hRAWplusJESplusMUONplusTAUMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJESplusMUONplusTAU_MEtPhi);
      hRAWplusJESplusMUONMEtPtCompGen_->Fill((RAWplusJESplusMUON_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()));
      hRAWplusJESplusMUONMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJESplusMUON_MEtPt);
      hRAWplusJESplusMUONMEtPhiCompGen_->Fill(RAWplusJESplusMUON_MEtPhi - theEventMET.genMET()->phi());
      hRAWplusJESplusMUONMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJESplusMUON_MEtPhi);
      hRAWplusJESMEtPtCompGen_->Fill((RAWplusJES_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()));
      hRAWplusJESMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAWplusJES_MEtPt);
      hRAWplusJESMEtPhiCompGen_->Fill(RAWplusJES_MEtPhi - theEventMET.genMET()->phi());
      hRAWplusJESMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAWplusJES_MEtPhi);
      hRAWMEtPtCompGen_->Fill((RAW_MEtPt - theEventMET.genMET()->pt())/TMath::Sqrt(theEventMET.genMET()->pt()));
      hRAWMEtPtRecVsGen_->Fill(theEventMET.genMET()->pt(), RAW_MEtPt);
      hRAWMEtPhiCompGen_->Fill(RAW_MEtPhi - theEventMET.genMET()->phi());
      hRAWMEtPhiRecVsGen_->Fill(theEventMET.genMET()->phi(), RAW_MEtPhi);

      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt_->Fill(RAWplusJESplusMUONplusTAU_MEtPt - theEventMET.genMET()->pt());
      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px_->Fill(RAWplusJESplusMUONplusTAU_MEtPx - theEventMET.genMET()->px());
      hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py_->Fill(RAWplusJESplusMUONplusTAU_MEtPy - theEventMET.genMET()->py());
      hGenMEtDeltaRAWplusJESplusMUONMEt_Pt_->Fill(RAWplusJESplusMUON_MEtPt - theEventMET.genMET()->pt());
      hGenMEtDeltaRAWplusJESplusMUONMEt_Px_->Fill(RAWplusJESplusMUON_MEtPx - theEventMET.genMET()->px());
      hGenMEtDeltaRAWplusJESplusMUONMEt_Py_->Fill(RAWplusJESplusMUON_MEtPy - theEventMET.genMET()->py());
      hGenMEtDeltaRAWplusJESMEt_Pt_->Fill(RAWplusJES_MEtPt - theEventMET.genMET()->pt());
      hGenMEtDeltaRAWplusJESMEt_Px_->Fill(RAWplusJES_MEtPx - theEventMET.genMET()->px());
      hGenMEtDeltaRAWplusJESMEt_Py_->Fill(RAWplusJES_MEtPy - theEventMET.genMET()->py());
      hGenMEtDeltaRAWMEt_Pt_->Fill(RAW_MEtPt - theEventMET.genMET()->pt());
      hGenMEtDeltaRAWMEt_Px_->Fill(RAW_MEtPx - theEventMET.genMET()->px());
      hGenMEtDeltaRAWMEt_Py_->Fill(RAW_MEtPy - theEventMET.genMET()->py());

      hGenMEt_Pt_->Fill(theEventMET.genMET()->pt());
      hGenMEt_Phi_->Fill(theEventMET.genMET()->phi());
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

