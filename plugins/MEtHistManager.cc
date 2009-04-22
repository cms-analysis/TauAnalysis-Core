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

    hRAWplusJESplusMUONplusTAU_MEtPt_ = dqmStore.book1D("hRAWplusJESplusMUONplusTAU_MEtPt", "hRAWplusJESplusMUONplusTAU_MEtPt", 75, 0., 150.);
    hRAWplusJESplusMUONplusTAU_MEtPhi_ = dqmStore.book1D("hRAWplusJESplusMUONplusTAU_MEtPhi", "hRAWplusJESplusMUONplusTAU_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hRAWplusJESplusMUONplusTAU_MEtPx_ = dqmStore.book1D("hRAWplusJESplusMUONplusTAU_MEtPx", "hRAWplusJESplusMUONplusTAU_MEtPx", 150, -150., 150.);
    hRAWplusJESplusMUONplusTAU_MEtPy_ = dqmStore.book1D("hRAWplusJESplusMUONplusTAU_MEtPy", "hRAWplusJESplusMUONplusTAU_MEtPy", 150, -150., 150.);

    hRAWplusJESplusMUON_MEtPt_ = dqmStore.book1D("hRAWplusJESplusMUON_MEtPt", "hRAWplusJESplusMUON_MEtPt", 75, 0., 150.);
    hRAWplusJESplusMUON_MEtPhi_ = dqmStore.book1D("hRAWplusJESplusMUON_MEtPhi", "hRAWplusJESplusMUON_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hRAWplusJESplusMUON_MEtPx_ = dqmStore.book1D("hRAWplusJESplusMUON_MEtPx", "hRAWplusJESplusMUON_MEtPx", 150, -150., 150.);
    hRAWplusJESplusMUON_MEtPy_ = dqmStore.book1D("hRAWplusJESplusMUON_MEtPy", "hRAWplusJESplusMUON_MEtPy", 150, -150., 150.);

    hRAWplusMUONplusTAU_MEtPt_ = dqmStore.book1D("hRAWplusMUONplusTAU_MEtPt", "hRAWplusMUONplusTAU_MEtPt", 75, 0., 150.);
    hRAWplusMUONplusTAU_MEtPhi_ = dqmStore.book1D("hRAWplusMUONplusTAU_MEtPhi", "hRAWplusMUONplusTAU_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hRAWplusMUONplusTAU_MEtPx_ = dqmStore.book1D("hRAWplusMUONplusTAU_MEtPx", "hRAWplusMUONplusTAU_MEtPx", 150, -150., 150.);
    hRAWplusMUONplusTAU_MEtPy_ = dqmStore.book1D("hRAWplusMUONplusTAU_MEtPy", "hRAWplusMUONplusTAU_MEtPy", 150, -150., 150.);

    hRAW_MEtPt_ = dqmStore.book1D("hRAW_MEtPt", "hRAW_MEtPt", 75, 0., 150.);
    hRAW_MEtPhi_ = dqmStore.book1D("hRAW_MEtPhi", "hRAW_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hRAW_MEtPx_ = dqmStore.book1D("hRAW_MEtPx", "hRAW_MEtPx", 150, -150., 150.);
    hRAW_MEtPy_ = dqmStore.book1D("hRAW_MEtPy", "hRAW_MEtPy", 150, -150., 150.);

    hRAWplusJES_MEtPt_ = dqmStore.book1D("hRAWplusJES_MEtPt", "hRAWplusJES_MEtPt", 75, 0., 150.);
    hRAWplusJES_MEtPhi_ = dqmStore.book1D("hRAWplusJES_MEtPhi", "hRAWplusJES_MEtPhi", 36, -TMath::Pi(), +TMath::Pi());
    hRAWplusJES_MEtPx_ = dqmStore.book1D("hRAWplusJES_MEtPx", "hRAWplusJES_MEtPx", 150, -150., 150.);
    hRAWplusJES_MEtPy_ = dqmStore.book1D("hRAWplusJES_MEtPy", "hRAWplusJES_MEtPy", 150, -150., 150.);

    hMUON_MExCorrection_ = dqmStore.book1D("hMUON_MExCorrection", "hMUON_MExCorrection", 150, -150., 150.);
    hMUON_MEyCorrection_ = dqmStore.book1D("hMUON_MEyCorrection", "hMUON_MEyCorrection", 150, -150., 150.);
    hTAU_MExCorrection_ = dqmStore.book1D("hTAU_MExCorrection", "hTAU_MExCorrection", 150, -150., 150.);
    hTAU_MEyCorrection_ = dqmStore.book1D("hTAU_MEyCorrection", "hTAU_MEyCorrection", 150, -150., 150.);
    hJES_MExCorrection_ = dqmStore.book1D("hJES_MExCorrection", "hJES_MExCorrection", 150, -150., 150.);
    hJES_MEyCorrection_ = dqmStore.book1D("hJES_MEyCorrection", "hJES_MEyCorrection", 150, -150., 150.);

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

    hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", "hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Pt", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", "hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Px", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", "hGenMEtDeltaRAWplusJESplusMUONplusTAUMEt_Py", 150, -150.0, 150.0);

    hGenMEtDeltaRAWplusJESplusMUONMEt_Pt_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONMEt_Pt", "hGenMEtDeltaRAWplusJESplusMUONMEt_Pt", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESplusMUONMEt_Px_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONMEt_Px", "hGenMEtDeltaRAWplusJESplusMUONMEt_Px", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESplusMUONMEt_Py_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESplusMUONMEt_Py", "hGenMEtDeltaRAWplusJESplusMUONMEt_Py", 150, -150.0, 150.0);

    hGenMEtDeltaRAWplusJESMEt_Pt_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESMEt_Pt", "hGenMEtDeltaRAWplusJESMEt_Pt", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESMEt_Px_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESMEt_Px", "hGenMEtDeltaRAWplusJESMEt_Px", 150, -150.0, 150.0);
    hGenMEtDeltaRAWplusJESMEt_Py_ = dqmStore.book1D("hGenMEtDeltaRAWplusJESMEt_Py", "hGenMEtDeltaRAWplusJESMEt_Py", 150, -150.0, 150.0);

    hGenMEtDeltaRAWMEt_Pt_ = dqmStore.book1D("hGenMEtDeltaRAWMEt_Pt", "hGenMEtDeltaRAWMEt_Pt", 150, -150.0, 150.0);
    hGenMEtDeltaRAWMEt_Px_ = dqmStore.book1D("hGenMEtDeltaRAWMEt_Px", "hGenMEtDeltaRAWMEt_Px", 150, -150.0, 150.0);
    hGenMEtDeltaRAWMEt_Py_ = dqmStore.book1D("hGenMEtDeltaRAWMEt_Py", "hGenMEtDeltaRAWMEt_Py", 150, -150.0, 150.0);

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

//    MUONplusTAU_MEtPhi = (MUONplusTAU_MEtPhi>=0) ? MUONplusTAU_MEtPhi : (MUONplusTAU_MEtPhi + (2.0 * TMath::Pi()));

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

