#ifndef TauAnalysis_Core_CompositePtrCandidateT1T2MEtSVfitHistManager_h  
#define TauAnalysis_Core_CompositePtrCandidateT1T2MEtSVfitHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "TauAnalysis/Core/interface/FakeRateJetWeightExtractor.h"
#include "TauAnalysis/RecoTools/interface/PATLeptonTrackExtractor.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <vector>
#include <string>

template<typename T1, typename T2> class SVfitHistManagerEntryTemplateSpecific;

template<typename T1, typename T2>
class CompositePtrCandidateT1T2MEtSVfitHistManager : public HistManagerBase 
{
 public:  
  explicit CompositePtrCandidateT1T2MEtSVfitHistManager(const edm::ParameterSet&);
  ~CompositePtrCandidateT1T2MEtSVfitHistManager();
  
 private:
//--- histogram booking and filling functions 
//    inherited from HistManagerBase class
  void bookHistogramsImp();
  void fillHistogramsImp(const edm::Event&, const edm::EventSetup&, double);

//--- auxiliary functions
  double getDiTauCandidateWeight(const CompositePtrCandidateT1T2MEt<T1,T2>&);

//--- configuration parameters
  edm::InputTag diTauCandidateSrc_;
  edm::InputTag genParticleSrc_;
  edm::InputTag vertexSrc_;

  bool requireGenMatch_;

  typedef std::vector<int> vint;
  vint pdgIdsElectron_;
  vint pdgIdsMuon_;
  vint pdgIdsPhoton_;
  vint pdgIdsJet_;

  typedef std::vector<std::string> vstring;
  vstring algorithmNames_;

  typedef std::vector<double> vdouble;
  std::map<std::string, vdouble> massHypotheses_; // key = algorithmName, value = vector of mass hypotheses

  std::map<std::string, vstring> polarizationHypotheses_; // key = algorithmName, value = vector of polarization hypotheses 
                                                          // to be checked for each mass hypothesis

//--- "helper" class for accessing weight values
//    associated to tau decay products
//    (efficiency/fake-rate with which the tau-jet passes the tau id. criteria)
  std::vector<FakeRateJetWeightExtractor<T1>*> diTauLeg1WeightExtractors_;
  std::vector<FakeRateJetWeightExtractor<T2>*> diTauLeg2WeightExtractors_;

//--- "helper" classes for accessing the tracks 
//    of the two tau decay products
  PATLeptonTrackExtractor<T1> trackExtractorLeg1_;
  PATLeptonTrackExtractor<T2> trackExtractorLeg2_;

  typedef std::vector<SVfitHistManagerEntryTemplateSpecific<T1,T2>*> SVfitHistManagerEntryCollection;
  SVfitHistManagerEntryCollection svFitAlgorithmHistManagers_;

  std::map<std::string, MonitorElement*> hMassLRvsRLbestLR_; // key = algorithmName
  std::map<std::string, MonitorElement*> hMassLRvsRLbestRL_; // key = algorithmName
  std::map<std::string, MonitorElement*> hMassLLvsRRbestLL_; // key = algorithmName
  std::map<std::string, MonitorElement*> hMassLLvsRRbestRR_; // key = algorithmName

  struct massHypothesisEntryType
  {
    MonitorElement* hMass_; 
    MonitorElement* hMassGenLeg2Electron_;
    MonitorElement* hMassGenLeg2Muon_;
    MonitorElement* hMassGenLeg2Photon_;
    MonitorElement* hMassGenLeg2Jet_;
    std::vector<MonitorElement*> hMassVsNumVertices_; 
    MonitorElement* hPolarizationHypothesis_;
    MonitorElement* hX1res_;
    MonitorElement* hX2res_;
  };

  std::map<std::string, std::vector<massHypothesisEntryType*> > massHypothesisEntries_; // key = algorithmName, 
                                                                                        // value = vector of massHypothesisEntries
  vdouble vertexPtThresholds_;

  MonitorElement* hDiTauCandidateWeightPosLog_;
  MonitorElement* hDiTauCandidateWeightNegLog_;
  MonitorElement* hDiTauCandidateWeightZero_;
  MonitorElement* hDiTauCandidateWeightLinear_;
};

#endif  


