#ifndef TauAnalysis_Core_CompositePtrCandidateT1T2MEtSVfitHistManager_h  
#define TauAnalysis_Core_CompositePtrCandidateT1T2MEtSVfitHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
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

  bool requireGenMatch_;

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

  struct massHypothesisEntry
  {
    massHypothesisEntry(MonitorElement* hMass, MonitorElement* hPolarizationHypothesis, MonitorElement* hX1res, MonitorElement* hX2res)
      : hMass_(hMass),
	hPolarizationHypothesis_(hPolarizationHypothesis),
	hX1res_(hX1res),
	hX2res_(hX2res)
    {}
    MonitorElement* hMass_;   
    MonitorElement* hPolarizationHypothesis_;
    MonitorElement* hX1res_;
    MonitorElement* hX2res_;
  };

  std::map<std::string, std::vector<massHypothesisEntry*> > massHypothesisEntries_; // key = algorithmName, 
                                                                                    // value = vector of massHypothesisEntries

  MonitorElement* hDiTauCandidateWeightPosLog_;
  MonitorElement* hDiTauCandidateWeightNegLog_;
  MonitorElement* hDiTauCandidateWeightZero_;
  MonitorElement* hDiTauCandidateWeightLinear_;
};

#endif  


