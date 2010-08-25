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

//--- "helper" class for accessing weight values
//    associated to tau decay products
//    (efficiency/fake-rate with which the tau-jet passes the tau id. criteria)
  std::vector<FakeRateJetWeightExtractor<T1>*> diTauLeg1WeightExtractors_;
  std::vector<FakeRateJetWeightExtractor<T2>*> diTauLeg2WeightExtractors_;

//--- "helper" classes for accessing the tracks 
//    of the two tau decay products
  PATLeptonTrackExtractor<T1> trackExtractorLeg1_;
  PATLeptonTrackExtractor<T2> trackExtractorLeg2_;

//--- histograms
  MonitorElement* hX1_;
  MonitorElement* hX2_;
/*
  MonitorElement* hX1vsGenX1_;
  MonitorElement* hX1vsGenX1Profile_;
  MonitorElement* hX2vsGenX2_;
  MonitorElement* hX2vsGenX2Profile_;
 */
  MonitorElement* hMass_; 
  MonitorElement* hMass1stSolution_;
  MonitorElement* hGenLeg1RecLeg2Mass1stSolution_;
  MonitorElement* hRecLeg1GenLeg2Mass1stSolution_;
  MonitorElement* hMass2ndSolution_;
  MonitorElement* hGenLeg1RecLeg2Mass2ndSolution_;
  MonitorElement* hRecLeg1GenLeg2Mass2ndSolution_;
  MonitorElement* hMass3rdSolution_;
  MonitorElement* hGenLeg1RecLeg2Mass3rdSolution_;
  MonitorElement* hRecLeg1GenLeg2Mass3rdSolution_;
  MonitorElement* hMass4thSolution_;  
  MonitorElement* hGenLeg1RecLeg2Mass4thSolution_;
  MonitorElement* hRecLeg1GenLeg2Mass4thSolution_;
  MonitorElement* hMassAverage_; 
  MonitorElement* hMassNumSolutionsAveraged_;
  MonitorElement* hMassBestMatch_; 
  MonitorElement* hMassVsLogLikelihood_; 
  MonitorElement* hLogLikelihood_;
  MonitorElement* hDecayTimeLeg1_;
  MonitorElement* hDecayTimeLeg2_;
  MonitorElement* hSVfitStatus_; 

  MonitorElement* hDiTauCandidateWeightPosLog_;
  MonitorElement* hDiTauCandidateWeightNegLog_;
  MonitorElement* hDiTauCandidateWeightZero_;
  MonitorElement* hDiTauCandidateWeightLinear_;
};

#endif  

