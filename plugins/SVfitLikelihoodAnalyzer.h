#ifndef TauAnalysis_Core_SVfitLikelihoodAnalyzer_h
#define TauAnalysis_Core_SVfitLikelihoodAnalyzer_h

/** \class SVfitLikelihoodAnalyzer
 *
 * Auxiliary class for recomputing likelihood values
 * for different solutions of SVfit algorithm
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.9 $
 *
 * $Id: SVfitLikelihoodAnalyzer.h,v 1.9 2010/02/28 16:03:09 veelken Exp $
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"

#include "TauAnalysis/Core/interface/AnalyzerPluginBase.h"

#include <vector>
#include <string>

template<typename T1, typename T2>
class SVfitLikelihoodAnalyzer : public AnalyzerPluginBase
{
 public:
  // constructor 
  explicit SVfitLikelihoodAnalyzer(const edm::ParameterSet&);
  
  // destructor
  virtual ~SVfitLikelihoodAnalyzer();

  void beginJob() {}
  void analyze(const edm::Event&, const edm::EventSetup&, double);
  void endJob() {}
				       
 private:
//--- configuration parameters
  edm::InputTag diTauCandidateSrc_;

  typedef std::vector<std::string> vstring;
  struct svFitAlgorithmType
  {
    std::string algorithmName_;
    vstring polarizationHypotheses_;
  };
  std::vector<svFitAlgorithmType> svFitAlgorithms_;

//--- likelihood functions  
  std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> svFitLikelihoodFunctions_;
};

#endif  

