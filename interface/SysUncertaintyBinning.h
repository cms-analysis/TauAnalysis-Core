#ifndef TauAnalysis_Core_SysUncertaintyBinning_h  
#define TauAnalysis_Core_SysUncertaintyBinning_h

/** \class SysUncertaintyBinning
 *
 * Store number of events passing selection in different bins of SysUncertaintyBinning object.
 *
 * Class is intended to be used for estimating systematic uncertainties on Monte Carlo level, 
 * but is restricted to use reconstruction level information only
 * (cannot compute acceptance, purity, stability).
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SysUncertaintyBinning.h,v 1.1 2010/01/07 13:22:06 veelken Exp $
 *
 */

#include "TauAnalysis/Core/interface/BinningBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <string>

class SysUncertaintyBinning : public BinningBase
{
 public:
  SysUncertaintyBinning();
  SysUncertaintyBinning(const edm::ParameterSet&);
  ~SysUncertaintyBinning();

  const std::string& name() const { return name_; }

  const BinGrid* binGrid() const { return binGrid_; }

  void bin(const std::vector<double>&, double = 1.);

  void print(std::ostream&) const;

 protected:
  virtual std::vector<std::string> encodeStringRep() const;
  virtual void decodeStringRep(std::vector<std::string>&);

  typedef std::vector<std::string> vstring;
  vstring systematics_;

  struct binEntryType
  {
    double binContent_;
    double binSumw2_;
  };

  typedef std::map<std::string, binEntryType> binEntryMapType;
  std::vector<binEntryMapType> binEntries_;
  unsigned numBins_;

  void printBinEntries(std::ostream&, const binEntryMapType&) const;
};

#endif  


