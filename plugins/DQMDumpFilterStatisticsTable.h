#ifndef TauAnalysis_Core_DQMDumpFilterStatisticsTable_h
#define TauAnalysis_Core_DQMDumpFilterStatisticsTable_h

/** \class DQMDumpFilterStatisticsTable
 *  
 *  Class to print-out cut-flow information contained in FilterStatisticsTables
 *
 *  $Date: 2009/03/04 12:14:19 $
 *  $Revision: 1.1 $
 *  \author Christian Veelken, UC Davis
 */

// framework & common header files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/Core/interface/FilterStatisticsService.h"

#include <string>
#include <vector>

class DQMDumpFilterStatisticsTable : public edm::EDAnalyzer
{
 public:
  explicit DQMDumpFilterStatisticsTable(const edm::ParameterSet&);
  virtual ~DQMDumpFilterStatisticsTable();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();  

private:
  typedef std::vector<std::string> vstring;
  vstring dqmDirectories_;

  vstring columnsSummaryTable_;

  FilterStatisticsService* filterStatisticsService_;
  std::vector<FilterStatisticsTable*> filterStatisticsTables_;

  int cfgError_;
};

#endif


