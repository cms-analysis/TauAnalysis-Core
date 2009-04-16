#include "TauAnalysis/Core/plugins/DQMDumpFilterStatisticsTable.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

// framework & common header files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//DQM services
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TPRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <map>

DQMDumpFilterStatisticsTable::DQMDumpFilterStatisticsTable(const edm::ParameterSet& cfg)
{
  //std::cout << "<DQMDumpFilterStatisticsTable::DQMDumpFilterStatisticsTable>:" << std::endl;

  cfgError_ = 0;

  filterStatisticsService_ = new FilterStatisticsService();

  dqmDirectories_ = cfg.getParameter<vstring>("dqmDirectories");
  if ( dqmDirectories_.size() == 0 ) {
    edm::LogError("DQMDumpFilterStatisticsTable") << " Configuration Parameter dqmDirectories = " << format_vstring(dqmDirectories_)
						  << " contains no Entries --> skipping !!";
    cfgError_ = 1;
  }

  columnsSummaryTable_ = ( cfg.exists("columnsSummaryTable") ) ? cfg.getParameter<vstring>("columnsSummaryTable") : vstring();
}

DQMDumpFilterStatisticsTable::~DQMDumpFilterStatisticsTable() 
{
  delete filterStatisticsService_;

  for ( std::vector<FilterStatisticsTable*>::iterator it = filterStatisticsTables_.begin();
	it != filterStatisticsTables_.end(); ++it ) {
    delete (*it);
  }
}

void DQMDumpFilterStatisticsTable::analyze(const edm::Event&, const edm::EventSetup&)
{
//--- nothing to be done yet
}

typedef std::map<int, double> row_type;
typedef std::map<int, row_type> table_type;

void printSummaryTable(std::ostream& stream, unsigned widthNameColumn, unsigned widthNumberColumns,
		       const std::string& summaryTableType, 
		       const std::vector<std::string>& columnLabels, const std::vector<std::string>& filterTitles,
		       table_type& table, size_t numFilters, size_t numProcesses)
{
  std::cout << "Summary Table for " << summaryTableType << "Selection:" << std::endl;
  std::cout << std::endl;
  for ( std::vector<std::string>::const_iterator columnLabel = columnLabels.begin();
	columnLabel != columnLabels.end(); ++columnLabel ) {
    if ( columnLabel == columnLabels.begin() ) {
      stream << std::setw(widthNameColumn) << std::left << (*columnLabel);
    } else {
      stream << " "; 
      for ( unsigned iCharacter = 0; iCharacter < (widthNumberColumns - columnLabel->length()); ++iCharacter ) {
	stream << " ";
      }
      stream << " " << std::setw(columnLabel->length()) << std::left << (*columnLabel);
    }
  }
  stream << std::endl;
  for ( unsigned iCharacter = 0; iCharacter < widthNameColumn + columnLabels.size()*(widthNumberColumns + 2); ++iCharacter ) {
    stream << "-";
  }
  stream << std::endl;
  for ( size_t iFilter = 0; iFilter < numFilters; ++iFilter ) {
    stream << std::setw(widthNameColumn) << std::left << filterTitles[iFilter];
    for ( size_t iProcess = 0; iProcess < numProcesses; ++iProcess ) {
      stream << " ";
      stream << std::setw(widthNumberColumns - 10) << std::setprecision(3) << std::right << table[iFilter][iProcess];
      for ( unsigned iCharacter = 0; iCharacter < 10; ++iCharacter ) stream << " ";
    }
  } 
  for ( unsigned iCharacter = 0; iCharacter < widthNameColumn + columnLabels.size()*(widthNumberColumns + 2); ++iCharacter ) {
    stream << "-";
  }
  stream << std::endl << std::endl;  
}

void DQMDumpFilterStatisticsTable::endJob()
{
  //std::cout << "<DQMDumpFilterStatisticsTable::endJob>:" << std::endl;

//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError ("endjob") << " Error in Configuration ParameterSet --> FilterStatisticsTables will NOT be printed-out !!";
    return;
  }

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("endJob") << " Failed to access dqmStore --> FilterStatisticsTables will NOT be printed-out !!";
    return;
  }

//--- load FilterStatisticsTables from DQM directories
  for ( vstring::const_iterator dqmDirectory = dqmDirectories_.begin();
	dqmDirectory != dqmDirectories_.end(); ++dqmDirectory ) {
    FilterStatisticsTable* filterStatisticsTable = filterStatisticsService_->loadFilterStatisticsTable(*dqmDirectory); 

    if ( filterStatisticsTable ) {
      filterStatisticsTables_.push_back(filterStatisticsTable);
    } else {
      edm::LogError ("DQMDumpFilterStatisticsTable") << " Failed to load FilterStatisticsTable"
						     << " from dqmDirectory = " << (*dqmDirectory) 
						     << " --> FilterStatisticsTables will NOT be printed-out !!";
      return;
    }
  }

//--- print FilterStatisticsTables
  for ( std::vector<FilterStatisticsTable*>::const_iterator filterStatisticsTable = filterStatisticsTables_.begin();
	filterStatisticsTable != filterStatisticsTables_.end(); ++filterStatisticsTable ) {
    (*filterStatisticsTable)->print(std::cout);
  }

//--- print filter statistics (cut-flow) summary tables
  for ( vstring::const_iterator columnSummaryTable = columnsSummaryTable_.begin();
	columnSummaryTable != columnsSummaryTable_.end(); ++columnSummaryTable ) {

    if ( filterStatisticsTables_.begin() == filterStatisticsTables_.end() ) continue;

    table_type table;

    std::vector<std::string> columnLabels;
    columnLabels.push_back(FilterStatisticsRow::columnLabels()[0]);

//--- check that number of rows in column
//    and labels of cuts match for all processes;
//    print error message, if not
    const FilterStatisticsTable* refFilterStatisticsTable = (*filterStatisticsTables_.begin());
    std::vector<std::string> refFilterTitleColumn = refFilterStatisticsTable->extractFilterTitleColumn();
    size_t numFilters = refFilterTitleColumn.size();

    size_t numProcesses = filterStatisticsTables_.size();
    for ( size_t iProcess = 0; iProcess < numProcesses; ++iProcess ) {
      FilterStatisticsTable* filterStatisticsTable = filterStatisticsTables_[iProcess];

      std::vector<std::string> filterTitleColumn = filterStatisticsTable->extractFilterTitleColumn();
      if ( filterTitleColumn.size() != numFilters ) {
	edm::LogError ("DQMDumpFilterStatisticsTable") << " Number of entries in Filter Title columns do not match"
						       << " --> FilterStatistics summary Tables will NOT be printed-out !!";
	return;
      } else {
	for ( size_t iFilter = 0; iFilter < numFilters; ++iFilter ) {
	  if ( filterTitleColumn[iFilter] != refFilterTitleColumn[iFilter] ) {
	    edm::LogError ("DQMDumpFilterStatisticsTable") << " Filter Title columns do not match"
							   << " --> FilterStatistics summary Tables will NOT be printed-out !!";
	    return;
	  }
	}
      } 

      std::vector<double> column = filterStatisticsTable->extractColumn(*columnSummaryTable, true);
      if ( column.size() != numFilters ) {
	edm::LogError ("DQMDumpFilterStatisticsTable") << " Number of entries in Title and Number columns do not match"
						       << " --> FilterStatistics summary Tables will NOT be printed-out !!";
	return;
      }

      for ( size_t iFilter = 0; iFilter < numFilters; ++iFilter ) {
	table[iFilter][iProcess] = column[iFilter];
      }

      columnLabels.push_back(filterStatisticsTable->name());
    }

    printSummaryTable(std::cout, 30, 20, *columnSummaryTable, columnLabels, refFilterTitleColumn, table, numFilters, numProcesses);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DQMDumpFilterStatisticsTable);
