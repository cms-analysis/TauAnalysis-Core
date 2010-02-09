#include "TauAnalysis/Core/interface/SysUncertaintyBinning.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/Core/interface/SysUncertaintyService.h"

#include "TauAnalysis/Core/interface/binningAuxFunctions.h"
#include "TauAnalysis/Core/interface/sysUncertaintyAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
 
const std::string meOptionsBinContent = std::string(meOptionsSeparator).append("a1").append(meOptionsSeparator).append("s1");
const std::string meOptionsBinError = std::string(meOptionsSeparator).append("a2").append(meOptionsSeparator).append("s1");

SysUncertaintyBinning::SysUncertaintyBinning()
{
  numBins_ = 0;
}

SysUncertaintyBinning::SysUncertaintyBinning(const edm::ParameterSet& cfg) 
  : BinningBase(cfg)
{
  //std::cout << "<SysUncertaintyBinning::SysUncertaintyBinning>:" << std::endl;

  vstring cfgSystematics = cfg.getParameter<vstring>("systematics");
  for ( vstring::const_iterator sysName = cfgSystematics.begin();
	sysName != cfgSystematics.end(); ++sysName ) {
    vstring sysNames_expanded = expandSysName(*sysName);
    systematics_.insert(systematics_.end(), sysNames_expanded.begin(), sysNames_expanded.end());
  }

  systematics_.insert(systematics_.begin(), SysUncertaintyService::getNameCentralValue());
  //std::cout << " systematics = " << format_vstring(systematics_) << std::endl;

  numBins_ = binGrid_->numBins();
  //std::cout << "numBins = " << numBins_ << std::endl;

  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    binEntryMapType binEntryMap;

    for ( vstring::const_iterator sysName = systematics_.begin();
	  sysName != systematics_.end(); ++sysName ) {
      binEntryMap.insert(std::pair<std::string, binEntryType>(*sysName, binEntryType()));
    }

    binEntries_.push_back(binEntryMap);
  }
}

SysUncertaintyBinning::~SysUncertaintyBinning()
{
//--- nothing to be done yet...
}

void SysUncertaintyBinning::bin(const std::vector<double>& x, double weight)
{
 if ( !edm::Service<SysUncertaintyService>().isAvailable() ) {
    edm::LogError ("bin") << " Failed to access SysUncertaintyService --> binning results will NOT be filled !!";
    return;
  }

  const SysUncertaintyService* sysUncertaintyService = &(*edm::Service<SysUncertaintyService>());

  const std::string& currentSystematic = sysUncertaintyService->getCurrentSystematic();

  int iBin = binGrid_->binNumber(x);
  if ( iBin >= 0 && iBin < (int)numBins_ ) {

//--- fill binning results corresponding to current systematic;
//    print error message if no histogram managers defined for current systematic
    std::map<std::string, binEntryType>::iterator binEntry = binEntries_[iBin].find(currentSystematic);
    if ( binEntry != binEntries_[iBin].end() ) {
      binEntry->second.binContent_ += weight;
      binEntry->second.binSumw2_ += weight*weight;
    } else {
      edm::LogError ("bin") << " No binning results defined for systematic = " << currentSystematic << " !!";
    }
  }
}

void SysUncertaintyBinning::print(std::ostream& stream) const
{
  stream << "<SysUncertaintyBinning::print>:" << std::endl;
  stream << " name = " << name_ << std::endl;
  
  if ( numBins_ >= 1 ) {
    const std::vector<std::string>& objVarNames = binGrid_->objVarNames();
    
    for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
      stream << " bin " << std::setw(2) << iBin << " (center: ";
      
      std::vector<double> binCenter = binGrid_->binCenter(iBin);
      if ( binCenter.size() != objVarNames.size() ) {
	edm::LogError ("SysUncertaintyBinning::print") << "Invalid dimension of bin-center vector !!";
	return;
      }
      
      unsigned numObjVarNames = objVarNames.size();
      for ( unsigned iObjVar = 0; iObjVar < numObjVarNames; ++iObjVar ) {
	stream << objVarNames[iObjVar] << " = " << std::setprecision(3) << std::fixed << binCenter[iObjVar];
	if ( iObjVar < (numObjVarNames - 1) ) stream << ", ";
      }
      
      stream << "): " << std::endl;
      
      printBinEntries(stream, binEntries_[iBin]);
    }
  }

  binEntryMapType binEntries_sum;  
  
  for ( vstring::const_iterator sysName = systematics_.begin();
	sysName != systematics_.end(); ++sysName ) {
    binEntryType& binEntrySum = binEntries_sum[*sysName];
    
    for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
      std::map<std::string, binEntryType>::const_iterator binEntry = binEntries_[iBin].find(*sysName);
      if ( binEntry == binEntries_[iBin].end() ) {
	edm::LogError ("print") << " No binning results defined for systematic = " << (*sysName) << " !!";
	return;
      }
      
      binEntrySum.binContent_ += binEntry->second.binContent_;
      binEntrySum.binSumw2_ += binEntry->second.binSumw2_;
    }
  }
  
  stream << " sum:" << std::endl;
  
  printBinEntries(stream, binEntries_sum);
}

bool containsSysName(const std::vector<std::string>& sysNames_skip, const std::string& sysName)
{
  for ( std::vector<std::string>::const_iterator sysName_skip = sysNames_skip.begin(); 
	sysName_skip != sysNames_skip.end(); ++sysName_skip ) {
    if ( sysName == (*sysName_skip) ) return true;
  }

  return false;
}

void SysUncertaintyBinning::printBinEntries(std::ostream& stream, const binEntryMapType& binEntries) const
{
//--- get central value (no systematic shifts/reweights applied)
//    as reference with respect to which effect of systematic shifts/reweights get computed
  std::map<std::string, binEntryType>::const_iterator binEntry_central = binEntries.find(SysUncertaintyService::getNameCentralValue());
  if ( binEntry_central == binEntries.end() ) {
    edm::LogError ("printBinEntries") << " No binning results defined for central value !!";
    return;
  }

  double binCentralValue = binEntry_central->second.binContent_;

  stream << "  central value:"
	 << " " << std::setprecision(3) << std::fixed << binCentralValue
	 << " +/- " << std::setprecision(3) << std::fixed << TMath::Sqrt(binEntry_central->second.binSumw2_) << " (stat.)" << std::endl;
  stream << "  systematic uncertainties:" << std::endl;
  
  TPRegexp regexpParser_bidirectional_entry("[[:alnum:]]+(Up|Down)");
  TPRegexp regexpParser_bidirectional_name("([[:alnum:]]+)(Up|Down)");
  TPRegexp regexpParser_array_entry("[[:alnum:]]+\\([[:digit:]]+\\)");
  TPRegexp regexpParser_array_name("([[:alnum:]]+)\\([[:digit:]]+\\)");
  
  vstring sysNames_skip;

  for ( vstring::const_iterator sysName = systematics_.begin();
	sysName != systematics_.end(); ++sysName ) {

    if ( (*sysName) == SysUncertaintyService::getNameCentralValue() ) continue;

    TString sysName_tstring = sysName->data();
    
    bool parseError = false;
    
    if ( regexpParser_bidirectional_entry.Match(sysName_tstring) >= 1 ) {
      TObjArray* subStrings = regexpParser_bidirectional_name.MatchS(sysName_tstring);
      
      if ( subStrings->GetEntries() == 3 ) {
	std::string sysName_bidirectional = ((TObjString*)subStrings->At(1))->GetString().Data();

	if ( containsSysName(sysNames_skip, sysName_bidirectional) ) continue;

	std::string sysName_up = std::string(sysName_bidirectional).append("Up");
	std::map<std::string, binEntryType>::const_iterator binEntry_up = binEntries.find(sysName_up);
	if ( binEntry_up == binEntries.end() ) {
	  edm::LogError ("printBinEntries") << " No binning results defined for systematic = " << sysName_up << " !!";
	  return;
	}
	double binShiftedValue_up = binEntry_up->second.binContent_;
	double sysShift_up = ( binCentralValue ) ? (binShiftedValue_up - binCentralValue)/binCentralValue : 0;
	
	std::string sysName_down = std::string(sysName_bidirectional).append("Down");
	std::map<std::string, binEntryType>::const_iterator binEntry_down = binEntries.find(sysName_down);
	if ( binEntry_down == binEntries.end() ) {
	  edm::LogError ("printBinEntries") << " No binning results defined for systematic = " << sysName_down << " !!";
	  continue;
	}
	double binShiftedValue_down = binEntry_down->second.binContent_;
	double sysShift_down = ( binCentralValue ) ? (binShiftedValue_down - binCentralValue)/binCentralValue : 0;
	
	stream << " " << std::setw(20) << sysName_bidirectional << ":"
	       << " up = " << std::setprecision(3) << std::fixed << sysShift_up*100. << "%,"
	       << " down = " << std::setprecision(3) << std::fixed << sysShift_down*100. << "%" << std::endl;
	
	sysNames_skip.push_back(sysName_bidirectional);
      } else {
	parseError = true;
      }
    } else if ( regexpParser_array_entry.Match(sysName_tstring) == 1 ){
      TObjArray* subStrings = regexpParser_array_name.MatchS(sysName_tstring);
      
      if ( subStrings->GetEntries() == 2 ) {
	std::string sysName_array = ((TObjString*)subStrings->At(1))->GetString().Data();
	
	if ( containsSysName(sysNames_skip, sysName_array) ) continue;
	
	double sysShiftsSum2 = 0.;
	
	TPRegexp regexpParser_element_entry(std::string(sysName_array).append("\\([[:digit:]]+\\)").data());
	
	for ( vstring::const_iterator sysName_element = sysName;
	      sysName_element != systematics_.end(); ++sysName_element ) {	  
	  if ( regexpParser_element_entry.Match(sysName_tstring) == 1 ) {
	    std::map<std::string, binEntryType>::const_iterator binEntry_element = binEntries.find(*sysName_element);
	    if ( binEntry_element == binEntries.end() ) {
	      edm::LogError ("printBinEntries") << " No binning results defined for systematic = " << (*sysName_element) << " !!";
	      continue;
	    }
	    double binShiftedValue_element = binEntry_element->second.binContent_;
	    double sysShift_element = ( binCentralValue ) ? (binShiftedValue_element - binCentralValue)/binCentralValue : 0;
	    
	    sysShiftsSum2 += (sysShift_element*sysShift_element);
	  }
	}
	
	stream << " " << std::setw(20) << sysName_array << ":"
	       << " " << std::setprecision(3) << std::fixed << TMath::Sqrt(sysShiftsSum2)*100. << "%" << std::endl;
	
	sysNames_skip.push_back(sysName_array);
      } else {
	parseError = true;
      }
    } else {
      std::map<std::string, binEntryType>::const_iterator binEntry = binEntries.find(*sysName);
      if ( binEntry == binEntries.end() ) {
	edm::LogError ("printBinEntries") << " No binning results defined for systematic = " << (*sysName) << " !!";
	continue;
      }
      double binShiftedValue = binEntry->second.binContent_;
      double sysShift = ( binCentralValue ) ? (binShiftedValue - binCentralValue)/binCentralValue : 0;
      
      stream << " " << std::setw(20) << (*sysName) << ":"
	     << " " << std::setprecision(3) << std::fixed << sysShift*100. << "%" << std::endl;
    }
    
    if ( parseError ) {
      edm::LogError ("printBinEntries") << " Failed to decode sysName = " << (*sysName) << " !!";
    }
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

std::vector<std::string> SysUncertaintyBinning::encodeStringRep() const
{
  std::vector<std::string> buffer = BinningBase::encodeStringRep();

  std::string systematics_string = encodeVStringStringRep(systematics_);
  std::string entry_systematics = encodeBinningStringRep("systematics", "string", systematics_string);

  for ( vstring::const_iterator sysName = systematics_.begin();
	sysName != systematics_.end(); ++sysName ) {
    for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
      std::map<std::string, binEntryType>::const_iterator binEntry = binEntries_[iBin].find(*sysName);
      if ( binEntry == binEntries_[iBin].end() ) {
	edm::LogError ("encodeStringRep") << " No binning results defined for systematic = " << (*sysName) << " !!";
	continue;
      }

      std::ostringstream meName_binContent;
      meName_binContent << "binContent_region" << (iBin + 1) << "_" << (*sysName) << meOptionsBinContent;
      std::ostringstream meValue_binContent;
      meValue_binContent << std::setprecision(3) << std::fixed << binEntry->second.binContent_;
      std::string entry_binContent = encodeBinningStringRep(meName_binContent.str(), "float", meValue_binContent.str());
      buffer.push_back(entry_binContent);
      
      std::ostringstream meName_binError;
      meName_binError << "binError_region" << (iBin + 1) << "_" << (*sysName) << meOptionsBinError;
      std::ostringstream meValue_binError;
      meValue_binError << std::setprecision(3) << std::fixed << TMath::Sqrt(binEntry->second.binSumw2_);
      std::string entry_binError = encodeBinningStringRep(meName_binError.str(), "float", meValue_binError.str());
      buffer.push_back(entry_binError);
    }
  }

  return buffer;
}

void SysUncertaintyBinning::decodeStringRep(std::vector<std::string>& buffer)
{
  BinningBase::decodeStringRep(buffer);

  TPRegexp regexpParser_numBins("numBins");
  TPRegexp regexpParser_systematics("systematics");
  TPRegexp regexpParser_binContents_entry(std::string("binContent_region[[:digit:]]+_[[:alnum:]]+").append(meOptionsBinContent));
  TPRegexp regexpParser_binContents_elements(std::string("binContent_region([[:digit:]]+)_([[:alnum:]]+)").append(meOptionsBinContent));
  TPRegexp regexpParser_binError_entry(std::string("binError_region[[:digit:]]+_[[:alnum:]]+").append(meOptionsBinError));
  TPRegexp regexpParser_binError_elements(std::string("binError_region([[:digit:]]+)_([[:alnum:]]+)").append(meOptionsBinError));

  bool numBins_initialized = false;
  std::vector<std::map<std::string, bool> > binContents_initialized;
  std::vector<std::map<std::string, bool> > binSumw2_initialized;
  
  for ( std::vector<std::string>::const_iterator entry = buffer.begin();
	entry != buffer.end(); ++entry ) {
    std::string meName, meType, meValue;
    int error = 0;
    decodeBinningStringRep(*entry, meName, meType, meValue, error);

    if ( error ) {
      edm::LogError ("SysUncertaintyBinning::decodeStringRep") << " Error in parsing string = " << (*entry) << " --> skipping !!";
      continue;
    }

    TString meName_tstring = meName.data();
    
    bool binNumber_error = false;

    if ( regexpParser_numBins.Match(meName_tstring) == 1 ) {
      numBins_ = (unsigned)atoi(meValue.data());
      binEntries_.resize(numBins_);
      binContents_initialized.resize(numBins_);
      binSumw2_initialized.resize(numBins_);
      numBins_initialized = true;
    } else if ( regexpParser_systematics.Match(meName_tstring) == 1 ) {
      int error = 0;
      std::vector<std::string> systematics = decodeVStringStringRep(meValue, error);
      
      if ( error ) {
	edm::LogError ("SysUncertaintyBinning::decodeStringRep") << " Error in parsing string = " << (*entry) << " --> skipping !!";
	continue;
      }

      systematics_ = systematics;
    } else if ( regexpParser_binContents_entry.Match(meName_tstring) == 1 ) {
      if ( !numBins_initialized ) {
	edm::LogError ("decodeStringRep") << " Need to initialize numBins before setting binContents !!";
	continue;
      }
      
      TObjArray* subStrings = regexpParser_binContents_elements.MatchS(meName_tstring);
      if ( subStrings->GetEntries() == 3 ) {
	int binNumber = (unsigned)atoi(((TObjString*)subStrings->At(1))->GetString().Data()) - 1;
	std::string sysName = ((TObjString*)subStrings->At(2))->GetString().Data();
	float binContent = atof(meValue.data());
	
	if ( binNumber >= 0 && binNumber < (int)numBins_ ) {
	  binEntryType& binEntry = binEntries_[binNumber][sysName];
	  binEntry.binContent_ = binContent;
	  binContents_initialized[binNumber][sysName] = true;
	} else {
	  edm::LogError ("decodeStringRep") << " Bin number = " << binNumber << " decoded from meName = " << meName
					    << " not within numBins = " << numBins_ << " range of binning object !!";
	  continue;
	}
      } else {
	binNumber_error = true;
      }
    } else if ( regexpParser_binError_entry.Match(meName_tstring) == 1 ) {
      if ( !numBins_initialized ) {
	edm::LogError ("SysUncertaintyBinning::decodeStringRep") << " Need to initialize numBins before setting binError !!";
	continue;
      }

      TObjArray* subStrings = regexpParser_binError_elements.MatchS(meName_tstring);
      if ( subStrings->GetEntries() == 3 ) {
	int binNumber = (unsigned)atoi(((TObjString*)subStrings->At(1))->GetString().Data()) - 1;
	std::string sysName = ((TObjString*)subStrings->At(2))->GetString().Data();
	float binError = atof(meValue.data());

	if ( binNumber >= 0 && binNumber < (int)numBins_ ) {
	  binEntryType& binEntry = binEntries_[binNumber][sysName];
	  binEntry.binSumw2_ = binError*binError;
	  binSumw2_initialized[binNumber][sysName] = true;
	} else {
	  edm::LogError ("decodeStringRep") << " Bin number = " << binNumber << " decoded from meName = " << meName
					    << " not within numBins = " << numBins_ << " range of binning object !!";
	  continue;
	}
      } else {
	binNumber_error = true;
      }
    }

    if ( binNumber_error ) {
      edm::LogError ("SysUncertaintyBinning::decodeStringRep") << " Failed to decode bin number from meName = " << meName << " !!";
      continue;
    }
  }

//--- check that all data-members of SysUncertaintyBinning object 
//    have been initialized
  if ( numBins_initialized ) {
    for ( vstring::const_iterator sysName = systematics_.begin();
	  sysName != systematics_.end(); ++sysName ) {
      for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
	if ( !binContents_initialized[iBin][*sysName] ) 
	  edm::LogError ("decodeStringRep") << " Failed to decode binContents[" << iBin << "] for sysName = " << (*sysName) << " !!";
	if ( !binSumw2_initialized[iBin][*sysName] ) 
	  edm::LogError ("decodeStringRep") << " Failed to decode binError[" << iBin << "] for sysName = " << (*sysName) << " !!";
      }
    }
  } else {
    edm::LogError ("SysUncertaintyBinning::decodeStringRep") << " Failed to decode numBins !!";
  }
}


