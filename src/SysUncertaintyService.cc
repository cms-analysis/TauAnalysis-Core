#include "TauAnalysis/Core/interface/SysUncertaintyService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

SysUncertaintyService* SysUncertaintyService::gSysUncertaintyService = 0;

const double defaultEvtReweight_error = -1.;

//
//-----------------------------------------------------------------------------------------------------------------------
//

SysUncertaintyService::WeightEntry::WeightEntry(const edm::InputTag& src)
  : WeightEntryBase(src)
{}

void SysUncertaintyService::WeightEntry::update(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<double> weight;
  evt.getByLabel(src_, weight);
  value_ = (*weight);
}

double SysUncertaintyService::WeightEntry::getWeight(const std::string&) const
{
  return value_;
}

SysUncertaintyService::WeightVectorEntry::WeightVectorEntry(const edm::InputTag& src, unsigned numValues)
  : WeightEntryBase(src),
    numValues_(numValues)
{
  values_.resize(numValues_);
}

void SysUncertaintyService::WeightVectorEntry::update(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<vdouble> weights;
  evt.getByLabel(src_, weights);
  values_ = (*weights);
  numValues_ = values_.size();
}

double SysUncertaintyService::WeightVectorEntry::getWeight(const std::string& sysName) const
{
  static TPRegexp regexpParser_array_entry("[[:alnum:]]+\\[[[:digit:]]+\\]");
  static TPRegexp regexpParser_array_index("[[:alnum:]]+\\[([[:digit:]]+)\\]");
  
  TString sysName_tstring = sysName.data();

  if ( regexpParser_array_entry.Match(sysName_tstring) == 1 ) {
    TObjArray* subStrings = regexpParser_array_index.MatchS(sysName_tstring);
      
    if ( subStrings->GetEntries() == 2 ) {
      unsigned index = (unsigned)atoi(((TObjString*)subStrings->At(1))->GetString().Data());
      
      if ( index < numValues_ ) {
	return values_[index];
      } else {
	edm::LogError ("getWeight") << "Invalid index = " << index << " (valid values = 0.." << (numValues_ - 1) << ") !!";
	return defaultEvtReweight_error;
      }
    } 
  }
  
//--- failure to parse sysName
  edm::LogError ("getWeight") << " Failed to decode sysName = " << sysName << " !!";
  return defaultEvtReweight_error;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

SysUncertaintyService::SysUncertaintyService(const edm::ParameterSet& cfg)
{
  std::cout << "<SysUncertaintyService::SysUncertaintyService>:" << std::endl;

//--- ensure that SysUncertaintyService is singleton
  assert(!gSysUncertaintyService);
  gSysUncertaintyService = this;

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgEvtReweightEntries = cfg.getParameter<vParameterSet>("config");
  for ( vParameterSet::const_iterator cfgEvtReweightEntry = cfgEvtReweightEntries.begin();
	cfgEvtReweightEntry != cfgEvtReweightEntries.end(); ++cfgEvtReweightEntry ) {
    std::string name = cfgEvtReweightEntry->getParameter<std::string>("name");
    std::cout << " name = " << name << std::endl;

    TPRegexp regexpParser_array_entry("[[:alnum:]]+\\[[[:digit:]]+\\]");
    TPRegexp regexpParser_array_elements("([[:alnum:]]+)\\[([[:digit:]]+)\\]");
    
    TString name_tstring = name.data();

    edm::InputTag src = cfgEvtReweightEntry->getParameter<edm::InputTag>("src");
    std::cout << " src = " << src << std::endl;

    std::string sysName = "undefined";
    WeightEntryBase* evtWeightEntry = 0;
    if ( regexpParser_array_entry.Match(name_tstring) == 1 ) {
      TObjArray* subStrings = regexpParser_array_elements.MatchS(name_tstring);
      
      if ( subStrings->GetEntries() == 3 ) {
	sysName = ((TObjString*)subStrings->At(1))->GetString().Data();
	unsigned numValues = (unsigned)atoi(((TObjString*)subStrings->At(2))->GetString().Data());
	evtWeightEntry = new WeightVectorEntry(src, numValues);
      } else {
	edm::LogError ("SysUncertaintyService") << " Failed to decode name = " << name << " Parameter !!";
	continue;
      }
    } else {
      sysName = name;
      evtWeightEntry = new WeightEntry(src);
    }

    evtReweightEntries_.insert(std::pair<std::string, WeightEntryBase*>(sysName, evtWeightEntry));
  }
}

SysUncertaintyService::~SysUncertaintyService()
{
  for ( std::map<std::string, WeightEntryBase*>::iterator it = evtReweightEntries_.begin();
	it != evtReweightEntries_.end(); ++it ) {
    delete it->second;
  }

  if ( gSysUncertaintyService == this ) gSysUncertaintyService = 0;
}

void SysUncertaintyService::update(const std::string& systematic, const edm::Event& evt, const edm::EventSetup& es)
{
  currentSystematic_ = systematic;

  for ( std::map<std::string, WeightEntryBase*>::iterator evtReweightEntry = evtReweightEntries_.begin();
	evtReweightEntry != evtReweightEntries_.end(); ++evtReweightEntry ) {
    evtReweightEntry->second->update(evt, es);
  }
}

const std::string& SysUncertaintyService::getCurrentSystematic() const
{
  return currentSystematic_;
}

double SysUncertaintyService::getWeight() const
{
  static TPRegexp regexpParser_array_entry("[[:alnum:]]+\\[[[:digit:]]+\\]");
  static TPRegexp regexpParser_array_name("([[:alnum:]]+)\\[[[:digit:]]+\\]");
  
  TString currentSystematic_tstring = currentSystematic_.data();

  std::string sysName = "undefined";
  if ( regexpParser_array_entry.Match(currentSystematic_tstring) == 1 ) {
    TObjArray* subStrings = regexpParser_array_name.MatchS(currentSystematic_tstring);
      
    if ( subStrings->GetEntries() == 2 ) {
      sysName = ((TObjString*)subStrings->At(1))->GetString().Data();
    } else {
      edm::LogError ("getReweight") << " Failed to decode name = " << currentSystematic_ << " of current systematic uncertainty !!";
      return defaultEvtReweight_error;
    }
  } else {
    sysName = currentSystematic_;
  }

  std::map<std::string, WeightEntryBase*>::const_iterator evtReweightEntry = evtReweightEntries_.find(sysName);
  if ( evtReweightEntry != evtReweightEntries_.end() ) {
    return evtReweightEntry->second->getWeight(currentSystematic_);
  } else {
//--- no reweight defined for current systematic
//    (e.g. current systematic is a shift of the muon Pt;
//     handled by repeating a cut on the shifted value, not by applying a reweight)
    return 1.;
  }
}

/*
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

DEFINE_ANOTHER_FWK_SERVICE_MAKER(SysUncertaintyService, edm::serviceregistry::ParameterSetMaker<SysUncertaintyService>);
 */
