#include "TauAnalysis/Core/plugins/ElecMuEventDump.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/Core/interface/eventDumpAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include <TMath.h>

ElecMuEventDump::ElecMuEventDump(const edm::ParameterSet& cfg)
  : GenericEventDump(cfg)
{
  //std::cout << "<ElecMuEventDump::ElecMuEventDump>:" << std::endl;
}

ElecMuEventDump::~ElecMuEventDump()
{
//--- nothing to be done yet...
}

void ElecMuEventDump::print(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
			    const std::vector<std::pair<std::string, bool> >& filterResults_cumulative, 
			    const std::vector<std::pair<std::string, bool> >& filterResults_individual, 
			    double eventWeight) const
{  
  //std::cout << "<ElecMuEventDump::print>:" << std::endl; 

  printEventHeaderInfo(iEvent, eventWeight);

  printEventTriggerInfo(iEvent);

  printEventSelectionInfo(filterResults_cumulative, filterResults_individual, outputStream_);

  *outputStream_ << ">>GENERATOR LEVEL INFORMATION<<" << std::endl;

  edm::Handle<edm::View<reco::GenParticle> > genParticleCollection;
  iEvent.getByLabel(genParticleSource_, genParticleCollection);
  edm::Handle<edm::View<reco::GenJet> > genTauJetCollection;
  iEvent.getByLabel(genTauJetSource_, genTauJetCollection);
  printGenParticleInfo(genParticleCollection, genTauJetCollection, outputStream_);

  *outputStream_ << ">>RECONSTRUCTION LEVEL INFORMATION<<" << std::endl;

  printElectronInfo(iEvent);
  printMuonInfo(iEvent);

  printDiTauCandidateInfo(iEvent);

  printMissingEtInfo(iEvent);

  printJetInfo(iEvent);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(EventDumpPluginFactory, ElecMuEventDump, "ElecMuEventDump");
