// -*- C++ -*-
//
// Package:    UserCode/llvv_fwk
// Class:      BaseNTuplizer
//
/**\class BaseNTuplizer BaseNTuplizer.cc UserCode/llvv_fwk/plugins/BaseNTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  ccaputo
//         Created:  Fri, 19 Jun 2018 16:44:21 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "UserCode/llvv_fwk/interface/SummuryTuple.h"
#include "UserCode/llvv_fwk/interface/NTupleUtils.h"

// root includes

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

using namespace std;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BaseNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit BaseNTuplizer(const edm::ParameterSet&);
      ~BaseNTuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void Initialize(const edm::Event& iEvent);

      // ----------member data ---------------------------


      /// file service
      edm::Service<TFileService> fs;

      llttSummuryTree*  ntuple;
      // configure the process

      bool isMC;
      unsigned int numEvents;
      double numEventsWeight;
      double sumOfSWeights;


    private:
      edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
      edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
      edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;

      //edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrixTAG_;
      // edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      // edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      // edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

      std::string sampleType;

    protected:
      edm::Handle            < GenEventInfoProduct > genEventInfoHandle;
      edm::Handle < std::vector<PileupSummaryInfo> > puInfoH;
      edm::Handle                < LHEEventProduct > lheEPHandle;

};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BaseNTuplizer::BaseNTuplizer(const edm::ParameterSet& iConfig):
  ntuple(new llttSummuryTree( &fs->file(), "summuryllttTree" )),
  isMC(iConfig.getParameter<bool>("isMC")),
  genEventInfoToken_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  // triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  // triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  numEvents = 0;
  numEventsWeight = 0;
  sumOfSWeights = 0;
}


BaseNTuplizer::~BaseNTuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void BaseNTuplizer::Initialize(const edm::Event& iEvent){
  if (isMC){
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    iEvent.getByToken(PUInfoToken_, puInfoH);
    iEvent.getByToken(lheEventProductToken_, lheEPHandle);
  }
}

//
// member functions
//

// ------------ method called for each event  ------------
void
BaseNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 using namespace edm;

 //load all the objects we will need to access

  Initialize(iEvent);
  numEvents++;
  numEventsWeight+=genEventInfoHandle->weight();
  sumOfSWeights+= (genEventInfoHandle->weight() * genEventInfoHandle->weight());

}


// ------------ method called once each job just before starting event loop  ------------
void
BaseNTuplizer::beginJob()
{
  TH1F::SetDefaultSumw2(kTRUE);
}

// ------------ method called once each job just after ending the event loop  ------------
void
BaseNTuplizer::endJob()
{

  ntuple->fillEventInfo(numEvents,numEventsWeight,sumOfSWeights);
  // ntuple->GetTuple().totalNumberOfEvents = numEvents;
  // ntuple->GetTuple().sumOfGenWeights  = numEventsWeight;
  // ntuple->GetTuple().sumOfGenSquaredWeights   = sumOfSWeights;
  ntuple->Fill();
  //fs->cd();
  ntuple->Write();
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BaseNTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BaseNTuplizer);
