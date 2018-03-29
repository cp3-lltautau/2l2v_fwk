// -*- C++ -*-
//
// Package:    UserCode/llvv_fwk
// Class:      ZHTauTauAnalyzer
//
/**\class ZHTauTauAnalyzer ZHTauTauAnalyzer.cc UserCode/llvv_fwk/plugins/ZHTauTauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  mmusich
//         Created:  Fri, 10 Nov 2017 10:44:21 GMT
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

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ClassicSVfit
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/rochcor2016.h"
#include "UserCode/llvv_fwk/interface/RoccoR_Moriond17.h"
#include "UserCode/llvv_fwk/interface/muresolution_run2.h"
#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/METFilter.h"
#include "UserCode/llvv_fwk/interface/TrigUtils.h"
#include "UserCode/llvv_fwk/interface/EwkCorrections.h"
#include "UserCode/llvv_fwk/interface/ZZatNNLO.h"
#include "UserCode/llvv_fwk/interface/FRWeights.h"

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

///   #### UGLY define this bunch of stuff at global scope

//b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
//the scale factors are taken as average numbers from the pT dependent curves see:
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
//for Loose WP  //sfb is from page 7 https://indico.cern.ch/event/557018/contributions/2246312/attachments/1310986/1961665/csvSF_rwt_July18th_2016.pdf
float beff(0.827), sfb(0.980), sfbunc(0.015);
float leff(0.132), sfl(1.05), sflunc(0.12);

//double btagLoose = 0.605; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
double btagLoose = 0.5426;  //Moriond17 recommendation Loose
double btagMedium = 0.8484;


// Auxilliary enums

enum CRTypes {
  CR10,
  CR01,
  CR11,
  DEFAULT
};

enum HiggsFinalStates {
  elel=11*11,
  elmu=11*13,
  elha=11*15,
  mumu=13*13,
  muha=13*15,
  haha=15*15,
  none=999
};

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ZHTauTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZHTauTauAnalyzer(const edm::ParameterSet&);
      ~ZHTauTauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void Initialize(const edm::Event& iEvent);

      // auxilliary methods
      double getSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2);
      LorentzVector getClassicSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2);

      CRTypes checkBkgCR(std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2, double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut, reco::VertexCollection vtx);
      bool passHiggsCuts(std::vector<patUtils::GenericLepton> selLeptons,
			 int higgsCandL1, int higgsCandL2,
			 double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut, bool requireId, reco::VertexCollection vtx);
      double closestJet(const LorentzVector& obj, pat::JetCollection& selJets, int& closestJetIndex);
      float getTheFRWeight(std::vector<patUtils::GenericLepton>& selLeptons, pat::JetCollection& selJets,
			   int higgsCandL1, int higgsCandL2, FRWeights theFRWeightTool,double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut,
			   CRTypes CRType);
      std::vector<patUtils::GenericLepton> getLepVariations(  std::vector<patUtils::GenericLepton>& selLeptons, float factor);
      std::vector<patUtils::GenericLepton> getTauVariations( std::vector<patUtils::GenericLepton>& selLeptons,float factor);
      int tauDecayMode(const reco::GenParticle *genParticle);

      // ----------member data ---------------------------


      /// file service
      edm::Service<TFileService> fs;

      SmartSelectionMonitor mon;
      // configure the process

      bool isMC;
      double xsec;
      double xsecWeight;
      double PUNorm[3] = {1,1,1};
      int mctruthmode;
      TString dtag;
      TString suffix;
      //FIXME std::vector<std::string> urls;
      TString outUrl;
      TString jecDir;
      TString rocChorPath;

      bool filterOnlyEE, filterOnlyMUMU, filterOnlyEMU, filterOnlyPhoton, filterOnlyE, filterOnlyMU;

      bool isV0JetsMC,isWGmc,isZGmc,isMC_ZZ,isMC_ZZ2l2nu,isMC_WZ,isMC_WZ3lnu,isMC_QCD,isMC_GJet,is2015data,is2015MC,is2016data,is2016MC;

      //tree info
      TString dirname;

      //systematics
      bool runSystematics;
      bool runSVfit;

      std::vector<TString> varNames;
      std::vector<string> jetVarNames = {"", "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"};
      size_t nvarsToInclude;

      //ELECTROWEAK CORRECTION WEIGHTS
      std::vector<std::vector<float>> ewkTable, ZZ_NNLOTable;

      std::vector<double> dataPileupDistributionDouble;

      // corrections

      std::shared_ptr<TRandom3> rgenMuon_;
      std::shared_ptr<TRandom3> rgenEle_;

      pat::MET::METCorrectionLevel metcor;
      FactorizedJetCorrector *jesCor;
      JetCorrectionUncertainty *totalJESUnc;
      RoccoR_Moriond17  *muCorMoriond17;
      PhotonEnergyCalibratorRun2 PhotonEnCorrector;
      EnergyScaleCorrection_class eScaler;
      LeptonEfficiencySF lepEff;
      FRWeights theFRWeightTool;
      metf::MetFilter metFilter;
      edm::LumiReWeighting* LumiWeights;
      utils::cmssw::PuShifter_t PuShifters;

      BTagSFUtil btsfutil;
      BTagCalibration btagCalib;
      BTagCalibrationReader80X btagCal80X;

      // vector of cuts for optimization

      std::vector<float>     optim_Cuts_sumPt;
      std::vector<int>       optim_Cuts_taIso;
      std::vector<double>    optim_Cuts_muIso;
      std::vector<double>    optim_Cuts_elIso;
      std::vector<std::string> tauIDiso;

    private:
      edm::EDGetTokenT         < reco::VertexCollection > vtxMiniAODToken_;
      edm::EDGetTokenT                         < double > rhoToken_;
      edm::EDGetTokenT            < pat::MuonCollection >     muonsMiniAODToken_;
      edm::EDGetTokenT        < pat::ElectronCollection > electronsMiniAODToken_;
      //edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_,eleMediumIdMapToken_;
      edm::EDGetTokenT           < EcalRecHitCollection > recHitEBToken_;
      edm::EDGetTokenT           < EcalRecHitCollection > recHitEEToken_;
      edm::EDGetTokenT             < pat::TauCollection > tausMiniAODToken_;
      edm::EDGetTokenT < pat::PackedCandidateCollection > packedCandidateToken_;
      edm::EDGetTokenT             < pat::JetCollection > jetsMiniAODToken_;
      edm::EDGetTokenT             < pat::METCollection > pfMETAODToken_;
      edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;
      edm::EDGetTokenT            < edm::TriggerResults > metFilterResultsToken_;
      edm::EDGetTokenT    < reco::GenParticleCollection > genParticleToken_;
      edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
      edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
      edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;

      //edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrixTAG_;
      // edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      // edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      // edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

      std::string sampleType;

    protected:
      edm::Handle         < reco::VertexCollection > vtxHandle;
      edm::Handle                         < double > rhoHandle;
      edm::Handle            < pat::MuonCollection > muonsHandle;
      edm::Handle        < pat::ElectronCollection > electronsHandle;
      //edm::Handle<edm::ValueMap<bool> > tight_id_decisions, medium_id_decisions;
      edm::Handle           < EcalRecHitCollection > recHitCollectionEBHandle;
      edm::Handle           < EcalRecHitCollection > recHitCollectionEEHandle;
      edm::Handle             < pat::TauCollection > tausHandle;
      edm::Handle < pat::PackedCandidateCollection > pfCandidatesHandle;
      edm::Handle             < pat::JetCollection > jetsHandle;
      edm::Handle             < pat::METCollection > metsHandle;
      edm::Handle            < edm::TriggerResults > triggerResultsHandle;
      edm::Handle            < edm::TriggerResults > metFilterResultsHandle;
      /* Only for MC */
      edm::Handle    < reco::GenParticleCollection > genHandle;
      edm::Handle            < GenEventInfoProduct > genEventInfoHandle;
      edm::Handle < std::vector<PileupSummaryInfo> > puInfoH;
      edm::Handle                < LHEEventProduct > lheEPHandle;

      // edm::TriggerNames const& triggerNames;
      //
      // edm::TriggerNames const& getTriggerResults() const {return triggerNames;}

      //edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrix;
      //
      // edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
      // edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
      // const edm::Handle<edm::View<reco::Vertex> > GetVertexCollection() const {return vertices;}
      // const edm::Handle<edm::ValueMap<bool> >     GetTElectronIDValueMap() const {return tight_id_decisions;}
      // const edm::Handle<edm::ValueMap<bool> >     GetMElectronIDValueMap() const {return medium_id_decisions;}
      // const edm::Handle<edm::View<pat::MET> >     GetPFMet()			const {return pfMETs;}
      // const edm::Handle<edm::TriggerResults>      GetTriggerBits()      const {return triggerBits;}
      // const edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> GetMETCovMatrix() {return metCovMatrix;}
      // const std::string							  GetSampleType()		const {return sampleType;}
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
ZHTauTauAnalyzer::ZHTauTauAnalyzer(const edm::ParameterSet& iConfig):
  isMC(iConfig.getParameter<bool>("isMC")),
  xsec(iConfig.getParameter<double>("xsec")),
  mctruthmode(iConfig.getParameter<int>("mctruthmode")),
  dtag(iConfig.getParameter<std::string>("dtag")),
  suffix(iConfig.getParameter<std::string>("suffix")),
  //FIXME urls(iConfig.getUntrackedParameter<std::vector<std::string> >("input")),
  outUrl(iConfig.getParameter<std::string>("outfile")),
  jecDir(iConfig.getParameter<std::string>("jecDir")),
  rocChorPath(iConfig.getUntrackedParameter<std::string>("rocChorPath",std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/rcdata.2016.v3")),
  dirname(iConfig.getParameter<std::string>("dirName")),
  runSystematics(iConfig.getParameter<bool>("runSystematics")),
  runSVfit(iConfig.getParameter<bool>("runSystematics")),
  dataPileupDistributionDouble(iConfig.getParameter< std::vector<double> >("datapileup")),
  vtxMiniAODToken_(mayConsume<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
  rhoToken_(mayConsume<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonsMiniAODToken_(mayConsume<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronsMiniAODToken_(mayConsume<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  recHitEBToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEBSrc"))),
  recHitEEToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEESrc"))),
  // eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  // eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  tausMiniAODToken_(mayConsume< pat::TauCollection >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
  packedCandidateToken_(mayConsume< pat::PackedCandidateCollection >(iConfig.getParameter<edm::InputTag>("packCandSrc"))),
  jetsMiniAODToken_(mayConsume< pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  pfMETAODToken_(mayConsume<pat::METCollection>(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultSrc"))),
  metFilterResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterResultSrc"))),
//metCovMatrixTAG_(consumes<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >>(iConfig.getParameter<edm::InputTag>("metCov"))),
  genParticleToken_(mayConsume<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  genEventInfoToken_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  // triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  // triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  if(!isMC){
    if(dtag.Contains("DoubleElectron")) filterOnlyEE=true;
    if(dtag.Contains("DoubleMuon"))     filterOnlyMUMU=true;
    if(dtag.Contains("MuEG"))           filterOnlyEMU=true;
    if(dtag.Contains("SinglePhoton"))   filterOnlyPhoton=true;
    if(dtag.Contains("SingleMuon"))     filterOnlyMU=true;
    if(dtag.Contains("SingleElectron")) filterOnlyE=true;
  }

  isV0JetsMC   = false;//isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets")));  #FIXME should be reactivated as soon as we have exclusive jet samples
  isWGmc       = (isMC && dtag.Contains("WG"));
  isZGmc       = (isMC && dtag.Contains("ZG"));
  isMC_ZZ      = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  isMC_ZZ2l2nu = isMC && ( string(dtag.Data()).find("TeV_ZZ2l2nu")  != string::npos);
  isMC_WZ      = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  isMC_WZ3lnu  = isMC && ( string(dtag.Data()).find("TeV_WZ3lnu")  != string::npos);
  isMC_QCD     = (isMC && dtag.Contains("QCD"));
  isMC_GJet    = (isMC && dtag.Contains("GJet"));
  is2015data   = (!isMC && dtag.Contains("2015"));
  is2015MC     = (isMC && dtag.Contains("2015"));
  is2016data   = (!isMC && dtag.Contains("2016"));
  is2016MC     = (isMC);

  // push back here the empty string for no systematics
  varNames.push_back("");

  if(runSystematics){
    if(true){
      varNames.push_back("_scale_umetup"); varNames.push_back("_scale_umetdown");    //unclustered met
      varNames.push_back("_res_jup");      varNames.push_back("_res_jdown");    //jet energy resolution
      varNames.push_back("_scale_jup");    varNames.push_back("_scale_jdown");  //jet energy scale
      varNames.push_back("_scale_mup");    varNames.push_back("_scale_mdown");  //muon energy scale
      varNames.push_back("_scale_eup");    varNames.push_back("_scale_edown");  //electron energy scale
      varNames.push_back("_puup");         varNames.push_back("_pudown");      //pileup uncertainty
      varNames.push_back("_eff_bup");      varNames.push_back("_eff_bdown");    //btag veto
      varNames.push_back("_lepveto");                                           //3rd lepton veto
      varNames.push_back("_th_factup");    varNames.push_back("_th_factdown"); //factorization and renormalization scales
      varNames.push_back("_th_pdf");                                           //pdf
      varNames.push_back("_th_alphas");                                         //alpha_s (QCD)
    }
    if(isMC_ZZ || isMC_WZ){
      varNames.push_back("_th_ewkup"); varNames.push_back("_th_ewkdown"); //EWK+QCD corrections
    }
  }

  nvarsToInclude=varNames.size();

  // ewk correction tables
  if(isMC_ZZ2l2nu){
    ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
    ZZ_NNLOTable = ZZatNNLO::readFile_and_loadTable(dtag);
  }

  if(isMC_WZ3lnu){
    ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
  }

}


ZHTauTauAnalyzer::~ZHTauTauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void ZHTauTauAnalyzer::Initialize(const edm::Event& iEvent){
  iEvent.getByToken(vtxMiniAODToken_, vtxHandle);
  iEvent.getByToken(rhoToken_, rhoHandle);
  iEvent.getByToken(muonsMiniAODToken_, muonsHandle);
  iEvent.getByToken(electronsMiniAODToken_, electronsHandle);
  // iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
  // iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(recHitEBToken_, recHitCollectionEBHandle);
  iEvent.getByToken(recHitEEToken_, recHitCollectionEEHandle);
  iEvent.getByToken(tausMiniAODToken_, tausHandle);
  iEvent.getByToken(packedCandidateToken_, pfCandidatesHandle);
  iEvent.getByToken(jetsMiniAODToken_, jetsHandle);
  iEvent.getByToken(pfMETAODToken_, metsHandle);
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle);
  iEvent.getByToken(metFilterResultsToken_, metFilterResultsHandle);


  if (isMC){
    iEvent.getByToken(genParticleToken_, genHandle);
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
ZHTauTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 using namespace edm;

 //load all the objects we will need to access

  Initialize(iEvent);

  // vector<pat::TriggerObjectStandAlone> triggerObjects;
  // edm::Handle< vector<pat::TriggerObjectStandAlone> > triggerObjectsHandle;
  // iEvent.getByLabel("selectedPatTrigger",triggerObjectsHandle);

  reco::VertexCollection vtx;
  if(vtxHandle.isValid()){ vtx = *vtxHandle;}

  double rho = 0;
  if(rhoHandle.isValid()){ rho = *rhoHandle;}

  pat::MuonCollection muons;
  if(muonsHandle.isValid()){ muons = *muonsHandle;}

  pat::ElectronCollection electrons;
  if(electronsHandle.isValid()){ electrons = *electronsHandle;}

  EcalRecHitCollection recHitCollectionEB;
  EcalRecHitCollection recHitCollectionEE;
  if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}
  if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

  pat::TauCollection taus;
  if(tausHandle.isValid()){ taus = *tausHandle;}

  pat::PackedCandidateCollection pfCandidates;
  if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }


  pat::JetCollection jets;
  if(jetsHandle.isValid()){ jets = *jetsHandle;}

  pat::METCollection mets;
  if(metsHandle.isValid()){ mets = *metsHandle;}
  pat::MET met = mets[0];

  // pat::METCollection puppimets;
  // edm::Handle< pat::METCollection > puppimetsHandle;
  // iEvent.getByLabel("slimmedMETsPuppi",puppimetsHandle);
  // if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
  // // LorentzVector puppimet = puppimets[0].p4();

  float weight = xsecWeight;
  //float shapeWeight = 1.0;
  double puWeightUp = 1.0;
  double puWeightDown = 1.0;
  float puWeight(1.0);
  float weightNoLepSF(1.0);

  //##############################################   EVENT LOOP STARTS   ##############################################
  //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

  //Skip bad lumi
  //if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;
  // Muon SF 2016 Era dependent
  //  patUtils::CutVersion::CutSet muonSFs;
  //     if ( is2016data ) muonSFs = (dtag.Contains("2016H") || dtag.Contains("2016G")) ? patUtils::CutVersion::Moriond17Cut_GH : patUtils::CutVersion::Moriond17Cut_BCDEF;
  //     if ( !is2016data && is2016MC ) {
  //        TRandom3 *rgen = new TRandom3(0);
  //        int uniformValue = (int) rgen->Uniform(0, 100);
  //        //std::cout<<" Random:  "<<uniformValue<<std::endl;
  //        if ( uniformValue < 65 ) muonSFs = patUtils::CutVersion::Moriond17Cut_BCDEF;
  //  else muonSFs = patUtils::CutVersion::Moriond17Cut_GH;
  //     }

  reco::GenParticleCollection gen;
  GenEventInfoProduct eventInfo;
  int decayType = 0;
  int ZbosonType = -1;
  if(isMC){

    if(genHandle.isValid()){ gen = *genHandle;}

    if(genEventInfoHandle.isValid()){ eventInfo = *genEventInfoHandle;}

    //WEIGHT for NLO negative interference
    weight *= eventInfo.weight();

    //WEIGHT for Pileup
    int ngenITpu = 0;
    for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
      if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); } //getPU_NumInteractions();
    }

    if (ngenITpu == 0) return; // It prevents to fill vtxraw with -nan values

    puWeight          = 1; //FIXME LumiWeights->weight(ngenITpu) * PUNorm[0];
    // if ( puWeight == 0 ){
    //   std::cout<<"  Some Problem in with ngenITpu= " << ngenITpu << "vtx size= "<<vtx.size() <<std::endl;
    // }
    //FIXME puWeightUp  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
    //FIXME puWeightDown = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
    weight *= puWeight;

    //Z and Higgs GEN Level flavour
    if(isMC){

      for( auto& genParticle : gen){
	 if( abs( genParticle.pdgId() )  == 23 && genParticle.status() == 62){
	   const reco::GenParticleRefVector& daughterRefs = genParticle.daughterRefVector();
	   for(auto& daughter: daughterRefs) {
	     const reco::GenParticle *lepton(daughter.get());
	     int lpdgId = abs(lepton->pdgId());
	     ZbosonType = 3;
	     if( lpdgId == 13 ) ZbosonType = 0;
	     if( lpdgId == 11 ) ZbosonType = 1;
	     if( lpdgId == 15 ) ZbosonType = 2;
	     //if( lpdgId == 13 || lpdgId == 11) cout<<"##GEN##  Z Lepton:  pt = "<<lepton->pt()<<"  eta = "<<lepton->eta()<<"  phi = "<<lepton->phi()<<endl;
	   }
	 }
	 if( abs( genParticle.pdgId() )  == 25 && genParticle.status() == 62){

	   const reco::GenParticleRefVector& daughterRefs = genParticle.daughterRefVector();
	   decayType = 1;
	   for(auto& daughter: daughterRefs) {
	     const reco::GenParticle *tau(daughter.get());
	     decayType *= tauDecayMode(tau);
	     // GeneratorTau tauGEN = (*tau);
	     // if (tauGEN.computeDecayMode(tau)==0) decayType *= 2; //electron decay
	     // if (tauGEN.computeDecayMode(tau)==1) decayType *= 1; //muon decay
	     // if (tauGEN.computeDecayMode(tau)>1) decayType *= 3;  //hadron decay
	     // if (tauGEN.computeDecayMode(tau)==8) decayType = 0;
	     //     cout<<"    - Daughter "<<daughter->pdgId()<<"    "<< tauGEN.computeDecayMode(tau) <<endl;
	   }
	 }
      }
    }

    //GEN LEVEL FILTERING
    if(isMC && (mctruthmode==15 || mctruthmode==1113)){// && (string(dtag.Data()).find("Z#rightarrow")==0 || isMC_ZZ2l2nu))
      int prodId = 1;
      for( unsigned int k=0; k<gen.size(); ++k){
	 if( gen[k].isHardProcess() && ( abs( gen[k].pdgId() ) == 11 || abs( gen[k].pdgId() ) == 13 || abs( gen[k].pdgId() )==15 ) ) prodId*=gen[k].pdgId();
      }
       if(mctruthmode==15   && abs(prodId)!=225)return; //skip not tautau
       if(mctruthmode==1113 && abs(prodId)==225)return; //skip tautau
    }
   }

  weightNoLepSF = weight;


  //apply trigger and require compatibilitiy of the event with the PD
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*triggerResultsHandle);
  edm::TriggerResultsByName tr = edm::TriggerResultsByName(triggerResultsHandle.product(),&triggerNames);
  //if(!tr.isValid()) return;   // # FIXME TO BE CHECKED

  bool mumuTrigger(true); bool muTrigger(true);	bool eeTrigger(true); bool eTrigger(true); bool emuTrigger(true);

  int metFilterValue = 0;

  bool filterbadPFMuon = true;
  bool filterbadChCandidate = true;
  bool filterbadMuonHIP = true;
  bool filterduplicateMuonHIP = true;
  std::unique_ptr<std::vector<reco::Muon*>> outbadMuon(new std::vector<reco::Muon*>());
  std::unique_ptr<std::vector<reco::Muon*>> outduplicateMuon(new std::vector<reco::Muon*>());

  if(is2016data || is2016MC){
    mumuTrigger        = utils::passTriggerPatterns(tr,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    muTrigger          = utils::passTriggerPatterns(tr,"HLT_IsoMu22_v*","HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
    eeTrigger          = utils::passTriggerPatterns(tr,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //,"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",);
    eTrigger           = utils::passTriggerPatterns(tr,"HLT_Ele27_eta2p1_WPLoose_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*");
    emuTrigger         = utils::passTriggerPatterns(tr,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");

    metFilterValue = metFilter.passMetFilterInt( iEvent, is2016data );
    // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
    filterbadChCandidate = metFilter.passBadChargedCandidateFilter(muons,pfCandidates); if (!filterbadChCandidate) {  metFilterValue=9; }
    filterbadPFMuon = metFilter.passBadPFMuonFilter(muons,pfCandidates); if (!filterbadPFMuon) { metFilterValue=8; }
    filterbadMuonHIP = metFilter.BadGlobalMuonTaggerFilter(vtx,muons,outbadMuon,false); if (!filterbadMuonHIP) { metFilterValue=10; }
    filterduplicateMuonHIP = metFilter.BadGlobalMuonTaggerFilter(vtx,muons,outduplicateMuon,true); if (!filterduplicateMuonHIP) { metFilterValue=11; }
  }

  bool passTrigger        = mumuTrigger||muTrigger||eeTrigger||eTrigger;//||emuTrigger;

  if(  mumuTrigger)mon.fillHisto("trigger", "raw", 0 , weight);
  if(    muTrigger)mon.fillHisto("trigger", "raw", 1 , weight);
  if(    eeTrigger)mon.fillHisto("trigger", "raw", 2 , weight);
  if(     eTrigger)mon.fillHisto("trigger", "raw", 3 , weight);
  if(   emuTrigger)mon.fillHisto("trigger", "raw", 4 , weight);

  if(!isMC && passTrigger){ //avoid double counting of events from different PD
    if(filterOnlyMUMU)     { passTrigger = mumuTrigger;}
    if(filterOnlyMU)       { passTrigger = muTrigger     && !mumuTrigger;}
    if(filterOnlyEE)       { passTrigger = eeTrigger     && !muTrigger  && !mumuTrigger;}
    if(filterOnlyE)        { passTrigger = eTrigger      && !eeTrigger  && !muTrigger && !mumuTrigger; }
    if(filterOnlyEMU)      { passTrigger = emuTrigger    && !eTrigger   && !eeTrigger && !muTrigger && !mumuTrigger; }
  }

  if(passTrigger){
    if(  mumuTrigger)mon.fillHisto("trigger", "cleaned", 0 , weight);
    if(    muTrigger)mon.fillHisto("trigger", "cleaned", 1 , weight);
    if(    eeTrigger)mon.fillHisto("trigger", "cleaned", 2 , weight);
    if(     eTrigger)mon.fillHisto("trigger", "cleaned", 3 , weight);
    if(   emuTrigger)mon.fillHisto("trigger", "cleaned", 4 , weight);
  }

  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
  mon.fillHisto("eventflow"           , "all", 0, weight);
  mon.fillHisto("eventflowNoWeights"  , "all", 0, 1);
  mon.fillHisto("eventflowNoLepSF"    , "all", 0, weightNoLepSF);
  mon.fillHisto("eventflow_hMC"    , "all", 0, decayType, weightNoLepSF);
  mon.fillHisto("eventflow_ZMC"    , "all", 0, ZbosonType, weightNoLepSF);

  if(!passTrigger) return;



  mon.fillHisto("eventflow"           , "all", 1, weight);
  mon.fillHisto("eventflowNoWeights"  , "all", 1, 1);
  mon.fillHisto("eventflowNoLepSF"    , "all", 1, weightNoLepSF);
  mon.fillHisto("eventflow_hMC"    , "all", 1, decayType, weightNoLepSF);
  mon.fillHisto("eventflow_ZMC"    , "all", 1, ZbosonType, weightNoLepSF);
  //##############################################   EVENT PASSED THE TRIGGER   ######################################
  if (metFilterValue==10 || metFilterValue==11) { metFilterValue=0; }
  if( metFilterValue!=0 ) return;	 //Note this must also be applied on MC
  mon.fillHisto("eventflow"           , "all", 2, weight);
  mon.fillHisto("eventflowNoWeights"  , "all", 2, 1);
  mon.fillHisto("eventflowNoLepSF"    , "all", 2, weightNoLepSF);
  mon.fillHisto("eventflow_hMC"    , "all", 2, decayType, weightNoLepSF);
  mon.fillHisto("eventflow_ZMC"    , "all", 2, ZbosonType, weightNoLepSF);

  // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
  //	  if (!filterbadPFMuon || !filterbadChCandidate) continue;
  //##############################################   EVENT PASSED MET FILTER   #######################################

  if(isV0JetsMC){
    if(lheEPHandle.isValid()){
      mon.fillHisto("nup","",lheEPHandle->hepeup().NUP,1);
      if(lheEPHandle->hepeup().NUP>5) return;
      mon.fillHisto("nupfilt","",lheEPHandle->hepeup().NUP,1);
    }else{
      printf("Handle to externalLHEProducer is invalid --> Can not ignore V0+Jet events from inclusive samples\n");
    }
  }

  //Electroweak corrections to ZZ and WZ simulations
  double ewkCorrectionsWeight = 1.;
  double ewkCorrections_error = 0.;
  if(isMC_ZZ2l2nu || isMC_WZ3lnu) ewkCorrectionsWeight = EwkCorrections::getEwkCorrections(dtag, gen, ewkTable, eventInfo, ewkCorrections_error);
  double ewkCorrections_up = (ewkCorrectionsWeight + ewkCorrections_error)/ewkCorrectionsWeight;
  double ewkCorrections_down = (ewkCorrectionsWeight - ewkCorrections_error)/ewkCorrectionsWeight;

  //final event weight
  weight *= ewkCorrectionsWeight;

  //NNLO corrections on ZZ2l2nu
  double ZZ_NNLOcorrectionsWeight =1.;
  double mzz = - 404; // will be filled by getNNLOCorrections
  if(isMC_ZZ2l2nu) ZZ_NNLOcorrectionsWeight = ZZatNNLO::getNNLOCorrections(dtag, gen, ZZ_NNLOTable, mzz);
  if(isMC_ZZ2l2nu) mon.fillHisto("mzz", "qqZZ_atNLO", mzz, weight);
  weight *= ZZ_NNLOcorrectionsWeight;
  if(isMC_ZZ2l2nu) mon.fillHisto("mzz", "qqZZ_atNNLO", mzz, weight);

  //
  //
  // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
  //
  //

  //
  // PHOTON ANALYSIS
  //
  pat::PhotonCollection selPhotons;
  //int nPho55=0; int nPho100=0;
  //for(size_t ipho=0; ipho<photons.size(); ipho++){
  //   pat::Photon photon = photons[ipho];
  //   mon.fillHisto("phopt", "trg", photon.pt(), weight);
  //   mon.fillHisto("phoeta", "trg", photon.eta(), weight);

  //   //calibrate photon energy
  //   PhotonEnCorrector.calibrate(photon, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID());

  //   if(photon.pt()<55)continue;
  //   if(fabs(photon.superCluster()->eta())>1.4442 ) continue;
  //   if(!patUtils::passId(photon, rho, patUtils::llvvPhotonId::Tight)) continue;

  //   selPhotons.push_back(photon);
  //   if(photon.pt()>55)nPho55++;
  //   if(photon.pt()>100)nPho100++;
  //}

  //
  // LEPTON ANALYSIS
  //

  //start by merging electrons and muons
  std::vector<patUtils::GenericLepton> leptons;
  for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}
  for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}
  std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

  std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
  LorentzVector muDiff(0,0,0,0);
  LorentzVector elDiff(0,0,0,0);
  for(size_t ilep=0; ilep<leptons.size(); ilep++){
    bool passKin(true); //passId(true),passIso(true);
    bool passIsoWPforFakeRate(true);
    bool passVeryLooseLepton(true), passLooseLepton(true), passSoftMuon(true);
    int lid=leptons[ilep].pdgId();

    //no need for charge info any longer
    lid=abs(lid);
    TString lepStr( lid==13 ? "mu" : "e");

    //veto nearby photon (loose electrons are many times photons...)
    double minDRlg(9999.);
    for(size_t ipho=0; ipho<selPhotons.size(); ipho++){
      minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
    }
    if(minDRlg<0.1) continue;

    //veto leptons overlaping with other lep
    bool overlapWithLepton=false;
    for(int l1=0; l1<(int)selLeptons.size();++l1){
      if(deltaR(leptons[ilep].p4(), selLeptons[l1])<0.1){overlapWithLepton=true; break;}
    }if(overlapWithLepton)continue;

    //Cut based identification
    // passId
    //passId = lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
    //patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
    // passLooseLepton
    passLooseLepton &= lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
      patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
    passVeryLooseLepton &= passLooseLepton;
    // passSoftMuon
    passSoftMuon &= lid==11 ? false : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Soft, patUtils::CutVersion::CutSet::ICHEP16Cut);

    //isolation
    //  passIso
    // passIso = lid==11 ? patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
    // patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut);
    // passLooseLepton
    passLooseLepton &= lid==11 ? patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut) :
      patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Loose, patUtils::CutVersion::CutSet::Moriond17Cut);

    // passVeryLooseLepton
    passVeryLooseLepton &= lid==11 ?  patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::VeryLoose, patUtils::CutVersion::CutSet::ICHEP16Cut) :
      patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::VeryLoose, patUtils::CutVersion::CutSet::ICHEP16Cut);

    passIsoWPforFakeRate = lid==11 ?  patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut) :
      patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut);

    //apply muon corrections
    if(abs(lid)==13 && passVeryLooseLepton){
      //if(abs(lid)==13 && passIsoWPforFakeRate){
      passSoftMuon=false;
      if(is2016MC || is2016data){
	 if(muCorMoriond17){

	   muDiff -= leptons[ilep].p4();

	   double pt  = leptons[ilep].pt();
	   double eta = leptons[ilep].eta();
	   double phi = leptons[ilep].phi();
	   int charge = leptons[ilep].charge();
	   TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
	   // cout<<"PT Befor Correction: "<< p4.Pt() << endl;
	   int ntrk = leptons[ilep].mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();

	   if(is2016MC){

	     double u1 = rgenMuon_->Uniform();
	     double u2 = rgenMuon_->Uniform();
	     double mcSF = muCorMoriond17->kScaleAndSmearMC(charge, pt, eta, phi, ntrk, u1, u2, 0, 0);

	     leptons[ilep].mu.setP4(LorentzVector(p4.Px()*mcSF,p4.Py()*mcSF,p4.Pz()*mcSF,p4.E()*mcSF ) );
	     leptons[ilep] = patUtils::GenericLepton(leptons[ilep].mu);

	   }else if (is2016data){

	     double dataSF = muCorMoriond17->kScaleDT(charge, pt, eta, phi, 0, 0);

	     leptons[ilep].mu.setP4(LorentzVector(p4.Px()*dataSF,p4.Py()*dataSF,p4.Pz()*dataSF,p4.E()*dataSF ) );
	     leptons[ilep] = patUtils::GenericLepton(leptons[ilep].mu);
	   }

	   //  leptons[ilep].mu.setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
	   //  leptons[ilep] = patUtils::GenericLepton(leptons[ilep].mu); //recreate the generic lepton to be sure that the p4 is ok
	   muDiff += leptons[ilep].p4();
	 }
      }
    }// end muons correction

     //apply electron corrections
    if(abs(lid)==11  && passVeryLooseLepton){
      //if(abs(lid)==11 && passIsoWPforFakeRate){
      //std::cout<<"START ---- "<<std::endl;
       elDiff -= leptons[ilep].p4();
       //const EcalRecHitCollection* recHits = (leptons[ilep].el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product();
       unsigned int gainSeed = patUtils::GainSeed(leptons[ilep].el, (leptons[ilep].el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product() );

       if(!isMC){

         double scale_corr=eScaler.ScaleCorrection(iEvent.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed);
         //At this point, the new data energy will be:
         // E_new=E_old*(scale_corr);
         TLorentzVector p4(leptons[ilep].el.px(),leptons[ilep].el.py(),leptons[ilep].el.pz(),leptons[ilep].el.energy());
         leptons[ilep].el.setP4(LorentzVector(p4.Px()*scale_corr,p4.Py()*scale_corr,p4.Pz()*scale_corr,p4.E()*scale_corr ) );
         leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
       }
       if(isMC){
         //std::cout<<"Before  pt  ---- "<<leptons[ilep].el.p4()<<std::endl;
         double sigma=eScaler.getSmearingSigma(iEvent.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed,0,0);
         //Put the last two inputs at 0,0 for the nominal value of sigma
         //Now smear the MC energy
         //TRandom3 *rgen_ = new TRandom3(0);
         double smearValue = rgenEle_->Gaus(1, sigma) ;
         //std::cout<<"smearing  ---- "<<smearValue<<std::endl;
         TLorentzVector p4(leptons[ilep].el.px(),leptons[ilep].el.py(),leptons[ilep].el.pz(),leptons[ilep].el.energy());
         leptons[ilep].el.setP4(LorentzVector(p4.Px()*smearValue,p4.Py()*smearValue,p4.Pz()*smearValue,p4.E()*smearValue ) );
         //std::cout<<"After  pt  ---- "<<leptons[ilep].el.p4()<<std::endl;
         // std::cout<<"\n";
         leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
       }
        elDiff += leptons[ilep].p4();
    }

    // Compute relIso after corrections
    leptons[ilep].addUserFloat("relIso",  patUtils::relIso(leptons[ilep], rho) ); //compute it once for all

    mon.fillHisto( lid == 11 ? "eleiso" : "muiso" ,  "controlPlots" , leptons[ilep].userFloat("relIso"), 1);
    //kinematics
    float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
    if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
    if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
    passVeryLooseLepton &= passKin;
    passSoftMuon    &= passKin;
    if(lid==13){
      if(leptons[ilep].pt()<10) passVeryLooseLepton=false;
      if(leptons[ilep].pt()<3)  passSoftMuon=false;
      if(leptons[ilep].pt()<10) passKin=false;
    }else if(lid==11){
      if(leptons[ilep].pt()<10) passVeryLooseLepton=false;
      if(leptons[ilep].pt()<10) passKin=false;
    }
    //if(leptons[ilep].pt()<25) passKin=false;

    //if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]);
    if(passVeryLooseLepton && passKin)            selLeptons.push_back(leptons[ilep]); //we need loose lepton for FR
    if(passIsoWPforFakeRate && passKin)                 extraLeptons.push_back(leptons[ilep]);
  }

  std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
  std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);

  //update the met for lepton energy scales
  met.setP4(met.p4() - muDiff - elDiff); //note this also propagates to all MET uncertainties
  met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
  met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
  met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
  met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction


  //
  //TAU ANALYSIS
  //

  pat::TauCollection selTaus;
  int ntaus(0);
  for(size_t itau=0; itau<taus.size(); ++itau){
    pat::Tau& tau = taus[itau];
    if(tau.pt()<20. || fabs(tau.eta()) >2.3) continue;

    bool overlapWithLepton(false);
    for(int l1=0; l1<(int)selLeptons.size();++l1){
      if(deltaR(tau, selLeptons[l1])< (selLeptons[l1].pdgId() == 15 ? 0.5 : 0.3) ){overlapWithLepton=true; break;}
    }
    if(overlapWithLepton) continue;

    //	if(!tau.isPFTau()) continue; // Only PFTaus
    //	if(tau.emFraction() >=2.) continue;

    // we need to apply a very loose selection here (Lucia's suggestion)
    if(!tau.tauID("againstElectronLooseMVA6")) continue;
    if(!tau.tauID("againstMuonLoose3")) continue;
    if(!tau.tauID("decayModeFinding")) continue;

    selTaus.push_back(tau);
    selLeptons.push_back(tau);
    extraLeptons.push_back(tau);
    ntaus++;
  }
  std::sort(selTaus.begin(),     selTaus.end(),      utils::sort_CandidatesByPt);
  std::sort(selLeptons.begin(),  selLeptons.end(),   utils::sort_CandidatesByPt);
  std::sort(extraLeptons.begin(),extraLeptons.end(), utils::sort_CandidatesByPt);

  //
  //JET/MET ANALYSIS
  //

  //add scale/resolution uncertainties and propagate to the MET
  //utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC);

  //select the jets
  std::map<string, pat::JetCollection> selJetsVar;
  std::map<string, int   > njetsVar;
  std::map<string, int   > nbtagsVar;
  std::map<string, double> mindphijmetVar;
  for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){mindphijmetVar[jetVarNames[ivar]] = 9999.0;}  //initialize

  for(size_t ijet=0; ijet<jets.size(); ijet++){
    pat::Jet jet = jets[ijet]; //copy the jet, such that we can update it

    if(jet.pt()<20 || fabs(jet.eta())>4.7 ) continue;

    //mc truth for this jet
    //const reco::GenJet* genJet=jet.genJet();
    TString jetType( jet.genJet() && (jet.genJet())->pt()>0 ? "truejetsid" : "pujetsid" );

    //cross-clean with selected leptons and photons  (DISABLED AS WE NEED THOSE FOR FR STUDY)
    //double minDRlj(9999.); for(size_t ilep=0; ilep<selLeptons.size(); ilep++){if(abs(selLeptons[ilep].pdgId())>13){continue;}  minDRlj = TMath::Min( minDRlj, deltaR(jet,selLeptons[ilep]) );}  //ignore taus for the cross-cleaning
    //double minDRlg(9999.); for(size_t ipho=0; ipho<selPhotons.size(); ipho++)  minDRlg = TMath::Min( minDRlg, deltaR(jet,selPhotons[ipho]) );
    //if(minDRlj<0.4 || minDRlg<0.4) continue;

    //jet id
    bool passPFloose = patUtils::passPFJetID("Loose", jet);
    bool passLooseSimplePuId = patUtils::passPUJetID(jet); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
    //bool passLooseSimplePuId = jet.userInt("pileupJetId:fullId") & (1 << 2);
    if(jet.pt()>30){
      mon.fillHisto(jetType,"",fabs(jet.eta()),0);
      if(passPFloose)                        mon.fillHisto("jetId", jetType,fabs(jet.eta()),1);
      if(passLooseSimplePuId)                mon.fillHisto("jetId", jetType,fabs(jet.eta()),2);
      if(passPFloose && passLooseSimplePuId) mon.fillHisto("jetId", jetType,fabs(jet.eta()),3);
    }
    if(!passPFloose || !passLooseSimplePuId) continue;

    //check for btagging
    bool overlapWithTau(false);
    for(int l1=0; l1<(int)selTaus.size();++l1){
      if(deltaR(jet, selTaus[l1])< 0.5  ){overlapWithTau=true; break;}
    }

    if(!overlapWithTau && jet.pt()>30 && fabs(jet.eta())<2.5){
      bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>btagMedium);
      bool hasCSVtagUp = hasCSVtag;
      bool hasCSVtagDown = hasCSVtag;

      //update according to the SF measured by BTV
      if(isMC){
	 int flavId=jet.partonFlavour();  double eta=jet.eta();
	 btsfutil.SetSeed(iEvent.eventAuxiliary().event()*10 + ijet*10000);
	 if(abs(flavId)==5){
	   //  74X recommendation
					                 //     btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	   //  80X recommendation
	   btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCal80X.eval_auto_bounds("up", BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCal80X.eval_auto_bounds("down", BTagEntry::FLAV_B   , eta, jet.pt()), beff);
	 }else if(abs(flavId)==4){
	   //  74X recommendation
	   //     btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	   //  80X recommendation
	   btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCal80X.eval_auto_bounds("up", BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCal80X.eval_auto_bounds("down", BTagEntry::FLAV_C   , eta, jet.pt()), beff);
	 }else{
	   //  74X recommendation
	   //     btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCalL  .eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
	   //     btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
	   //  80X recommendation
	   btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_UDSG   , eta, jet.pt()), leff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCal80X.eval_auto_bounds("up", BTagEntry::FLAV_UDSG   , eta, jet.pt()), leff);
	   btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCal80X.eval_auto_bounds("down", BTagEntry::FLAV_UDSG   , eta, jet.pt()), leff);
	 }
      }

      if(hasCSVtag    )nbtagsVar[""          ]++;
      if(hasCSVtagUp  )nbtagsVar["_eff_bup"  ]++;
      if(hasCSVtagDown)nbtagsVar["_eff_bdown"]++;
    }

    for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){
      if(!isMC && ivar>0) continue;
      pat::Jet varJet = jet;
      if(ivar!=0) varJet.setP4(jet.p4());// * jet.userFloat(jetVarNames[ivar]));
      selJetsVar[jetVarNames[ivar]].push_back(varJet);

      if(varJet.pt()>30){
	 njetsVar[jetVarNames[ivar]]++;

	 float dphijmet=fabs(deltaPhi(met.corP4(metcor).phi(), varJet.phi()));
	 if(dphijmet<mindphijmetVar[jetVarNames[ivar]]) mindphijmetVar[jetVarNames[ivar]]=dphijmet;
      }
    }
  }
  //sort all jet collection by pT
  for(auto jetCollIt = selJetsVar.begin(); jetCollIt!=selJetsVar.end(); jetCollIt++){
    std::sort(jetCollIt->second.begin(), jetCollIt->second.end(), utils::sort_CandidatesByPt);
  }


  // LOOP ON SYSTEMATIC VARIATION FOR THE STATISTICAL ANALYSIS
  double initialWeight = weight;           //save weight
  //compute scale uncertainty once and for all
  std::pair<double, double> scaleUncVar = patUtils::scaleVariationCMSSW(iEvent);  //compute it only once

  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
    if(!isMC && ivar>0 ) continue; //loop on variation only for MC samples

    //start from a nominal
    float weight = initialWeight;

    //Theoretical Uncertanties: PDF, Alpha and Scale
    if(varNames[ivar]=="_th_factup")     weight *= std::max(0.9, std::min(scaleUncVar.first , 1.1));
    if(varNames[ivar]=="_th_factdown")   weight *= std::max(0.9, std::min(scaleUncVar.second, 1.1));
    if(varNames[ivar]=="_th_alphas")     weight *= patUtils::alphaVariationCMSSW(iEvent);
    if(varNames[ivar]=="_th_pdf")        weight *= patUtils::pdfVariationCMSSW(iEvent);

    //EwkCorrections variation
    if ( varNames[ivar]=="_th_ewkup")    weight *= ewkCorrections_up;
    if ( varNames[ivar]=="_th_ewkdown")  weight *= ewkCorrections_down;

    //pileup variations
    if(varNames[ivar]=="_puup")          weight *= puWeightUp;
    if(varNames[ivar]=="_pudown")        weight *= puWeightDown;

    //recompute MET with variation
    LorentzVector imet = met.corP4(metcor);
    if(varNames[ivar]=="_scale_jup")      imet = met.shiftedP4(pat::MET::METUncertainty::JetEnUp           , metcor);
    if(varNames[ivar]=="_scale_jdown")    imet = met.shiftedP4(pat::MET::METUncertainty::JetEnDown         , metcor);
    if(varNames[ivar]=="_res_jup")        imet = met.shiftedP4(pat::MET::METUncertainty::JetResUp          , metcor);
    if(varNames[ivar]=="_res_jdown")      imet = met.shiftedP4(pat::MET::METUncertainty::JetResDown        , metcor);
    if(varNames[ivar]=="_scale_umetup")   imet = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp   , metcor);
    if(varNames[ivar]=="_scale_umetdown") imet = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown , metcor);
    if(varNames[ivar]=="_scale_mup")      imet = met.shiftedP4(pat::MET::METUncertainty::MuonEnUp          , metcor);
    if(varNames[ivar]=="_scale_mdown")    imet = met.shiftedP4(pat::MET::METUncertainty::MuonEnDown        , metcor);
    if(varNames[ivar]=="_scale_eup")      imet = met.shiftedP4(pat::MET::METUncertainty::ElectronEnUp      , metcor);
    if(varNames[ivar]=="_scale_edown")    imet = met.shiftedP4(pat::MET::METUncertainty::ElectronEnDown    , metcor);

    //to be implemented
    //	       if(varNames[ivar]=="_tesup")   selLeptons=getTauVariations(selLeptons,1.03);
    //	       if(varNames[ivar]=="_tesdown") selLeptons=getTauVariations(selLeptons,0.97);


    auto& selJets      = selJetsVar[""];        if(selJetsVar    .find(varNames[ivar].Data())!=selJetsVar    .end())selJets     = selJetsVar    [varNames[ivar].Data()];
    auto& njets        = njetsVar [""];         if(njetsVar      .find(varNames[ivar].Data())!=njetsVar      .end())njets       = njetsVar      [varNames[ivar].Data()];
    auto& nbtags       = nbtagsVar[""];         if(nbtagsVar     .find(varNames[ivar].Data())!=nbtagsVar     .end())nbtags      = nbtagsVar     [varNames[ivar].Data()];
    auto& mindphijmet  = mindphijmetVar[""];    if(mindphijmetVar.find(varNames[ivar].Data())!=mindphijmetVar.end())mindphijmet = mindphijmetVar[varNames[ivar].Data()];

    //
    // ASSIGN CHANNEL
    //

    std::vector<TString> chTags;
    TString evCat;
    int dilId(1);
    int dilLep1, dilLep2;
    LorentzVector leadingLep, trailerLep, zll, zlltmp;
    //get the Z candidate
    dilLep1=-1; dilLep2=-1; dilId=-1;
    zll = LorentzVector(0.,0.,0.,0.);

    for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
      if(abs(selLeptons[l1].pdgId())==15)continue;

      double leadPtCutValue  = abs(selLeptons[l1].pdgId())==11 ? 24.0 : 18.0;
      if( selLeptons[l1].pt()< leadPtCutValue ) continue;
      // if(!( abs(selLeptons[l1].pdgId())==11 ? patUtils::passIso(selLeptons[l1].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
      //                                       patUtils::passIso(selLeptons[l1].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut)) ||
      //    !( abs(selLeptons[l1].pdgId())==11 ? patUtils::passId(selLeptons[l1].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
      //                                       patUtils::passId(selLeptons[l1].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut)) ) continue;

      for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
	 if(abs(selLeptons[l2].pdgId())==15)continue;

	 double trailPtCutValue = abs(selLeptons[l2].pdgId())==11 ? 13.0 : 10.0;
	 if( selLeptons[l2].pt() < trailPtCutValue ) continue;
	 // if(!( abs(selLeptons[l2].pdgId())==11 ? patUtils::passIso(selLeptons[l2].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
	 //                                       patUtils::passIso(selLeptons[l2].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut)) ||
	 //    !( abs(selLeptons[l2].pdgId())==11 ? patUtils::passId(selLeptons[l2].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
	 //                                       patUtils::passId(selLeptons[l2].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut)) ) continue;

	 if(abs(selLeptons[l1].pdgId())!=abs(selLeptons[l2].pdgId())) continue; 				 //SAME FLAVOUR PAIR
	 if(selLeptons[l1].pdgId()*selLeptons[l2].pdgId()>=0) continue;					 //OPPOSITE SIGN

	 zlltmp = (selLeptons[l1].p4()+selLeptons[l2].p4());
	 if( fabs(zlltmp.mass() - 91.2) < fabs(zll.mass()-91.2) ){    //BEST MASS [76.2,106.2]
	   dilLep1 = l1;
	   dilLep2 = l2;
	   zll=zlltmp;
	   leadingLep=selLeptons[l1].p4();
	   trailerLep=selLeptons[l2].p4();
	   dilId = selLeptons[l1].pdgId() * selLeptons[l2].pdgId();
	 }
      }
    }
    //get the Z candiate (end)

    //The "all" tags considered all cases, also events without a Z
    chTags.push_back("all");

    bool isDileptonCandidate = false;
    if(dilId!=-1){
      //check the channel
      if( abs(dilId)==121){ chTags.push_back("ll"); chTags.push_back("ee");   isDileptonCandidate=true; }
      if( abs(dilId)==169){ chTags.push_back("ll"); chTags.push_back("mumu"); isDileptonCandidate=true; }

      if(is2016MC) {
	 if(abs(dilId)==121){
	   weight *= lepEff.getRecoEfficiency( selLeptons[dilLep1].el.superCluster()->eta(), abs(selLeptons[dilLep1].pdgId())).first; //Reconstruction eff
	   weight *= lepEff.getRecoEfficiency( selLeptons[dilLep2].el.superCluster()->eta(), abs(selLeptons[dilLep2].pdgId())).first; //Reconstruction eff
	 }
	 else if(abs(dilId)==169){
	   weight *= lepEff.getTrackingEfficiency( selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].pdgId())).first; //Tracking eff
	   weight *= lepEff.getTrackingEfficiency( selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].pdgId())).first; //Tracking eff

	   weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ISO w.r.t ID
	   weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;
	 }
	 weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].pdgId()),  abs(selLeptons[dilLep1].pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID
	 weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].pdgId()),  abs(selLeptons[dilLep2].pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID

	 // Trigger Eff
	 if(isMC && abs(dilId)==169)weight *= lepEff.getTriggerEfficiencySF(selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), dilId,is2016MC).first;
	 if(isMC && abs(dilId)==121)weight *= lepEff.getTriggerEfficiencySF(selLeptons[dilLep1].pt(), selLeptons[dilLep1].el.superCluster()->eta(), selLeptons[dilLep2].pt(), selLeptons[dilLep2].el.superCluster()->eta(), dilId,is2016MC).first;
      }
    }


    if(!isDileptonCandidate) continue;
    /************************* EVENT HAS a Z-like candidate ***************************/

    // cout<<"  ##RECO##  Z Lepton 1:  pt = "<<selLeptons[dilLep1].pt()<<"  eta = "<<selLeptons[dilLep1].eta()<<"  phi = "<<selLeptons[dilLep1].phi()<<endl;
    // if ( selLeptons[dilLep1].genParticle() ) cout<<"    ##RECO (GEN Match)##  Z Lepton 1:  pt = "<<(selLeptons[dilLep1].genParticle())->pt()<<"  eta = "<<(selLeptons[dilLep1].genParticle())->eta()<<"  phi = "<<(selLeptons[dilLep1].genParticle())->phi()<<endl;
    // cout<<"  ##RECO##  Z Lepton 2:  pt = "<<selLeptons[dilLep2].pt()<<"  eta = "<<selLeptons[dilLep2].eta()<<"  phi = "<<selLeptons[dilLep2].phi()<<endl;
    // if ( selLeptons[dilLep2].genParticle() ) cout<<"    ##RECO (GEN Match)##  Z Lepton 2:  pt = "<<(selLeptons[dilLep2].genParticle())->pt()<<"  eta = "<<(selLeptons[dilLep2].genParticle())->eta()<<"  phi = "<<(selLeptons[dilLep2].genParticle())->phi()<<endl;
    //
    // int triggerType = abs(dilId)==121 ? 82 : 83;
    // std::cout << "\n TRIGGER OBJECTS " << std::endl;
    // for (pat::TriggerObjectStandAlone obj : *triggerObjectsHandle) { // note: not "const &" since we want to call unpackPathNames
    //     obj.unpackPathNames(names);
    //
    //     bool typeMatched = false;
    //     for (unsigned h = 0; h < obj.filterIds().size(); ++h) {
    //       typeMatched |= (obj.filterIds()[h] == triggerType) ;
    //     }
    //     if (!typeMatched) continue;
    //
    //     std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << ", pdgId "<< obj.pdgId() << std::endl;
    //     // Print trigger object collection and type
    //     std::cout << "\t   Collection: " << obj.collection() << std::endl;
    //     std::cout << "\t   Type IDs:   ";
    //     for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
    //     std::cout << std::endl;
    //     // Print associated trigger filters
    //     std::cout << "\t   Filters:    ";
    //     for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
    //     std::cout << std::endl;
    //     std::vector< std::string > pathNamesAll = obj.pathNames(false);
    //     std::vector< std::string > pathNamesLast = obj.pathNames(true);
    //     // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    //     // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    //     // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    //     std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    //     for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
    //         bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
    //         bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
    //         bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
    //         bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
    //         std::cout << "   " << pathNamesAll[h];
    //         if (isBoth) std::cout << "(L,3)";
    //         if (isL3 && !isBoth) std::cout << "(*,3)";
    //         if (isLF && !isBoth) std::cout << "(L,*)";
    //         if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    if (ivar == 0 ){
      mon.fillHisto("eventflow"           , chTags, 3, weight);
      mon.fillHisto("eventflowNoWeights"  , chTags, 3, 1);
      mon.fillHisto("eventflowNoLepSF"    , chTags, 3, weightNoLepSF);
      mon.fillHisto("eventflow_hMC"       , chTags, 3, decayType, weightNoLepSF);
      mon.fillHisto("eventflow_ZMC"       , chTags, 3, ZbosonType, weightNoLepSF);
      mon.fillHisto("zllmass","controlPlots",zll.mass(),weight);
    }


    bool passZmass = (fabs(zll.mass()-91.2)<30.0);
    bool passZpt   = (zll.pt()>20);
    bool passBJetVetoMain = (nbtags ==0);
    bool passLepVetoMain = true;

    int higgsCandL1=-1, higgsCandL2=-1;
    LorentzVector higgsCand;
    int HiggsShortId=-1, higgsCandId;
    std::vector<TString> chTagsMain=chTags;

    int NCleanedJetMain = 0;
    bool passDPhiCut    = 0;
    bool passHiggsLoose = 0;
    bool passHiggsMain  = 0;
    double higgsCand_SVFitMass = 0;
    double classiSVFit_mass = 0;
    LorentzVector higgsCand_ClassicSVFit;
    LorentzVector higgsCand_SVFit;
    LorentzVector higgsCandH;
    LorentzVector higgsCandH_SVFit;


    //LEPTON FAKE RATE ANALYSIS Z+1jets  (no systematics taken into account here)
    if(ivar==0 && passZmass && (int) extraLeptons.size()==3){  //Request exactly one Z + 1 additional lepton
      bool IdentifiedThirdLepton=false;
      double tmass=-999;
      for(int i=0   ;i<(int)extraLeptons.size() && !IdentifiedThirdLepton;i++){
	 // if((i==dilLep1) || (i==dilLep2)) continue;
	 if(deltaR(extraLeptons[i],  selLeptons[dilLep1])<0.1 || deltaR(extraLeptons[i],  selLeptons[dilLep2])<0.1)continue;
	 if(abs(extraLeptons[i].pdgId())==11||abs(extraLeptons[i].pdgId())==13||abs(extraLeptons[i].pdgId())==15){
	   tmass = TMath::Sqrt(2*extraLeptons[i].pt()*met.pt()*(1-TMath::Cos(deltaPhi(met.phi(), extraLeptons[i].phi()))));
	 }
	 if(abs(extraLeptons[i].pdgId())==11 || abs(extraLeptons[i].pdgId())==13 || abs(extraLeptons[i].pdgId())==15){
	   int closestJetIndexL1=-1; double pTL1=-1; double etaL1=-1;
	   double dRminL1 = closestJet(extraLeptons[i].p4(), selJets, closestJetIndexL1);
	   if(closestJetIndexL1>=0 && dRminL1<0.5){pTL1=selJets[closestJetIndexL1].pt(); etaL1=abs(selJets[closestJetIndexL1].eta());}
	   else{pTL1=extraLeptons[i].pt(); etaL1=abs(extraLeptons[i].eta());}

	   TString PartName = "FR_";//+chTags.at(1)+"_";
	   if     (abs(extraLeptons[i].pdgId())==11)PartName += "El";
	   else if(abs(extraLeptons[i].pdgId())==13)PartName += "Mu";
	   else if(abs(extraLeptons[i].pdgId())==15)PartName += "Ta";
	   else PartName+= abs(selLeptons[i].pdgId());

	   std::vector<TString> TagsFR;

	   if(abs(extraLeptons[i].pdgId())==11 || abs(extraLeptons[i].pdgId())==13){
	     bool passId = false;
	     if(abs(extraLeptons[i].pdgId())==11) passId = patUtils::passId(extraLeptons[i].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
	     if(abs(extraLeptons[i].pdgId())==13) passId = patUtils::passId(extraLeptons[i].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
	     float relIso = patUtils::relIso(extraLeptons[i], rho);

	     if(true                 )TagsFR.push_back(PartName);
	     if(passId && relIso<=0.1)TagsFR.push_back(PartName+("_Id_Iso01"));
	     if(passId && relIso<=0.2)TagsFR.push_back(PartName+("_Id_Iso02"));
	     if(passId && relIso<=0.3)TagsFR.push_back(PartName+("_Id_Iso03"));

	     if(passId && relIso<=0.3)IdentifiedThirdLepton=true;
	   }else{

	     bool IdL         = extraLeptons[i].tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	     bool IdM         = extraLeptons[i].tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
	     bool IdL_MVA     = extraLeptons[i].tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
	     bool IdM_MVA     = extraLeptons[i].tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
	     bool IdL_MVA_R03 = extraLeptons[i].tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
	     bool IdM_MVA_R03 = extraLeptons[i].tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");

	     if(true                 )TagsFR.push_back(PartName);
	     if(IdL                  )TagsFR.push_back(PartName+("_Id_IsoLo"));
	     if(IdM                  )TagsFR.push_back(PartName+("_Id_IsoMe"));
	     if(IdL_MVA              )TagsFR.push_back(PartName+("_Id_IsoLo_MVA"));
	     if(IdM_MVA              )TagsFR.push_back(PartName+("_Id_IsoMe_MVA"));
	     if(IdL_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoLo_MVAR03"));
	     if(IdM_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoMe_MVAR03"));
	   }

	   if(tmass<30){
	     unsigned int NTags = TagsFR.size();
	     for(unsigned int iTags=0;iTags<NTags;iTags++){
	       TagsFR.push_back(TagsFR[iTags] + TString("_TMCut"));
	     }
	   }

	   auto TagsFRJet = TagsFR;
	   auto TagsFRLep = TagsFR;

	   for(unsigned int iTags=0;iTags<TagsFR.size();iTags++){
	     TagsFRJet.push_back(TagsFR[iTags] + (etaL1<1.4                   ?TString("_B"):TString("_E")));
	     TagsFRLep.push_back(TagsFR[iTags] + (abs(extraLeptons[i].eta())<1.4?TString("_B"):TString("_E")));
	   }

	   mon.fillHisto("wrtJetPt", TagsFRJet, pTL1              , weight);
	   if(closestJetIndexL1>=0 && dRminL1<0.5) mon.fillHisto("wrtJetPt_v2", TagsFRJet, pTL1              , weight);
	   mon.fillHisto("wrtLepPt", TagsFRLep, extraLeptons[i].pt(), weight);
	 }
      }//close loop on leptons

    }//close FR study Zmass

    int higgsMuonCand = -1;
    int higgsEleCand  = -1;
    double collinearMass = 0;
    //SIGNAL ANALYSIS Z+2Leptons  (no systematics taken into account here)
    if(passZmass && passZpt && (int)selLeptons.size()>=4){  //Request at least 4 leptons
      //printf("%30s %2i --> ", "BEFORE", -1); for(int l=0   ;l<(int)selLeptons.size();l++){ printf("%i ", selLeptons[l].pdgId());}printf("\n");

      //Get the Higgs candidate
      higgsCandL1=-1;
      higgsCandL2=-1;
      higgsCand = LorentzVector(0.,0.,0.,0.);
      HiggsShortId=-1;
      higgsCandId=0;

      for(int l=0   ;l<(int)selLeptons.size();l++){
	 if(l==dilLep1 || l==dilLep2)continue;
	 if(higgsCandL1<0){higgsCandL1=l;continue;}
	 if(higgsCandL2<0){higgsCandL2=l;break;}//ordered in pT, so all done
      }

      string ChannelName = "none";   string signName = "";
      if(higgsCandL1>=0 && higgsCandL2>=0){
	 higgsCandId=selLeptons[higgsCandL1].pdgId()*selLeptons[higgsCandL2].pdgId();
	 higgsCand = LorentzVector(selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());
	 if(higgsCandId<0){signName="_OS";}else{signName="_SS";}
	 if(higgsCandId<0){HiggsShortId = 0;}else{HiggsShortId =12;}
	 if(abs(selLeptons[dilLep1].pdgId())==11){HiggsShortId += 0;}else{HiggsShortId += 6;}
	 switch(abs(higgsCandId)){
	 case HiggsFinalStates::elel :  ChannelName  = "elel";  HiggsShortId+= 0; break;
	 case HiggsFinalStates::mumu :  ChannelName  = "mumu";  HiggsShortId+= 1; break;
	 case HiggsFinalStates::elmu :  ChannelName  = "elmu";  HiggsShortId+= 2; break;
	 case HiggsFinalStates::elha :  ChannelName  = "elha";  HiggsShortId+= 3; break;
	 case HiggsFinalStates::muha :  ChannelName  = "muha";  HiggsShortId+= 4; break;
	 case HiggsFinalStates::haha :  ChannelName  = "haha";  HiggsShortId+= 5; break;
	 default:     ChannelName  = "none";  HiggsShortId =-1; break;
	 }
	 // if (abs(higgsCandId) == 195) {
	 if(selLeptons[higgsCandL2].pdgId() == 13) higgsMuonCand =  higgsCandL2;
	 if(selLeptons[higgsCandL1].pdgId() == 13) higgsMuonCand =  higgsCandL1;

	 if(selLeptons[higgsCandL2].pdgId() == 11) higgsEleCand =  higgsCandL2;
	 if(selLeptons[higgsCandL1].pdgId() == 11) higgsEleCand =  higgsCandL1;
	 // }

	 double tauCrossTau_z  = selLeptons[higgsCandL1].px() * selLeptons[higgsCandL2].py() - selLeptons[higgsCandL1].py()*selLeptons[higgsCandL2].px();
	 double metCrossTau1_z = met.px() * selLeptons[higgsCandL1].py() - met.py() * selLeptons[higgsCandL1].px();
	 double metCrossTau2_z = met.px() * selLeptons[higgsCandL2].py() - met.py() * selLeptons[higgsCandL2].px();

	 double tau1Efraction  = tauCrossTau_z / (tauCrossTau_z + metCrossTau2_z);
	 double tau2Efraction  = tauCrossTau_z / (tauCrossTau_z + metCrossTau1_z);

	 double Den = tau1Efraction * tau2Efraction;

	 if (Den > 0) collinearMass = higgsCand.mass() / ( TMath::Sqrt(tau1Efraction*tau2Efraction) );

	 // cout << " \n---------   Testing Collinear Mass  ----------- "<<endl;
	 // cout << "  Decay products pdg id:  "<< ChannelName<<"   Event id: "<<ev.eventAuxiliary().event()<< endl;
	 // cout << " ---------------------------------------------- "<< endl;
	 // cout << "     mass value = "<< higgsCand.mass()
	 //      << "\n   chi_1 = "<< tau1Efraction<<"  chi_2 = "<< tau2Efraction
	 //      << "\n Num = "<< tauCrossTau_z << " Met x Tau_1 = " << metCrossTau1_z << "Met x Tau_2 = " << metCrossTau2_z << endl;
	 // cout << " \t ------->  Collinear mass = "<< collinearMass << "  <---------" <<endl;
      }

      chTagsMain.push_back(chTagsMain[chTagsMain.size()-1] + signName + ChannelName);
      //Get the Higgs candidate (end)
      //printf("%30s %2i --> %i %i %i %i\n", (chTagsMain[chTagsMain.size()-1]).Data(), HiggsShortId, selLeptons[dilLep1].pdgId(), selLeptons[dilLep2].pdgId(), selLeptons[higgsCandL1].pdgId(), selLeptons[higgsCandL2].pdgId());

      //reweight the event to account for lept eff.
      if(isMC && higgsCandL1>=0 && abs(selLeptons[higgsCandL1].pdgId())<15){

	 int id( abs(selLeptons[higgsCandL1].pdgId()) );

	 if(id==11)weight *= lepEff.getRecoEfficiency( selLeptons[higgsCandL1].el.superCluster()->eta(), id).first; //Reconstruction eff
	 else if(id==13)weight *= lepEff.getTrackingEfficiency( selLeptons[higgsCandL1].eta(), id).first; //Tracking eff
	 weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), id,  id ==11 ? "loose"    : "loose"   ,patUtils::CutVersion::Moriond17Cut).first : 1.0; //ID
	 if(id==13){ weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), id, "looseiso_looseid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;} //ISO w.r.t ID
      }

      if(isMC && higgsCandL2>=0 && abs(selLeptons[higgsCandL2].pdgId())<15){

	 int id( abs(selLeptons[higgsCandL2].pdgId()) );

	 if(id==11)weight *= lepEff.getRecoEfficiency( selLeptons[higgsCandL2].el.superCluster()->eta(), id).first; //Reconstruction eff
	 else if(id==13)weight *= lepEff.getTrackingEfficiency( selLeptons[higgsCandL2].eta(), id).first; //Tracking eff
	 weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), id,  id ==11 ? "loose"    : "loose"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID
	 if(id==13){ weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), id, "looseiso_looseid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;} //ISO w.r.t ID
      }


      //check how many additional light jets are present
      NCleanedJetMain = 0;
      for(int j1=0;j1<(int)selJets.size();j1++){
	 if(dilLep1    !=-1 && deltaR(selJets[j1]   , selLeptons[dilLep1 ])<0.4) continue;
	 if(dilLep2    !=-1 && deltaR(selJets[j1]   , selLeptons[dilLep2 ])<0.4) continue;
	 if(higgsCandL1         !=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandL1      ])<0.4) continue;
	 if(higgsCandL2         !=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandL2      ])<0.4) continue;
	 NCleanedJetMain++;
      }

      passDPhiCut    =  (fabs(deltaPhi(zll.phi(), met.phi()))>1.5);
      passDPhiCut    = true;
      passHiggsLoose = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.3, 0.3, "decayModeFinding", 0., false, vtx);
      passHiggsMain  = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.1, 0.15, "byMediumIsolationMVArun2v1DBoldDMwLT", 0., true, vtx);

      //SVFIT MASS
      higgsCand_SVFit = higgsCand;

      //FIXME gives a lot of warning currently
      if(runSVfit && passZmass && passZpt && selLeptons.size()>=4 && passLepVetoMain && passBJetVetoMain && passDPhiCut){
        // std::cout<<"START SVFIT\n";
        // cout<<"============================================================="<<endl;
        // cout<<"    Higgs mass = "<<higgsCand.M()<<"   Higgs pt = "<<higgsCand.Pt()<<endl;
        higgsCand_SVFitMass = getSVFit(met, selLeptons, higgsCandL1, higgsCandL2);  //compute svfit mass in a smart way
        //  higgsCand_ClassicSVFit = getClassicSVFit(met, selLeptons, higgsCandL1, higgsCandL2);
	 //  higgsCand_SVFit = higgsCand_ClassicSVFit;
	 //  classiSVFit_mass = higgsCand_ClassicSVFit.mass();
        // cout<<"============================================================="<<endl;
        // std::cout<<"END SVFIT\n";
      }

      //build the higgs candH
      higgsCandH = zll + higgsCand;
      higgsCandH_SVFit = zll + higgsCand_SVFit;
    }

    // bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
    // bool passBtags(nbtags==0);
    // bool passMinDphijmet( njets==0 || mindphijmet>0.5);
    // bool removeDump(false);

    //
    // NOW FOR THE CONTROL PLOTS
    //

    if(ivar==0){//fill plots only for nominal
      // mon.fillHisto("eventflow"           , chTagsMain, 4, weight);
      // mon.fillHisto("eventflowNoWeights"  , chTagsMain, 4, 1);
      // mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 4, weightNoLepSF);
      if(selLeptons.size()>=2){
	 mon.fillHisto("nlep"           ,   chTags, selLeptons.size(), weight);
	 mon.fillHisto("NLep_vs_TauDecay"           ,   chTags, selLeptons.size(), decayType,weight);
	 mon.fillHisto("NTau_vs_TauDecay"           ,   chTags, selTaus.size(), decayType,weight);
	 mon.fillHisto("eventflow"           , chTagsMain, 4, weight);
	 mon.fillHisto("eventflowNoWeights"  , chTagsMain, 4, 1);
	 mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 4, weightNoLepSF);
	 mon.fillHisto("eventflow_hMC"       , chTagsMain, 4, decayType, weightNoLepSF);
	 mon.fillHisto("eventflow_ZMC"       , chTagsMain, 4, ZbosonType, weightNoLepSF);
	 mon.fillHisto("zllmass"          ,   chTagsMain, zll.mass(),    weight);
	 if(passZmass){
	   mon.fillHisto("eventflow"           , chTagsMain, 5, weight);
	   mon.fillHisto("eventflowNoWeights"  , chTagsMain, 5, 1);
	   mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 5, weightNoLepSF);
	   mon.fillHisto("eventflow_hMC"       , chTagsMain, 5, decayType, weightNoLepSF);
	   mon.fillHisto("eventflow_ZMC"       , chTagsMain, 5, ZbosonType, weightNoLepSF);
	   //pu control
	   mon.fillHisto("nvtx"        ,   chTagsMain, vtx.size(),      weight);
	   mon.fillHisto("nvtxraw"     ,   chTagsMain, vtx.size(),      weight/puWeight);
	   mon.fillHisto("nvtxpuweight"     ,   chTagsMain, vtx.size(),      weightNoLepSF);
	   mon.fillHisto("rho"         ,   chTagsMain, rho,       weight);

	   //Z kinematics control
	   mon.fillHisto("leadpt"      ,   chTagsMain, leadingLep.pt(), weight);
	   mon.fillHisto("leadeta"     ,   chTagsMain, leadingLep.eta(), weight);
	   mon.fillHisto("trailerpt"   ,   chTagsMain, trailerLep.pt(), weight);
	   mon.fillHisto("trailereta"  ,   chTagsMain, trailerLep.eta(), weight);
	   mon.fillHisto("leppt"       ,   chTagsMain, leadingLep.pt(), weight);
	   mon.fillHisto("leppt"       ,   chTagsMain, trailerLep.pt(), weight);
	   mon.fillHisto("lepeta"      ,   chTagsMain, leadingLep.eta(), weight);
	   mon.fillHisto("lepeta"      ,   chTagsMain, trailerLep.eta(), weight);

	   //analyze dilepton kinematics
	   mon.fillHisto("zllpt"         ,   chTagsMain, zll.pt(),      weight);
	   mon.fillHisto("zlleta"        ,   chTagsMain, zll.eta(),     weight);
	   mon.fillHisto("zlly"          ,   chTagsMain, zll.Rapidity(),weight);

	   if(passZpt){
	     mon.fillHisto("eventflow"           , chTagsMain, 6, weight);
	     mon.fillHisto("eventflowNoWeights"  , chTagsMain, 6, 1);
	     mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 6, weightNoLepSF);
	     mon.fillHisto("eventflow_hMC"       , chTagsMain, 6, decayType, weightNoLepSF);
	     mon.fillHisto("eventflow_ZMC"       , chTagsMain, 6, ZbosonType, weightNoLepSF);

	     mon.fillHisto("ntaus"           ,  chTags, selTaus.size(), weight);
	     mon.fillHisto("tauleadpt"       ,  chTagsMain,   selTaus.size()>0?selTaus[0].pt():-1,  weight);
	     mon.fillHisto("tauleadeta"      ,  chTagsMain,   selTaus.size()>0?selTaus[0].eta():-10, weight);
	     mon.fillHisto("tautrailerpt"    ,  chTagsMain,   selTaus.size()>1?selTaus[1].pt():-1,  weight);
	     mon.fillHisto("tautrailereta"   ,  chTagsMain,   selTaus.size()>1?selTaus[1].eta():-10, weight);
	     mon.fillHisto("taupt"           ,  chTags, selTaus.size()>0?selTaus[0].pt():-1, weight);
	     mon.fillHisto("taupt"           ,  chTags, selTaus.size()>1?selTaus[1].pt():-1, weight);
	     mon.fillHisto("taueta"          ,  chTagsMain,   selTaus.size()>0?selTaus[0].eta():-10, weight);
	     mon.fillHisto("taueta"          ,  chTagsMain,   selTaus.size()>1?selTaus[1].eta():-10, weight);

	     if(selLeptons.size()>=4){
	       mon.fillHisto("eventflow"           , chTagsMain, 7, weight);
	       mon.fillHisto("eventflowNoWeights"  , chTagsMain, 7, 1);
	       mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 7, weightNoLepSF);
	       mon.fillHisto("eventflow_hMC"       , chTagsMain, 7, decayType, weightNoLepSF);
	       mon.fillHisto("eventflow_ZMC"       , chTagsMain, 7, ZbosonType, weightNoLepSF);
	       mon.fillHisto("yields"          ,"all_4Lep",                HiggsShortId, weight);
	       mon.fillHisto("yields"          ,"all_4Lep_NoWeigth",                HiggsShortId, 1);

	       if (higgsMuonCand > -1)
		 {
		   mon.fillHisto("higgsMuonpt",        chTagsMain, selLeptons[higgsMuonCand].pt(), weight);
		   mon.fillHisto("higgsMuoneta",       chTagsMain, selLeptons[higgsMuonCand].eta(),weight);
		   mon.fillHisto("higgsMuoniso",       chTagsMain, selLeptons[higgsMuonCand].userFloat("relIso"),weight);
		   for(int l1=0;l1<(int)selLeptons.size();l1++){
		     mon.fillHisto("higgsMuonDeltaRLep", chTagsMain,deltaR(selLeptons[l1], selLeptons[higgsMuonCand]),weight);
		   }
		   for(int j1=0;j1<(int)selJets.size();j1++){
		     mon.fillHisto("higgsMuonDeltaRJets",chTagsMain, deltaR(selJets[j1], selLeptons[higgsMuonCand]) ,weight);
		   }
		 }
	       if (higgsEleCand > -1)
		 {
		   mon.fillHisto("higgsElept",        chTagsMain, selLeptons[higgsEleCand].pt(), weight);
		   mon.fillHisto("higgsEleeta",       chTagsMain, selLeptons[higgsEleCand].eta(),weight);
		   mon.fillHisto("higgsEleiso",       chTagsMain, selLeptons[higgsEleCand].userFloat("relIso"),weight);

		   for(int l1=0;l1<(int)selLeptons.size();l1++){
		     mon.fillHisto("higgsEleDeltaRLep", chTagsMain,deltaR(selLeptons[l1], selLeptons[higgsEleCand]),weight);
                   }
		   for(int j1=0;j1<(int)selJets.size();j1++){
		     mon.fillHisto("higgsEleDeltaRJets",chTagsMain, deltaR(selJets[j1], selLeptons[higgsEleCand]) ,weight);
                   }
		 }

	       if(passLepVetoMain){
		 mon.fillHisto("eventflow"           , chTagsMain, 8, weight);
		 mon.fillHisto("eventflowNoWeights"  , chTagsMain, 8, 1);
		 mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 8, weightNoLepSF);
		 mon.fillHisto("eventflow_hMC"       , chTagsMain, 8, decayType, weightNoLepSF);
		 mon.fillHisto("eventflow_ZMC"       , chTagsMain, 8, ZbosonType, weightNoLepSF);

		 mon.fillHisto("nbtags"    , chTags, nbtags,  weight);
		 mon.fillHisto("njets"     , chTags, njets,   weight);

		 if(passBJetVetoMain){
		   mon.fillHisto("eventflow"           , chTagsMain, 9, weight);
		   mon.fillHisto("eventflowNoWeights"  , chTagsMain, 9, 1);
		   mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 9, weightNoLepSF);
		   mon.fillHisto("eventflow_hMC"       , chTagsMain, 9, decayType, weightNoLepSF);
		   mon.fillHisto("eventflow_ZMC"       , chTagsMain, 9, ZbosonType, weightNoLepSF);

		   mon.fillHisto("dPhi_AZ"    , chTagsMain, deltaPhi(higgsCand.phi(), zll.phi()),    weight);
		   mon.fillHisto("dPhi_AMet"  , chTagsMain, deltaPhi(higgsCand.phi(), met.phi()),    weight);
		   mon.fillHisto("dPhi_ZMet"  , chTagsMain, deltaPhi(zll.phi(), met.phi()),    weight);
		   mon.fillHisto("met"      	, chTagsMain, met.pt()         , weight);

		   if(passDPhiCut){
		     mon.fillHisto("eventflow"           , chTagsMain, 10, weight);
		     mon.fillHisto("eventflowNoWeights"  , chTagsMain, 10, 1);
		     mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 10, weightNoLepSF);
		     mon.fillHisto("eventflow_hMC"       , chTagsMain, 10, decayType, weightNoLepSF);
		     mon.fillHisto("eventflow_ZMC"       , chTagsMain, 10, ZbosonType, weightNoLepSF);
		     if(passHiggsLoose){
		       mon.fillHisto("sumpt",   chTagsMain, selLeptons[higgsCandL1].pt()+selLeptons[higgsCandL2].pt(), weight);
		       if(passHiggsMain){
			 mon.fillHisto("eventflow"           , chTagsMain, 11, weight);
			 mon.fillHisto("eventflowNoWeights"  , chTagsMain, 11, 1);
			 mon.fillHisto("eventflowNoLepSF"    , chTagsMain, 11, weightNoLepSF);
			 mon.fillHisto("eventflow_hMC"       , chTagsMain, 11, decayType, weightNoLepSF);
			 mon.fillHisto("eventflow_ZMC"       , chTagsMain, 11, ZbosonType, weightNoLepSF);

			 mon.fillHisto("yields"          ,chTagsMain,                HiggsShortId, weight);
			 mon.fillHisto("yieldsOS"     ,chTagsMain,                HiggsShortId, weight);

			 mon.fillHisto("Apt"       	, chTagsMain, higgsCand.pt(),    weight);
			 mon.fillHisto("AmassFine"           , chTagsMain, higgsCand.mass(),  weight);
			 mon.fillHisto("AmassFineSVFit"           , chTagsMain, higgsCand_SVFitMass,  weight);
			 mon.fillHisto("AmassFineClassicSVFit"           , chTagsMain, classiSVFit_mass,  weight);
			 mon.fillHisto("AmassFineCollinear"           , chTagsMain, collinearMass,  weight);
			 // cout<< " h mass = "<<higgsCand.mass()<<" -  h coll_mass = "<<collinearMass<<" -  weight = "<<weight<<endl;
			 mon.fillHisto("Amass"           , chTagsMain, higgsCand.mass(),  weight);
			 mon.fillHisto("Amasssvfit"      , chTagsMain, higgsCand_SVFitMass,  weight);
			 mon.fillHisto("AmassClassicsvfit"      , chTagsMain, classiSVFit_mass,  weight);
			 mon.fillHisto("Hmass"           , chTagsMain, higgsCandH.mass(),  weight);
			 mon.fillHisto("Hpt"             , chTagsMain, higgsCandH.pt(),  weight);
			 mon.fillHisto("Hmasssvfit"   , chTagsMain, higgsCandH_SVFit.mass(),  weight);

			 mon.fillHisto("Anjets"    	, chTagsMain, NCleanedJetMain      , weight);
			 mon.fillHisto("Amet"      	, chTagsMain, met.pt()         , weight);
		       } // HiggsCuts
		     }   // Loose HiggsCuts
		   }     // DPhi Cut
		 }
	       }
	     }
	   }
	 }
      }
    } // end filling plots for nominal


    if( passZmass && passZpt && selLeptons.size()>=4 && passLepVetoMain && passBJetVetoMain && passDPhiCut && passHiggsLoose){
      for(unsigned int index=0; index<optim_Cuts_sumPt.size();index++){
	 bool passHiggs = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],true,vtx);
	 if(passHiggs){
	   mon.fillHisto(TString("Hsvfit_shapes")+varNames[ivar],chTagsMain,index,higgsCandH_SVFit.mass(),weight);
	   mon.fillHisto(TString("Asvfit_shapes")+varNames[ivar],chTagsMain,index,higgsCand_SVFitMass,weight);

	 } else {   //if ( runSystematics ){

	   // Control regions

	   CRTypes theCR =  checkBkgCR(selLeptons, higgsCandL1, higgsCandL2, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],vtx);

	   float theFRWeight=1;
	   mon.fillHisto("CRCounts","controlPlots",theCR,1);

	   if(theCR==CRTypes::CR10){
	     // CR10
	     theFRWeight*=getTheFRWeight(selLeptons, selJets, higgsCandL1, higgsCandL2, theFRWeightTool, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],CRTypes::CR10);

	     mon.fillHisto(TString("Hsvfit_shapes_CR10")+varNames[ivar],chTagsMain,index,higgsCandH_SVFit.mass(),weight*theFRWeight);
	     mon.fillHisto(TString("Asvfit_shapes_CR10")+varNames[ivar],chTagsMain,index,higgsCand_SVFitMass,weight*theFRWeight);
	     // cout<<" CRTypes::CR10 - FR weight value = "<<theFRWeight<<endl;

	   } else if (theCR==CRTypes::CR01) {
	     // CR01
	     theFRWeight*=getTheFRWeight(selLeptons, selJets, higgsCandL1, higgsCandL2, theFRWeightTool, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],CRTypes::CR01);

	     mon.fillHisto(TString("Hsvfit_shapes_CR01")+varNames[ivar],chTagsMain,index,higgsCandH_SVFit.mass(),weight*theFRWeight);
	     mon.fillHisto(TString("Asvfit_shapes_CR01")+varNames[ivar],chTagsMain,index,higgsCand_SVFitMass,weight*theFRWeight);
		// cout<<" CRTypes::CR01 - FR weight value = "<<theFRWeight<<endl;
	   } else {

	     // CR11
	     theFRWeight*=getTheFRWeight(selLeptons, selJets, higgsCandL1, higgsCandL2, theFRWeightTool, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],CRTypes::CR11);
	     mon.fillHisto(TString("Hsvfit_shapes_CR11")+varNames[ivar],chTagsMain,index,higgsCandH_SVFit.mass(),weight*theFRWeight);
	     mon.fillHisto(TString("Asvfit_shapes_CR11")+varNames[ivar],chTagsMain,index,higgsCand_SVFitMass,weight*theFRWeight);
	     // cout<<" CRTypes::CR11 - FR weight value = "<<theFRWeight<<endl;
	   }
	   // cout<<" FR weight value = "<<theFRWeight<<endl;
	 }

	 if(index==0 && selLeptons.size()>=2 && passZmass && passZpt && selLeptons.size()>=4 && passLepVetoMain && passBJetVetoMain ){
	   mon.fillHisto(TString("metsys")+varNames[ivar], chTagsMain, imet.pt(), weight);
	 }
      }//end of the loop on cutIndex
    }
  }//END SYSTEMATIC LOOP
}


// ------------ method called once each job just before starting event loop  ------------
void
ZHTauTauAnalyzer::beginJob()
{

  TH1F::SetDefaultSumw2(kTRUE);
  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  printf("Definition of plots");

  //event selection
  TH1 *h1=mon.addHistogram(fs->make<TH1F>("eventflow", ";;Events", 13,0,13) );
  h1->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1->GetXaxis()->SetBinLabel(2,"Trigger");
  h1->GetXaxis()->SetBinLabel(3,"METFilter");
  h1->GetXaxis()->SetBinLabel(4,"Zcandidate");
  h1->GetXaxis()->SetBinLabel(5,"Nlep#geq2");
  h1->GetXaxis()->SetBinLabel(6,"Zmass");
  h1->GetXaxis()->SetBinLabel(7,"Zkin");
  h1->GetXaxis()->SetBinLabel(8,"Nlep+Ntau#geq4");
  h1->GetXaxis()->SetBinLabel(9,"Lep Veto");
  h1->GetXaxis()->SetBinLabel(10,"Btag Veto");
  h1->GetXaxis()->SetBinLabel(11,"#Delta #phi Z-MET");
  h1->GetXaxis()->SetBinLabel(12,"di-#tau Cand");

  TH2 *eventflow_hMC= (TH2*) mon.addHistogram( fs->make<TH2F> ("eventflow_hMC", ";;Events", 13,0,13,10,0,10) );
  eventflow_hMC->GetXaxis()->SetBinLabel(1,"InitialEv");
  eventflow_hMC->GetXaxis()->SetBinLabel(2,"Trigger");
  eventflow_hMC->GetXaxis()->SetBinLabel(3,"METFilter");
  eventflow_hMC->GetXaxis()->SetBinLabel(4,"Zcandidate");
  eventflow_hMC->GetXaxis()->SetBinLabel(5,"Nlep#geq2");
  eventflow_hMC->GetXaxis()->SetBinLabel(6,"Zmass");
  eventflow_hMC->GetXaxis()->SetBinLabel(7,"Zkin");
  eventflow_hMC->GetXaxis()->SetBinLabel(8,"Nlep+Ntau#geq4");
  eventflow_hMC->GetXaxis()->SetBinLabel(9,"Lep Veto");
  eventflow_hMC->GetXaxis()->SetBinLabel(10,"Btag Veto");
  eventflow_hMC->GetXaxis()->SetBinLabel(11,"#Delta #phi Z-MET");
  eventflow_hMC->GetXaxis()->SetBinLabel(12,"di-#tau Cand");
  eventflow_hMC->GetYaxis()->SetBinLabel(1,"Unknown");
  eventflow_hMC->GetYaxis()->SetBinLabel(2,"#tau_{#mu} #tau_{#mu}");
  eventflow_hMC->GetYaxis()->SetBinLabel(3,"#tau_{#mu} #tau_{e}");
  eventflow_hMC->GetYaxis()->SetBinLabel(4,"#tau_{#mu} #tau_{h}");
  eventflow_hMC->GetYaxis()->SetBinLabel(5,"#tau_{e} #tau_{e}");
  eventflow_hMC->GetYaxis()->SetBinLabel(6,"");
  eventflow_hMC->GetYaxis()->SetBinLabel(7,"#tau_{e} #tau_{h}");
  eventflow_hMC->GetYaxis()->SetBinLabel(8,"");
  eventflow_hMC->GetYaxis()->SetBinLabel(10,"#tau_{h} #tau_{h}");

  TH2 *eventflow_ZMC= (TH2*) mon.addHistogram( fs->make<TH2F> ("eventflow_ZMC", ";;Events", 13,0,13,4,0,4) );
  eventflow_ZMC->GetXaxis()->SetBinLabel(1,"InitialEv");
  eventflow_ZMC->GetXaxis()->SetBinLabel(2,"Trigger");
  eventflow_ZMC->GetXaxis()->SetBinLabel(3,"METFilter");
  eventflow_ZMC->GetXaxis()->SetBinLabel(4,"Zcandidate");
  eventflow_ZMC->GetXaxis()->SetBinLabel(5,"Nlep#geq2");
  eventflow_ZMC->GetXaxis()->SetBinLabel(6,"Zmass");
  eventflow_ZMC->GetXaxis()->SetBinLabel(7,"Zkin");
  eventflow_ZMC->GetXaxis()->SetBinLabel(8,"Nlep+Ntau#geq4");
  eventflow_ZMC->GetXaxis()->SetBinLabel(9,"Lep Veto");
  eventflow_ZMC->GetXaxis()->SetBinLabel(10,"Btag Veto");
  eventflow_ZMC->GetXaxis()->SetBinLabel(11,"#Delta #phi Z-MET");
  eventflow_ZMC->GetXaxis()->SetBinLabel(12,"di-#tau Cand");
  eventflow_ZMC->GetYaxis()->SetBinLabel(1,"Z #rightarrow #mu#mu");
  eventflow_ZMC->GetYaxis()->SetBinLabel(2,"Z #rightarrow ee");
  eventflow_ZMC->GetYaxis()->SetBinLabel(3,"Z #rightarrow #tau#tau");
  eventflow_ZMC->GetYaxis()->SetBinLabel(4,"Z #rightarrow #nu or quarks");

  TH1 *h1_NoW=mon.addHistogram( fs->make<TH1F> ("eventflowNoWeights", ";;Events", 13,0,13) );
  h1_NoW->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1_NoW->GetXaxis()->SetBinLabel(2,"Trigger");
  h1_NoW->GetXaxis()->SetBinLabel(3,"METFilter");
  h1_NoW->GetXaxis()->SetBinLabel(4,"Zcandidate");
  h1_NoW->GetXaxis()->SetBinLabel(5,"Nlep#geq2");
  h1_NoW->GetXaxis()->SetBinLabel(6,"Zmass");
  h1_NoW->GetXaxis()->SetBinLabel(7,"Zkin");
  h1_NoW->GetXaxis()->SetBinLabel(8,"Nlep+Ntau#geq4");
  h1_NoW->GetXaxis()->SetBinLabel(9,"Lep Veto");
  h1_NoW->GetXaxis()->SetBinLabel(10,"Btag Veto");
  h1_NoW->GetXaxis()->SetBinLabel(11,"#Delta #phi Z-MET");
  h1_NoW->GetXaxis()->SetBinLabel(12,"di-#tau Cand");

  TH1 *h1_NoLepSF=mon.addHistogram( fs->make<TH1F> ("eventflowNoLepSF", ";;Events", 13,0,13) );
  h1_NoLepSF->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1_NoLepSF->GetXaxis()->SetBinLabel(2,"Trigger");
  h1_NoLepSF->GetXaxis()->SetBinLabel(3,"METFilter");
  h1_NoLepSF->GetXaxis()->SetBinLabel(4,"Zcandidate");
  h1_NoLepSF->GetXaxis()->SetBinLabel(5,"Nlep#geq2");
  h1_NoLepSF->GetXaxis()->SetBinLabel(6,"Zmass");
  h1_NoLepSF->GetXaxis()->SetBinLabel(7,"Zkin");
  h1_NoLepSF->GetXaxis()->SetBinLabel(8,"Nlep+Ntau#geq4");
  h1_NoLepSF->GetXaxis()->SetBinLabel(9,"Lep Veto");
  h1_NoLepSF->GetXaxis()->SetBinLabel(10,"Btag Veto");
  h1_NoLepSF->GetXaxis()->SetBinLabel(11,"#Delta #phi Z-MET");
  h1_NoLepSF->GetXaxis()->SetBinLabel(12,"di-#tau Cand");

  TH1 *h_tr= mon.addHistogram( fs->make<TH1F> ("trigger", ";;Events", 10,0,10) );
  h_tr->GetXaxis()->SetBinLabel(1,"#mu#mu");
  h_tr->GetXaxis()->SetBinLabel(2,"#mu");
  h_tr->GetXaxis()->SetBinLabel(3,"ee");
  h_tr->GetXaxis()->SetBinLabel(4,"e");
  h_tr->GetXaxis()->SetBinLabel(5,"e#mu");
  h_tr->GetXaxis()->SetBinLabel(6,"#gamma");

  TH1 *h2=mon.addHistogram( fs->make<TH1F> ("yields", ";;Events", 25,0,25) );
  h2->GetXaxis()->SetBinLabel(1,"OS eeee");
  h2->GetXaxis()->SetBinLabel(2,"OS ee#mu#mu");
  h2->GetXaxis()->SetBinLabel(3,"OS eee#mu");
  h2->GetXaxis()->SetBinLabel(4,"OS eee#tau");
  h2->GetXaxis()->SetBinLabel(5,"OS ee#mu#tau");
  h2->GetXaxis()->SetBinLabel(6,"OS ee#tau#tau");
  h2->GetXaxis()->SetBinLabel(7,"OS #mu#muee");
  h2->GetXaxis()->SetBinLabel(8,"OS #mu#mu#mu#mu");
  h2->GetXaxis()->SetBinLabel(9,"OS #mu#mue#mu");
  h2->GetXaxis()->SetBinLabel(10,"OS #mu#mue#tau");
  h2->GetXaxis()->SetBinLabel(11,"OS #mu#mu#mu#tau");
  h2->GetXaxis()->SetBinLabel(12,"OS #mu#mu#tau#tau");
  h2->GetXaxis()->SetBinLabel(13,"SS eeee");
  h2->GetXaxis()->SetBinLabel(14,"SS ee#mu#mu");
  h2->GetXaxis()->SetBinLabel(15,"SS eee#mu");
  h2->GetXaxis()->SetBinLabel(16,"SS eee#tau");
  h2->GetXaxis()->SetBinLabel(17,"SS ee#mu#tau");
  h2->GetXaxis()->SetBinLabel(18,"SS ee#tau#tau");
  h2->GetXaxis()->SetBinLabel(19,"SS #mu#muee");
  h2->GetXaxis()->SetBinLabel(20,"SS #mu#mu#mu#mu");
  h2->GetXaxis()->SetBinLabel(21,"SS #mu#mue#mu");
  h2->GetXaxis()->SetBinLabel(22,"SS #mu#mue#tau");
  h2->GetXaxis()->SetBinLabel(23,"SS #mu#mu#mu#tau");
  h2->GetXaxis()->SetBinLabel(24,"SS #mu#mu#tau#tau");

  TH1 *h3=mon.addHistogram( fs->make<TH1F> ("yieldsOS", ";;Events", 12,0,12) );
  h3->GetXaxis()->SetBinLabel(1,"OS eeee");
  h3->GetXaxis()->SetBinLabel(2,"OS ee#mu#mu");
  h3->GetXaxis()->SetBinLabel(3,"OS eee#mu");
  h3->GetXaxis()->SetBinLabel(4,"OS eee#tau");
  h3->GetXaxis()->SetBinLabel(5,"OS ee#mu#tau");
  h3->GetXaxis()->SetBinLabel(6,"OS ee#tau#tau");
  h3->GetXaxis()->SetBinLabel(7,"OS #mu#muee");
  h3->GetXaxis()->SetBinLabel(8,"OS #mu#mu#mu#mu");
  h3->GetXaxis()->SetBinLabel(9,"OS #mu#mue#mu");
  h3->GetXaxis()->SetBinLabel(10,"OS #mu#mue#tau");
  h3->GetXaxis()->SetBinLabel(11,"OS #mu#mu#mu#tau");
  h3->GetXaxis()->SetBinLabel(12,"OS #mu#mu#tau#tau");

  TH1 *h_jetID = mon.addHistogram( fs->make<TH1F> ("jetId",";;Events",3,0,3));
  h_jetID->GetXaxis()->SetBinLabel(1,"PFloose");
  h_jetID->GetXaxis()->SetBinLabel(2,"LooseSimplePUId");
  h_jetID->GetXaxis()->SetBinLabel(3,"PFloose & LooseSimplePUId");

  mon.addHistogram( fs->make<TH1F>( "muiso"     ,  ";I_{#mu};Events", 100,0.,1.) );
  mon.addHistogram( fs->make<TH1F>( "eleiso"     ,  ";I_{ele};Events", 100,0.,1.) );
  //Z leptons kinematics control
  mon.addHistogram( fs->make<TH1F>( "leadpt"      ,  ";p_{T}^{lead} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "leadeta"     ,  ";#eta_{lead};Events", 50,-2.6,2.6) );
  // mon.addHistogram( fs->make<TH1F>( "leadiso"     ,  ";#I_{lead};Events", 20,0.,1.) );
  mon.addHistogram( fs->make<TH1F>( "trailerpt"   ,  ";p_{T}^{trail} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "trailereta"  ,  ";#eta_{trail};Events", 50,-2.6,2.6) );
  // mon.addHistogram( fs->make<TH1F>( "traileriso"     ,  ";#I_{trail};Events", 20,0.,1.) );
  mon.addHistogram( fs->make<TH1F>( "leppt"       ,  ";p_{T}^{lepton} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "lepeta"      ,  ";#eta_{lepton};Events", 50,-2.6,2.6) );

  // zll control
  mon.addHistogram( fs->make<TH1F>( "zlly",      		";y_{ll};Events", 50,-6,6) );
  mon.addHistogram( fs->make<TH1F>( "zlleta",    		";#eta_{ll};Events", 50,-10,10) );
  mon.addHistogram( fs->make<TH1F>( "zllpt",     		";p_{T}^{ll} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "zllmass",   		";M_{ll} (GeV);Events/2 GeV", 80,20,180) );
  mon.addHistogram( fs->make<TH1F>( "nlep",      	";Number of Leptons;Events", 10,0,10) );

  mon.addHistogram( fs->make<TH1F>( "sumpt",            ";L_{T} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "dPhi_AZ",          ";#DeltaPhi(#tau#tau,ll);Events",50,-3,3));
  mon.addHistogram( fs->make<TH1F>( "dPhi_AMet",        ";#Delta#phi(#tau#tau,#slash{E}_{T});Events",50,-3,3));
  mon.addHistogram( fs->make<TH1F>( "met",             ";#slash{E}_{T} (GeV);Events/10 GeV",50,0,500));

  mon.addHistogram( fs->make<TH1F>( "Amet",             ";#slash{E}_{T} (GeV);Events/10 GeV",50,0,500));
  mon.addHistogram( fs->make<TH1F>( "Anjets",           ";Number of Jets;Events",10,-0.5,9.5));
  mon.addHistogram( fs->make<TH1F>( "Apt",              ";p_{T}^{#tau#tau} (GeV);Events/10 GeV",50,0,500));
  mon.addHistogram( fs->make<TH1F>( "Hpt",              ";p_{T}^{ll#tau#tau} (GeV);Events/10 GeV",50,0,500));

  double bins[]={5, 30,70,110,190,300,550,1800};
  int nbins=sizeof(bins)/sizeof(double) - 1;
  mon.addHistogram( fs->make<TH1F>( "AmassFine",        ";M_{#tau#tau} (GeV);Events",100,0,500));
  mon.addHistogram( fs->make<TH1F>( "AmassFineCollinear",        ";M_{#tau#tau} (GeV);Events",100,0,500));
  mon.addHistogram( fs->make<TH1F>( "AmassFineSVFit",        ";SVFit M_{#tau#tau} (GeV);Events",100,0,500));
  mon.addHistogram( fs->make<TH1F>( "AmassFineClassicSVFit",        ";ClassicSVFit M_{#tau#tau} (GeV);Events",100,0,500));
  mon.addHistogram( fs->make<TH1F>( "Amass",            ";M_{#tau#tau} (GeV);Events",nbins,bins));
  mon.addHistogram( fs->make<TH1F>( "Hmass",            ";M_{ll#tau#tau} (GeV);Events",nbins,bins));
  mon.addHistogram( fs->make<TH1F>( "Amasssvfit",       ";SVFit M_{#tau#tau} (GeV);Events",nbins,bins));
  mon.addHistogram( fs->make<TH1F>( "AmassClassicsvfit",       ";ClassicSVFit M_{#tau#tau} (GeV);Events",nbins,bins));
  mon.addHistogram( fs->make<TH1F>( "Hmasssvfit",       ";SVFit M_{ll#tau#tau} (GeV);Events",nbins,bins));

  //pu control
  mon.addHistogram( fs->make<TH1F>( "nvtx",";Vertices;Events",50,0,50) );
  mon.addHistogram( fs->make<TH1F>( "nvtxraw",";Vertices;Events",50,0,50) );
  mon.addHistogram( fs->make<TH1F>( "nvtxpuweight",";Vertices;Events",50,0,50) );
  mon.addHistogram( fs->make<TH1F>( "rho",";#rho;Events",50,0,25) );

  //tau control
  mon.addHistogram( fs->make<TH1F>( "ntaus",      	";Number of Taus;Events", 10,0,10) );
  mon.addHistogram( fs->make<TH1F>( "tauleadpt",  	";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "tauleadeta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( fs->make<TH1F>( "tautrailerpt",  	";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "tautrailereta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( fs->make<TH1F>( "taupt",  		";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "taueta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );

  //Higgs leptons control
  mon.addHistogram( fs->make<TH1F>( "higgsMuonpt"      ,  ";p_{T}^{#mu} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "higgsMuoneta"     ,  ";#eta_{#mu};Events", 50,-2.6,2.6) );
  mon.addHistogram( fs->make<TH1F>( "higgsMuoniso"     ,  ";I_{#mu};Events", 20,0.,1.) );
  mon.addHistogram( fs->make<TH1F>( "higgsMuonDeltaRLep"     ,  ";#Delta R(#mu,lepton);Events", 40,0.,2.) );
  mon.addHistogram( fs->make<TH1F>( "higgsMuonDeltaRJets"     ,  ";#Delta R(#mu,jet);Events", 40,0.,2.) );

  mon.addHistogram( fs->make<TH1F>( "higgsElept"      ,  ";p_{T}^{e} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( fs->make<TH1F>( "higgsEleeta"     ,  ";#eta_{e};Events", 50,-2.6,2.6) );
  mon.addHistogram( fs->make<TH1F>( "higgsEleiso"     ,  ";I_{e};Events", 50,0.,1.) );
  mon.addHistogram( fs->make<TH1F>( "higgsEleDeltaRLep"     ,  ";#Delta R(e,lepton);Events", 40,0.,2.) );
  mon.addHistogram( fs->make<TH1F>( "higgsEleDeltaRJets"     ,  ";#Delta R(e,jet);Events", 40,0.,2.) );

  //extra leptons in the event

  TH1 *hbtags=mon.addHistogram( fs->make<TH1F>("nbtags",   ";b-tag multiplicity;Events",5,0,5) );
  TH1 *hjets=mon.addHistogram( fs->make<TH1F>("njets",  ";Jet multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=hjets->GetXaxis()->GetNbins(); ibin++){
    TString label("");
    if(ibin==hjets->GetXaxis()->GetNbins()) label +="#geq";
    else                                    label +="=";
    label += (ibin-1);
    hjets->GetXaxis()->SetBinLabel(ibin,label);
    hbtags->GetXaxis()->SetBinLabel(ibin,label);
  }

  //fake rate histograms
  mon.addHistogram( fs->make<TH1F>( "CRCounts"  ,  ";CR type;Events", 4,-0.5,3.5));

  float ptbinsJets[] = {10, 20, 30, 40, 60, 80, 100, 125, 150, 175,250};
  int ptbinsJetsN = sizeof(ptbinsJets)/sizeof(float)-1;
  mon.addHistogram( fs->make<TH1F>( "wrtJetPt",  ";Jet p_{T} (GeV);Events",ptbinsJetsN,ptbinsJets));
  mon.addHistogram( fs->make<TH1F>( "wrtLepPt",  ";Lep p_{T} (GeV);Events",ptbinsJetsN,ptbinsJets));

  //
  // HISTOGRAMS FOR OPTIMIZATION and STATISTICAL ANALYSIS
  //

  std::vector<double> eleIsoValues = {0.3, 0.2, 0.1};
  std::vector<double> muIsoValues  = {0.3, 0.2, 0.1};
  tauIDiso = {"byLooseIsolationMVArun2v1DBoldDMwLT","byMediumIsolationMVArun2v1DBoldDMwLT"};
  // std::vector<std::string> tauIDiso = {"byLooseCombinedIsolationDeltaBetaCorr3Hits","byLooseIsolationMVArun2v1DBoldDMwLT",
  //                                       "byLooseIsolationMVArun2v1DBdR03oldDMwLT"};


  for(unsigned int elIso=0; elIso<eleIsoValues.size(); elIso++){
    for(unsigned int muIso=0; muIso<muIsoValues.size(); muIso++){
      for(unsigned int taIso=0; taIso<tauIDiso.size(); taIso++){
        for(float sumPt=0;sumPt<=200;sumPt+=20){
          optim_Cuts_elIso.push_back( eleIsoValues.at(elIso) );
          optim_Cuts_muIso.push_back( muIsoValues.at(muIso) );
          optim_Cuts_taIso.push_back(taIso);
          optim_Cuts_sumPt.push_back(sumPt);
        }
      }
    }
  }

  mon.addHistogram( fs->make<TH2F>("NLep_vs_TauDecay",";NLep; #tau decay mode",8,0,8,10,0,10));
  mon.addHistogram( fs->make<TH2F>("NTau_vs_TauDecay",";NTau; #tau decay mode",4,0,4,10,0,10));

  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(), 4, 0, 4)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "eIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(2, "muIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(3, "tauIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(4, "sumPt>");

  for(unsigned int index=0;index<optim_Cuts_sumPt.size();index++){
    Hoptim_cuts->Fill(index,0.0,optim_Cuts_elIso[index]);
    Hoptim_cuts->Fill(index,1.0,optim_Cuts_muIso[index]);
    Hoptim_cuts->Fill(index,2.0,optim_Cuts_taIso[index]);
    Hoptim_cuts->Fill(index,3.0,optim_Cuts_sumPt[index]);
  }

  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( fs->make<TH1F> ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
    Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    mon.addHistogram( fs->make<TH2F> (TString("Hsvfit_shapes")+varNames[ivar],";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Hsvfit_shapes_CR01")+varNames[ivar],";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Hsvfit_shapes_CR10")+varNames[ivar],";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Hsvfit_shapes_CR11")+varNames[ivar],";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Asvfit_shapes")+varNames[ivar],";cut index;M_{#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Asvfit_shapes_CR01")+varNames[ivar],";cut index;M_{#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Asvfit_shapes_CR10")+varNames[ivar],";cut index;M_{#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH2F> (TString("Asvfit_shapes_CR11")+varNames[ivar],";cut index;M_{#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
    mon.addHistogram( fs->make<TH1F>(TString("metsys")+varNames[ivar],                   ";#slash{E}_{T} (GeV);Events/10 GeV",50,0,500));
  }

  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  //MC normalization (to 1/pb)
  // if(isMC) xsecWeight=xsec/utils::getTotalNumberOfEvents(urls, false, true);//need to use the slow method in order to take NLO negative events into account
  xsecWeight = 1;

  //MET CORRection level
  metcor = pat::MET::METCorrectionLevel::Type1XY;

  //jet energy scale and uncertainties
  if (is2016data) {
    if(dtag.Contains("2016H")){
      jecDir+="Moriond17_80X/PromptReco_DATA/";
    } else {
      jecDir+="Moriond17_80X/2016SeptRepro_DATA/";
    }
  } else {
    jecDir+="Moriond17_80X/Summer16_re-digi-reco_MC/";
  }

  jesCor = utils::cmssw::getJetCorrector(jecDir,isMC);
  TString pf(isMC ? "MC" : "DATA");
  totalJESUnc = new JetCorrectionUncertainty((jecDir+"/"+pf+"_Uncertainty_AK4PFchs.txt").Data());

  //muon energy scale and uncertainties
  muCorMoriond17 = new RoccoR_Moriond17(std::string(rocChorPath.Data()));

  //--- photonESC initializes an EnergyScaleCorrection_class that provides the full file path
  std::string photonESC = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_74x_pho";
  PhotonEnCorrector = PhotonEnergyCalibratorRun2(isMC, false, photonESC);
  PhotonEnCorrector.initPrivateRng(new TRandom(1234));

  eScaler = EnergyScaleCorrection_class("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele");
  eScaler.doScale=true;
  eScaler.doSmearings=true;

  //Fake rate tool
  theFRWeightTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/test/zhtautau/FR_Weights.root").c_str());

  btagCalib = BTagCalibration("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_Moriond17_B_H.csv");

  btagCal80X = BTagCalibrationReader80X(BTagEntry::OP_LOOSE, "central", {"up", "down"});

  // setup calibration readers 80X
  btagCal80X.load(btagCalib, BTagEntry::FLAV_B, "comb");
  btagCal80X.load(btagCalib, BTagEntry::FLAV_C, "comb");
  btagCal80X.load(btagCalib, BTagEntry::FLAV_UDSG, "incl");

  //pileup weighting
  LumiWeights = NULL;

  //FIXME if(isMC){
  //   std::vector<float> dataPileupDistribution;
  //   for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){
  //     dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
  //   }
  //   std::vector<float> mcPileupDistribution;
  //
  //   utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
  //   while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  //   while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
  //   LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
  //   PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
  //   utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  // }

  if(!isMC) {
    if(is2015data){
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_csc2015.txt");
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_ecalscn1043093.txt");
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_csc2015.txt");
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_ecalscn1043093.txt");
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_csc2015.txt");
      metFilter.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_ecalscn1043093.txt");
    }
    if(is2016data){}
  }

  // initialize the TRandom
  rgenMuon_ = std::shared_ptr<TRandom3>(new TRandom3(0));
  rgenEle_  = std::shared_ptr<TRandom3>(new TRandom3(1));

}

// ------------ method called once each job just after ending the event loop  ------------
void
ZHTauTauAnalyzer::endJob()
{
  fs->cd();
  mon.Write();
}


// here go all the auxilliary methods

//***********************************************************************************************//
double ZHTauTauAnalyzer::getSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2)
//***********************************************************************************************//
{
  if(higgsCandL1<0 || higgsCandL2<0) return 0.;

  TMatrixD covMET(2, 2); // PFMET significance matrix

  covMET[0][0] = met.getSignificanceMatrix()(0,0);
  covMET[0][1] = met.getSignificanceMatrix()(0,1);
  covMET[1][0] = met.getSignificanceMatrix()(1,0);
  covMET[1][1] = met.getSignificanceMatrix()(1,1);

  // std::cout<<"MET MATRIX: " << covMET[0][0] << " " << covMET[0][1] << " " << covMET[1][0] << " " << covMET[1][1] << "\n";

  int dlid = abs( selLeptons[higgsCandL1].pdgId() * selLeptons[higgsCandL2].pdgId() );

  //std::cout<<"\n"<<selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << " Di-Tau ID ------------> " << dlid << std::endl;

  if ( dlid == HiggsFinalStates::elha || dlid == HiggsFinalStates::muha){
    if (abs(selLeptons[higgsCandL1].pdgId()) == 15) {
      // Switching leptons in semileptonic pairs:
      // e and mu should be passed as first MesuredTauLepton
      int temp = higgsCandL1;
      higgsCandL1  = higgsCandL2;
      higgsCandL2  = temp;
    }
  }

  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;

  if ( dlid == HiggsFinalStates::elha ){
    //std::cout<< " ETau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), svFitStandalone::electronMass) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if( dlid == HiggsFinalStates::muha ){
    //std::cout<< " MuTau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), svFitStandalone::muonMass) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if ( dlid == HiggsFinalStates::haha ){
    //std::cout<< " TauTau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), selLeptons[higgsCandL1].mass(), selLeptons[higgsCandL1].tau.decayMode()) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if (dlid == HiggsFinalStates::elmu ){
    //std::cout<< " EMu Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), svFitStandalone::electronMass) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), svFitStandalone::muonMass) );
  }
  else if (dlid == HiggsFinalStates::mumu){
    //std::cout<< " EE Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
                    selLeptons[higgsCandL1].phi(), svFitStandalone::muonMass) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
                    selLeptons[higgsCandL2].phi(), svFitStandalone::muonMass) );
  }
  else if (dlid == HiggsFinalStates::elel){
    //std::cout<< " EE Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
                selLeptons[higgsCandL1].phi(), svFitStandalone::electronMass) );
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
                selLeptons[higgsCandL2].phi(), svFitStandalone::electronMass) );
  }
  else return (selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4()).mass();


  SVfitStandaloneAlgorithm algo(measuredTauLeptons, met.px(), met.py() , covMET, 0);
  algo.addLogM(false);
  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);
  TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  algo.shiftVisPt(true, inputFile_visPtResolution);

  algo.integrateMarkovChain();

  double mass = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getMass(); // full mass of tau lepton pair in units of GeV

  //double mass = algo.getMass(); // Full SVFit mass - return value is in units of GeV
  //double transverse_mass = algo.getTransverseMass(); // Transverse SVFit mass

  // if ( algo.isValidSolution() ) {
  //   std::cout << "\t-------- StandaloneSVfit --------"<<endl;
  //   std::cout << "found mass = " << mass << std::endl;
  // } else {
  //   std::cout << "sorry -- status of NLL is not valid [" << algo.fitStatus() << "]" << std::endl;
  // }

  delete inputFile_visPtResolution;

  //return LorentzVector(selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());
  return mass;
}


//***********************************************************************************************//
LorentzVector ZHTauTauAnalyzer::getClassicSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2)
//***********************************************************************************************//
{
  using namespace classic_svFit;

  LorentzVector pEmpty(0,0,0,0);
  if(higgsCandL1<0 || higgsCandL2<0) return pEmpty;

  TMatrixD covMET(2, 2); // PFMET significance matrix

  covMET[0][0] = met.getSignificanceMatrix()(0,0);
  covMET[0][1] = met.getSignificanceMatrix()(0,1);
  covMET[1][0] = met.getSignificanceMatrix()(1,0);
  covMET[1][1] = met.getSignificanceMatrix()(1,1);

  // std::cout<<"MET MATRIX: " << covMET[0][0] << " " << covMET[0][1] << " " << covMET[1][0] << " " << covMET[1][1] << "\n";

  int dlid = abs( selLeptons[higgsCandL1].pdgId() * selLeptons[higgsCandL2].pdgId() );

  //std::cout<<"\n"<<selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << " Di-Tau ID ------------> " << dlid << std::endl;

  if ( dlid == HiggsFinalStates::elha || dlid == HiggsFinalStates::muha){
    if (abs(selLeptons[higgsCandL1].pdgId()) == 15) {
      // Switching leptons in semileptonic pairs:
      // e and mu should be passed as first MesuredTauLepton
      int temp = higgsCandL1;
      higgsCandL1  = higgsCandL2;
      higgsCandL2  = temp;
    }
  }

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

  if ( dlid == HiggsFinalStates::elha ){
    //std::cout<< " ETau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if( dlid == HiggsFinalStates::muha ){
    //std::cout<< " MuTau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), classic_svFit::muonMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if ( dlid == HiggsFinalStates::haha ){
    //std::cout<< " TauTau Pair --- > "<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), selLeptons[higgsCandL1].mass(), selLeptons[higgsCandL1].tau.decayMode()) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass(), selLeptons[higgsCandL2].tau.decayMode()) );
  }
  else if (dlid == HiggsFinalStates::elmu ){
    //std::cout<< " EMu Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), classic_svFit::muonMass) );
  }
  else if (dlid == HiggsFinalStates::mumu){
    //std::cout<< " EE Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
        measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), classic_svFit::muonMass) );
	measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), classic_svFit::muonMass) );
  }
  else if (dlid == HiggsFinalStates::elel){
    //std::cout<< " EE Pair  --->"<< selLeptons[higgsCandL1].pdgId() << "  " << selLeptons[higgsCandL2].pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(),
								    selLeptons[higgsCandL1].phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(),
								    selLeptons[higgsCandL2].phi(), classic_svFit::electronMass) );
  }
  else return (selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());


  int verbosity = 0;
  ClassicSVfit *svFitAlgo = new ClassicSVfit(verbosity);
#ifdef USE_SVFITTF
  // HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  // svFitAlgo->setHadTauTF(hadTauTF);
  // svFitAlgo->enableHadTauTF();
#endif
  //svFitAlgo.addLogM_fixed(false);
  svFitAlgo->addLogM_fixed(true, 6.);
  //svFitAlgo.addLogM_dynamic(true, "(m/1000.)*15.");
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  svFitAlgo->integrate(measuredTauLeptons, met.px(), met.py() , covMET);
  // bool isValidSolution = svFitAlgo->isValidSolution();

  DiTauSystemHistogramAdapter* DiTauSystemPtr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo->getHistogramAdapter());

  double mass = DiTauSystemPtr->getMass();
  //double massErr = DiTauSystemPtr->getMassErr();
  //double transverseMass = DiTauSystemPtr->getTransverseMass();
  //double transverseMassErr = DiTauSystemPtr->getTransverseMassErr();
  double pt = DiTauSystemPtr->getPt();
  //double ptErr = DiTauSystemPtr->getPtErr();
  double eta = DiTauSystemPtr->getEta();
  double phi = DiTauSystemPtr->getPhi();

  double p  = pt*TMath::CosH(eta);
  double px = pt*TMath::Cos(phi);
  double py = pt*TMath::Sin(phi);
  double pz = pt*TMath::SinH(eta);
  double energy = TMath::Sqrt(p*p + mass*mass);
  LorentzVector p4(px, py, pz, energy);
  // if ( isValidSolution ) {
  //   std::cout << "\t-------- ClassicSVfit --------"<<endl;
  //   std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " \n,"
  //             << " transverse mass = " << transverseMass << " +/- " << transverseMassErr << " \n,"
  //             << " pt = " << pt << " +/- " << ptErr << std::endl;
  // } else {
  //   std::cout << "sorry, failed to find valid solution !!" << std::endl;
  // }
  delete svFitAlgo;
#ifdef USE_SVFITTF
  // delete hadTauTF;
#endif
  return p4;
}

//**********************************************************************************************//
CRTypes ZHTauTauAnalyzer::checkBkgCR(std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2, double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut, reco::VertexCollection vtx){
//**********************************************************************************************//

  using namespace patUtils; //<---- needed for the cut version

  CRTypes theCR = CRTypes::DEFAULT;
  std::vector<patUtils::GenericLepton*> HiggsLegs = {&(selLeptons[higgsCandL1]), &(selLeptons[higgsCandL2])};

  std::vector<bool> passId;
  std::vector<bool> passIso;



  for(auto lepIt=HiggsLegs.begin();lepIt!=HiggsLegs.end();lepIt++){
    patUtils::GenericLepton* lep = (*lepIt);
    if(abs(lep->pdgId())==11){
      // cout << " Electorn - Pass ID: "<< patUtils::passId(lep->el, vtx[0], llvvElecId::Loose, CutVersion::CutSet::ICHEP16Cut, true)
      //      << "\n Iso: "<< lep->userFloat("relIso")<<endl;
      passId.push_back(patUtils::passId(lep->el, vtx[0], patUtils::llvvElecId::Loose, CutVersion::CutSet::ICHEP16Cut, true));
      passIso.push_back((lep->userFloat("relIso") <= isoElCut));
    }else if(abs(lep->pdgId())==13){
      passId.push_back(patUtils::passId(lep->mu, vtx[0], patUtils::llvvMuonId::Loose, CutVersion::CutSet::ICHEP16Cut));
      passIso.push_back((lep->userFloat("relIso") <= isoMuCut));
    }else if(abs(lep->pdgId())==15){
      passId.push_back(lep->tau.tauID("againstElectronTightMVA6") && lep->tau.tauID("againstMuonLoose3"));
      passIso.push_back(bool(lep->tau.tauID(isoHaCut)));
    }
  }

  passId.resize(2);
  passIso.resize(2);

  if( (!passId[0] || !passIso[0]) && (passId[1]&&passIso[1]) ) {
    theCR=CRTypes::CR10;
  } else if (  (passId[0] || passIso[0]) && (!passId[1] || !passIso[1] ) ) {
    theCR=CRTypes::CR01;
  } else if (  (!passId[0] || !passIso[0]) && (!passId[1] || !passIso[1])  ) {
    theCR=CRTypes::CR11;
  }

  return theCR;

}

//**********************************************************************************************//
bool ZHTauTauAnalyzer::passHiggsCuts(std::vector<patUtils::GenericLepton> selLeptons,
				     int higgsCandL1, int higgsCandL2,
				     double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut, bool requireId, reco::VertexCollection vtx)
//**********************************************************************************************//
{
  if(higgsCandL1<0 || higgsCandL2<0)return false;
  std::vector<patUtils::GenericLepton*> HiggsLegs = {&(selLeptons[higgsCandL1]), &(selLeptons[higgsCandL2])};
  // std::vector< std::shared_ptr<patUtils::GenericLepton> > HiggsLegs;
  // auto lep1 = std::make_shared<patUtils::GenericLepton>(&(selLeptons[higgsCandL1]));
  // auto lep2 = std::make_shared<patUtils::GenericLepton>(&(selLeptons[higgsCandL2]));
  // HiggsLegs.push_back(lep1); HiggsLegs.push_back(lep2);

  bool passId=true;
  bool passIso=true;
  float sumpt = 0;
  for(auto& lep: HiggsLegs){
    //patUtils::GenericLepton* lep = (*lepIt);

    using namespace patUtils; //<----

    if(abs(lep->pdgId())==11){
      //        passId  &= patUtils::passId(lep->el, vtx[0], patUtils::llvvElecId::Loose);
      passId  &= patUtils::passId(lep->el, vtx[0], patUtils::llvvElecId::Loose, CutVersion::CutSet::ICHEP16Cut, true);
      passIso &= (lep->userFloat("relIso") <= isoElCut);
    }else if(abs(lep->pdgId())==13){
      //        passId  &= patUtils::passId(lep->mu, vtx[0], patUtils::llvvMuonId::Loose);
      passId  &= patUtils::passId(lep->mu, vtx[0], patUtils::llvvMuonId::Loose, CutVersion::CutSet::ICHEP16Cut);
      passIso &= (lep->userFloat("relIso") <= isoMuCut);
    }else if(abs(lep->pdgId())==15){
      passId  &= lep->tau.tauID("againstElectronTightMVA6") && lep->tau.tauID("againstMuonLoose3");
      passIso &= bool(lep->tau.tauID(isoHaCut));
    }
    sumpt += lep->pt();
  }
  return sumpt>sumPtCut && passIso && (passId || !requireId);
}



//**********************************************************************************************//
double ZHTauTauAnalyzer::closestJet(const LorentzVector& obj, pat::JetCollection& selJets, int& closestJetIndex)
//**********************************************************************************************//
{
  double dRMin = 1E100;  closestJetIndex = -1;
  for(int j=0;j<(int)selJets.size();j++){
    double dR = deltaR(selJets[j].p4(), obj);
    if(dR<dRMin){dRMin=dR; closestJetIndex=j;}
  }
  return dRMin;
}

 //**********************************************************************************************//
float ZHTauTauAnalyzer::getTheFRWeight(std::vector<patUtils::GenericLepton>& selLeptons, pat::JetCollection& selJets,
				       int higgsCandL1, int higgsCandL2, FRWeights theFRWeightTool,double isoElCut, double isoMuCut, std::string isoHaCut, float sumPtCut,
				       CRTypes CRType)
//**********************************************************************************************//
{
  float theFinalWeight=1;

  std::vector<patUtils::GenericLepton*> HiggsLegs = {&(selLeptons[higgsCandL1]), &(selLeptons[higgsCandL2])};
  std::vector<float> theWeights;

  for(auto& lep: HiggsLegs){
    //patUtils::GenericLepton* lep = (*lepIt);

    int closestJetIndex=-1; double pTclosestJet=-1; double etaclosestJet=-1;
    std::string etabin = "";
    std::string isobin = "";

    double dRmin = closestJet(lep->p4(), selJets, closestJetIndex);
    if(closestJetIndex>=0 && dRmin<0.5){pTclosestJet=selJets[closestJetIndex].pt(); etaclosestJet=abs(selJets[closestJetIndex].eta());}

    if(sumPtCut<30){
      etabin = etaclosestJet <1.4 ? "_B" : "_E";
    } else {
      etabin = etaclosestJet <1.4 ? "_TMCut_B" : "_TMCut_E";
    }

    //cout<<" isoElCut: "<< isoElCut << (isoElCut<=0.3 ? ":Pass ":":FALSE ")<< " - isoMuCut: " << isoMuCut << " - isoHaCut: "<< isoHaCut << endl;
    //cout << std::setprecision(20) << fabs(isoElCut - 0.3) << endl;
    if(abs(lep->pdgId())==11){

      if(isoElCut<=0.3) isobin= "_Id_Iso03weight";
      if(isoElCut<=0.2) isobin= "_Id_Iso02weight";
      if(isoElCut<=0.1) isobin= "_Id_Iso01weight";

      theWeights.push_back(theFRWeightTool.getWeight("FR_El",etabin,isobin,"_wrtJetPt",pTclosestJet));

    } else if (abs(lep->pdgId())==13){

      if(isoMuCut<=0.3) isobin= "_Id_Iso03weight";
      if(isoMuCut<=0.2) isobin= "_Id_Iso02weight";
      if(isoMuCut<=0.1) isobin= "_Id_Iso01weight";

      theWeights.push_back(theFRWeightTool.getWeight("FR_Mu",etabin,isobin,"_wrtJetPt",pTclosestJet));

     } else if(abs(lep->pdgId())==15){

      if (isoHaCut.find("CombinedIsolationDeltaBetaCorr3Hits") != std::string::npos){
        if(isoHaCut.find("byLooseCombinedIsolationDeltaBetaCorr3Hits") != std::string::npos)    isobin = "_Id_IsoLoweight";
        if(isoHaCut.find("byMediumCombinedIsolationDeltaBetaCorr3Hits")!= std::string::npos) isobin = "_Id_IsoMeweight";
      }
      if (isoHaCut.find("IsolationMVArun2v1DBoldDMwLT") != std::string::npos){
        if(isoHaCut.find("byLooseIsolationMVArun2v1DBoldDMwLT") != std::string::npos)    isobin = "_Id_IsoLo_MVAweight";
        if(isoHaCut.find("byMediumIsolationMVArun2v1DBoldDMwLT")!= std::string::npos) isobin = "_Id_IsoMe_MVAweight";
      }
      theWeights.push_back(theFRWeightTool.getWeight("FR_Ta",etabin,isobin,"_wrtJetPt",pTclosestJet));

    }
  }

  if(CRType==CRTypes::CR10){
    theFinalWeight = theWeights[0] / (1 - theWeights[0]);
  } else if (CRType==CRTypes::CR01){
    theFinalWeight = theWeights[1] / (1 - theWeights[1]);
  } else if (CRType==CRTypes::CR11){
    theFinalWeight = - (theWeights[0]*theWeights[1]) / (1 + theWeights[0]*theWeights[1] - theWeights[0] -theWeights[1] );
  } else {
    cout << "unrecognized Control Region" << endl;
  }

  return theFinalWeight;

}

//**********************************************************************************************//
std::vector<patUtils::GenericLepton> ZHTauTauAnalyzer::getLepVariations(  std::vector<patUtils::GenericLepton>& selLeptons, float factor)
//**********************************************************************************************//
{
  std::vector<patUtils::GenericLepton> selLeptonsNew;
  for(size_t ilep=0; ilep<selLeptons.size(); ilep++){
    if(abs(selLeptons[ilep].pdgId())!=15){
      patUtils::GenericLepton selLeptonNew = selLeptons[ilep];
      selLeptonNew.setP4(selLeptons[ilep].p4() * factor);
      selLeptonsNew.push_back(selLeptonNew);
    }else{
      selLeptonsNew.push_back(selLeptons[ilep]);
    }
  }
  std::sort(selLeptonsNew.begin(), selLeptonsNew.end(), utils::sort_CandidatesByPt);
  return selLeptonsNew;
}

//**********************************************************************************************//
std::vector<patUtils::GenericLepton> ZHTauTauAnalyzer::getTauVariations( std::vector<patUtils::GenericLepton>& selLeptons,float factor)
//**********************************************************************************************//
{
  std::vector<patUtils::GenericLepton> selLeptonsNew;
  for(size_t ilep=0; ilep<selLeptons.size(); ilep++){
    if(abs(selLeptons[ilep].pdgId())==15){
      patUtils::GenericLepton selLeptonNew = selLeptons[ilep];
      selLeptonNew.setP4(selLeptons[ilep].p4() * factor);
      selLeptonsNew.push_back(selLeptonNew);
    }else{
      selLeptonsNew.push_back(selLeptons[ilep]);
    }
  }
  std::sort(selLeptonsNew.begin(), selLeptonsNew.end(), utils::sort_CandidatesByPt);
  return selLeptonsNew;
}

//**********************************************************************************************//
int ZHTauTauAnalyzer::tauDecayMode(const reco::GenParticle *genParticle)
//**********************************************************************************************//
{
  int decayType = 0;
  const reco::GenParticleRefVector& daughterRefs = (*genParticle).daughterRefVector();

  for(auto& daughterRef: daughterRefs) {
     const reco::GenParticle *daughter(daughterRef.get());
     int dpdgId = fabs(daughter->pdgId());
     if (dpdgId == 13 or dpdgId == 11 or dpdgId == 111 or dpdgId == 211 or dpdgId == 311 or dpdgId == 321){
       if (dpdgId == 13) decayType = 1;
       if (dpdgId == 11) decayType = 2;
       if (dpdgId == 111 or dpdgId == 211 or dpdgId == 311 or dpdgId == 321) decayType = 3;
       return decayType;
     }
     else if (dpdgId == 15){
       return tauDecayMode(daughter);
     }
   }
  return decayType; // Unknown hadronic decay mode
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZHTauTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZHTauTauAnalyzer);
