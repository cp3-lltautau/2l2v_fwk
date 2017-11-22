// -*- C++ -*-
//
// Package:    UserCode/ZHTauTauAnalyzer
// Class:      ZHTauTauAnalyzer
// 
/**\class ZHTauTauAnalyzer ZHTauTauAnalyzer.cc UserCode/ZHTauTauAnalyzer/plugins/ZHTauTauAnalyzer.cc

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
      SmartSelectionMonitor mon;
  
      /// file service 
      edm::Service<TFileService> fs;
    
      // configure the process

      bool isMC;        
      double xsec;      
      int mctruthmode;  
      TString dtag;     
      TString suffix;   
      std::vector<std::string> urls; 
      TString outUrl;  
      TString jecDir; 
      TString rocChorPath;
  
      bool filterOnlyEE, filterOnlyMUMU, filterOnlyEMU, filterOnlyPhoton, filterOnlyE, filterOnlyMU;

      bool isV0JetsMC,isWGmc,isZGmc,isMC_ZZ,isMC_ZZ2l2nu,isMC_WZ,isMC_WZ3lnu,isMC_QCD,isMC_GJet,is2015data,is2015MC,is2016data,is2016MC;     

      //tree info
      TString dirname; 

      //systematics
      bool runSystematics;                        

      std::vector<TString> varNames;
      std::vector<string> jetVarNames = {"", "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"};
      size_t nvarsToInclude;

      //ELECTROWEAK CORRECTION WEIGHTS
      std::vector<std::vector<float>> ewkTable, ZZ_NNLOTable;

      std::vector<double> dataPileupDistributionDouble; 

      // corrections
      pat::MET::METCorrectionLevel metcor;
      FactorizedJetCorrector *jesCor;
      JetCorrectionUncertainty *totalJESUnc;
      RoccoR_Moriond17  *muCorMoriond17; 
      PhotonEnergyCalibratorRun2 PhotonEnCorrector;
      EnergyScaleCorrection_class eScaler;
      LeptonEfficiencySF lepEff;
      FRWeights theFRWeightTool;
      patUtils::MetFilter metFilter;
      edm::LumiReWeighting* LumiWeights; 
      utils::cmssw::PuShifter_t PuShifters;
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
  urls(iConfig.getUntrackedParameter<std::vector<std::string> >("input")),	     
  outUrl(iConfig.getParameter<std::string>("outfile")),                         
  jecDir(iConfig.getParameter<std::string>("jecDir")),
  rocChorPath(iConfig.getUntrackedParameter<std::string>("rocChorPath",std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/rcdata.2016.v3")),
  dirname(iConfig.getParameter<std::string>("dirName")),
  runSystematics(iConfig.getParameter<bool>("runSystematics")),
  dataPileupDistributionDouble(iConfig.getParameter< std::vector<double> >("datapileup"))
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


//
// member functions
//

// ------------ method called for each event  ------------
void
ZHTauTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //load all the objects we will need to access

   vector<pat::TriggerObjectStandAlone> triggerObjects;
   edm::Handle< vector<pat::TriggerObjectStandAlone> > triggerObjectsHandle;
   iEvent.getByLabel("selectedPatTrigger",triggerObjectsHandle);
   
   /*
   edm::Handle< edm::TriggerResults > triggerBitsHandle;
   triggerBitsHandle.getByLabel(ev,"TriggerResults","","HLT");
    
   const edm::TriggerNames &names = ev.triggerNames(*triggerBitsHandle);
   // auto names = tr.triggerNames();
   // for (auto& name: names){
   //  if (tr.accept(name))  cout<<" Trigger: "<<name <<"  "<< (tr.accept(name) ? "PASS" : "fail (or not run)")
   //             << std::endl;
   // }
   
   */

   reco::VertexCollection vtx;
   edm::Handle< reco::VertexCollection > vtxHandle;
   iEvent.getByLabel("offlineSlimmedPrimaryVertices",vtxHandle);
   if(vtxHandle.isValid()){ vtx = *vtxHandle;}

   //   double rho = 0;
   edm::Handle< double > rhoHandle;
   iEvent.getByLabel("fixedGridRhoFastjetAll",rhoHandle);
   //if(rhoHandle.isValid()){ rho = *rhoHandle;}
   
   pat::MuonCollection muons;
   edm::Handle< pat::MuonCollection > muonsHandle;
   iEvent.getByLabel("slimmedMuons",muonsHandle);
   if(muonsHandle.isValid()){ muons = *muonsHandle;}

   pat::ElectronCollection electrons;
   edm::Handle< pat::ElectronCollection > electronsHandle;
   iEvent.getByLabel("slimmedElectrons",electronsHandle);
   if(electronsHandle.isValid()){ electrons = *electronsHandle;}
   
   EcalRecHitCollection recHitCollectionEB;
   EcalRecHitCollection recHitCollectionEE;
   edm::Handle<EcalRecHitCollection> recHitCollectionEBHandle;
   edm::Handle<EcalRecHitCollection> recHitCollectionEEHandle;
   iEvent.getByLabel("reducedEgamma","reducedEBRecHits",recHitCollectionEBHandle );
   iEvent.getByLabel("reducedEgamma","reducedEERecHits",recHitCollectionEEHandle );
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}
    
   pat::JetCollection jets;
   edm::Handle< pat::JetCollection > jetsHandle;
   iEvent.getByLabel("slimmedJets",jetsHandle);
   if(jetsHandle.isValid()){ jets = *jetsHandle;}
   
   pat::METCollection mets;
   edm::Handle< pat::METCollection > metsHandle;
   iEvent.getByLabel("slimmedMETs",   metsHandle);
   if(metsHandle.isValid()){ mets = *metsHandle;}
   pat::MET met = mets[0];
   
   pat::METCollection puppimets;
   edm::Handle< pat::METCollection > puppimetsHandle;
   iEvent.getByLabel("slimmedMETsPuppi",puppimetsHandle);
   if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
   // LorentzVector puppimet = puppimets[0].p4();
   
   pat::TauCollection taus;
   edm::Handle< pat::TauCollection > tausHandle;
   iEvent.getByLabel("slimmedTaus",tausHandle);
   if(tausHandle.isValid()){ taus = *tausHandle;}
  

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
  mon.addHistogram( fs->make<TH1F>( "wrtJetPt",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( fs->make<TH1F>( "wrtLepPt",  ";Lep p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));

  //
  // HISTOGRAMS FOR OPTIMIZATION and STATISTICAL ANALYSIS
  //

  std::vector<double> eleIsoValues = {0.3, 0.2, 0.1};
  std::vector<double> muIsoValues  = {0.3, 0.2, 0.1};
  std::vector<std::string> tauIDiso = {"byLooseIsolationMVArun2v1DBoldDMwLT","byMediumIsolationMVArun2v1DBoldDMwLT"};
  // std::vector<std::string> tauIDiso = {"byLooseCombinedIsolationDeltaBetaCorr3Hits","byLooseIsolationMVArun2v1DBoldDMwLT",
  //                                       "byLooseIsolationMVArun2v1DBdR03oldDMwLT"};

  std::vector<float>    optim_Cuts_sumPt;
  std::vector<int>      optim_Cuts_taIso;
  std::vector<double>    optim_Cuts_muIso;
  std::vector<double>    optim_Cuts_elIso;

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
  double xsecWeight = 1.0;
  if(isMC) xsecWeight=xsec/utils::getTotalNumberOfEvents(urls, false, true);//need to use the slow method in order to take NLO negative events into account

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

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  float beff(0.68), sfb(0.99), sfbunc(0.015);
  float leff(0.13), sfl(1.05), sflunc(0.12);

  //double btagLoose = 0.605; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
  double btagLoose = 0.5426;  //Moriond17 recommendation Loose
  double btagMedium = 0.8484;

  // setup calibration readers 80X
  BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_Moriond17_B_H.csv");

  BTagCalibrationReader80X btagCal80X   (BTagEntry::OP_LOOSE, "central", {"up", "down"});
  btagCal80X.load(btagCalib, BTagEntry::FLAV_B, "comb");
  btagCal80X.load(btagCalib, BTagEntry::FLAV_C, "comb");
  btagCal80X.load(btagCalib, BTagEntry::FLAV_UDSG, "incl");

  beff = 0.827; sfb = 0.980; //for Loose WP  //sfb is from page 7 https://indico.cern.ch/event/557018/contributions/2246312/attachments/1310986/1961665/csvSF_rwt_July18th_2016.pdf
  leff = 0.132;

  //pileup weighting
  LumiWeights = NULL;
  double PUNorm[] = {1,1,1};
  if(isMC){
    std::vector<float> dataPileupDistribution; 
    for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){
      dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    }
    std::vector<float> mcPileupDistribution;

    utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
    LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }
    
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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZHTauTauAnalyzer::endJob() 
{
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
  bool isValidSolution = svFitAlgo->isValidSolution();

  DiTauSystemHistogramAdapter* DiTauSystemPtr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo->getHistogramAdapter());

  double mass = DiTauSystemPtr->getMass();
  double massErr = DiTauSystemPtr->getMassErr();
  double transverseMass = DiTauSystemPtr->getTransverseMass();
  double transverseMassErr = DiTauSystemPtr->getTransverseMassErr();
  double pt = DiTauSystemPtr->getPt();
  double ptErr = DiTauSystemPtr->getPtErr();
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
