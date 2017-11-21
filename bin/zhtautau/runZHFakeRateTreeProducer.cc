#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
//#include "DataFormats/FWLite/interface/ChainEvent.h"
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

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

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


//**********************************************************************************************//
double closestJet(const LorentzVector& obj, pat::JetCollection& selJets, int& closestJetIndex)
//**********************************************************************************************//
{
  double dRMin = 1E100;  closestJetIndex = -1;
  for(int j=0;j<(int)selJets.size();j++){
    double dR = deltaR(selJets[j].p4(), obj);
    if(dR<dRMin){dRMin=dR; closestJetIndex=j;}
  }
  return dRMin;
}




int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
  TString outUrl = runProcess.getParameter<std::string>("outfile");

  //good lumi MASK
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false), filterOnlyPhoton(false), filterOnlyE(false), filterOnlyMU(false);
  if(!isMC){
    if(dtag.Contains("DoubleElectron"))   filterOnlyEE=true;
    if(dtag.Contains("DoubleMuon"))    filterOnlyMUMU=true;
    if(dtag.Contains("MuEG"))        filterOnlyEMU=true;
    if(dtag.Contains("SinglePhoton"))filterOnlyPhoton=true;
    if(dtag.Contains("SingleMuon"))    filterOnlyMU=true;
    if(dtag.Contains("SingleElectron"))filterOnlyE=true;
  }
  bool isV0JetsMC(false);//isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets")));  #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_ZZ2l2nu  = isMC && ( string(dtag.Data()).find("TeV_ZZ2l2nu")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  bool isMC_WZ3lnu  = isMC && ( string(dtag.Data()).find("TeV_WZ3lnu")  != string::npos);
  bool isMC_QCD = (isMC && dtag.Contains("QCD"));
  bool isMC_GJet = (isMC && dtag.Contains("GJet"));
  bool is2015data = (!isMC && dtag.Contains("2015"));
  bool is2015MC = (isMC && dtag.Contains("2015"));
  bool is2016data = (!isMC && dtag.Contains("2016"));
  bool is2016MC = (isMC);


  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");

  std::vector<string> jetVarNames = {"", "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"};

  size_t nvarsToInclude=varNames.size();


  //ELECTROWEAK CORRECTION WEIGHTS
  std::vector<std::vector<float>> ewkTable, ZZ_NNLOTable;
  if(isMC_ZZ2l2nu){
        ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
        ZZ_NNLOTable = ZZatNNLO::readFile_and_loadTable(dtag);
  }
  if(isMC_WZ3lnu){
        ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
  }


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  printf("Definition of plots");

  //event selection


  float ptbinsJets[] = {10, 20, 30, 40, 60, 80, 100, 125, 150, 175,250,350,1000};
  int ptbinsJetsN = sizeof(ptbinsJets)/sizeof(float)-1;
  mon.addHistogram( new TH1F( "wrtJetPt",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtJetPt_v2",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtLepPt",  ";Lep p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));


  //create a tree and related variables to save higgs candidate info for each cutIndex values
  unsigned int treeEventId;
  unsigned int treeLumiId;
  unsigned int treeRunId;
  int          treeZId;
  int          treeLeg1Id;
  float        treeLeg1ClosestJetPt;
  float        treeLeg1ClosestJetEta;
  float        treeLeg1ClosestJetPt_LepClean;
  float        treeLeg1ClosestJetEta_LepClean;
  float        treeLeg1Pt;
  float        treeLeg1Eta;
  float        treeLeg1Iso;
  int          treeLeg1LepID_loose;
  int          treeLeg1LepID_medium;
  int          treeLeg1LepID_tight;
  float        treeLeg1DeltaR;
  int          treeLeg1TauIsoLoose3;
  int          treeLeg1TauIsoMedium3;
  int          treeLeg1TauIsoLooseMVA;
  int          treeLeg1TauIsoMediumMVA;
  float        treeWeight;
  float        treeWeightNoLepSF;
  float        treeLepEff1;

  TFile *ofile=TFile::Open("out.root", "recreate");

  TTree* tree = new TTree("CandTree","CandTree");
  tree->Branch("eventId", &treeEventId , string("eventId/i" ).c_str());
  tree->Branch("lumiId" , &treeLumiId  , string("lumiId/i"  ).c_str());
  tree->Branch("runId"  , &treeRunId   , string("runId/i"   ).c_str());
  tree->Branch("zId",     &treeZId     , string("zId/I" ).c_str());
  tree->Branch("leg1Id" , &treeLeg1Id  , string("leg1Id/I"  ).c_str());
  tree->Branch("leg1ClosestJetPt" , &treeLeg1ClosestJetPt  , string("leg1ClosestJetPt/F"  ).c_str());
  tree->Branch("leg1ClosestJetEta" , &treeLeg1ClosestJetEta  , string("leg1ClosestJetEta/F"  ).c_str());
  tree->Branch("treeLeg1ClosestJetPt_LepClean" , &treeLeg1ClosestJetPt_LepClean  , string("treeLeg1ClosestJetPt_LepClean/F"  ).c_str());
  tree->Branch("treeLeg1ClosestJetEta_LepClean" , &treeLeg1ClosestJetEta_LepClean  , string("treeLeg1ClosestJetEta_LepClean/F"  ).c_str());
  tree->Branch("leg1Pt" , &treeLeg1Pt  , string("leg1Pt/F"  ).c_str());
  tree->Branch("leg1Eta", &treeLeg1Eta , string("leg1Eta/F" ).c_str());
  tree->Branch("leg1Iso", &treeLeg1Iso , string("leg1Iso/F" ).c_str());
  tree->Branch("leg1LepIDloose", &treeLeg1LepID_loose , string("leg1LepIDloose/I" ).c_str());
  tree->Branch("leg1LepIDmedium", &treeLeg1LepID_medium , string("leg1LepIDmedium/I" ).c_str());
  tree->Branch("leg1LepIDtight", &treeLeg1LepID_tight , string("leg1LepIDtight/I" ).c_str());
  tree->Branch("leg1DeltaR", &treeLeg1DeltaR , string("leg1DeltaR/F" ).c_str());
  tree->Branch("tau1Loose3", &treeLeg1TauIsoLoose3 , string("tau1Loose3/i" ).c_str());
  tree->Branch("tau1Medium3", &treeLeg1TauIsoMedium3 , string("tau1Medium3/i" ).c_str());
  tree->Branch("tau1LooseMVA", &treeLeg1TauIsoLooseMVA , string("tau1LooseMVA/i" ).c_str());
  tree->Branch("tau1MediumMVA", &treeLeg1TauIsoMediumMVA , string("tau1MediumMVA/i" ).c_str());
  tree->Branch("weight" , &treeWeight  , string("weight/F"  ).c_str());
  tree->Branch("weightNoLepSF" , &treeWeightNoLepSF  , string("weightNoLepSF/F"  ).c_str());
  tree->Branch("lepEff1" , &treeLepEff1  , string("lepEff1/F"  ).c_str());


  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  //MC normalization (to 1/pb)
  double xsecWeight = 1.0;
  if(isMC) xsecWeight=xsec/utils::getTotalNumberOfEvents(urls, false, true);//need to use the slow method in order to take NLO negative events into account

  //MET CORRection level
  pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;

  //jet energy scale and uncertainties
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  if (is2016data) {
      if(dtag.Contains("2016H")) jecDir+="Moriond17_80X/PromptReco_DATA/";
      else jecDir+="Moriond17_80X/2016SeptRepro_DATA/";
  } else jecDir+="Moriond17_80X/Summer16_re-digi-reco_MC/";

  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  TString pf(isMC ? "MC" : "DATA");
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/"+pf+"_Uncertainty_AK4PFchs.txt").Data());

  //muon energy scale and uncertainties
  TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
  gSystem->ExpandPathName(muscleDir);
  std::string roccorPath = string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/rcdata.2016.v3";
  RoccoR_Moriond17  *muCorMoriond17 = new RoccoR_Moriond17(roccorPath);

  std::string photonESC = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_74x_pho";//string(std::getenv("CMSSW_BASE"))+"/src/EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_74x_pho";
  // std::string photonSmearings = "../../../../EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_74x_pho_smearings.dat";
  // std::string correctionFile = isMC ? photonSmearings : photonScales;

  //--- photonESC initializes an EnergyScaleCorrection_class that provides the full file path
  PhotonEnergyCalibratorRun2 PhotonEnCorrector(isMC, false, photonESC);
  PhotonEnCorrector.initPrivateRng(new TRandom(1234));

  EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele");
  eScaler.doScale=true;
  eScaler.doSmearings=true;


  //lepton efficiencies
  LeptonEfficiencySF lepEff;

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


  // from Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
  // beff = 0.747; sfb = 0.899; //for Loose WP  //sfb is not actually used as it's taken from btagCal
  //
  //beff = 0.836; sfb = 0.920; //for Loose WP  //sfb is from page 7 https://indico.cern.ch/event/557018/contributions/2246312/attachments/1310986/1961665/csvSF_rwt_July18th_2016.pdf
  //leff = 0.139;
  beff = 0.827; sfb = 0.980; //for Loose WP  //sfb is from page 7 https://indico.cern.ch/event/557018/contributions/2246312/attachments/1310986/1961665/csvSF_rwt_July18th_2016.pdf
  leff = 0.132;

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
    std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
    std::vector<float> mcPileupDistribution;

    utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
    LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning

  patUtils::MetFilter metFilter;
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
  string debugText = "";

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0)
  std::shared_ptr<TRandom3> rgenMuon_(new TRandom3(0));
  std::shared_ptr<TRandom3> rgenEle_(new TRandom3(1));

  printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");
  for(unsigned int f=0;f<urls.size();f++){
    TFile* file = TFile::Open(urls[f].c_str() );
    fwlite::Event ev(file);
    printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)urls.size());
    int iev=0;
    int treeStep(ev.size()/50);
    for(ev.toBegin(); !ev.atEnd(); ++ev){ iev++;
      if(iev%treeStep==0){printf(".");fflush(stdout);}
      float weight = xsecWeight;
      float shapeWeight = 1.0;
      double puWeightUp = 1.0;
      double puWeightDown = 1.0;
      float puWeight(1.0);
      float weightNoLepSF(1.0);

      //##############################################   EVENT LOOP STARTS   ##############################################
      //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      //Skip bad lumi
      if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;
      treeEventId=0;
  		treeLumiId=0;
  		treeRunId=0;
  		treeZId=0;
  		treeLeg1Id=0;
  		treeLeg1ClosestJetPt=0;
  		treeLeg1ClosestJetEta=0;
  		treeLeg1Pt=0;
  		treeLeg1Eta=0;
  		treeLeg1Iso=0;
  		treeLeg1LepID_loose=0;
      treeLeg1LepID_medium=0;
  		treeLeg1LepID_tight=0;
  		treeLeg1DeltaR=0;
  		treeLeg1TauIsoLoose3=0;
  		treeLeg1TauIsoMedium3=0;
  		treeLeg1TauIsoLooseMVA=0;
  		treeLeg1TauIsoMediumMVA=0;
  		treeWeight=0;
  		treeLepEff1=0;



      reco::GenParticleCollection gen;
      GenEventInfoProduct eventInfo;
      int decayType = 0;
      int ZbosonType = -1;
      if(isMC){
        fwlite::Handle< reco::GenParticleCollection > genHandle;
        genHandle.getByLabel(ev, "prunedGenParticles");
        if(genHandle.isValid()){ gen = *genHandle;}

        fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
        genEventInfoHandle.getByLabel(ev, "generator");
        if(genEventInfoHandle.isValid()){ eventInfo = *genEventInfoHandle;}

        //WEIGHT for NLO negative interference
        weight *= eventInfo.weight();


        //WEIGHT for Pileup
        int ngenITpu = 0;
        fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
        puInfoH.getByLabel(ev, "slimmedAddPileupInfo");
        for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
          if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); } //getPU_NumInteractions();
        }

        if (ngenITpu == 0) continue; // It prevents to fill vtxraw with -nan values

        puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
        // if ( puWeight == 0 ){
        //   std::cout<<"  Some Problem in with ngenITpu= " << ngenITpu << "vtx size= "<<vtx.size() <<std::endl;
        // }
        puWeightUp  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
        puWeightDown = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
        weight *= puWeight;
      }

    weightNoLepSF = weight;

    //apply trigger and require compatibilitiy of the event with the PD
    edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
    if(!tr.isValid())return false;

    bool mumuTrigger(true); bool muTrigger(true);	bool eeTrigger(true); bool eTrigger(true); bool emuTrigger(true);

    int metFilterValue = 0;

	  bool filterbadPFMuon = true;
	  bool filterbadChCandidate = true;
	  bool filterbadMuonHIP = true;
	  bool filterduplicateMuonHIP = true;
	  std::unique_ptr<std::vector<reco::Muon*>> outbadMuon(new std::vector<reco::Muon*>());
	  std::unique_ptr<std::vector<reco::Muon*>> outduplicateMuon(new std::vector<reco::Muon*>());

    if(is2016data || is2016MC){

      mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
      muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu22_v*","HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
      eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //,"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",);
      eTrigger           = utils::passTriggerPatterns(tr, "HLT_Ele27_eta2p1_WPLoose_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*");
      emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");

      metFilterValue = metFilter.passMetFilterInt( ev, is2016data );

	    // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
	    filterbadChCandidate = metFilter.passBadChargedCandidateFilter(ev); if (!filterbadChCandidate) {  metFilterValue=9; }
	    filterbadPFMuon = metFilter.passBadPFMuonFilter(ev); if (!filterbadPFMuon) { metFilterValue=8; }
	    filterbadMuonHIP = metFilter.BadGlobalMuonTaggerFilter(ev,outbadMuon,false); if (!filterbadMuonHIP) { metFilterValue=10; }
	    filterduplicateMuonHIP = metFilter.BadGlobalMuonTaggerFilter(ev,outduplicateMuon,true); if (!filterduplicateMuonHIP) { metFilterValue=11; }
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
    if(!passTrigger)continue;
    //##############################################   EVENT PASSED THE TRIGGER   ######################################

	  if (metFilterValue==10 || metFilterValue==11) { metFilterValue=0; }
          if( metFilterValue!=0 ) continue;	 //Note this must also be applied on MC

	  // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
	  //	  if (!filterbadPFMuon || !filterbadChCandidate) continue;
    //##############################################   EVENT PASSED MET FILTER   #######################################


    vector<pat::TriggerObjectStandAlone> triggerObjects;
    fwlite::Handle< vector<pat::TriggerObjectStandAlone> > triggerObjectsHandle;
    triggerObjectsHandle.getByLabel(ev,"selectedPatTrigger");

    fwlite::Handle< edm::TriggerResults > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev,"TriggerResults","","HLT");

    const edm::TriggerNames &names = ev.triggerNames(*triggerBitsHandle);

    //load all the objects we will need to access
    reco::VertexCollection vtx;
    fwlite::Handle< reco::VertexCollection > vtxHandle;
    vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
    if(vtxHandle.isValid()){ vtx = *vtxHandle;}

    double rho = 0;
    fwlite::Handle< double > rhoHandle;
    rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
    if(rhoHandle.isValid()){ rho = *rhoHandle;}

    pat::MuonCollection muons;
    fwlite::Handle< pat::MuonCollection > muonsHandle;
    muonsHandle.getByLabel(ev, "slimmedMuons");
    if(muonsHandle.isValid()){ muons = *muonsHandle;}

    pat::ElectronCollection electrons;
    fwlite::Handle< pat::ElectronCollection > electronsHandle;
    electronsHandle.getByLabel(ev, "slimmedElectrons");
    if(electronsHandle.isValid()){ electrons = *electronsHandle;}

    EcalRecHitCollection recHitCollectionEB;
    EcalRecHitCollection recHitCollectionEE;
    fwlite::Handle<EcalRecHitCollection> recHitCollectionEBHandle;
    fwlite::Handle<EcalRecHitCollection> recHitCollectionEEHandle;
    recHitCollectionEBHandle.getByLabel(ev, "reducedEgamma","reducedEBRecHits" );
    recHitCollectionEEHandle.getByLabel(ev, "reducedEgamma","reducedEERecHits" );
    if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}
    if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

    pat::JetCollection jets;
    fwlite::Handle< pat::JetCollection > jetsHandle;
    jetsHandle.getByLabel(ev, "slimmedJets");
    if(jetsHandle.isValid()){ jets = *jetsHandle;}

    pat::METCollection mets;
    fwlite::Handle< pat::METCollection > metsHandle;
    metsHandle.getByLabel(ev, "slimmedMETs");
    if(metsHandle.isValid()){ mets = *metsHandle;}
    pat::MET met = mets[0];

    pat::TauCollection taus;
    fwlite::Handle< pat::TauCollection > tausHandle;
    tausHandle.getByLabel(ev, "slimmedTaus");
    if(tausHandle.isValid()){ taus = *tausHandle;}

    if(isV0JetsMC){
      fwlite::Handle< LHEEventProduct > lheEPHandle;
      lheEPHandle.getByLabel(ev, "externalLHEProducer");
      if(lheEPHandle.isValid()){
        mon.fillHisto("nup","",lheEPHandle->hepeup().NUP,1);
        if(lheEPHandle->hepeup().NUP>5) continue;
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
    // PHOTON ANALYSIS
    //
    pat::PhotonCollection selPhotons;

    //
    // LEPTON ANALYSIS
    //

    //start by merging electrons and muons
    std::vector<patUtils::GenericLepton> leptons;
    for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}
    for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}
    std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

    std::vector<patUtils::GenericLepton> selLeptons;
    LorentzVector muDiff(0,0,0,0);
    LorentzVector elDiff(0,0,0,0);
    for(size_t ilep=0; ilep<leptons.size(); ilep++){
      bool passKin(true);
      bool passIsoWPforFakeRate(true), passLooseLepton(true);
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
      //isolation
      // passLooseLepton
      passLooseLepton &= lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
                                   patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIsoWPforFakeRate &= lid==11 ?  patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut) :
                                        patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut);

  //     //apply muon corrections
  //     //if(abs(lid)==13 && passIso && passId)
      if(abs(lid)==13 && passIsoWPforFakeRate && passLooseLepton ){
          if(is2016MC || is2016data){
	           if(muCorMoriond17){
	              muDiff -= leptons[ilep].p4();

	               float qter =1.0;
                 double pt  = leptons[ilep].pt();
                 double eta = leptons[ilep].eta();
                 double phi = leptons[ilep].phi();
                 int charge = leptons[ilep].charge();
                 TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
                 // cout<<"PT Befor Correction: "<< p4.Pt() << endl;
                 int ntrk = leptons[ilep].mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
                 if(is2016MC){
                   //muCor2016->momcor_mc  (p4, lid<0 ? -1 :1, ntrk, qter);
                   //TRandom3 *rgen_ = new TRandom3(0);
                   double u1 = rgenMuon_->Uniform();
                   double u2 = rgenMuon_->Uniform();

                   double mcSF = muCorMoriond17->kScaleAndSmearMC(charge, pt, eta, phi, ntrk, u1, u2, 0, 0);

                   leptons[ilep].mu.setP4(LorentzVector(p4.Px()*mcSF,p4.Py()*mcSF,p4.Pz()*mcSF,p4.E()*mcSF ) );
                   leptons[ilep] = patUtils::GenericLepton(leptons[ilep].mu);

                   //cout<<"\t PT After Correction (2016): "<< p4.Pt() << " -- (Moriond17): "<< p4_2.Pt()<< endl;
                 }else if (is2016data){
                   // muCor2016->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
                   double dataSF = muCorMoriond17->kScaleDT(charge, pt, eta, phi, 0, 0);

                   leptons[ilep].mu.setP4(LorentzVector(p4.Px()*dataSF,p4.Py()*dataSF,p4.Pz()*dataSF,p4.E()*dataSF ) );
                   leptons[ilep] = patUtils::GenericLepton(leptons[ilep].mu);
                    // cout<<"\t PT After Correction: "<< leptons[ilep].pt() << endl;
                  }
	      muDiff += leptons[ilep].p4();
	    }
	  }
  }

      //apply electron corrections
      //if(abs(lid)==11  && passIso && passId)
      if(abs(lid)==11  && passIsoWPforFakeRate && passLooseLepton){
        //std::cout<<"START ---- "<<std::endl;
        elDiff -= leptons[ilep].p4();
        //const EcalRecHitCollection* recHits = (leptons[ilep].el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product();
        unsigned int gainSeed = patUtils::GainSeed(leptons[ilep].el, (leptons[ilep].el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product() );

        if(!isMC){

          double scale_corr=eScaler.ScaleCorrection(ev.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed);
          //At this point, the new data energy will be:
          // E_new=E_old*(scale_corr);
          TLorentzVector p4(leptons[ilep].el.px(),leptons[ilep].el.py(),leptons[ilep].el.pz(),leptons[ilep].el.energy());
          leptons[ilep].el.setP4(LorentzVector(p4.Px()*scale_corr,p4.Py()*scale_corr,p4.Pz()*scale_corr,p4.E()*scale_corr ) );
          leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
        }
        if(isMC){
          //std::cout<<"Before  pt  ---- "<<leptons[ilep].el.p4()<<std::endl;
          double sigma=eScaler.getSmearingSigma(ev.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed,0,0);
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

      //kinematics
      float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
      if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
      if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
      if(lid==13 && leptons[ilep].pt()<10) passKin=false;
      else if(lid==11 && leptons[ilep].pt()<10) passKin=false;

      if(passIsoWPforFakeRate && passKin)            selLeptons.push_back(leptons[ilep]); //we need loose lepton for FR
      }
      std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);

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

        // we need to apply a very loose selection here (Lucia's suggestion)
        if(!tau.tauID("againstElectronLooseMVA6")) continue;
        if(!tau.tauID("againstMuonLoose3")) continue;
        if(!tau.tauID("decayModeFinding")) continue;

        selTaus.push_back(tau);
        selLeptons.push_back(tau);
        ntaus++;
      }
      std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
      std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);


      if ( selLeptons.size() != 3 ) continue;
      /*******************************************************************************/
      /***    Cut on number of selLeptons     **/



      //
      //JET/MET ANALYSIS
      //

      //add scale/resolution uncertainties and propagate to the MET
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC);

      //select the jets
      pat::JetCollection selJets;
      std::map<string, int   > nbtagsVar;

      for(size_t ijet=0; ijet<jets.size(); ijet++){
        pat::Jet jet = jets[ijet]; //copy the jet, such that we can update it

        if(jet.pt()<15 || fabs(jet.eta())>4.7 ) continue;

        //jet id
        bool passPFloose = patUtils::passPFJetID("Loose", jet);
        bool passLooseSimplePuId = patUtils::passPUJetID(jet); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
        //bool passLooseSimplePuId = jet.userInt("pileupJetId:fullId") & (1 << 2);
        if(!passPFloose || !passLooseSimplePuId) continue;

        selJets.push_back(jet);
      }
      //sort all jet collection by pT
      std::sort(selJets.begin(), selJets.end(), utils::sort_CandidatesByPt);


        //
        // ASSIGN CHANNEL
        //

        std::vector<TString> chTags;
        TString evCat;
        int dilId(1);
        patUtils::GenericLepton dilLep1;
        patUtils::GenericLepton dilLep2;
        double BestMass;
        LorentzVector leadingLep, trailerLep, zll, zlltmp;
        //get the Z candidate
        dilId=-1;
        BestMass=0;
        zll = LorentzVector(0.,0.,0.,0.);

        for(auto& lep1: selLeptons){
          if(abs(lep1.pdgId())==15)continue;

          double leadPtCutValue  = abs(lep1.pdgId())==11 ? 24.0 : 18.0;
          if( lep1.pt()< leadPtCutValue ) continue;
          if(!( abs(lep1.pdgId())==11 ? patUtils::passIso(lep1.el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
                                                patUtils::passIso(lep1.mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut)) ||
             !( abs(lep1.pdgId())==11 ? patUtils::passId(lep1.el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
                                        patUtils::passId(lep1.mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut)) ) continue;

          for(auto& lep2: selLeptons){
            if ( &lep1 == &lep2 || abs(lep2.pdgId())==15)continue;

            double trailPtCutValue = abs(lep2.pdgId())==11 ? 13.0 : 10.0;
            if( lep2.pt() < trailPtCutValue ) continue;
            if(!( abs(lep2.pdgId())==11 ? patUtils::passIso(lep2.el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
                                          patUtils::passIso(lep2.mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut)) ||
               !( abs(lep2.pdgId())==11 ? patUtils::passId(lep2.el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
                                          patUtils::passId(lep2.mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut)) ) continue;


            if(abs(lep1.pdgId())!=abs(lep2.pdgId())) continue; 				 //SAME FLAVOUR PAIR
            if(lep1.pdgId()*lep2.pdgId()>=0) continue;					 //OPPOSITE SIGN

            zlltmp = (lep1.p4()+lep2.p4());
            if( fabs(zlltmp.mass() - 91.2) < fabs(zll.mass()-91.2) ){    //BEST MASS [76.2,106.2]
              dilLep1 = patUtils::GenericLepton(lep1);
              dilLep2 = patUtils::GenericLepton(lep2);
              zll=zlltmp;
              leadingLep=lep1.p4();
              trailerLep=lep2.p4();
              dilId = lep1.pdgId() * lep2.pdgId();
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
              weight *= lepEff.getRecoEfficiency( dilLep1.el.superCluster()->eta(), abs(dilLep1.pdgId())).first; //Reconstruction eff
              weight *= lepEff.getRecoEfficiency( dilLep2.el.superCluster()->eta(), abs(dilLep2.pdgId())).first; //Reconstruction eff
            }
            else if(abs(dilId)==169){
              weight *= lepEff.getTrackingEfficiency( dilLep1.eta(), abs( dilLep1.pdgId())).first; //Tracking eff
              weight *= lepEff.getTrackingEfficiency( dilLep2.eta(), abs( dilLep2.pdgId())).first; //Tracking eff

              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep1.pt(), dilLep1.eta(), abs(dilLep1.pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ISO w.r.t ID
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep2.pt(), dilLep2.eta(), abs(dilLep2.pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;
            }
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep1.pt(), dilLep1.eta(), abs( dilLep1.pdgId()),  abs( dilLep1.pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep2.pt(), dilLep2.eta(), abs( dilLep2.pdgId()),  abs( dilLep2.pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID

            // Trigger Eff
            if(isMC && abs(dilId)==169)weight *= lepEff.getTriggerEfficiencySF( dilLep1.pt(), dilLep1.eta(), dilLep2.pt(), dilLep2.eta(), dilId,is2016MC).first;
       		  if(isMC && abs(dilId)==121)weight *= lepEff.getTriggerEfficiencySF( dilLep1.pt(), dilLep1.el.superCluster()->eta(), dilLep2.pt(), dilLep2.el.superCluster()->eta(), dilId,is2016MC).first;
         }

        }


        bool passZmass = (fabs(zll.mass()-91.2)<30.0);
        bool passZpt   = (zll.pt()>20);

        if( !isDileptonCandidate && !passZmass ) continue;
        /************************* EVENT HAS a Z-like candidate ***************************/


        //LEPTON FAKE RATE ANALYSIS Z+1jets  (no systematics taken into account here)
        if( passZmass && (int)selLeptons.size()==3){  //Request exactly one Z + 1 additional lepton
          double tmass=-999;
          treeEventId   = ev.eventAuxiliary().event();
    			treeLumiId    = ev.eventAuxiliary().luminosityBlock();
    			treeRunId     = ev.eventAuxiliary().run();
          treeZId       = dilId;
          for(auto& fakeLepton: selLeptons){
            if( ( &fakeLepton == &dilLep1 ) || ( &fakeLepton == &dilLep2 ) ) continue;
            if( deltaR( fakeLepton,  dilLep1 )<0.1 || deltaR( fakeLepton, dilLep2)<0.1 ) continue;
            if(abs( fakeLepton.pdgId() ) == 11 || abs( fakeLepton.pdgId() ) == 13 || abs( fakeLepton.pdgId() ) == 15 ){
              tmass = TMath::Sqrt(2* fakeLepton.pt()*met.pt()*(1-TMath::Cos(deltaPhi(met.phi(), fakeLepton.phi()))));
            }
            if(abs(fakeLepton.pdgId())==11 || abs(fakeLepton.pdgId())==13 || abs(fakeLepton.pdgId())==15){
              int closestJetIndexL1=-1; double pTL1=-1; double etaL1=-1;
              double dRminL1 = closestJet(fakeLepton.p4(), selJets, closestJetIndexL1);
              bool jetFound = closestJetIndexL1>=0 && dRminL1<0.5;
              if(jetFound){pTL1=selJets[closestJetIndexL1].pt(); etaL1=abs(selJets[closestJetIndexL1].eta());}
              else{pTL1=fakeLepton.pt(); etaL1=abs(fakeLepton.eta());}
              treeLeg1ClosestJetPt   = pTL1;
        			treeLeg1ClosestJetEta  = etaL1;
              treeLeg1ClosestJetPt_LepClean   = jetFound ? pTL1  : -999;
        			treeLeg1ClosestJetEta_LepClean  = jetFound ? etaL1 : -999;
        			treeLeg1Pt             = fakeLepton.pt();
        			treeLeg1DeltaR         = dRminL1;

              treeLeg1Eta  = fakeLepton.eta();
        			treeLeg1Id   = fakeLepton.pdgId();
        			treeLeg1Iso  = (abs(treeLeg1Id)!=15)?patUtils::relIso(fakeLepton,rho):-999.;

              TString PartName = "FR_";//+chTags.at(1)+"_";
              if     (abs(fakeLepton.pdgId())==11)PartName += "El";
              else if(abs(fakeLepton.pdgId())==13)PartName += "Mu";
              else if(abs(fakeLepton.pdgId())==15)PartName += "Ta";
              else PartName+= abs(fakeLepton.pdgId());


              std::vector<TString> TagsFR;

              if(abs(fakeLepton.pdgId())==11 || abs(fakeLepton.pdgId())==13){
                bool passId = false;
                if(abs(fakeLepton.pdgId())==11) passId = patUtils::passId(fakeLepton.el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                if(abs(fakeLepton.pdgId())==13) passId = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                float relIso = patUtils::relIso(fakeLepton, rho);

                if(true                 )TagsFR.push_back(PartName);
                if(passId && relIso<=0.1){
                  TagsFR.push_back(PartName+("_Id_Iso01"));
                  treeLeg1LepID_tight = 1;
                }
                if(passId && relIso<=0.2){
                  TagsFR.push_back(PartName+("_Id_Iso02"));
                  treeLeg1LepID_medium = 1;
                }
                if(passId && relIso<=0.3){
                  TagsFR.push_back(PartName+("_Id_Iso03"));
                  treeLeg1LepID_loose = 1;
                }
              }else{
                bool IdL         = fakeLepton.tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
                bool IdM         = fakeLepton.tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
                bool IdL_MVA     = fakeLepton.tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
                bool IdM_MVA     = fakeLepton.tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
                bool IdL_MVA_R03 = fakeLepton.tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
                bool IdM_MVA_R03 = fakeLepton.tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");

                if(true                 )TagsFR.push_back(PartName);
                if(IdL                  )TagsFR.push_back(PartName+("_Id_IsoLo"));
                if(IdM                  )TagsFR.push_back(PartName+("_Id_IsoMe"));
                if(IdL_MVA              ){
                  TagsFR.push_back(PartName+("_Id_IsoLo_MVA"));
                  treeLeg1TauIsoLooseMVA = 1;
                }
                if(IdM_MVA              ){
                  TagsFR.push_back(PartName+("_Id_IsoMe_MVA"));
                  treeLeg1TauIsoMediumMVA = 1;
                }
                if(IdL_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoLo_MVAR03"));
                if(IdM_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoMe_MVAR03"));
              }

              if(tmass<30){
                int NTags = TagsFR.size();
                for(unsigned int iTags=0;iTags<NTags;iTags++){
                  TagsFR.push_back(TagsFR[iTags] + TString("_TMCut"));
                }
              }

              auto TagsFRJet = TagsFR;
              auto TagsFRLep = TagsFR;

              for(unsigned int iTags=0;iTags<TagsFR.size();iTags++){
                TagsFRJet.push_back(TagsFR[iTags] + (etaL1<1.4                   ?TString("_B"):TString("_E")));
                TagsFRLep.push_back(TagsFR[iTags] + (abs(fakeLepton.eta())<1.4?TString("_B"):TString("_E")));
              }

              mon.fillHisto("wrtJetPt", TagsFRJet, pTL1              , weight);
              if(closestJetIndexL1>=0 && dRminL1<0.5) mon.fillHisto("wrtJetPt_v2", TagsFRJet, pTL1              , weight);
              mon.fillHisto("wrtLepPt", TagsFRLep, fakeLepton.pt(), weight);
            }
          }//close loop on leptons
          treeWeightNoLepSF = weightNoLepSF;
          treeWeight = weight;
          if(tree) tree->Fill();
        }//close FR study Zmass

    }
    printf("\n");
    delete file;
  }

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  TString terminationCmd = "";
  //save control plots to file
  printf("Results save in local directory and moved to %s\n", outUrl.Data());

  //save all to the file
  terminationCmd += TString("mv out.root ") + outUrl + ";";
  //TFile
  ofile->cd();
  mon.Write();
  if(tree){tree->SetDirectory(ofile); tree->Write();}
  ofile->Close();

  if(!isMC && debugText!=""){
    TString outTxtUrl= outUrl + ".txt";
    terminationCmd += TString("mv out.txt ") + outTxtUrl + ";";
    FILE* outTxtFile = fopen("out.txt", "w");
    fprintf(outTxtFile, "%s", debugText.c_str());
    printf("TextFile URL = %s\n",outTxtUrl.Data());
    if(outTxtFile)fclose(outTxtFile);
  }

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
    terminationCmd += TString("mv out.json ") + ((outUrl.ReplaceAll(".root",""))+".json") + ";";
    goodLumiFilter.FindLumiInFiles(urls);
    goodLumiFilter.DumpToJson("out.json");
  }

  system(terminationCmd.Data());

}
