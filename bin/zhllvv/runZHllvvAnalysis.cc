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

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"
#include "UserCode/llvv_fwk/interface/muresolution_run2.h"
#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/TrigUtils.h"
#include "UserCode/llvv_fwk/interface/EwkCorrections.h"
#include "UserCode/llvv_fwk/interface/ZZatNNLO.h"

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

// Additional functions


LorentzVector getSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2){
  if(higgsCandL1<0 || higgsCandL2<0) return LorentzVector(0,0,0,0);

  TMatrixD covMET(2, 2); // PFMET significance matrix
// FIXME in MINIAODv2 74X, covariance is always 0000
//  covMET[0][0] = met.getSignificanceMatrix()(0,0);
//  covMET[0][1] = met.getSignificanceMatrix()(0,1);
//  covMET[1][0] = met.getSignificanceMatrix()(1,0);
//  covMET[1][1] = met.getSignificanceMatrix()(1,1);
//  std::cout<<"MET MATRIX: " << covMET[0][0] << " " << covMET[0][1] << " " << covMET[1][0] << " " << covMET[1][1] << "\n";

  covMET[0][0] = 0.95;  covMET[0][1] = 0.05; covMET[1][0] = 0.05; covMET[1][1] = 0.95;

  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL1].pdgId())==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL1].pdgId())==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay,selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), selLeptons[higgsCandL1].phi(), selLeptons[higgsCandL1].mass() ));
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL2].pdgId())==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL2].pdgId())==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay, selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), selLeptons[higgsCandL2].phi(), selLeptons[higgsCandL2].mass() ));

  SVfitStandaloneAlgorithm algo(measuredTauLeptons, met.px(), met.py() , covMET, 0);
  algo.addLogM(false);
  algo.fit();
  if(algo.isValidSolution()){
 //   return algo.fittedDiTauSystem();
  }
  return LorentzVector(selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());
}

//**********************************************************************************************//
bool passHiggsCuts(std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2, float isoElCut, float isoMuCut, const char* isoHaCut, float sumPtCut, bool requireId, reco::VertexCollection vtx)
//**********************************************************************************************//
{
  if(higgsCandL1<0 || higgsCandL2<0)return false;
  std::vector<patUtils::GenericLepton*> HiggsLegs = {&(selLeptons[higgsCandL1]), &(selLeptons[higgsCandL2])};

  bool passId=true;
  bool passIso=true;
  float sumpt = 0;
  for(auto lepIt=HiggsLegs.begin();lepIt!=HiggsLegs.end();lepIt++){
     patUtils::GenericLepton* lep = (*lepIt);
     if(abs(lep->pdgId())==11){
       passId  &= patUtils::passId(lep->el, vtx[0], patUtils::llvvElecId::Loose,patUtils::CutVersion::ICHEP16Cut);
       passIso &= (lep->userFloat("relIso") <= isoElCut);
     }else if(abs(lep->pdgId())==13){
       passId  &= patUtils::passId(lep->mu, vtx[0], patUtils::llvvMuonId::Loose,patUtils::CutVersion::ICHEP16Cut);
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

//**********************************************************************************************//
std::vector<patUtils::GenericLepton> getLepVariations(  std::vector<patUtils::GenericLepton>& selLeptons, float factor)
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
std::vector<patUtils::GenericLepton> getTauVariations( std::vector<patUtils::GenericLepton>& selLeptons,float factor)
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

int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

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
      if(dtag.Contains("DoubleEle"))   filterOnlyEE=true;
      if(dtag.Contains("DoubleMu"))    filterOnlyMUMU=true;
      if(dtag.Contains("MuEG"))        filterOnlyEMU=true;
      if(dtag.Contains("SinglePhoton"))filterOnlyPhoton=true;
      if(dtag.Contains("SingleMu"))    filterOnlyMU=true;
      if(dtag.Contains("SingleElectron"))filterOnlyE=true;
  }
  bool isV0JetsMC(false);//isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets")));  #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_ZZ2l2nu  = isMC && ( string(dtag.Data()).find("TeV_ZZ2l2nu")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  bool isMC_QCD = (isMC && dtag.Contains("QCD"));
  bool isMC_GJet = (isMC && dtag.Contains("GJet"));
  bool is2015data = (!isMC && dtag.Contains("2015"));
  bool is2015MC = (isMC && dtag.Contains("2015"));
  bool is2016data = (!isMC && dtag.Contains("2016"));
  bool is2016MC = (isMC && dtag.Contains("2016"));

  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");

  std::vector<string> jetVarNames = {""}; //, "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"};

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
     if(isMC_ZZ){
        varNames.push_back("_th_ewkup"); varNames.push_back("_th_ewkdown"); //EWK+QCD corrections
     }
  }
  size_t nvarsToInclude=varNames.size();


  //ELECTROWEAK CORRECTION WEIGHTS
  std::vector<std::vector<float>> ewkTable, ZZ_NNLOTable;
  if(isMC_ZZ2l2nu){
     ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
     ZZ_NNLOTable = ZZatNNLO::readFile_and_loadTable(dtag);
  }

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  printf("Definition of plots");

  // zhllvv Event selection
  TH1 *h1zllvv=mon.addHistogram( new TH1F ("eventflow_zhllvv", ";;Events", 10,0,10) );
  h1zllvv->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1zllvv->GetXaxis()->SetBinLabel(2,"Nlep#geq2");
  h1zllvv->GetXaxis()->SetBinLabel(3,"ZmassVeto");
  h1zllvv->GetXaxis()->SetBinLabel(4,"jetCount");
  h1zllvv->GetXaxis()->SetBinLabel(5,"TopVeto");
  h1zllvv->GetXaxis()->SetBinLabel(6,"EtMiss");
  h1zllvv->GetXaxis()->SetBinLabel(7,"#Delta #phi ll-MET");
  h1zllvv->GetXaxis()->SetBinLabel(8,"pTunbalance");

  TH1 *h_tr= mon.addHistogram( new TH1F ("trigger", ";;Events", 10,0,10) );
  h_tr->GetXaxis()->SetBinLabel(1,"#mu#mu");
  h_tr->GetXaxis()->SetBinLabel(2,"#mu");
  h_tr->GetXaxis()->SetBinLabel(3,"ee");
  h_tr->GetXaxis()->SetBinLabel(4,"e");
  h_tr->GetXaxis()->SetBinLabel(5,"e#mu");
  h_tr->GetXaxis()->SetBinLabel(6,"#gamma");

  //pu control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) );
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) );
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) );

  //lepton control
  mon.addHistogram( new TH1F( "leadpt",  	";p_{T}^{#ell} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta", 	";#eta_{#ell};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "trailerpt",  	";p_{T}^{#ell} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", 	";#eta_{#ell};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "leppt",  	";p_{T}^{#ell} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "lepeta", 	";#eta_{#ell};Events", 50,-2.6,2.6) );

  //tau control
  mon.addHistogram( new TH1F( "ntaus",      	";Number of Taus;Events", 10,0,10) );
  mon.addHistogram( new TH1F( "tauleadpt",  	";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "tauleadeta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "tautrailerpt",  	";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "tautrailereta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "taupt",  	";p_{T}^{#tau} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "taueta", 	";#eta_{#tau};Events", 50,-2.6,2.6) );

  //extra leptons in the event

  TH1 *hbtags=mon.addHistogram( new TH1F("nbtags",   ";b-tag multiplicity;Events",5,0,5) );
  TH1 *hjets=mon.addHistogram( new TH1F("njets",  ";Jet multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=hjets->GetXaxis()->GetNbins(); ibin++){
      TString label("");
      if(ibin==hjets->GetXaxis()->GetNbins()) label +="#geq";
      else                                    label +="=";
      label += (ibin-1);
      hjets->GetXaxis()->SetBinLabel(ibin,label);
      hbtags->GetXaxis()->SetBinLabel(ibin,label);
    }

  //fake rate histograms
  float ptbinsJets[] = {10, 20, 30, 40, 60, 80, 100, 125, 150, 175,250};
  int ptbinsJetsN = sizeof(ptbinsJets)/sizeof(float)-1;
  mon.addHistogram( new TH1F( "wrtJetPt",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtLepPt",  ";Lep p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));


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
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  TString pf(isMC ? "MC" : "DATA");
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/"+pf+"_Uncertainty_AK4PFchs.txt").Data());

  //muon energy scale and uncertainties
  TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
  gSystem->ExpandPathName(muscleDir);
  rochcor2015* muCor = new rochcor2015();  //replace the MuScleFitCorrector we used at run1

  /*
    New energy scale comment out
  */
  // ElectronEnergyCalibratorRun2 ElectronEnCorrector;
  // string EGammaEnergyCorrectionFile = "";
  // PhotonEnergyCalibratorRun2 PhotonEnCorrector;
  // EpCombinationTool theEpCombinationTool;
  // if(is2015MC || is2015data){
  // 	EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015";
  // 	PhotonEnCorrector = PhotonEnergyCalibratorRun2(isMC, false, EGammaEnergyCorrectionFile);
  // 	PhotonEnCorrector.initPrivateRng(new TRandom(1234));
  // 	theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC
  // 	ElectronEnCorrector = ElectronEnergyCalibratorRun2(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
  // 	ElectronEnCorrector.initPrivateRng(new TRandom(1234));
  // }
  // if(is2016MC || is2016data){
  // 	string EleEnergyCorrectionFile = "UserCode/llvv_fwk/data/jec/80X_ichepV1_2016_ele";
  // 	string PhoEnergyCorrectionFile = "UserCode/llvv_fwk/data/jec/80X_ichepV1_2016_pho";
  // 	PhotonEnCorrector = PhotonEnergyCalibratorRun2(isMC, false, PhoEnergyCorrectionFile);
  // 	PhotonEnCorrector.initPrivateRng(new TRandom(1234));
  // 	theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC
  // 	ElectronEnCorrector = ElectronEnergyCalibratorRun2(theEpCombinationTool, isMC, false, EleEnergyCorrectionFile);
  // 	ElectronEnCorrector.initPrivateRng(new TRandom(1234));
  // }

  EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/ScalesSmearings//Moriond17_23Jan_ele");
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

  double btagLoose = 0.605; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
  // setup calibration readers
  BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
  BTagCalibrationReader btagCal   (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalUp (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "up"     );  // sys up
  BTagCalibrationReader btagCalDn (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "down"   );  // sys down
  BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_LOOSE, "comb", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_LOOSE, "comb", "up"     );  // sys up
  BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_LOOSE, "comb", "down"   );  // sys down

  // from Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
  beff = 0.747; sfb = 0.899; //for Loose WP  //sfb is not actually used as it's taken from btagCal

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

          //##############################################   EVENT LOOP STARTS   ##############################################
          //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

         //Skip bad lumi
         if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;

         reco::GenParticleCollection gen;
         GenEventInfoProduct eventInfo;
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
            puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
            puWeightUp  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
            puWeightDown = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
            weight *= puWeight;

            //GEN LEVEL FILTERING
            if(isMC && (mctruthmode==15 || mctruthmode==1113)){// && (string(dtag.Data()).find("Z#rightarrow")==0 || isMC_ZZ2l2nu))
                int prodId = 1;
                for( unsigned int k=0; k<gen.size(); ++k){
                        if( gen[k].isHardProcess() && ( abs( gen[k].pdgId() ) == 11 || abs( gen[k].pdgId() ) == 13 || abs( gen[k].pdgId() )==15 ) ) prodId*=gen[k].pdgId();
                }
                if(mctruthmode==15   && abs(prodId)!=225)continue; //skip not tautau
                if(mctruthmode==1113 && abs(prodId)==225)continue; //skip tautau
            } else if (isMC && mctruthmode==121416 ){
              int prodId = 1;
              for( unsigned int k=0; k<gen.size(); ++k){
                if( gen[k].isHardProcess() && ( abs( gen[k].pdgId() ) == 12 || abs( gen[k].pdgId() ) == 14 || abs( gen[k].pdgId() )==16 ) ) prodId*=gen[k].pdgId();
              }
              if(abs(prodId)!=144  && abs(prodId)!=196  && abs(prodId)!=256 )continue; //skip not nunu
            }
          }

          //apply trigger and require compatibilitiy of the event with the PD
          edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
          if(!tr.isValid())return false;

	  bool mumuTrigger(true); 
	  bool muTrigger(true);	
	  bool eeTrigger(true); 
	  bool eTrigger(true); 
	  bool emuTrigger(true);

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

          pat::METCollection puppimets;
          fwlite::Handle< pat::METCollection > puppimetsHandle;
          puppimetsHandle.getByLabel(ev, "slimmedMETsPuppi");
          if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
          LorentzVector puppimet = puppimets[0].p4();

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

  	 //Electroweak corrections to ZZ and WZ(soon) simulations
     	 double ewkCorrectionsWeight = 1.;
     	 double ewkCorrections_error = 0.;
     	 if(isMC_ZZ2l2nu) ewkCorrectionsWeight = EwkCorrections::getEwkCorrections(dtag, gen, ewkTable, eventInfo, ewkCorrections_error);
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
             bool passKin(true),passId(true),passIso(true);
             bool passLooseLepton(true), passSoftMuon(true), passSoftElectron(true), passVetoElectron(true);
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
             passId           = lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight,patUtils::CutVersion::ICHEP16Cut);
             passLooseLepton &= lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose,patUtils::CutVersion::ICHEP16Cut) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose,patUtils::CutVersion::ICHEP16Cut);
             passSoftMuon &= lid==11? false : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Soft,patUtils::CutVersion::ICHEP16Cut);

             //isolation
             passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight,patUtils::CutVersion::ICHEP16Cut,0) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight,patUtils::CutVersion::ICHEP16Cut);
             passLooseLepton &= lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Loose,patUtils::CutVersion::ICHEP16Cut,0) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Loose,patUtils::CutVersion::ICHEP16Cut);

             //apply muon corrections
             if(abs(lid)==13 && passIso && passId){
                 passSoftMuon=false;
                 if(muCor){
                   float qter;
                   TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
                   if(isMC){muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
                   }else{   muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
                   }

                   muDiff -= leptons[ilep].p4();
                   leptons[ilep].setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
                   muDiff += leptons[ilep].p4();
                 }
               }

             //apply electron corrections
	     if(abs(lid)==11  && passIso && passId){
	       //std::cout<<"START ---- "<<std::endl;
	       elDiff -= leptons[ilep].p4();
	       const EcalRecHitCollection* recHits = (leptons[ilep].el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product();
	       unsigned int gainSeed = patUtils::GainSeed(leptons[ilep].el,recHits);

	       if(!isMC){

		 double scale_corr=eScaler.ScaleCorrection(ev.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed);
		 //At this point, the new data energy will be:
		 // E_new=E_old*(scale_corr);
		 TLorentzVector p4(leptons[ilep].el.px(),leptons[ilep].el.py(),leptons[ilep].el.pz(),leptons[ilep].el.energy());
		 leptons[ilep].el.setP4(LorentzVector(p4.Px()*scale_corr,p4.Py()*scale_corr,p4.Pz()*scale_corr,p4.E()*scale_corr ) );
		 //leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
	       }
	       if(isMC){
		 //std::cout<<"Before  pt  ---- "<<leptons[ilep].el.p4()<<std::endl;
		 double sigma=eScaler.getSmearingSigma(ev.eventAuxiliary().run(),leptons[ilep].el.isEB(),leptons[ilep].el.r9(), leptons[ilep].el.superCluster()->eta(), leptons[ilep].el.et(),gainSeed,0,0);
		 //Put the last two inputs at 0,0 for the nominal value of sigma
		 //Now smear the MC energy
		 TRandom3 *rgen_ = new TRandom3(0);
		 double smearValue = rgen_->Gaus(1, sigma) ;
		 //std::cout<<"smearing  ---- "<<smearValue<<std::endl;
		 TLorentzVector p4(leptons[ilep].el.px(),leptons[ilep].el.py(),leptons[ilep].el.pz(),leptons[ilep].el.energy());
		 leptons[ilep].el.setP4(LorentzVector(p4.Px()*smearValue,p4.Py()*smearValue,p4.Pz()*smearValue,p4.E()*smearValue ) );
		 //std::cout<<"After  pt  ---- "<<leptons[ilep].el.p4()<<std::endl;
		 // std::cout<<"\n";
		 //leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok

	       }
         elDiff += leptons[ilep].p4();
	     }

	     // conput the relative isolation on the momentum-corrected leptons
	     leptons[ilep].addUserFloat("relIso",  patUtils::relIso(leptons[ilep], rho) ); //compute it once for all

              //kinematics
             float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
             if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
             if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
             passLooseLepton &= passKin;
             passSoftMuon    &= passKin;
             if(lid==13){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
               if(leptons[ilep].pt()<3)  passSoftMuon=false;
             }else if(lid==11){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
             }
             if(leptons[ilep].pt()<25) passKin=false;

             //if(passId && passIso && passKin)        selLeptons.push_back(leptons[ilep]);
             if(passLooseLepton && passKin)            selLeptons.push_back(leptons[ilep]); //we need loose lepton for FR
             else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);
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
              if(deltaR(tau, selLeptons[l1])<0.1){overlapWithLepton=true; break;}
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
           ntaus++;
         }
         std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
         std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);


         //
         //JET/MET ANALYSIS
         //

         //add scale/resolution uncertainties and propagate to the MET
         utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC);

         //select the jets
         std::map<string, pat::JetCollection> selJetsVar;
         std::map<string, int   > njetsVar;
         std::map<string, int   > nbtagsVar;
         std::map<string, double> mindphijmetVar;
         for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){mindphijmetVar[jetVarNames[ivar]] = 9999.0;}  //initialize

         for(size_t ijet=0; ijet<jets.size(); ijet++){
             pat::Jet jet = jets[ijet]; //copy the jet, such that we can update it

             if(jet.pt()<15 || fabs(jet.eta())>4.7 ) continue;

             //mc truth for this jet
             const reco::GenJet* genJet=jet.genJet();
             TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );

             //cross-clean with selected leptons and photons  (DISABLED AS WE NEED THOSE FOR FR STUDY)
             //double minDRlj(9999.); for(size_t ilep=0; ilep<selLeptons.size(); ilep++){if(abs(selLeptons[ilep].pdgId())>13){continue;}  minDRlj = TMath::Min( minDRlj, deltaR(jet,selLeptons[ilep]) );}  //ignore taus for the cross-cleaning
             //double minDRlg(9999.); for(size_t ipho=0; ipho<selPhotons.size(); ipho++)  minDRlg = TMath::Min( minDRlg, deltaR(jet,selPhotons[ipho]) );
             //if(minDRlj<0.4 || minDRlg<0.4) continue;

             //jet id
             bool passPFloose = patUtils::passPFJetID("Loose", jet);
             bool passLooseSimplePuId = true; //patUtils::passPUJetID(jet); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
             if(jet.pt()>30){
                 mon.fillHisto(jetType,"",fabs(jet.eta()),0);
                 if(passPFloose)                        mon.fillHisto("jetId", jetType,fabs(jet.eta()),1);
                 if(passLooseSimplePuId)                mon.fillHisto("jetId", jetType,fabs(jet.eta()),2);
                 if(passPFloose && passLooseSimplePuId) mon.fillHisto("jetId", jetType,fabs(jet.eta()),3);
             }
             if(!passPFloose || !passLooseSimplePuId) continue;


            //check for btagging
            if(jet.pt()>30 && fabs(jet.eta())<2.5){
              bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>btagLoose);
              bool hasCSVtagUp = hasCSVtag;
              bool hasCSVtagDown = hasCSVtag;

              //update according to the SF measured by BTV
              if(isMC){
                  int flavId=jet.partonFlavour();  double eta=jet.eta();
                  if      (abs(flavId)==5){  btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
                  }else if(abs(flavId)==4){  btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
                  }else{                     btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCalL  .eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
                                             btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
                  }
              }

              if(hasCSVtag    )nbtagsVar[""          ]++;
              if(hasCSVtagUp  )nbtagsVar["_eff_bup"  ]++;
              if(hasCSVtagDown)nbtagsVar["_eff_bdown"]++;
            }


            for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){
               pat::Jet varJet = jet;
               if(ivar!=0) varJet.setP4(jet.p4() * jet.userFloat(jetVarNames[ivar]));
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
         std::pair<double, double> scaleUncVar = patUtils::scaleVariation(ev);  //compute it only once

         for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
          if(!isMC && ivar>0 ) continue; //loop on variation only for MC samples

          //start from a nominal
          float weight = initialWeight;

           //Theoretical Uncertanties: PDF, Alpha and Scale
           if(varNames[ivar]=="_th_factup")     weight *= std::max(0.9, std::min(scaleUncVar.first , 1.1));
           if(varNames[ivar]=="_th_factdown")   weight *= std::max(0.9, std::min(scaleUncVar.second, 1.1));
           if(varNames[ivar]=="_th_alphas")     weight *= patUtils::alphaVariation(ev);
           if(varNames[ivar]=="_th_pdf")        weight *= patUtils::pdfVariation(ev);

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
            double BestMass;
            LorentzVector leadingLep, trailerLep, zll, zlltmp;
            //get the Z candidate
              dilLep1=-1; dilLep2=-1; dilId=-1;
              BestMass=0;
              zll = LorentzVector(0.,0.,0.,0.);

              for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
                if(abs(selLeptons[l1].pdgId())==15)continue;
                if( selLeptons[l1].pt()<20 ) continue;
                if(!(abs(selLeptons[l1].pdgId())==11?patUtils::passIso(selLeptons[l1].el,  patUtils::llvvElecIso::Tight,patUtils::CutVersion::ICHEP16Cut, 0) : patUtils::passIso(selLeptons[l1].mu,  patUtils::llvvMuonIso::Tight,patUtils::CutVersion::ICHEP16Cut))) continue;

                for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
                  if(abs(selLeptons[l2].pdgId())==15)continue;
                  if( selLeptons[l2].pt()<20 ) continue;
                  if(!(abs(selLeptons[l2].pdgId())==11?patUtils::passIso(selLeptons[l2].el,  patUtils::llvvElecIso::Tight,patUtils::CutVersion::ICHEP16Cut,0) : patUtils::passIso(selLeptons[l2].mu,  patUtils::llvvMuonIso::Tight,patUtils::CutVersion::ICHEP16Cut))) continue;

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

            bool isDileptonCandidate = false;
            if(dilId!=-1){
                //check the channel
                if( abs(dilId)==121){ chTags.push_back("ll"); chTags.push_back("ee");   isDileptonCandidate=true; }
                if( abs(dilId)==169){ chTags.push_back("ll"); chTags.push_back("mumu"); isDileptonCandidate=true; }

                weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].pdgId()),  abs(selLeptons[dilLep1].pdgId()) ==11 ? "tight"    : "tight"    , patUtils::CutVersion::ICHEP16Cut).first : 1.0; //ID
                weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].pdgId()),  abs(selLeptons[dilLep1].pdgId()) ==11 ? "tightiso" : "tightiso",patUtils::CutVersion::ICHEP16Cut).first : 1.0; //ISO w.r.t ID
                weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].pdgId()),  abs(selLeptons[dilLep2].pdgId()) ==11 ? "tight"    : "tight",patUtils::CutVersion::ICHEP16Cut).first : 1.0; //ID
                weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].pdgId()),  abs(selLeptons[dilLep2].pdgId()) ==11 ? "tightiso" : "tightiso",patUtils::CutVersion::ICHEP16Cut).first : 1.0; //ISO w.r.t ID
            }

            std::vector<TString> tags(1,"all");
            for(size_t ich=0; ich<chTags.size(); ich++){
              tags.push_back( chTags[ich] );
            }

            if(!isDileptonCandidate) continue;
            bool passZmass = (fabs(zll.mass()-91.2)<15);
	    bool passZpt   = (zll.pt()>20);

	    bool passZVeto  = (fabs(zll.mass()-91.2)>15);
	    bool EtMiss = (met.corP4(metcor).E()>120);
	    bool jetCount = (njets<2);

	    LorentzVector diLept = (selLeptons[0].p4()+selLeptons[1].p4());
	    float dphillmet=fabs(deltaPhi(met.corP4(metcor).phi(),diLept.phi()));
	    float pTU = fabs( met.corP4(metcor).E() -  diLept.pt() ) / diLept.pt();

	    bool DeltaPhillMet  = (dphillmet > 2.7);
	    bool pTunbalance = (pTU < 0.25);

            //bool passLLpt   = (zll.pt()>50);

            bool passMass = passZmass;
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
            LorentzVector higgsCand_SVFit;
            LorentzVector higgsCandH;
            LorentzVector higgsCandH_SVFit;

            //LEPTON FAKE RATE ANALYSIS Z+1jets  (no systematics taken into account here)
            if(ivar==0 && passZmass && (int)selLeptons.size()==3){  //Request exactly one Z + 1 additional lepton
               bool IdentifiedThirdLepton=false;
               double tmass=-999;
                for(int i=0   ;i<(int)selLeptons.size() && !IdentifiedThirdLepton;i++){
                  if((i==dilLep1) || (i==dilLep2)) continue;
                  if(deltaR(selLeptons[i],  selLeptons[dilLep1])<0.1 || deltaR(selLeptons[i],  selLeptons[dilLep2])<0.1)continue;
                  if(abs(selLeptons[i].pdgId())==11||abs(selLeptons[i].pdgId())==13||abs(selLeptons[i].pdgId())==15){
                    tmass = TMath::Sqrt(2*selLeptons[i].pt()*met.pt()*(1-TMath::Cos(deltaPhi(met.phi(), selLeptons[i].phi()))));
                  }
                  if(abs(selLeptons[i].pdgId())==11 || abs(selLeptons[i].pdgId())==13 || abs(selLeptons[i].pdgId())==15){
                    int closestJetIndexL1=-1; double pTL1=-1; double etaL1=-1;
                    double dRminL1 = closestJet(selLeptons[i].p4(), selJets, closestJetIndexL1);
                    if(closestJetIndexL1>=0 && dRminL1<0.5){pTL1=selJets[closestJetIndexL1].pt(); etaL1=abs(selJets[closestJetIndexL1].eta());}
                    else{pTL1=selLeptons[i].pt(); etaL1=abs(selLeptons[i].eta());}


                    TString PartName = "FR_";
                    if     (abs(selLeptons[i].pdgId())==11)PartName += "El";
                    else if(abs(selLeptons[i].pdgId())==13)PartName += "Mu";
                    else if(abs(selLeptons[i].pdgId())==15)PartName += "Ta";
                    else PartName+= abs(selLeptons[i].pdgId());


                    std::vector<TString> TagsFR;

                    if(abs(selLeptons[i].pdgId())==11 || abs(selLeptons[i].pdgId())==13){
                       bool passId = false;
                       if(abs(selLeptons[i].pdgId())==11) passId = patUtils::passId(selLeptons[i].el, vtx[0], patUtils::llvvElecId::Loose,patUtils::CutVersion::ICHEP16Cut);
                       if(abs(selLeptons[i].pdgId())==13) passId = patUtils::passId(selLeptons[i].mu, vtx[0], patUtils::llvvMuonId::Loose,patUtils::CutVersion::ICHEP16Cut);
                       float relIso = patUtils::relIso(selLeptons[i], rho);

                       if(true                 )TagsFR.push_back(PartName);
                       if(passId && relIso<=0.1)TagsFR.push_back(PartName+("_Id_Iso01"));
                       if(passId && relIso<=0.2)TagsFR.push_back(PartName+("_Id_Iso02"));
                       if(passId && relIso<=0.3)TagsFR.push_back(PartName+("_Id_Iso03"));

                        if(passId && relIso<=0.3)IdentifiedThirdLepton=true;
                    }else{
                       bool IdL = selLeptons[i].tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
                       bool IdM = selLeptons[i].tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");

                       if(true                 )TagsFR.push_back(PartName);
                       if(IdL                  )TagsFR.push_back(PartName+("_Id_IsoLo"));
                       if(IdM                  )TagsFR.push_back(PartName+("_Id_IsoMe"));
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
                        TagsFRLep.push_back(TagsFR[iTags] + (abs(selLeptons[i].eta())<1.4?TString("_B"):TString("_E")));
                     }

                      mon.fillHisto("wrtJetPt", TagsFRJet, pTL1              , weight);
                      mon.fillHisto("wrtLepPt", TagsFRLep, selLeptons[i].pt(), weight);
                  }
                }//close loop on leptons

            }//close FR study Zmass


            //SIGNAL ANALYSIS Z+2Leptons  (no systematics taken into account here)
            if(passZmass && (int)selLeptons.size()>=4){  //Request at least 4 leptons
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
                case 11*11:  ChannelName  = "elel";  HiggsShortId+= 0; break;
                case 11*13:  ChannelName  = "elmu";  HiggsShortId+= 1; break;
                case 11*15:  ChannelName  = "elha";  HiggsShortId+= 2; break;
                case 13*13:  ChannelName  = "mumu";  HiggsShortId+= 3; break;
                case 13*15:  ChannelName  = "muha";  HiggsShortId+= 4; break;
                case 15*15:  ChannelName  = "haha";  HiggsShortId+= 5; break;
                default:     ChannelName  = "none";  HiggsShortId =-1; break;
                }
               }
               chTagsMain.push_back(chTagsMain[chTagsMain.size()-1] + signName + ChannelName);
               //Get the Higgs candidate (end)

               //printf("%30s %2i --> %i %i %i %i\n", (chTagsMain[chTagsMain.size()-1]).Data(), HiggsShortId, selLeptons[dilLep1].pdgId(), selLeptons[dilLep2].pdgId(), selLeptons[higgsCandL1].pdgId(), selLeptons[higgsCandL2].pdgId());

               //reweight the event to account for lept eff.
               if(isMC && higgsCandL1>=0 && abs(selLeptons[higgsCandL1].pdgId())<15){
		 weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), abs(selLeptons[higgsCandL1].pdgId()), abs(selLeptons[higgsCandL1].pdgId()) ==11 ? "tight" : "tight",patUtils::CutVersion::ICHEP16Cut ).first;
		 weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), abs(selLeptons[higgsCandL1].pdgId()), abs(selLeptons[higgsCandL1].pdgId()) ==11 ? "tightiso" : "tightiso",patUtils::CutVersion::ICHEP16Cut ).first;
               }
               if(isMC && higgsCandL2>=0 && abs(selLeptons[higgsCandL2].pdgId())<15){
		 weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), abs(selLeptons[higgsCandL2].pdgId()), abs(selLeptons[higgsCandL2].pdgId()) ==11 ? "tight" : "tight",patUtils::CutVersion::ICHEP16Cut ).first;
		 weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), abs(selLeptons[higgsCandL2].pdgId()), abs(selLeptons[higgsCandL2].pdgId()) ==11 ? "tightiso" : "tightiso", patUtils::CutVersion::ICHEP16Cut).first;
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
               passHiggsLoose = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.5, 0.5, "decayModeFinding", 0., false, vtx);
               passHiggsMain  = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.3, 0.3, "byLooseCombinedIsolationDeltaBetaCorr3Hits", 20., true, vtx);

               //SVFIT MASS
               higgsCand_SVFit = higgsCand;

               //build the higgs candH
               higgsCandH = zll + higgsCand;
               higgsCandH_SVFit = zll + higgsCand_SVFit;
            }


            bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
            bool passBtags(nbtags==0);
            bool passMinDphijmet( njets==0 || mindphijmet>0.5);
            bool removeDump(false);

            //
            // NOW FOR THE CONTROL PLOTS
            //

	    if(ivar==0){//fill plots only for nominal
	      mon.fillHisto("eventflow_zhllvv", chTagsMain,0, weight);
	      if(selLeptons.size()>=2){
		mon.fillHisto("eventflow_zhllvv",chTagsMain,1, weight);
		mon.fillHisto("leadpt"      ,   chTagsMain, leadingLep.pt(), weight);
		mon.fillHisto("leadeta"     ,   chTagsMain, leadingLep.eta(), weight);
		mon.fillHisto("trailerpt"   ,   chTagsMain, trailerLep.pt(), weight);
		mon.fillHisto("trailereta"  ,   chTagsMain, trailerLep.eta(), weight);
		mon.fillHisto("leppt"       ,   chTagsMain, leadingLep.pt(), weight);
		mon.fillHisto("leppt"       ,   chTagsMain, trailerLep.pt(), weight);
		mon.fillHisto("lepeta"      ,   chTagsMain, leadingLep.eta(), weight);
		mon.fillHisto("lepeta"      ,   chTagsMain, trailerLep.eta(), weight);
		if(passZVeto){
		  mon.fillHisto("eventflow_zhllvv",chTagsMain,2, weight);
		  mon.fillHisto("njets"     , chTags, njets,   weight);
		  mon.fillHisto("nvtx"        ,   chTagsMain, vtx.size(),      weight);
		  mon.fillHisto("nvtxraw"     ,   chTagsMain, vtx.size(),      weight/puWeight);
		  mon.fillHisto("rho"         ,   chTagsMain, rho,       weight);
		  if(jetCount){
		    mon.fillHisto("eventflow_zhllvv", chTagsMain,3, weight);
		    mon.fillHisto("nbtags"    , chTags, nbtags,  weight);
		    if(passBJetVetoMain){
		      mon.fillHisto("eventflow_zhllvv", chTagsMain,4, weight);
		      mon.fillHisto("met"      	, chTagsMain, met.pt()         , weight);
		      if(EtMiss){
			mon.fillHisto("eventflow_zhllvv", chTagsMain,5, weight);
			mon.fillHisto("dPhi_llMET", chTagsMain, DeltaPhillMet,  weight);
			if(DeltaPhillMet){
			  mon.fillHisto("eventflow_zhllvv", chTagsMain,6, weight);
			  mon.fillHisto("llMETUnbalance", chTagsMain, pTU,  weight);
			  if(pTunbalance){
			    mon.fillHisto("eventflow_zhllvv", chTagsMain,7, weight);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }

	 }
	 //END SYSTEMATIC LOOP

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
  TFile *ofile=TFile::Open("out.root", "recreate");
  mon.Write();
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
