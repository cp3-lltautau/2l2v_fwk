#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/regex.hpp>
#include <boost/format.hpp>

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
#include "UserCode/llvv_fwk/interface/NTupleUtils.h"
#include "UserCode/llvv_fwk/interface/TrigUtils.h"
#include "UserCode/llvv_fwk/interface/EwkCorrections.h"
#include "UserCode/llvv_fwk/interface/ZZatNNLO.h"
#include "UserCode/llvv_fwk/interface/FRWeights.h"
//ClassicSVfit
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"


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
#include <bits/stdc++.h>

using namespace std;
//***********************************************************************************************//
LorentzVector getClassicSVFit(pat::MET met, patUtils::GenericLepton firstLepton, patUtils::GenericLepton secondLepton)
//***********************************************************************************************//
{
  using namespace classic_svFit;

  LorentzVector pEmpty(0,0,0,0);

  TMatrixD covMET(2, 2); // PFMET significance matrix

  covMET[0][0] = met.getSignificanceMatrix()(0,0);
  covMET[0][1] = met.getSignificanceMatrix()(0,1);
  covMET[1][0] = met.getSignificanceMatrix()(1,0);
  covMET[1][1] = met.getSignificanceMatrix()(1,1);

  // std::cout<<"MET MATRIX: " << covMET[0][0] << " " << covMET[0][1] << " " << covMET[1][0] << " " << covMET[1][1] << "\n";

  int dlid = abs( firstLepton.pdgId() * secondLepton.pdgId() );

  //std::cout<<"\n"<<firstLepton.pdgId() << "  " << secondLepton.pdgId() << " Di-Tau ID ------------> " << dlid << std::endl;

  if ( dlid == 11*15 || dlid == 13*15){
    if (abs(firstLepton.pdgId()) == 15) {
      // Switching leptons in semileptonic pairs:
      // e and mu should be passed as first MesuredTauLepton
      patUtils::GenericLepton temp = firstLepton;
      firstLepton  = secondLepton;
      secondLepton = temp;
    }
  }

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

  if ( dlid == 11*15 ){
    //std::cout<< " ETau Pair --- > "<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), secondLepton.mass(), secondLepton.tau.decayMode()) );
  }
  else if( dlid == 13*15 ){
    //std::cout<< " MuTau Pair --- > "<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), classic_svFit::muonMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), secondLepton.mass(), secondLepton.tau.decayMode()) );
  }
  else if ( dlid == 15*15 ){
    //std::cout<< " TauTau Pair --- > "<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), firstLepton.mass(), firstLepton.tau.decayMode()) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), secondLepton.mass(), secondLepton.tau.decayMode()) );
  }
  else if (dlid == 13*11 ){
    //std::cout<< " EMu Pair  --->"<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), classic_svFit::muonMass) );
  }
  else if (dlid == 13*13){
    //std::cout<< " EE Pair  --->"<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
        measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), classic_svFit::muonMass) );
	measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), classic_svFit::muonMass) );
  }
  else if (dlid == 11*11){
    //std::cout<< " EE Pair  --->"<< firstLepton.pdgId() << "  " << secondLepton.pdgId() << std::endl;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, firstLepton.pt(), firstLepton.eta(),
								    firstLepton.phi(), classic_svFit::electronMass) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, secondLepton.pt(), secondLepton.eta(),
								    secondLepton.phi(), classic_svFit::electronMass) );
  }
  else return (firstLepton.p4()+secondLepton.p4());


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

pat::JetCollection  skimJetsCollection (const std::vector<patUtils::GenericLepton> diLeptons, pat::JetCollection& selJets){

  pat::JetCollection  jetCollection;
  auto firstLepton  = diLeptons.at(0);
  auto secondLepton = diLeptons.at(1);

  for(auto& jet: selJets){
    if (deltaR( jet,  firstLepton )<0.4 || deltaR( jet, secondLepton )<0.4 ) continue;
    else jetCollection.push_back(jet);
  }
  return jetCollection;
}

//**************************************************************************************************************************//
std::pair<int,pat::JetCollection>  skimJetsCollectionInfo (const std::vector<patUtils::GenericLepton> diLeptons, pat::JetCollection& selJets, pat::TauCollection selTaus)
//**************************************************************************************************************************//
{
  std::pair<int,pat::JetCollection> jetCollentionInfo;
  //pat::JetCollection  jetCollection;
  auto firstLepton  = diLeptons.at(0);
  auto secondLepton = diLeptons.at(1);

  LorentzVector zll = (firstLepton.p4()+secondLepton.p4());

  for(auto& jetCorr: selJets){

    // cout << "  Jet :  "<< jetCorr << endl;
    const pat::Jet jet = jetCorr.correctedJet("Uncorrected");
    //cout << "  Un Jet :  "<< jet << endl;
//    if (deltaR( jet,  firstLepton )<0.1 || deltaR( jet, secondLepton )<0.1 ) {
    LorentzVector zl1_jet = firstLepton.p4()+jet.p4();
    LorentzVector zl2_jet = secondLepton.p4()+jet.p4();
    if (abs( zl1_jet.mass() - zll.mass() )< 10. || abs( zl2_jet.mass() - zll.mass() )< 10. ) {
      continue;
    }
    else jetCollentionInfo.second.push_back(jet);

    bool overlapWithTau(false);
    for(int l1=0; l1<(int)selTaus.size();++l1){
      if(deltaR(jet, selTaus[l1])< 0.5  ){overlapWithTau=true; break;}
    }

    if(!overlapWithTau && jet.pt()>30 && fabs(jet.eta())<2.5){
      bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.8484);
      if (hasCSVtag) jetCollentionInfo.first++;
    }
  }
  return jetCollentionInfo;
}


//**************************************************************************************************************************//
unsigned char leptonIDmap (patUtils::GenericLepton lepton, reco::Vertex& vtx)
//**************************************************************************************************************************//
{
  unsigned char lepIDmap = 0xF; // 1111
  if(abs(lepton.pdgId())==11 || abs(lepton.pdgId())==13){
    bool passIDloose(false),  passIDtight(false);
    bool passIDvloose(false), passIDsoft(false); // Muon WP
    bool passIDmedium(false); // Ele WP

    if(abs(lepton.pdgId())==11) {
      passIDloose  = patUtils::passId(lepton.el, vtx, patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIDmedium = patUtils::passId(lepton.el, vtx, patUtils::llvvElecId::Medium, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIDtight  = patUtils::passId(lepton.el, vtx, patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
      lepIDmap &= (passIDloose? 1 << 1: 0 << 1) | (passIDmedium? 1 << 2 : 0 << 2) | (passIDtight ? 1 << 3 : 0 << 3);
    }
    if(abs(lepton.pdgId())==13) {
      passIDvloose = patUtils::passId(lepton.mu, vtx, patUtils::llvvMuonId::FRLoose, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIDloose  = patUtils::passId(lepton.mu, vtx, patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIDsoft   = patUtils::passId(lepton.mu, vtx, patUtils::llvvMuonId::Soft, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIDtight  = patUtils::passId(lepton.mu, vtx, patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
      lepIDmap &= (passIDvloose? 1 : 0) | (passIDloose? 1 << 1 : 0 << 1) | (passIDsoft ? 1 << 2 : 0 << 2) | (passIDtight? 1 << 3 : 0 << 3) ;
    //cout << " (After checking) Muon ID: "<< (int) ( (passIDloose? 1 : 0) | (passIDsoft? 1 << 1 : 0 << 1) | (passIDtight ? 1 << 2 : 0 << 2) ) << endl;
    }
  }
  else if(abs(lepton.pdgId())==15){
    bool IdL         = lepton.tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
    bool IdM         = lepton.tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
    bool IdVL_MVA    = lepton.tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
    bool IdL_MVA     = lepton.tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    bool IdM_MVA     = lepton.tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    bool IdT_MVA     = lepton.tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    bool IdL_MVA_R03 = lepton.tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
    bool IdM_MVA_R03 = lepton.tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");

    lepIDmap &= (IdVL_MVA? 1 : 0) | (IdL_MVA ? 1 << 1 : 0 << 1)| (IdM_MVA? 1 << 2 : 0 << 2) | (IdT_MVA ? 1 << 3 : 0 << 3);
  }

return lepIDmap;
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
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_ZZ2l2nu  = isMC && ( string(dtag.Data()).find("TeV_ZZ2l2nu")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  bool isMC_WZ3lnu  = isMC && ( string(dtag.Data()).find("TeV_WZ3lnu")  != string::npos);
  bool is2015data = (!isMC && dtag.Contains("2015"));
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
  //Z leptons kinematics control
  mon.addHistogram( new TH1F( "leadpt"      ,  ";p_{T}^{lead} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta"     ,  ";#eta_{lead};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "leadiso"     ,  ";#I_{lead};Events", 20,0.,1.) );
  mon.addHistogram( new TH1F( "trailerpt"   ,  ";p_{T}^{trail} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta"  ,  ";#eta_{trail};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "traileriso"     ,  ";#I_{trail};Events", 20,0.,1.) );
  mon.addHistogram( new TH1F( "leppt"       ,  ";p_{T}^{lepton} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "lepeta"      ,  ";#eta_{lepton};Events", 50,-2.6,2.6) );

  // zll control
  mon.addHistogram( new TH1F( "zlly",      		";y_{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zlleta",    		";#eta_{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zllpt",     		";p_{T}^{ll} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "zllmass",   		";M_{ll} (GeV);Events/2 GeV", 80,20,180) );
  mon.addHistogram( new TH1F( "nlep",      	";Number of Leptons;Events", 10,0,10) );

  float ptbinsJets[] = {10, 20, 30, 40, 60, 80, 100, 125, 150, 175,250,350,1000};
  int ptbinsJetsN = sizeof(ptbinsJets)/sizeof(float)-1;
  mon.addHistogram( new TH1F( "wrtJetPt",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtJetPt_v2",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtLepPt",  ";Lep p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));

  LorentzVector null_p4(0., 0., 0., 0.);

  //create a tree and related variables to save higgs candidate info for each cutIndex values
  ntupleutils::Weights *eventWeights = nullptr;
  ntupleutils::analysis::gen_truth::LheSummary *lheSummary = nullptr;
  unsigned int  treeEventId;
  unsigned int  treeLumiId;
  unsigned int  treeRunId;
  float         treeNumEvents;
  bitset<4>     treeTriggerFiredBits;
  int           treeDoubleTriggerMatch;
  int           treeLepIsMatched;
  int           treeZId;
  int           treeHiggsPairCharge;
  int           treeHiggsPairFlavour;
  float         treeWeight;
  float         treeWeightNoLepSF;
  float         treeLep3Eff;
  float         treeLep4Eff;
  int           treeNlep;
  int           treeNbtag;
  LorentzVector treeJetCandidate;
  LorentzVector treeLep1Candidate;
  int           treeLep1GenMatch;
  int           treeLep1PdgId;
  LorentzVector treeLep2Candidate;
  int           treeLep2GenMatch;
  int           treeLep2PdgId;
  LorentzVector treeLep3Candidate;
  int           treeLep3GenMatch;
  int           treeLep3PdgId;
  unsigned char treeLep3ID;
  float         treeLep3Iso;
  float         treeTMass;
  int           treeLeg1LepID_loose;
  int           treeLeg1LepID_medium;
  int           treeLeg1LepID_tight;
  int           treeLeg1TauIsoLooseMVA;
  int           treeLeg1TauIsoMediumMVA;
  LorentzVector treeLep4Candidate;
  int           treeLep4GenMatch;
  int           treeLep4PdgId;
  unsigned char treeLep4ID;
  float         treeLep4Iso;
  LorentzVector treeDiTauSVfitCandidate;

  TFile *ofile=TFile::Open("out.root", "recreate");

  TTree* tree = new TTree("CandTree","CandTree");
  tree->Branch("EventWeights",&eventWeights);
  tree->Branch("LHESummary",&lheSummary);
  tree->Branch("eventId", &treeEventId , string("eventId/i" ).c_str());
  tree->Branch("lumiId" , &treeLumiId  , string("lumiId/i"  ).c_str());
  tree->Branch("runId"  , &treeRunId   , string("runId/i"   ).c_str());
  tree->Branch("numEvents" , &treeNumEvents  , string("numEvents/F"  ).c_str());
  tree->Branch("triggerFiredBits",     &treeTriggerFiredBits    , string("triggerFiredBits/b" ).c_str());
  tree->Branch("doubleTriggerMatch",     &treeDoubleTriggerMatch    , string("doubleTriggerMatch/I" ).c_str());
  tree->Branch("lepIsMatched" , &treeLepIsMatched , string("lepIsMatched/I").c_str());
  tree->Branch("zId",     &treeZId     , string("zId/I" ).c_str());
  tree->Branch("higgsPairCharge",     &treeHiggsPairCharge     , string("higgsPairCharge/I" ).c_str());
  tree->Branch("higgsPairFlavour",     &treeHiggsPairFlavour     , string("higgsPairFlavour/I" ).c_str());
  tree->Branch("weight" , &treeWeight  , string("weight/F"  ).c_str());
  tree->Branch("weightNoLepSF" , &treeWeightNoLepSF  , string("weightNoLepSF/F"  ).c_str());
  tree->Branch("lep3Eff" , &treeLep3Eff  , string("lep3Eff/F"  ).c_str());
  tree->Branch("lep4Eff" , &treeLep4Eff  , string("lep4Eff/F"  ).c_str());
  tree->Branch("nlep" , &treeNlep  , string("nlep/I"  ).c_str());
  tree->Branch("nbtag" , &treeNbtag  , string("nbtag/I"  ).c_str());

  tree->Branch("jetMatched",     &treeJetCandidate );
  /* Lep 1 */
  tree->Branch("lep1Candidate",     &treeLep1Candidate );
  tree->Branch("lep1GenMatch",     &treeLep1GenMatch     , string("lep1GenMatch/I" ).c_str());
  tree->Branch("lep1PdgId",     &treeLep1PdgId     , string("lep1PdgId/I" ).c_str());
  /* Lep 2 */
  tree->Branch("lep2Candidate",     &treeLep2Candidate );
  tree->Branch("lep2GenMatch",     &treeLep2GenMatch     , string("lep2GenMatch/I" ).c_str());
  tree->Branch("lep2PdgId",     &treeLep2PdgId     , string("lep2PdgId/I" ).c_str());
  /* Lep 3 */
  tree->Branch("lep3Candidate",     &treeLep3Candidate );
  tree->Branch("lep3GenMatch",     &treeLep3GenMatch     , string("lep3GenMatch/I" ).c_str());
  tree->Branch("lep3PdgId" , &treeLep3PdgId  , string("lep3PdgId/I"  ).c_str());
  tree->Branch("lep3ID", &treeLep3ID , string("lep3ID/b" ).c_str());
  tree->Branch("lep3Iso", &treeLep3Iso , string("lep3Iso/F" ).c_str());
  tree->Branch("TMass" , &treeTMass  , string("TMass/F"  ).c_str());
  tree->Branch("leg1LepIDloose", &treeLeg1LepID_loose , string("leg1LepIDloose/I" ).c_str());
  tree->Branch("leg1LepIDmedium", &treeLeg1LepID_medium , string("leg1LepIDmedium/I" ).c_str());
  tree->Branch("leg1LepIDtight", &treeLeg1LepID_tight , string("leg1LepIDtight/I" ).c_str());
  tree->Branch("tau1LooseMVA", &treeLeg1TauIsoLooseMVA , string("tau1LooseMVA/i" ).c_str());
  tree->Branch("tau1MediumMVA", &treeLeg1TauIsoMediumMVA , string("tau1MediumMVA/i" ).c_str());
  /* Lep 4 */
  tree->Branch("lep4Candidate",     &treeLep4Candidate );
  tree->Branch("lep4GenMatch",     &treeLep4GenMatch     , string("lep4GenMatch/I" ).c_str());
  tree->Branch("lep4PdgId" , &treeLep4PdgId  , string("lep4PdgId/I"  ).c_str());
  tree->Branch("lep4ID", &treeLep4ID , string("lep4ID/b" ).c_str());
  tree->Branch("lep4Iso", &treeLep4Iso , string("lep4Iso/F" ).c_str());
  /* di-tau system */
  tree->Branch("diTauSVfitCandidate",     &treeDiTauSVfitCandidate );

  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  //MC normalization (to 1/pb)
  treeNumEvents = 0;
  double xsecWeight = 1.0;
  treeNumEvents = utils::getTotalNumberOfEvents(urls, false, true);
  if(isMC) xsecWeight=xsec/treeNumEvents; //utils::getTotalNumberOfEvents(urls, false, true);//need to use the slow method in order to take NLO negative events into account

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
      eventWeights = new ntupleutils::Weights();
      eventWeights->SetXSecWeight(xsecWeight);
      float shapeWeight = 1.0;
      double puWeightUp = 1.0;
      double puWeightDown = 1.0;
      float puWeight(1.0);
      float weightNoLepSF(1.0);

      //##############################################   EVENT LOOP STARTS   ##############################################
      //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      //Skip bad lumi
      if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;
      lheSummary = new ntupleutils::analysis::gen_truth::LheSummary();
      treeEventId=0;
  		treeLumiId=0;
  		treeRunId=0;
      treeDoubleTriggerMatch=0;
      treeLepIsMatched=0;
      treeTriggerFiredBits.reset();
  		treeZId=0;
      treeWeight = 0;
      treeWeightNoLepSF = 0;
      treeLep3Eff = 0;
      treeLep4Eff = 0;
      treeNlep = 0;
      treeNbtag = 0;
      treeHiggsPairCharge = -1;
      treeHiggsPairFlavour = -1;

      treeJetCandidate = null_p4;
      treeLep1Candidate = null_p4;
      treeLep1GenMatch = -1;
      treeLep1PdgId = 0;
      treeLep2Candidate = null_p4;
      treeLep2GenMatch = -1;
      treeLep2PdgId = 0;
      treeLep3Candidate = null_p4;
      treeLep3GenMatch = -1;
      treeLep3PdgId = 0;
      treeLep3ID = 0;
      treeLep3Iso = 0;
      treeTMass = 0;
      treeLeg1LepID_loose = 0;
      treeLeg1LepID_medium = 0;
      treeLeg1LepID_tight = 0;
      treeLeg1TauIsoLooseMVA = 0;
      treeLeg1TauIsoMediumMVA = 0;
      treeLep4Candidate = null_p4;
      treeLep4GenMatch = -1;
      treeLep4PdgId = 0;
      treeLep4ID = 0;
      treeLep4Iso = 0;
      treeDiTauSVfitCandidate = null_p4;

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

        eventWeights->SetMCWeight( eventInfo.weight() );


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

        eventWeights->SetPUWeights(puWeight,puWeightUp,puWeightDown);
      }



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

    if(  mumuTrigger){treeTriggerFiredBits.set(0);mon.fillHisto("trigger", "raw", 0 , weight);}
    if(    muTrigger){treeTriggerFiredBits.set(1);mon.fillHisto("trigger", "raw", 1 , weight);}
    if(    eeTrigger){treeTriggerFiredBits.set(2);mon.fillHisto("trigger", "raw", 2 , weight);}
    if(     eTrigger){treeTriggerFiredBits.set(3);mon.fillHisto("trigger", "raw", 3 , weight);}
    if(   emuTrigger)mon.fillHisto("trigger", "raw", 4 , weight);

    if(!isMC && passTrigger){ //avoid double counting of events from different PD
      treeTriggerFiredBits.reset();
      if(filterOnlyMUMU)     { passTrigger = mumuTrigger;}
      if(filterOnlyMU)       { passTrigger = muTrigger     && !mumuTrigger;}
      if(filterOnlyEE)       { passTrigger = eeTrigger     && !muTrigger  && !mumuTrigger;}
      if(filterOnlyE)        { passTrigger = eTrigger      && !eeTrigger  && !muTrigger && !mumuTrigger; }
      if(filterOnlyEMU)      { passTrigger = emuTrigger    && !eTrigger   && !eeTrigger && !muTrigger && !mumuTrigger; }

      if(filterOnlyMUMU && passTrigger) treeTriggerFiredBits.set(0);
      if(filterOnlyMU && passTrigger)   treeTriggerFiredBits.set(1);
      if(filterOnlyEE && passTrigger)   treeTriggerFiredBits.set(2);
      if(filterOnlyE && passTrigger)    treeTriggerFiredBits.set(3);
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


    if (isMC){
      fwlite::Handle< LHEEventProduct > lheEPHandle;
      lheEPHandle.getByLabel(ev, "externalLHEProducer");
      if(lheEPHandle.isValid()){
        auto& lheEventProduct = *lheEPHandle;
        auto lheSummaryValue = ntupleutils::analysis::gen_truth::ExtractLheSummary(lheEventProduct);
        lheSummary->n_partons   = lheSummaryValue.n_partons;
        lheSummary->n_b_partons = lheSummaryValue.n_b_partons;
        lheSummary->n_c_partons = lheSummaryValue.n_c_partons;
        lheSummary->HT          = lheSummaryValue.HT;
      }
    }
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

    eventWeights->SetEwkMCWeight(ewkCorrectionsWeight);
    eventWeights->SetZZMCWeight(ZZ_NNLOcorrectionsWeight);

    weightNoLepSF = weight;

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
        if(deltaR(leptons[ilep].p4(), selLeptons[l1])<0.3){overlapWithLepton=true; break;}
      }if(overlapWithLepton)continue;

      //Cut based identification
      //isolation
      // passLooseLepton
      passLooseLepton &= lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
                                   patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
      passIsoWPforFakeRate &= lid==11 ?  patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut) :
                                        patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::FakeRateWP, patUtils::CutVersion::CutSet::ICHEP16Cut);
      bool applyCorr = true;
  //     //apply muon corrections
  //     //if(abs(lid)==13 && passIso && passId)
      if(applyCorr && abs(lid)==13 && passIsoWPforFakeRate && passLooseLepton && leptons[ilep].pt()>20 ){
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
      if(applyCorr && abs(lid)==11  && passIsoWPforFakeRate && passLooseLepton){
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
        pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
        if(!(fabs(packedLeadTauCand->dz()) < 0.2)) continue;

        selTaus.push_back(tau);
        selLeptons.push_back(tau);
        ntaus++;
      }
      std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
      std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);

      mon.fillHisto("nlep"           ,   "controlPlots", selLeptons.size(), weight);
      mon.fillHisto("nlep"           ,   "controlPlots_Int", (int)selLeptons.size(), weight);

      if ( selLeptons.size() < 3 ) continue;
      /*******************************************************************************/
      /***    Cut on number of selLeptons     **/
      mon.fillHisto("nlep"           ,   "controlPlots_AfterCut", selLeptons.size(), weight);
      mon.fillHisto("nlep"           ,   "controlPlots_AfterCut", (int)selLeptons.size(), weight);


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

        if(jet.pt()<20 || fabs(jet.eta())>4.7 ) continue;

        //TString jetType( jet.genJet() && (jet.genJet())->pt()>0 ? "truejetsid" : "pujetsid" );

        //jet id
        bool passPFloose = patUtils::passPFJetID("Loose", jet);
        bool passLooseSimplePuId = patUtils::passPUJetID(jet); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
        bool passLooseSimplePuId_new = jet.userInt("pileupJetId:fullId") & (1 << 2);
        if(passLooseSimplePuId)          mon.fillHisto( "jetId", "controlPlots", 0, jet.pt(), weight);
        if(passLooseSimplePuId_new)      mon.fillHisto( "jetId", "controlPlots", 1, jet.pt(), weight);
        // cout << "passLooseSimplePuId: "<< passLooseSimplePuId <<endl;
        // cout << "Passes loose: " << bool(jet.userInt("pileupJetId:fullId") & (1 << 2)) << " - medium: " << bool(jet.userInt("pileupJetId:fullId") & (1 << 1)) <<" - tight: "<< bool(jet.userInt("pileupJetId:fullId") & (1 << 0)) << endl;
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
        patUtils::GenericLepton *dilLep1=nullptr;
        patUtils::GenericLepton *dilLep2=nullptr;
        double BestMass;
        LorentzVector leadingLep, trailerLep, zll, zlltmp;
        //get the Z candidate
        dilId=-1;
        BestMass=0;
        zll = LorentzVector(0.,0.,0.,0.);

        for(auto& lep1: selLeptons){
          if(abs(lep1.pdgId())==15)continue;

          double leadPtCutValue  = abs(lep1.pdgId())==11 ? 27.0 : 19.0;
          if( lep1.pt()< leadPtCutValue ) continue;
          if(!( abs(lep1.pdgId())==11 ? patUtils::passIso(lep1.el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
                                                patUtils::passIso(lep1.mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut)) ||
             !( abs(lep1.pdgId())==11 ? patUtils::passId(lep1.el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
                                        patUtils::passId(lep1.mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut)) ) continue;

          for(auto& lep2: selLeptons){
            //cout << " lep1 : " << &lep1 << " lep2 : " << &lep2 << " - check: " << (&lep1 == &lep2 ? "True" : "False") << endl;
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
              dilLep1 = &lep1;
              dilLep2 = &lep2;
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
              weight *= lepEff.getRecoEfficiency( dilLep1->el.superCluster()->eta(), abs(dilLep1->pdgId())).first; //Reconstruction eff
              weight *= lepEff.getRecoEfficiency( dilLep2->el.superCluster()->eta(), abs(dilLep2->pdgId())).first; //Reconstruction eff
            }
            else if(abs(dilId)==169){
              weight *= lepEff.getTrackingEfficiency( dilLep1->eta(), abs( dilLep1->pdgId())).first; //Tracking eff
              weight *= lepEff.getTrackingEfficiency( dilLep2->eta(), abs( dilLep2->pdgId())).first; //Tracking eff

              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep1->pt(), dilLep1->eta(), abs(dilLep1->pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ISO w.r.t ID
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep2->pt(), dilLep2->eta(), abs(dilLep2->pdgId()), "tightiso_tightid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;
            }
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep1->pt(), dilLep1->eta(), abs( dilLep1->pdgId()),  abs( dilLep1->pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID
              weight *= isMC ? lepEff.getLeptonEfficiency( dilLep2->pt(), dilLep2->eta(), abs( dilLep2->pdgId()),  abs( dilLep2->pdgId()) ==11 ? "tight"    : "tight"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID

            // Trigger Eff
            if(isMC && abs(dilId)==169)weight *= lepEff.getTriggerEfficiencySF( dilLep1->pt(), dilLep1->eta(), dilLep2->pt(), dilLep2->eta(), dilId,is2016MC).first;
       		  if(isMC && abs(dilId)==121)weight *= lepEff.getTriggerEfficiencySF( dilLep1->pt(), dilLep1->el.superCluster()->eta(), dilLep2->pt(), dilLep2->el.superCluster()->eta(), dilId,is2016MC).first;
         }

        }


        bool passZmass = (fabs(zll.mass()-91.2)<30.0);
        bool passZpt   = (zll.pt()>20);

        if( !isDileptonCandidate || !passZmass ) continue;
        /************************* EVENT HAS a Z-like candidate ***************************/
    	      mon.fillHisto("leadpt"      ,   chTags, leadingLep.pt(), weight);
              mon.fillHisto("leadeta"     ,   chTags, leadingLep.eta(), weight);
              mon.fillHisto("trailerpt"   ,   chTags, trailerLep.pt(), weight);
              mon.fillHisto("trailereta"  ,   chTags, trailerLep.eta(), weight);
              mon.fillHisto("leppt"       ,   chTags, leadingLep.pt(), weight);
              mon.fillHisto("leppt"       ,   chTags, trailerLep.pt(), weight);
              mon.fillHisto("lepeta"      ,   chTags, leadingLep.eta(), weight);
              mon.fillHisto("lepeta"      ,   chTags, trailerLep.eta(), weight);

              //analyze dilepton kinematics
              mon.fillHisto("zllpt"         ,   chTags, zll.pt(),      weight);
              mon.fillHisto("zlleta"        ,   chTags, zll.eta(),     weight);
              mon.fillHisto("zlly"          ,   chTags, zll.Rapidity(),weight);
	mon.fillHisto("zllmass"          ,   chTags, zll.mass(),    weight);
// cout<<"  ##RECO##  Z Lepton 1:  pt = "<<(*dilLep1).pt()<<"  eta = "<<(*dilLep1).eta()<<"  phi = "<<(*dilLep1).phi()<<endl;
        // // if ( selLeptons[dilLep1].genParticle() ) cout<<"    ##RECO (GEN Match)##  Z Lepton 1:  pt = "<<(selLeptons[dilLep1].genParticle())->pt()<<"  eta = "<<(selLeptons[dilLep1].genParticle())->eta()<<"  phi = "<<(selLeptons[dilLep1].genParticle())->phi()<<endl;
        // cout<<"  ##RECO##  Z Lepton 2:  pt = "<<(*dilLep2).pt()<<"  eta = "<<(*dilLep2).eta()<<"  phi = "<<(*dilLep2).phi()<<endl;
        // if ( selLeptons[dilLep2].genParticle() ) cout<<"    ##RECO (GEN Match)##  Z Lepton 2:  pt = "<<(selLeptons[dilLep2].genParticle())->pt()<<"  eta = "<<(selLeptons[dilLep2].genParticle())->eta()<<"  phi = "<<(selLeptons[dilLep2].genParticle())->phi()<<endl;
        //
        // int triggerType = abs(dilId)==121 ? 82 : 83;
        // std::cout << "\n TRIGGER OBJECTS " << std::endl;
        std::vector<patUtils::GenericLepton> diLeptons;
        diLeptons.push_back(*dilLep1);
        diLeptons.push_back(*dilLep2);

        // Removing jets that match with leptons from Z boson
        // auto selJetsSkimmed = skimJetsCollection(diLeptons,selJets);
        auto selJetsSkimmedInfo = skimJetsCollectionInfo(diLeptons,selJets,selTaus);
        auto selJetsSkimmed     = selJetsSkimmedInfo.second;
        auto njetSkimmed        = selJetsSkimmed.size();
        auto nbtagSkimmed       = selJetsSkimmedInfo.first;


        treeEventId   = ev.eventAuxiliary().event();
        treeLumiId    = ev.eventAuxiliary().luminosityBlock();
        treeRunId     = ev.eventAuxiliary().run();
        treeZId       = dilId;
        treeLep1Candidate = dilLep1->p4();
        treeLep1PdgId     = dilLep1->pdgId();
        if (isMC) treeLep1GenMatch = (int) ntupleutils::analysis::gen_truth::LeptonGenMatch(dilLep1->p4(), gen).first;
        treeLep2Candidate = dilLep2->p4();
        treeLep2PdgId     = dilLep2->pdgId();
        if (isMC) treeLep2GenMatch = (int) ntupleutils::analysis::gen_truth::LeptonGenMatch(dilLep2->p4(), gen).first;
        treeNbtag               = nbtagSkimmed;
        treeNlep                = selLeptons.size();

        std::map <std::string,std::vector<std::string> > muonTriggerToFilter = { {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",{"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2"}},
                                                                                  {"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", {"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"}}
                                                                                };
        std::map <std::string,std::vector<std::string> > eleTriggerToFilter  = { {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",{"hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter"}}
                                                                                };
        auto triggerToFilterMap = abs(dilId)==121? eleTriggerToFilter: muonTriggerToFilter;
        std::vector<pat::TriggerObjectStandAlone> triggerObjectMatchedCollection;

        for (auto lepton: diLeptons){
        for (pat::TriggerObjectStandAlone obj : *triggerObjectsHandle) { // note: not "const &" since we want to call unpackPathNames
            obj.unpackPathNames(names);

            // bool typeMatched = false;
            // for (unsigned h = 0; h < obj.filterIds().size(); ++h) {
            //   typeMatched |= (obj.filterIds()[h] == triggerType) ;
            // }
            // if (!typeMatched) continue;
            bool filterPass = true;
            if(deltaR( obj,  lepton )>0.5 ) continue;
            for(std::map<std::string,std::vector<std::string> >::iterator iter = triggerToFilterMap.begin(); iter != triggerToFilterMap.end(); ++iter)
            {
              auto k =  iter->first;
              for (auto& filter : iter->second){
                if(!obj.hasFilterLabel(filter)) filterPass &= false;
              }
            if (!filterPass) continue;
            // std::cout << "---- " << k << endl;

            // // Print trigger object collection and type
            // std::cout << "\t   Collection: " << obj.collection() << std::endl;
            // std::cout << "\t   Type IDs:   ";
            // for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
            // std::cout << std::endl;
            // // Print associated trigger filters
            // std::cout << "\t   Filters:    ";
            // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
            // std::cout << std::endl;
            std::vector< std::string > pathNamesAll = obj.pathNames(false);
            std::vector< std::string > pathNamesLast = obj.pathNames(true);
            // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
            // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
            // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
            // std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
              for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

                 bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
                 if (!isBoth) continue;
                 static const std::string regex_format = "^%1%[0-9]+$";
                 const std::string regex_str = boost::str(boost::format(regex_format) % k);
                 boost::regex regex = boost::regex(regex_str);
                 bool matched = boost::regex_match(pathNamesAll[h], regex);
                 if (matched) {
                   triggerObjectMatchedCollection.push_back(obj);
                   // cout << " ------>>  Trigger Path: "<<pathNamesAll[h]<<" (L,3) "<<endl;
                   // std::cout << std::endl;
                   // std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << ", pdgId "<< obj.pdgId() << std::endl;
                 }
            //     bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
            //     bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
            //     bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
            //     std::cout << "   " << pathNamesAll[h];
            //     if (isBoth) std::cout << "(L,3)";
            //     if (isL3 && !isBoth) std::cout << "(*,3)";
            //     if (isLF && !isBoth) std::cout << "(L,*)";
            //     if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
             }
           }
          // std::cout << std::endl;
        }
      }

      if( triggerObjectMatchedCollection.size()==2 ) treeDoubleTriggerMatch =1;
        // std::cout << std::endl;

        //LEPTON FAKE RATE ANALYSIS Z+1jets  (no systematics taken into account here)
        if( passZmass && treeNlep == 3 ){ //Request exactly one Z + 1 additional lepton
          double tmass=-999;
         // cout<< " EventiD : "<< treeEventId << endl;
          for(auto& fakeLepton: selLeptons){
            if( ( &fakeLepton == dilLep1 ) || ( &fakeLepton == dilLep2 ) ) continue;
           // cout << "   -  fakeLepton:  "<< &fakeLepton << " pt: "<< fakeLepton.pt()<<endl;
            if( deltaR( fakeLepton,  *dilLep1 )<0.1 || deltaR( fakeLepton, *dilLep2)<0.1 ) continue;
           // cout << "      -  fakeLepton:  "<< &fakeLepton << " pt: "<< fakeLepton.pt()<<endl;
            if(abs( fakeLepton.pdgId() ) == 11 || abs( fakeLepton.pdgId() ) == 13 || abs( fakeLepton.pdgId() ) == 15 ){
              tmass = TMath::Sqrt(2* fakeLepton.pt()*met.pt()*(1-TMath::Cos(deltaPhi(met.phi(), fakeLepton.phi()))));
            }
            if(abs(fakeLepton.pdgId())==11 || abs(fakeLepton.pdgId())==13 || abs(fakeLepton.pdgId())==15){
              int closestJetIndexL1=-1; double pTL1=-1; double etaL1=-1;
              double dRminL1 = closestJet(fakeLepton.p4(), selJetsSkimmed, closestJetIndexL1);
              bool jetFound = closestJetIndexL1>=0 && dRminL1<0.5;
              if(jetFound){
		              pTL1=selJetsSkimmed[closestJetIndexL1].pt();
		              etaL1=abs(selJetsSkimmed[closestJetIndexL1].eta());
                  treeJetCandidate = selJetsSkimmed[closestJetIndexL1].p4();
	            }
              else{pTL1=fakeLepton.pt(); etaL1=abs(fakeLepton.eta());}

              if (isMC) treeLep3GenMatch= (int) ntupleutils::analysis::gen_truth::LeptonGenMatch(fakeLepton.p4(), gen).first;

              treeLepIsMatched       = jetFound;

              treeLep3Candidate      = fakeLepton.p4();
              treeTMass              = tmass;
        			treeLep3PdgId  = fakeLepton.pdgId();

              treeLep3Iso  = (abs(treeLep3PdgId)!=15)? patUtils::relIso(fakeLepton,rho): fakeLepton.tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

              TString PartName = "FR_";//+chTags.at(1)+"_";
              if     (abs(fakeLepton.pdgId())==11)PartName += "El";
              else if(abs(fakeLepton.pdgId())==13)PartName += "Mu";
              else if(abs(fakeLepton.pdgId())==15)PartName += "Ta";
              else PartName+= abs(fakeLepton.pdgId());


              std::vector<TString> TagsFR;

              unsigned char lepIDmap = 0xF; // 1111
              if(abs(fakeLepton.pdgId())==11 || abs(fakeLepton.pdgId())==13){
                bool passId = false;
                bool passIDloose(false),  passIDtight(false);
                bool passIDvloose(false), passIDsoft(false); // Muon WP
                bool passIDmedium(false); // Ele WP

                if(abs(fakeLepton.pdgId())==11) {
                  passId = patUtils::passId(fakeLepton.el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDloose  = patUtils::passId(fakeLepton.el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDmedium = patUtils::passId(fakeLepton.el, vtx[0], patUtils::llvvElecId::Medium, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDtight  = patUtils::passId(fakeLepton.el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  lepIDmap &= (passIDloose? 1 << 1: 0 << 1) | (passIDmedium? 1 << 2 : 0 << 2) | (passIDtight ? 1 << 3 : 0 << 3);
                }
                if(abs(fakeLepton.pdgId())==13) {
                  passId = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDvloose = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::FRLoose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDloose = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDsoft  = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::Soft, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  passIDtight = patUtils::passId(fakeLepton.mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
                  lepIDmap &= (passIDvloose? 1 : 0) | (passIDloose? 1 << 1 : 0 << 1) | (passIDsoft ? 1 << 2 : 0 << 2) | (passIDtight? 1 << 3 : 0 << 3) ;
                  //cout << " (After checking) Muon ID: "<< (int) ( (passIDloose? 1 : 0) | (passIDsoft? 1 << 1 : 0 << 1) | (passIDtight ? 1 << 2 : 0 << 2) ) << endl;
                }
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
                bool IdVL_MVA    = fakeLepton.tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
                bool IdL_MVA     = fakeLepton.tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
                bool IdM_MVA     = fakeLepton.tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
                bool IdT_MVA     = fakeLepton.tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
                bool IdL_MVA_R03 = fakeLepton.tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
                bool IdM_MVA_R03 = fakeLepton.tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");

                lepIDmap &= (IdVL_MVA? 1 : 0) | (IdL_MVA ? 1 << 1 : 0 << 1)| (IdM_MVA? 1 << 2 : 0 << 2) | (IdT_MVA ? 1 << 3 : 0 << 3);

                if(true                 )TagsFR.push_back(PartName);
                // if(IdL                  )TagsFR.push_back(PartName+("_Id_IsoLo"));
                // if(IdM                  )TagsFR.push_back(PartName+("_Id_IsoMe"));
                if(IdL_MVA              ){
                  TagsFR.push_back(PartName+("_Id_IsoLo_MVA"));
                  treeLeg1TauIsoLooseMVA = 1;
                }
                if(IdM_MVA              ){
                  TagsFR.push_back(PartName+("_Id_IsoMe_MVA"));
                  treeLeg1TauIsoMediumMVA = 1;
                }
                // if(IdL_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoLo_MVAR03"));
                // if(IdM_MVA_R03          )TagsFR.push_back(PartName+("_Id_IsoMe_MVAR03"));
              }
              treeLep3ID = lepIDmap;

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

        }//close FR study Zmass


        std::vector<patUtils::GenericLepton> higgsLikePair;
        //SIGNAL ANALYSIS Z+2Leptons  (no systematics taken into account here)
        if(passZmass && treeNlep>=4){  //Request at least 4 leptons

          using namespace ntupleutils;

          //Get the Higgs candidate
          int higgsCandId=0;

          for(auto& higgsLepton: selLeptons){
            if( ( &higgsLepton == dilLep1 ) || ( &higgsLepton == dilLep2 ) ) continue;
            if ( higgsLikePair.size() < 1 ) {
              higgsLikePair.push_back(higgsLepton);
              continue;
            }
            if ( higgsLikePair.size() < 2 ) {
              higgsLikePair.push_back(higgsLepton);
              break;
            }
          }

          if( higgsLikePair.size() == 2 ){
            auto firstLepton  = higgsLikePair.at(0);
            auto secondLepton = higgsLikePair.at(1);

            higgsCandId=firstLepton.pdgId()*secondLepton.pdgId();
            treeHiggsPairCharge = (int) (higgsCandId<0 ? category::tausPairCharge::OS : category::tausPairCharge::SS);

            switch(abs(higgsCandId)){
              case 11*11:  // ChannelName  = "elel";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::ee ;
                          break;
              case 13*13:  // ChannelName  = "mumu";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::mm;
                          break;
              case 11*13:  // ChannelName  = "elmu";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::em;
                          break;
              case 11*15:  // ChannelName  = "elha";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::et;
                          break;
              case 13*15:  // ChannelName  = "muha";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::mt;
                          break;
              case 15*15:  // ChannelName  = "haha";
                          treeHiggsPairFlavour = (int) category::tausPairFlavour::tt;
                          break;
              default:     // ChannelName  = "none";
                          break;
            }

            treeLep3Candidate = firstLepton.p4();
            treeLep3PdgId     = firstLepton.pdgId();
            if (isMC) treeLep3GenMatch= (int) ntupleutils::analysis::gen_truth::LeptonGenMatch(firstLepton.p4(), gen).first;
            treeLep3ID        = leptonIDmap(firstLepton, vtx[0]);
            treeLep3Iso       = (abs(treeLep3PdgId)!=15)? patUtils::relIso(firstLepton,rho): firstLepton.tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            treeLep3Eff       = 1;
            treeTMass         = TMath::Sqrt(2* firstLepton.pt()*met.pt()*(1-TMath::Cos(deltaPhi(met.phi(), firstLepton.phi()))));

            treeLep4Candidate = secondLepton.p4();
            treeLep4PdgId     = secondLepton.pdgId();
            if (isMC) treeLep4GenMatch= (int) ntupleutils::analysis::gen_truth::LeptonGenMatch(secondLepton.p4(), gen).first;
            treeLep4ID        = leptonIDmap(secondLepton, vtx[0]);
            treeLep4Iso       = (abs(treeLep4PdgId)!=15)? patUtils::relIso(secondLepton,rho): secondLepton.tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            treeLep4Eff       = 1;

            //reweight the event to account for lept eff.
            if(isMC && abs(firstLepton.pdgId())<15){

              int id( abs(firstLepton.pdgId()) );

              if(id==11)treeLep3Eff *= lepEff.getRecoEfficiency( firstLepton.el.superCluster()->eta(), id).first; //Reconstruction eff
              else if(id==13)treeLep3Eff *= lepEff.getTrackingEfficiency( firstLepton.eta(), id).first; //Tracking eff
              treeLep3Eff *= isMC ? lepEff.getLeptonEfficiency( firstLepton.pt(), firstLepton.eta(), id,  id ==11 ? "loose"    : "loose"   ,patUtils::CutVersion::Moriond17Cut).first : 1.0; //ID
              if(id==13){ treeLep3Eff *= isMC ? lepEff.getLeptonEfficiency( firstLepton.pt(), firstLepton.eta(), id, "looseiso_looseid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;} //ISO w.r.t ID
            }

            if(isMC && abs(secondLepton.pdgId())<15){

              int id( abs(secondLepton.pdgId()) );

              if(id==11)treeLep4Eff *= lepEff.getRecoEfficiency( secondLepton.el.superCluster()->eta(), id).first; //Reconstruction eff
              else if(id==13)treeLep4Eff *= lepEff.getTrackingEfficiency( secondLepton.eta(), id).first; //Tracking eff
              treeLep4Eff *= isMC ? lepEff.getLeptonEfficiency( secondLepton.pt(), secondLepton.eta(), id,  id ==11 ? "loose"    : "loose"   ,patUtils::CutVersion::Moriond17Cut ).first : 1.0; //ID
              if(id==13){ treeLep4Eff *= isMC ? lepEff.getLeptonEfficiency( secondLepton.pt(), secondLepton.eta(), id, "looseiso_looseid",patUtils::CutVersion::Moriond17Cut ).first : 1.0;} //ISO w.r.t ID
            }

           treeDiTauSVfitCandidate = getClassicSVFit(met, firstLepton, secondLepton);
	         }


        }

        treeWeightNoLepSF = weightNoLepSF;
        treeWeight = weight;

        eventWeights->SetLeptonsWeights(weight/weightNoLepSF);

        if(tree) tree->Fill();
        delete eventWeights;
        delete lheSummary;
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
