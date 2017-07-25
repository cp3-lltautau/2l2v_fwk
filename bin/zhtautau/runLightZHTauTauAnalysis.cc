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

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
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
#include "UserCode/llvv_fwk/interface/rochcor2016.h"
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
    if(dtag.Contains("DoubleMuon"))    filterOnlyMUMU=true;
  }

  bool is2016data = (!isMC && dtag.Contains("2016"));
  bool is2016MC = (isMC);

  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  // ---------> Look runZHTauTauAnalysis
  // -----------------------------------


  //ELECTROWEAK CORRECTION WEIGHTS
  // ---------> Look runZHTauTauAnalysis
  // -----------------------------------

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  printf("Definition of plots");

  //event selection
  TH1 *h1=mon.addHistogram( new TH1F ("eventflow", ";;Events", 10,0,10) );
  h1->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1->GetXaxis()->SetBinLabel(2,"Nlep#geq2");
  h1->GetXaxis()->SetBinLabel(3,"Zmass");
  h1->GetXaxis()->SetBinLabel(4,"Zkin");
  h1->GetXaxis()->SetBinLabel(5,"Nlep+Ntau#geq4");
  h1->GetXaxis()->SetBinLabel(6,"Lep Veto");
  h1->GetXaxis()->SetBinLabel(7,"Btag Veto");
  h1->GetXaxis()->SetBinLabel(8,"#Delta #phi Z-MET");
  h1->GetXaxis()->SetBinLabel(9,"di-#tau Cand");

  TH1 *h_tr= mon.addHistogram( new TH1F ("trigger", ";;Events", 10,0,10) );
  h_tr->GetXaxis()->SetBinLabel(1,"#mu#mu");
  h_tr->GetXaxis()->SetBinLabel(2,"#mu");
  h_tr->GetXaxis()->SetBinLabel(3,"ee");
  h_tr->GetXaxis()->SetBinLabel(4,"e");
  h_tr->GetXaxis()->SetBinLabel(5,"e#mu");
  h_tr->GetXaxis()->SetBinLabel(6,"#gamma");

  //Z leptons kinematics control
  mon.addHistogram( new TH1F( "leadpt"      ,  ";p_{T}^{lead} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta"     ,  ";#eta_{lead};Events", 50,-2.6,2.6) );
  // mon.addHistogram( new TH1F( "leadiso"     ,  ";#I_{lead};Events", 20,0.,1.) );
  mon.addHistogram( new TH1F( "trailerpt"   ,  ";p_{T}^{trail} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta"  ,  ";#eta_{trail};Events", 50,-2.6,2.6) );
  // mon.addHistogram( new TH1F( "traileriso"     ,  ";#I_{trail};Events", 20,0.,1.) );
  mon.addHistogram( new TH1F( "leppt"       ,  ";p_{T}^{lepton} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "lepeta"      ,  ";#eta_{lepton};Events", 50,-2.6,2.6) );

  // zll control
  mon.addHistogram( new TH1F( "zlly",      		";y_{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zlleta",    		";#eta_{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zllpt",     		";p_{T}^{ll} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "zllmass",   		";M_{ll} (GeV);Events/2 GeV", 80,20,180) );
  mon.addHistogram( new TH1F( "nlep",      	";Number of Leptons;Events", 10,0,10) );

  //pu control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) );
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) );
  mon.addHistogram( new TH1F( "nvtxpuweight",";Vertices;Events",50,0,50) );
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) );

  //fake rate histograms
  float ptbinsJets[] = {10, 20, 30, 40, 60, 80, 100, 125, 150, 175,250};
  int ptbinsJetsN = sizeof(ptbinsJets)/sizeof(float)-1;
  mon.addHistogram( new TH1F( "wrtJetPt",  ";Jet p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));
  mon.addHistogram( new TH1F( "wrtLepPt",  ";Lep p_{T} (GeV);Events",sizeof(ptbinsJets)/sizeof(float)-1,ptbinsJets));

  double bins[]={5, 30,70,110,190,300,550,1800};
  int nbins=sizeof(bins)/sizeof(double) - 1;

  mon.addHistogram( new TH1F( "Amass_fine",            ";M_{#tau#tau} (GeV);Events",200,0,1000));
  mon.addHistogram( new TH1F( "Hmass_fine",            ";M_{ll#tau#tau} (GeV);Events",200,0,1000));
  mon.addHistogram( new TH1F( "Amass",            ";M_{#tau#tau} (GeV);Events",nbins,bins));
  mon.addHistogram( new TH1F( "Hmass",            ";M_{ll#tau#tau} (GeV);Events",nbins,bins));

  //
  // HISTOGRAMS FOR OPTIMIZATION and STATISTICAL ANALYSIS
  //
  //

  std::vector<const char*> tauIDiso = {"byLooseCombinedIsolationDeltaBetaCorr3Hits","byLooseIsolationMVArun2v1DBoldDMwLT",
                                        "byLooseIsolationMVArun2v1DBdR03oldDMwLT"};
  std::vector<float>    optim_Cuts_sumPt;
  std::vector<int>      optim_Cuts_taIso;
  std::vector<float>    optim_Cuts_muIso;
  std::vector<float>    optim_Cuts_elIso;

  for(float elIso=0.3;elIso>=0.1;elIso-=0.1){
    for(float muIso=0.3;muIso>=0.1;muIso-=0.1){
      for(int taIso=0;taIso<tauIDiso.size();taIso++){
        for(float sumPt=0;sumPt<=200;sumPt+=20){
          optim_Cuts_elIso.push_back(elIso);
          optim_Cuts_muIso.push_back(muIso);
          optim_Cuts_taIso.push_back(taIso);
          optim_Cuts_sumPt.push_back(sumPt);
        }
      }
    }
  }

  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(), 4, 0, 4)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "eIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(2, "muIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(3, "tauIso<");
  Hoptim_cuts->GetYaxis()->SetBinLabel(4, "sumPt>");

  for(unsigned int index=0;index<optim_Cuts_sumPt.size();index++){
    Hoptim_cuts->Fill(index,0.0,optim_Cuts_elIso[index]); // mon.fillHisto("optim_cut","",)
    Hoptim_cuts->Fill(index,1.0,optim_Cuts_muIso[index]); // mon.fillHisto("optim_cut","",)
    Hoptim_cuts->Fill(index,2.0,optim_Cuts_taIso[index]); // mon.fillHisto("optim_cut","",)
    Hoptim_cuts->Fill(index,3.0,optim_Cuts_sumPt[index]); // mon.fillHisto("optim_cut","",)
  }

  mon.addHistogram( new TH2F (TString("Hsvfit_shapes"),";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );
  mon.addHistogram( new TH2F (TString("Asvfit_shapes"),";cut index;M_{ll#tau#tau};Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),nbins,bins) );


  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  //MC normalization (to 1/pb)
  double xsecWeight = 1.0;
  if(isMC) xsecWeight=xsec/utils::getTotalNumberOfEvents(urls, false, true);//need to use the slow method in order to take NLO negative events into account

  //MET CORRection level
  pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;

  //jet energy scale and uncertainties
  // ---------> Look runZHTauTauAnalysis
  // -----------------------------------

  //muon energy scale and uncertainties
  rochcor2016* muCor2016 = new rochcor2016();

  //electron energy scale and uncertainties
  EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele");
  eScaler.doScale=true;
  eScaler.doSmearings=true;

  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  // ---------> Look runZHTauTauAnalysis
  // -----------------------------------

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

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
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

      //Skip bad lumi
      if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock())) continue;

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

        if (ngenITpu == 0) continue; // It prevents to fill vtxraw with -nan values

        puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
        puWeightUp  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
        puWeightDown = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
        weight *= puWeight;

        //GEN LEVEL FILTERING
        // ---------> Look runZHTauTauAnalysis
        // -----------------------------------
      }// gen event weight and PU weight

      //apply trigger and require compatibilitiy of the event with the PD
      edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
      if(!tr.isValid())return false;

      bool mumuTrigger(true); bool muTrigger(true);	bool eeTrigger(true); bool eTrigger(true); bool emuTrigger(true);

      patUtils::MetFilter metFilter;

      int metFilterValue = 0;

  	  bool filterbadPFMuon = true;
  	  bool filterbadChCandidate = true;
  	  bool filterbadMuonHIP = true;
  	  bool filterduplicateMuonHIP = true;
  	  std::unique_ptr<std::vector<reco::Muon*>> outbadMuon(new std::vector<reco::Muon*>());
  	  std::unique_ptr<std::vector<reco::Muon*>> outduplicateMuon(new std::vector<reco::Muon*>());

      if(is2016data || is2016MC){
        // cout << " TRIGGER PATHs \n ";
        // for (auto &triggerName : tr.triggerNames() ){
        //   if ( tr.accept(triggerName) ) cout << triggerName << " --- " << tr.accept(triggerName) << endl;
        // }
        mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                                                        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu22_v*","HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
        eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //,"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",);
        eTrigger           = utils::passTriggerPatterns(tr, "HLT_Ele27_eta2p1_WPLoose_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*");
        metFilterValue = metFilter.passMetFilterInt( ev, is2016data );

  	    // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
  	    filterbadChCandidate = metFilter.passBadChargedCandidateFilter(ev); if (!filterbadChCandidate) {  metFilterValue=9; }
  	    filterbadPFMuon = metFilter.passBadPFMuonFilter(ev); if (!filterbadPFMuon) { metFilterValue=8; }
  	    filterbadMuonHIP = metFilter.BadGlobalMuonTaggerFilter(ev,outbadMuon,false); if (!filterbadMuonHIP) { metFilterValue=10; }
  	    filterduplicateMuonHIP = metFilter.BadGlobalMuonTaggerFilter(ev,outduplicateMuon,true); if (!filterduplicateMuonHIP) { metFilterValue=11; }
      }

      bool passTrigger        = mumuTrigger || eeTrigger;

      if(  mumuTrigger)mon.fillHisto("trigger", "raw", 0 , weight);
      if(    muTrigger)mon.fillHisto("trigger", "raw", 1 , weight);
      if(    eeTrigger)mon.fillHisto("trigger", "raw", 2 , weight);
      if(     eTrigger)mon.fillHisto("trigger", "raw", 3 , weight);

      if(!isMC && passTrigger){ //avoid double counting of events from different PD
        if(filterOnlyMUMU)     { passTrigger = mumuTrigger;}
      }

      if(passTrigger){
        if(  mumuTrigger)mon.fillHisto("trigger", "cleaned", 0 , weight);
        if(    muTrigger)mon.fillHisto("trigger", "cleaned", 1 , weight);
        if(    eeTrigger)mon.fillHisto("trigger", "cleaned", 2 , weight);
        if(     eTrigger)mon.fillHisto("trigger", "cleaned", 3 , weight);
      }

      //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
      if(!passTrigger)continue;
      //##############################################   EVENT PASSED THE TRIGGER   ######################################
  	  if (metFilterValue==10 || metFilterValue==11) { metFilterValue=0; }
      if( metFilterValue!=0 ) continue;	 //Note this must also be applied on MC
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

     pat::TauCollection taus;
     fwlite::Handle< pat::TauCollection > tausHandle;
     tausHandle.getByLabel(ev, "slimmedTaus");
     if(tausHandle.isValid()){ taus = *tausHandle;}

     EcalRecHitCollection recHitCollectionEB;
     EcalRecHitCollection recHitCollectionEE;
     fwlite::Handle<EcalRecHitCollection> recHitCollectionEBHandle;
     fwlite::Handle<EcalRecHitCollection> recHitCollectionEEHandle;
     recHitCollectionEBHandle.getByLabel(ev, "reducedEgamma","reducedEBRecHits" );
     recHitCollectionEEHandle.getByLabel(ev, "reducedEgamma","reducedEERecHits" );
     if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}
     if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

     //
     // LEPTON ANALYSIS
     //

     //start by merging electrons and muons
     std::vector<patUtils::GenericLepton> leptons;
     for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}
     for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}
     std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

     std::vector<patUtils::GenericLepton> selLeptons;

     for(size_t ilep=0; ilep<leptons.size(); ilep++){
       bool passKin(true),passId(true),passIso(true);
       bool passLooseLepton(true);

       int lid=abs( leptons[ilep].pdgId() );

       //veto leptons overlaping with other lep
       bool overlapWithLepton=false;
       for(int l1=0; l1<(int)selLeptons.size();++l1){
         if(deltaR(leptons[ilep].p4(), selLeptons[l1])<0.1){overlapWithLepton=true; break;}
       }if(overlapWithLepton)continue;

       //Cut based identification
       // passId
       passId = lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut, true) :
       patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut);
       //isolation
       //  passIso
       passIso = lid==11 ? patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::CutSet::ICHEP16Cut) :
       patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::CutSet::Moriond17Cut);

       //apply muon & electron corrections
       // ---------> Look runZHTauTauAnalysis
       // -----------------------------------

       // Compute relIso after corrections
       leptons[ilep].addUserFloat("relIso",  patUtils::relIso(leptons[ilep], rho) ); //compute it once for all

       //kinematics
       float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
       if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
       if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
       if(leptons[ilep].pt()<25) passKin=false;

       if(passId && passIso && passKin)  selLeptons.push_back(leptons[ilep]);
     }// end leptons selection

     std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);

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
     // ASSIGN CHANNEL
     //

     std::vector<TString> chTags;
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
       for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
         if(abs(selLeptons[l2].pdgId())==15)continue;
         if( selLeptons[l2].pt()<20 ) continue;

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

       if(isMC && abs(dilId)==169)weight *= lepEff.getTriggerEfficiencySF(selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), dilId,is2016MC).first;
		   if(isMC && abs(dilId)==121)weight *= lepEff.getTriggerEfficiencySF(selLeptons[dilLep1].pt(), selLeptons[dilLep1].el.superCluster()->eta(), selLeptons[dilLep2].pt(), selLeptons[dilLep2].el.superCluster()->eta(), dilId,is2016MC).first;
     }

     if(!isDileptonCandidate) continue;
     //##############################################   EVENT PASSED DiLepton CANDIDATE   #######################################
     bool passZmass = (fabs(zll.mass()-91.2)<15);
     bool passZpt   = (zll.pt()>20);
     std::vector<TString> chTagsMain=chTags;

     int higgsCandL1=-1, higgsCandL2=-1;
     LorentzVector higgsCand;
     int HiggsShortId=-1, higgsCandId;

     int NCleanedJetMain = 0;
     bool passDPhiCut    = 0;
     bool passHiggsLoose = 0;
     bool passHiggsMain  = 0;
     LorentzVector higgsCand_SVFit;
     LorentzVector higgsCandH;
     LorentzVector higgsCandH_SVFit;

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
         //std::cout << "SelLep num. "<<l<<":  "<< selLeptons[l].pdgId() << " - Momentum: "<< selLeptons[l].p4() << std::endl;
         if(l==dilLep1 || l==dilLep2)continue;
         if(higgsCandL1<0){higgsCandL1=l;continue;}
         if(higgsCandL2<0){higgsCandL2=l;break;}//ordered in pT, so all done
       }

       //std::cout << "---  Leading Higgs Lep: "<< selLeptons[higgsCandL1].pdgId() << " - Momentum: "<< selLeptons[higgsCandL1].p4() << std::endl;
       //std::cout << "--- Trailing Higgs Lep: "<< selLeptons[higgsCandL2].pdgId() << " - Momentum: "<< selLeptons[higgsCandL2].p4() << std::endl;

       string ChannelName = "none";   string signName = "";
       if(higgsCandL1>=0 && higgsCandL2>=0){
         higgsCandId=selLeptons[higgsCandL1].pdgId()*selLeptons[higgsCandL2].pdgId();
         higgsCand = LorentzVector(selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());
         if(higgsCandId<0){signName="_OS";}else{signName="_SS";}
         if(higgsCandId<0){HiggsShortId = 0;}else{HiggsShortId =12;}
         if(abs(selLeptons[dilLep1].pdgId())==11){HiggsShortId += 0;}else{HiggsShortId += 6;}
         switch(abs(higgsCandId)){
           case 11*11:  ChannelName  = "elel";  HiggsShortId+= 0; break;
           case 13*13:  ChannelName  = "mumu";  HiggsShortId+= 1; break;
           case 11*13:  ChannelName  = "elmu";  HiggsShortId+= 2; break;
           case 11*15:  ChannelName  = "elha";  HiggsShortId+= 3; break;
           case 13*15:  ChannelName  = "muha";  HiggsShortId+= 4; break;
           case 15*15:  ChannelName  = "haha";  HiggsShortId+= 5; break;
           default:     ChannelName  = "none";  HiggsShortId =-1; break;
         }
       }
       //std::cout << " ---------- Higgs flavour:  "<< signName + ChannelName << std::endl;
       chTagsMain.push_back(chTagsMain[chTagsMain.size()-1] + signName + ChannelName);
       //Get the Higgs candidate (end)

       //printf("%30s %2i --> %i %i %i %i\n", (chTagsMain[chTagsMain.size()-1]).Data(), HiggsShortId, selLeptons[dilLep1].pdgId(), selLeptons[dilLep2].pdgId(), selLeptons[higgsCandL1].pdgId(), selLeptons[higgsCandL2].pdgId());

       //reweight the event to account for lept eff.



       //check how many additional light jets are present


       // passDPhiCut    =  (fabs(deltaPhi(zll.phi(), met.phi()))>1.5);
       passHiggsLoose = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.3, 0.3, "decayModeFinding", 0., false, vtx);
       // passHiggsMain  = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, 0.3, 0.15, "byLooseCombinedIsolationDeltaBetaCorr3Hits", 20., true, vtx);

       //SVFIT MASS
       higgsCand_SVFit = higgsCand;

       //FIXME gives a lot of warning currently
    //  if(passZmass && passZpt && passDPhiCut && passHiggsLoose && passLepjetMain && passBJetVetoMain){
      // std::cout<<"START SVFIT\n";
    //   higgsCand_SVFit = getSVFit(met, selLeptons, higgsCandL1, higgsCandL2);  //compute svfit mass in a smart way
     //  std::cout<<"END SVFIT\n";
    //  }

       //build the higgs candH
       higgsCandH = zll + higgsCand;
       higgsCandH_SVFit = zll + higgsCand_SVFit;
     }


     mon.fillHisto("eventflow"       , chTagsMain,                 0, weight);
     if(selLeptons.size()>=2){
       mon.fillHisto("nlep"           ,   chTags, selLeptons.size(), weight);
       mon.fillHisto("eventflow"     ,   chTagsMain,                 1, weight);
       mon.fillHisto("zllmass"          ,   chTagsMain, zll.mass(),    weight);
       if(passZmass){
         mon.fillHisto("eventflow"   ,   chTagsMain,                 2, weight);
         //pu control
         mon.fillHisto("nvtx"        ,   chTagsMain, vtx.size(),      weight);
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
           mon.fillHisto("eventflow",   chTagsMain,                 3, weight);
           if(selLeptons.size()>=4){
             mon.fillHisto("eventflow",   chTagsMain,                 4, weight);
             if(passHiggsLoose){
               mon.fillHisto("Amass"           , chTagsMain, higgsCand.mass(),  weight);
               mon.fillHisto("Hmass"           , chTagsMain, higgsCandH.mass(),  weight);
               mon.fillHisto("Amass_fine"           , chTagsMain, higgsCand.mass(),  weight);
               mon.fillHisto("Hmass_fine"           , chTagsMain, higgsCandH.mass(),  weight);
             }
           }
         }
       }
     }

     if( selLeptons.size()>=2 && passZmass && passZpt && selLeptons.size()>=4 && passHiggsLoose){
       for(unsigned int index=0; index<optim_Cuts_sumPt.size();index++){
         bool passHiggs = passHiggsCuts(selLeptons, higgsCandL1, higgsCandL2, optim_Cuts_elIso[index], optim_Cuts_muIso[index], tauIDiso[optim_Cuts_taIso[index]], optim_Cuts_sumPt[index],true,vtx);
         std::cout<<" Leptons cut used for higgs selection optimization:\n"
                  <<" IsoEle = "<<optim_Cuts_elIso[index]<<" IsoMu = "<<optim_Cuts_muIso[index]
                  <<" IsoTau = "<<tauIDiso[optim_Cuts_taIso[index]] <<" SumPt = "<<optim_Cuts_sumPt[index]<<std::endl;
         std::cout<<" Higgs cut status: "<< (passHiggs ? "PASSED" : "NOT PASSED") << std::endl;
         std::cout<<" Masses values:   A = "<<higgsCand_SVFit.mass()<<" GeV --  H = "<<higgsCandH_SVFit.mass()<<" GeV"<<std::endl;
         if(passHiggs){
           mon.fillHisto(TString("Hsvfit_shapes"),chTagsMain,index,higgsCandH_SVFit.mass(),weight);
           mon.fillHisto(TString("Asvfit_shapes"),chTagsMain,index,higgsCand_SVFit.mass(),weight);
         }
       }
     }

    } // event loop closes here
    printf("\n");
    delete file;
  } // loop on files closes here

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

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
    terminationCmd += TString("mv out.json ") + ((outUrl.ReplaceAll(".root",""))+".json") + ";";
    goodLumiFilter.FindLumiInFiles(urls);
    goodLumiFilter.DumpToJson("out.json");
  }

  system(terminationCmd.Data());
} // main closes here
