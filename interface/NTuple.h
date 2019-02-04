#pragma once

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

#include "UserCode/llvv_fwk/interface/NTupleUtils.h"

#define TREE_STRUCTURE() \
  VAR(unsigned int , treeEventId) \
  VAR(unsigned int , treeLumiId) \
  VAR(unsigned int , treeRunId) \
  VAR(float        , treeNumEvents) \
  VAR(bitset<4>    , treeTriggerFiredBits) \
  VAR(bitset<7>    , treeUncertantiesBits) \
  VAR(int          , treeDoubleTriggerMatch) \
  VAR(int          , treeLepIsMatched) \
  VAR(int          , treeZId) \
  VAR(int          , treeHiggsPairCharge) \
  VAR(int          , treeHiggsPairFlavour) \
  VAR(float        , treeWeight) \
  VAR(float        , treeWeightNoLepSF) \
  VAR(float        , treeLep3Eff) \
  VAR(float        , treeLep4Eff) \
  VAR(int          , treeNlep) \
  VAR(int          , treeNbtag) \
  VAR(LorentzVector, treeJetCandidate) \
  VAR(LorentzVector, treeLep1Candidate) \
  VAR(int          , treeLep1GenMatch) \
  VAR(int          , treeLep1PdgId) \
  VAR(unsigned char, treeLep1ID) \
  VAR(float        , treeLep1Iso) \
  VAR(LorentzVector, treeLep2Candidate) \
  VAR(int          , treeLep2GenMatch) \
  VAR(int          , treeLep2PdgId) \
  VAR(unsigned char, treeLep2ID) \
  VAR(float        , treeLep2Iso) \
  VAR(LorentzVector, treeLep3Candidate) \
  VAR(int          , treeLep3GenMatch) \
  VAR(int          , treeLep3PdgId) \
  VAR(unsigned char, treeLep3ID) \
  VAR(float        , treeLep3Iso) \
  VAR(float        , treeTMass) \
  VAR(int          , treeLeg1LepID_loose) \
  VAR(int          , treeLeg1LepID_medium) \
  VAR(int          , treeLeg1LepID_tight) \
  VAR(int          , treeLeg1TauIsoLooseMVA) \
  VAR(int          , treeLeg1TauIsoMediumMVA) \
  VAR(LorentzVector, treeLep4Candidate) \
  VAR(int          , treeLep4GenMatch) \
  VAR(int          , treeLep4PdgId) \
  VAR(unsigned char, treeLep4ID) \
  VAR(float        , treeLep4Iso) \
  VAR(LorentzVector, treeMET) \
  VAR(float        , treeMETCov00) \
  VAR(float        , treeMETCov01) \
  VAR(float        , treeMETCov10) \
  VAR(float        , treeMETCov11) \
  VAR(LorentzVector, treeDiTauSVfitCandidate) \



  //**************************************************************************************************************************//
  inline unsigned char leptonIDmap (patUtils::GenericLepton lepton, reco::Vertex& vtx)
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


struct EventNTuple{
#define VAR(type, name) type name;
TREE_STRUCTURE()
#undef VAR


  void reset(){
    LorentzVector null_p4(0., 0., 0., 0.);

    treeEventId=0;
    treeLumiId=0;
    treeRunId=0;
    treeDoubleTriggerMatch=0;
    treeLepIsMatched=0;
    treeTriggerFiredBits.reset();
    treeUncertantiesBits.reset();
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
    treeLep1ID = 0;
    treeLep1Iso = 0;
    treeLep2Candidate = null_p4;
    treeLep2GenMatch = -1;
    treeLep2PdgId = 0;
    treeLep2ID = 0;
    treeLep2Iso = 0;
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
    treeMET = null_p4;
    treeDiTauSVfitCandidate = null_p4;
  }
};

class llttTree{

public:

  llttTree(TDirectory* _directory, const std::string& _name ):
  directory(_directory), name(_name){
    tree = new TTree(name, name);
    // cout<<"llttTree CREATED"<<endl;
    tree->SetDirectory(directory);
    cout<<"Initialising llttTree!!  "<<  tree <<  " name: "<<name<< endl;
    // tree->Branch("EventWeights",&ntuple.eventWeights);
    // tree->Branch("LHESummary",&ntuple.lheSummary);
    tree->Branch("eventId", &ntuple.treeEventId , string("eventId/i" ).c_str());
    tree->Branch("lumiId" , &ntuple.treeLumiId  , string("lumiId/i"  ).c_str());
    tree->Branch("runId"  , &ntuple.treeRunId   , string("runId/i"   ).c_str());
    tree->Branch("numEvents" , &ntuple.treeNumEvents  , string("numEvents/F"  ).c_str());
    tree->Branch("triggerFiredBits",     &ntuple.treeTriggerFiredBits    , string("triggerFiredBits/b" ).c_str());
    tree->Branch("uncertantiesBits",     &ntuple.treeUncertantiesBits    , string("uncertantiesBits/b" ).c_str());
    tree->Branch("doubleTriggerMatch",     &ntuple.treeDoubleTriggerMatch    , string("doubleTriggerMatch/I" ).c_str());
    tree->Branch("lepIsMatched" , &ntuple.treeLepIsMatched , string("lepIsMatched/I").c_str());
    tree->Branch("zId",     &ntuple.treeZId     , string("zId/I" ).c_str());
    tree->Branch("higgsPairCharge",     &ntuple.treeHiggsPairCharge     , string("higgsPairCharge/I" ).c_str());
    tree->Branch("higgsPairFlavour",     &ntuple.treeHiggsPairFlavour     , string("higgsPairFlavour/I" ).c_str());
    tree->Branch("weight" , &ntuple.treeWeight  , string("weight/F"  ).c_str());
    tree->Branch("weightNoLepSF" , &ntuple.treeWeightNoLepSF  , string("weightNoLepSF/F"  ).c_str());
    tree->Branch("lep3Eff" , &ntuple.treeLep3Eff  , string("lep3Eff/F"  ).c_str());
    tree->Branch("lep4Eff" , &ntuple.treeLep4Eff  , string("lep4Eff/F"  ).c_str());
    tree->Branch("nlep" , &ntuple.treeNlep  , string("nlep/I"  ).c_str());
    tree->Branch("nbtag" , &ntuple.treeNbtag  , string("nbtag/I"  ).c_str());

    tree->Branch("jetMatched",     &ntuple.treeJetCandidate );
    /* Lep 1 */
    tree->Branch("lep1Candidate",     &ntuple.treeLep1Candidate );
    tree->Branch("lep1GenMatch",     &ntuple.treeLep1GenMatch     , string("lep1GenMatch/I" ).c_str());
    tree->Branch("lep1PdgId",     &ntuple.treeLep1PdgId     , string("lep1PdgId/I" ).c_str());
    tree->Branch("lep1ID", &ntuple.treeLep1ID , string("lep1ID/b" ).c_str());
    tree->Branch("lep1Iso", &ntuple.treeLep1Iso , string("lep1Iso/F" ).c_str());
    /* Lep 2 */
    tree->Branch("lep2Candidate",     &ntuple.treeLep2Candidate );
    tree->Branch("lep2GenMatch",     &ntuple.treeLep2GenMatch     , string("lep2GenMatch/I" ).c_str());
    tree->Branch("lep2PdgId",     &ntuple.treeLep2PdgId     , string("lep2PdgId/I" ).c_str());
    tree->Branch("lep2ID", &ntuple.treeLep2ID , string("lep2ID/b" ).c_str());
    tree->Branch("lep2Iso", &ntuple.treeLep2Iso , string("lep2Iso/F" ).c_str());
    /* Lep 3 */
    tree->Branch("lep3Candidate",     &ntuple.treeLep3Candidate );
    tree->Branch("lep3GenMatch",     &ntuple.treeLep3GenMatch     , string("lep3GenMatch/I" ).c_str());
    tree->Branch("lep3PdgId" , &ntuple.treeLep3PdgId  , string("lep3PdgId/I"  ).c_str());
    tree->Branch("lep3ID", &ntuple.treeLep3ID , string("lep3ID/b" ).c_str());
    tree->Branch("lep3Iso", &ntuple.treeLep3Iso , string("lep3Iso/F" ).c_str());
    tree->Branch("TMass" , &ntuple.treeTMass  , string("TMass/F"  ).c_str());
    tree->Branch("leg1LepIDloose", &ntuple.treeLeg1LepID_loose , string("leg1LepIDloose/I" ).c_str());
    tree->Branch("leg1LepIDmedium", &ntuple.treeLeg1LepID_medium , string("leg1LepIDmedium/I" ).c_str());
    tree->Branch("leg1LepIDtight", &ntuple.treeLeg1LepID_tight , string("leg1LepIDtight/I" ).c_str());
    tree->Branch("tau1LooseMVA", &ntuple.treeLeg1TauIsoLooseMVA , string("tau1LooseMVA/i" ).c_str());
    tree->Branch("tau1MediumMVA", &ntuple.treeLeg1TauIsoMediumMVA , string("tau1MediumMVA/i" ).c_str());
    /* Lep 4 */
    tree->Branch("lep4Candidate",     &ntuple.treeLep4Candidate );
    tree->Branch("lep4GenMatch",     &ntuple.treeLep4GenMatch     , string("lep4GenMatch/I" ).c_str());
    tree->Branch("lep4PdgId" , &ntuple.treeLep4PdgId  , string("lep4PdgId/I"  ).c_str());
    tree->Branch("lep4ID", &ntuple.treeLep4ID , string("lep4ID/b" ).c_str());
    tree->Branch("lep4Iso", &ntuple.treeLep4Iso , string("lep4Iso/F" ).c_str());
    /* di-tau system */
    tree->Branch("met",     &ntuple.treeMET );
    tree->Branch("METCov00", &ntuple.treeMETCov00 , string("METCov00/F" ).c_str());
    tree->Branch("METCov01", &ntuple.treeMETCov01 , string("METCov01/F" ).c_str());
    tree->Branch("METCov10", &ntuple.treeMETCov10 , string("METCov10/F" ).c_str());
    tree->Branch("METCov11", &ntuple.treeMETCov11 , string("METCov11/F" ).c_str());
    tree->Branch("diTauSVfitCandidate",     &ntuple.treeDiTauSVfitCandidate );
  }

  ~llttTree(){
    if(directory) directory->Delete(name);
    else delete tree;
  }

  void fillEventInfo (unsigned int  _eventId, unsigned int  _lumiId, unsigned int  _runId){
    ntuple.treeEventId = _eventId;
    ntuple.treeLumiId  = _lumiId;
    ntuple.treeRunId   = _runId;
  }

  void setHiggsLegs(patUtils::GenericLepton firstLepton, patUtils::GenericLepton secondLepton, reco::Vertex& vtx){
    ntuple.treeLep3Candidate = firstLepton.p4();
    ntuple.treeLep3PdgId     = firstLepton.pdgId();
    ntuple.treeLep3ID        = leptonIDmap(firstLepton, vtx);
    // ntuple.treeLep3Iso       = (abs(firstLepton.pdgId())!=15)? patUtils::relIso(firstLepton,rho): firstLepton.tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    ntuple.treeLep3Eff       = 1;

    ntuple.treeLep4Candidate = secondLepton.p4();
    ntuple.treeLep4PdgId     = secondLepton.pdgId();
    ntuple.treeLep4ID        = leptonIDmap(secondLepton, vtx);
    // ntuple.treeLep4Iso       = (abs(secondLepton.pdgId())!=15)? patUtils::relIso(secondLepton,rho): secondLepton.tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    ntuple.treeLep4Eff       = 1;
  }

   void Fill() {
     tree->Fill();
     ntuple.reset();
   }

   EventNTuple GetTuple(){ return ntuple;}

   void Write()
    {
        if(directory)
            directory->WriteTObject(tree, tree->GetName(), "WriteDelete");
    }


private:
  TDirectory* directory;
  TString name;

  TTree* tree;
  EventNTuple ntuple;
};
