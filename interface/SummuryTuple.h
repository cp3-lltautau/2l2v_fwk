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
  VAR(unsigned int        , totalNumberOfEvents) \
  VAR(float        , sumOfGenWeights) \
  VAR(float        , sumOfGenSquaredWeights) \



struct SummuryTuple{
#define VAR(type, name) type name;
TREE_STRUCTURE()
#undef VAR


  void reset(){

    totalNumberOfEvents=0;
    sumOfGenWeights=0;
    sumOfGenSquaredWeights=0;
  }
};

class llttSummuryTree{

public:

  llttSummuryTree(TDirectory* _directory, const std::string& _name ):
  directory(_directory), name(_name){
    tree = new TTree(name, name);
    // cout<<"llttSummuryTree CREATED"<<endl;
    tree->SetDirectory(directory);
    cout<<"Initialising llttSummuryTree!!  "<<  tree <<  " name: "<<name<< endl;

    tree->Branch("totalNumberOfEvents", &ntuple.totalNumberOfEvents , string("totalNumberOfEvents/i" ).c_str());
    tree->Branch("sumOfGenWeights", &ntuple.sumOfGenWeights , string("sumOfGenWeights/F" ).c_str());
    tree->Branch("sumOfGenSquaredWeights", &ntuple.sumOfGenSquaredWeights , string("sumOfGenSquaredWeights/F" ).c_str());
  }

  ~llttSummuryTree(){
    if(directory) directory->Delete(name);
    else delete tree;
  }

  void fillEventInfo (unsigned int  _totalNumberOfEvents, double  _sumOfGenWeights, double _sumOfGenSquaredWeights){
    ntuple.totalNumberOfEvents = _totalNumberOfEvents;
    ntuple.sumOfGenWeights  = _sumOfGenWeights;
    ntuple.sumOfGenSquaredWeights   = _sumOfGenSquaredWeights;
  }


   void Fill() {
     tree->Fill();
     ntuple.reset();
   }

   SummuryTuple GetTuple(){ return ntuple;}

   void Write()
    {
        if(directory)
            directory->WriteTObject(tree, tree->GetName(), "WriteDelete");
    }


private:
  TDirectory* directory;
  TString name;

  TTree* tree;
  SummuryTuple ntuple;
};
