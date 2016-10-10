#include "UserCode/llvv_fwk/interface/FRWeights.h"
#include <TFile.h>
#include <TKey.h>
#include <TSystem.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

/*****************************************************************/
FRWeights::FRWeights()
/*****************************************************************/
{
}

/*****************************************************************/
FRWeights::~FRWeights()
/*****************************************************************/
{
}

/*****************************************************************/
bool FRWeights::init(const string& WeightsFileName)
/*****************************************************************/
{

  //std::cout<<"I AM INITIALIZING"<<std::endl;
  //std::cout<<"FRWeights::init() opening file"<<WeightsFileName.c_str()<<std::endl;

  WeightsFile = TFile::Open(WeightsFileName.c_str());
  if(!WeightsFile) {
    cout<<"FRWeights::init() ERROR: Cannot open weights file "<<WeightsFileName<<"\n";
    return false;
  } else {

    std::cout<<"FRWeights::init() retrieving weights from: "<<WeightsFileName.c_str()<<std::endl;
    
    std::vector<string> catL  = {"Fake e", "Fake #mu", "Fake #tau_{had}"};
    std::vector<string> binL  = {"Inc.", "Barrel", "Endcap", "Inc. (m_{T}>30)", "Barrel (m_{T}>30)", "Endcap (m_{T}>30)"};
    
    std::vector<string> cat   = {"FR_El", "FR_Mu", "FR_Ta"};
    std::vector<string> bin   = {"","_B","_E", "_TMCut", "_TMCut_B", "_TMCut_E"};
    std::vector<string> var   = {"_Id_Iso01weight", "_Id_Iso02weight", "_Id_Iso03weight", "_Id_IsoLoweight", "_Id_IsoMeweight"};
    std::vector<string> wrt   = {"_wrtJetPt", "_wrtLepPt"};
    
    Int_t theGraphs=0;

    for(unsigned int c=0;c<cat.size();c++){
      for(unsigned int b=0;b<bin.size();b++){
	for(unsigned int v=0;v<var.size();v++){
	  for(unsigned int w=0;w<wrt.size();w++){

	    //std::cout<<"Get:" << cat[c]+var[v]+bin[b]+wrt[w] << std::endl;

	    std::string theSearchString = cat[c]+var[v]+bin[b]+wrt[w];

	    TIter next((WeightsFile->GetListOfKeys()));
	    TKey *key;
	    while ((key = (TKey*)next())) {
	      // std::cout<<((TGraphErrors*)key->ReadObj())->GetName()<<std::endl;
	      if(((TGraphErrors*)key->ReadObj())->GetName()==("Graph_from_"+theSearchString)){
		//  std::cout<<((TGraphErrors*)key->ReadObj())->GetName()<<std::endl;
		if(FRWeightGraphs.count(theSearchString)==0){
		  TGraphErrors* theClone = (TGraphErrors*)(( key->ReadObj())->Clone() );
		  FRWeightGraphs.emplace(theSearchString,theClone);
		  theGraphs++;
		  //std::cout<<theGraphs<<std::endl;
		  //std::cout<<"Booking "<< theSearchString <<std::endl;
		} else {
		  std::cout<<" FRWeights::init() this key: "<< theSearchString<<" was already booked!"<<std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
    
    // std::cout<<" THE GRAPHS are:"<< theGraphs<<std::endl;
    // for (std::map<std::string,TGraphErrors*>::iterator it=FRWeightGraphs.begin(); it!= FRWeightGraphs.end(); ++it){
    //std::cout << it->first << " => " << (it->second)->GetName() << '\n';
    //std::cout<< (it->second)->Eval(10.) <<std::endl;
    //} 
    
    WeightsFile->Close();
    
    return true;
  } 
}

/*****************************************************************/
bool FRWeights::check(const bool& isverbose )
/*****************************************************************/
{
  if(isverbose){
    std::cout << "FRWeights::check() FRWeightGraphs.size() is " << FRWeightGraphs.size() << '\n';
    std::cout << "FRWeights::check() Available strings: "<<std::endl;
    for (std::map<std::string,TGraphErrors*>::iterator it=FRWeightGraphs.begin(); it!= FRWeightGraphs.end(); ++it){
      std::cout << it->first << " => " << (it->second)->GetName() << '\n';
      // std::cout << it->first << " => " << it->second <<std::endl;
      // std::cout<< (it->second)->Eval(10.) <<std::endl;
    }
  }
  
  if(FRWeightGraphs.size()!=0){
    return true;
  } else {
    return false;
  }
}


/*****************************************************************/
double FRWeights::getWeight(const std::string& cat ,const std::string& bin, const std::string& var,const std::string& wrt ,const double& pT)
/*****************************************************************/
{
 
  std::string theSearchString = cat+var+bin+wrt;
  //std::cout<<"Searching for string: "<<theSearchString<<std::endl;
   
  TGraphErrors* graph=NULL;
  
  for (std::map<std::string,TGraphErrors*>::iterator it=  FRWeightGraphs.begin(); it!=  FRWeightGraphs.end(); ++it){
    // std::cout<<"to Match:"<< it->first << " target: "<< theSearchString<<std::endl;
    if(it->first == theSearchString ){
      graph =it->second;
      std::cout<<"FRWeights::getWeight() retrieving weight from: " << graph->GetName()<<std::endl;
    } 
    // else {
    //   std::cout<<"couldn't match"<<std::endl;
    // }
  }

  double result=-1.;
  if(graph){
    result=(graph->Eval(pT));
    //  std::cout<<"got the graph"<<std::endl;
  } else {
    result = 1;
    std::cout<<"FRWeights::getWeight() didn't got the graph: falling back to no weight"<<std::endl;
  }

  //delete graph;
  return result;
}
