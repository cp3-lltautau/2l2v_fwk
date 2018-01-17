#ifndef ntupleutils_h
#define ntupleutils_h

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/PtrVector.h"

//need for the good lumi filter
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/LumiUtils.h"

// Electron ID
#include "RecoEgamma/ElectronIdentification/interface/VersionedPatElectronSelector.h"

#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include <Math/VectorUtil.h>

#define ENUM_CLASS(name, ...) \
namespace category{\
  enum class name { nodef, __VA_ARGS__ };\
}\

namespace ntupleutils{
  ENUM_CLASS(ZFlavour, ee, mm );
  ENUM_CLASS(tausPairFlavour, ee, em, et, mm, mt, tt );
  ENUM_CLASS(tausPairCharge, OS, SS );

  class EventInfo: public TObject{

  public:
    EventInfo():
     EventId( 0 ),
     LumiId( 0 ),
     RunId(0),
     zFlavour(category::ZFlavour::nodef),
     tausFlavour(category::tausPairFlavour::nodef),
     tausCharge(category::tausPairCharge::nodef) {
       weight = 0.0;
       weightNoLepSF = 0.0;
     }

    EventInfo(fwlite::Event& ev):
     EventId( ev.eventAuxiliary().event() ),
     LumiId( ev.eventAuxiliary().luminosityBlock() ),
     RunId(ev.eventAuxiliary().run()),
     zFlavour(category::ZFlavour::nodef),
     tausFlavour(category::tausPairFlavour::nodef),
     tausCharge(category::tausPairCharge::nodef) {
       weight = 0.0;
       weightNoLepSF = 0.0;
     }

     EventInfo(unsigned int _EventId, unsigned int _LumiId, unsigned int _RunId):
      EventId( _EventId ),
      LumiId( _LumiId ),
      RunId( _RunId ),
      zFlavour(category::ZFlavour::nodef),
      tausFlavour(category::tausPairFlavour::nodef),
      tausCharge(category::tausPairCharge::nodef) {
        weight = 0.0;
        weightNoLepSF = 0.0;
      }

    ~EventInfo(){}

    unsigned int GetEventID() const{return EventId;}
    unsigned int GetLumiID () const{return LumiId;}
    unsigned int GetRunID  () const{return RunId;}

    category::ZFlavour         GetZFlavour() const{return zFlavour;}
    category::tausPairFlavour  GetTausPairFlavour () const{return tausFlavour;}
    category::tausPairCharge   GetTausPairCharge () const{return tausCharge;}

    double GetEventWeight() const {return weight;}
    double GetEventWeightNoLepSF() const {return weightNoLepSF;}

    void SetZFlavour (category::ZFlavour zf){ zFlavour = zf; }
    void SetTausPairFlavour (category::tausPairFlavour tf){ tausFlavour = tf; }
    void SetTausPairCharge  (category::tausPairCharge tc){ tausCharge = tc; }
    void SetEventWeight (const float _weight) { weight = _weight; }
    void SetEventWeightNoLepSF (const float _weightNoLepSF ) { weightNoLepSF = _weightNoLepSF; }

  private:
    unsigned int    EventId;
    unsigned int    LumiId;
    unsigned int    RunId;
    category::ZFlavour        zFlavour;
    category::tausPairFlavour tausFlavour;
    category::tausPairCharge  tausCharge;

    float weight;
    float weightNoLepSF;

    ClassDef(EventInfo,1);
    // Add weights
  };

}

#endif
