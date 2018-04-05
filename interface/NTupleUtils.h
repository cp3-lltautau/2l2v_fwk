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


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
       is3Lep = false;
       is4Lep = false;
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

    bool eventHas3Leptons() const {return is3Lep;}
    bool eventHas4Leptons() const {return is4Lep;}

    void SetZFlavour (category::ZFlavour zf){ zFlavour = zf; }
    void SetTausPairFlavour (category::tausPairFlavour tf){ tausFlavour = tf; }
    void SetTausPairCharge  (category::tausPairCharge tc){ tausCharge = tc; }
    void SetEventWeight (const float _weight) { weight = _weight; }
    void SetEventWeightNoLepSF (const float _weightNoLepSF ) { weightNoLepSF = _weightNoLepSF; }
    void SetEventLeptonsNum (const int _nlep) {
      if( _nlep == 3 ) is3Lep  = true;
      else if( _nlep >= 4)  is4Lep  = true;
    }

  private:
    unsigned int    EventId;
    unsigned int    LumiId;
    unsigned int    RunId;
    category::ZFlavour        zFlavour;
    category::tausPairFlavour tausFlavour;
    category::tausPairCharge  tausCharge;
    bool is3Lep;
    bool is4Lep;

    float weight;
    float weightNoLepSF;

    ClassDef(EventInfo,1);
    // Add weights
  };


  class Weights: public TObject{

  public:

    Weights():
      xsecWeight(1),puWeight(1),
      puWeightUp(1),puWeightDown(1),
      mcWeight(1){

      ewkCorrectionsWeight     = 1;
      ZZ_NNLOcorrectionsWeight = 1;

      // lep1weight = 1;
      // lep2weight = 1;
      // lep3weight = 1;
      // lep4weight = 1;

      leptonsWeight = 1;
    }

    ~Weights(){}

    void SetXSecWeight (const double _xsecWeight) { xsecWeight = _xsecWeight; }
    void SetPUWeights (const double _puWeight, const double _puWeightUp = 1, const double _puWeightDown = 1) {
      puWeight     = _puWeight;
      puWeightUp   = _puWeightUp;
      puWeightDown = _puWeightDown;
    }
    void SetMCWeight    (const double _mcWeight)  { mcWeight = _mcWeight; }
    void SetEwkMCWeight (const double _ewkCorrectionsWeight) { ewkCorrectionsWeight = _ewkCorrectionsWeight; }
    void SetZZMCWeight  (const double _ZZ_NNLOcorrectionsWeight) { ZZ_NNLOcorrectionsWeight = _ZZ_NNLOcorrectionsWeight; }

    void SetLeptonsWeights (const double _leptonsWeight) { leptonsWeight = _leptonsWeight; }


    double GetPUWeight     () const { return puWeight; }
    double GetPUWeightUp   () const { return puWeightUp; }
    double GetPUWeightDown () const { return puWeightDown; }
    std::pair<double, double> GetPUWeightShifts () const { return std::make_pair(puWeightUp, puWeightDown); }

    double GetXSecWeight   () const { return xsecWeight; }
    double GetMCWeight     () const { return mcWeight; }
    double GetEwkMCWeight  () const { return ewkCorrectionsWeight; }
    double GetZZMCWeight   () const { return ZZ_NNLOcorrectionsWeight; }
    double GetAllMCWeights () const { return mcWeight*ewkCorrectionsWeight*ZZ_NNLOcorrectionsWeight*xsecWeight;}

    double GetLeptonsWeights () const { return leptonsWeight; }

    double GetEventWeight () const { return mcWeight*ewkCorrectionsWeight*ZZ_NNLOcorrectionsWeight*xsecWeight*leptonsWeight*puWeight; }

  private:
    double xsecWeight;
    double puWeight;
    double puWeightUp;
    double puWeightDown;

    double mcWeight;
    double ewkCorrectionsWeight;
    double ZZ_NNLOcorrectionsWeight;

    // double lep1weight;
    // double lep2weight;
    // double lep3weight;
    // double lep4weight;

    double leptonsWeight;

      ClassDef(Weights,1);
  };
  // class ParticleCandidate: public TObject{
  //
  // public:
  //   ParticleCandidate():
  //   p4(LorentzVector(0.,0.,0.,0.)),
  //   charge(-9), pdgId(0),
  //   isolation(-99),
  //   IDmap(0x0){}
  //
  //   ParticleCandidate(LorentzVector _p4, int _charge, int _pdgId, double _isolation, char _IDmap;):
  //   p4(_p4), charge(_charge), pdgId(_pdgId), isolation(_isolation), IDmap(_IDmap){}
  //
  //   ParticleCandidate()
  //
  // private:
  //   LorentzVector p4;
  //   int charge;
  //   int pdgId;
  //   double isolation;
  //   char IDmap;
  //
  //   ClassDef(ParticleCandidate,1);
  // };

  /*! Tools for working with MC generator truth.
This file is part of https://github.com/hh-italian-group/h-tautau. */


namespace analysis {

namespace gen_truth {

enum class GenMatch { Electron = 1, Muon = 2, TauElectron = 3,  TauMuon = 4, Tau = 5, NoMatch = 6 };

namespace tools {
template<typename Type>
std::set<Type> union_sets(std::initializer_list<std::set<Type>> sets)
{
    std::set<Type> result;
    for(const auto& set : sets)
        result.insert(set.begin(), set.end());
    return result;
}
}

using MatchResult = std::pair<GenMatch, const reco::GenParticle*>;

inline void FindFinalStateDaughters(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& daughters,
                                    const std::set<int>& pdg_to_exclude = {})
{
    if(!particle.daughterRefVector().size()) {
        const int abs_pdg = std::abs(particle.pdgId());
        if(!pdg_to_exclude.count(abs_pdg))
            daughters.push_back(&particle);
    } else {
        for(const auto& daughter : particle.daughterRefVector())
            FindFinalStateDaughters(*daughter, daughters, pdg_to_exclude);
    }
}

inline LorentzVector GetFinalStateMomentum(const reco::GenParticle& particle, bool excludeInvisible,
                                              bool excludeLightLeptons)
{
    using set = std::set<int>;
    using pair = std::pair<bool, bool>;
    static const set empty = {};
    static const set light_leptons = { 11, 13 };
    static const set invisible_particles = { 12, 14, 16 };
    static const set light_and_invisible = tools::union_sets({light_leptons, invisible_particles});

    static const std::map<pair, const set*> to_exclude {
        { pair(false, false), &empty }, { pair(true, false), &invisible_particles },
        { pair(false, true), &light_leptons }, { pair(true, true), &light_and_invisible },
    };

    std::vector<const reco::GenParticle*> daughters;
    FindFinalStateDaughters(particle, daughters, *to_exclude.at(pair(excludeInvisible, false)));

    LorentzVector p4;
    for(auto daughter : daughters){
	if(excludeLightLeptons && light_leptons.count(std::abs(daughter->pdgId())) && daughter->statusFlags().isDirectTauDecayProduct()) continue;
	p4 += daughter->p4();
    }
    return p4;
}

template<typename LVector>
MatchResult LeptonGenMatch(const LVector& p4, const std::vector<reco::GenParticle>& genParticles)
{
    static constexpr int electronPdgId = 11, muonPdgId = 13, tauPdgId = 15;
    static constexpr double dR2_threshold = std::pow(0.2, 2);

    static const std::map<int, double> pt_thresholds = {
        { electronPdgId, 8 }, { muonPdgId, 8 }, { tauPdgId, 15 }
    };

    using pair = std::pair<int, bool>;
    static const std::map<pair, GenMatch> genMatches = {
        { { electronPdgId, false }, GenMatch::Electron }, { { electronPdgId, true }, GenMatch::TauElectron },
        { { muonPdgId, false }, GenMatch::Muon }, { { muonPdgId, true }, GenMatch::TauMuon },
        { { tauPdgId, false }, GenMatch::Tau }, { { tauPdgId, true }, GenMatch::Tau }
    };

    MatchResult result(GenMatch::NoMatch, nullptr);
    double match_dr2 = dR2_threshold;



    for(const reco::GenParticle& particle : genParticles) {
        const bool isTauProduct = particle.statusFlags().isDirectPromptTauDecayProduct();
        if((!particle.statusFlags().isPrompt() && !isTauProduct) || !particle.statusFlags().isLastCopy()) continue;

        const int abs_pdg = std::abs(particle.pdgId());
        if(!pt_thresholds.count(abs_pdg)) continue;

        const auto particle_p4 = abs_pdg == tauPdgId ? GetFinalStateMomentum(particle, true, true) : particle.p4();

        const double dr2 = ROOT::Math::VectorUtil::DeltaR2(p4, particle_p4);
        if(dr2 >= match_dr2) continue;
        if(particle_p4.pt() <= pt_thresholds.at(abs_pdg)) continue;

        match_dr2 = dr2;
        result.first = genMatches.at(pair(abs_pdg, isTauProduct));
        result.second = &particle;
    }
    return result;
}

struct LheSummary {
    size_t n_partons = 0, n_b_partons = 0, n_c_partons = 0;
    double HT = 0.;
};

inline LheSummary ExtractLheSummary(const LHEEventProduct& lheEventProduct)
{
    static constexpr int c_quark = 4, b_quark = 5;
    static const std::set<int> quarks_and_gluons = { 1, 2, 3, 4, 5, 6, 21 };

    LheSummary summary;
    const lhef::HEPEUP& lheEvent = lheEventProduct.hepeup();
    const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
    // std::vector<analysis::LorentzVectorXYZ> h0_p4;
    for(size_t n = 0; n < lheParticles.size(); ++n) {
        const int absPdgId = std::abs(lheEvent.IDUP[n]);
        const int status = lheEvent.ISTUP[n];
        // if(absPdgId == H0) summary.m_H = lheParticles[n][4];
        // if(absPdgId == h0)
        //     h0_p4.push_back(analysis::LorentzVectorXYZ(lheParticles[n][0], lheParticles[n][1],
        //                                                       lheParticles[n][2], lheParticles[n][3]));
        if(status != 1 || !quarks_and_gluons.count(absPdgId)) continue;
        ++summary.n_partons;
        if(absPdgId == c_quark) ++summary.n_c_partons;
        if(absPdgId == b_quark) ++summary.n_b_partons;
        summary.HT += std::sqrt(std::pow(lheParticles[n][0], 2) + std::pow(lheParticles[n][1], 2));
    }
    // if(h0_p4.size() == 2) {
    //     const analysis::LorentzVectorXYZ H_p4 = h0_p4.at(0) + h0_p4.at(1);
    //     const auto boosted_h0 = ROOT::Math::VectorUtil::boost(h0_p4.at(0), H_p4.BoostToCM());
    //     summary.cosTheta_hh = ROOT::Math::VectorUtil::CosTheta(boosted_h0, ROOT::Math::Cartesian3D<>(0, 0, 1));
    //     summary.m_hh = H_p4.mass();
    // }
    return summary;
}

} // namespace gen_truth
} // namespace analysis

}

#endif
