#! /usr/bin/env python

import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle
from math import *
from array import array

def isAncestor(a,p) :
        if a == p :
                return True
        for i in xrange(0,p.numberOfMothers()) :
                if isAncestor(a,p.mother(i)) :
                         return True
        return False

class baseParticle:
    def __init__(self,prunedParticle):
        self.prunedParticle = prunedParticle
        self.pdgId = prunedParticle.pdgId()
        #Kinematics
        self.p4 = ROOT.TLorentzVector(prunedParticle.px(),prunedParticle.py(),prunedParticle.pz(),prunedParticle.energy())
        self.pt  = prunedParticle.pt()
        self.eta = prunedParticle.eta()
        self.phi = prunedParticle.phi()
        self.mass = prunedParticle.mass()
    ###
    # Overwrite comparison operators
    def __eq__(self,other):
        return self.pt == other.pt
    def __lt__(self,other):
        return self.pt < other.pt
    def __gt__(self,other):
        return self.pt > other.pt
    #Print output
    def __str__(self):
        return "PdgId : {0}  pt : {1}  eta : {2}  phi : {3}".format(self.pdgId,self.pt,self.eta,self.phi)

class genParticle(baseParticle):
    def __init__(self, prunedParticle):
        baseParticle.__init__(self, prunedParticle)
        self.numberOfDaughters = prunedParticle.numberOfDaughters()
        self.decayProducts = self.getDaughters(prunedParticle)
        self.leading  = self.getDaughters(prunedParticle)[0]
        self.trailing = self.getDaughters(prunedParticle)[1]
        if len(self.decayProducts)>1:
            mother = self.leading.p4 + self.trailing.p4
            self.motherMass = mother.M()
        self.decayMode  = abs(self.leading.pdgId)
        self.tausFinalState = self.getTau1Flavour() * self.getTau2Flavour()

    def getDaughters(self, p):
        fermions = []
        for i in range(0, p.numberOfDaughters()):
            daughter = p.daughter(i)
            fermion = baseParticle(daughter)
            fermions.append(fermion)
        return sorted(fermions, reverse=True) # Pt Descending ordering

    def daughtersDeltaR(self):
        deltaR = sqrt( pow((self.leading.eta - self.trailing.eta),2) + pow((self.leading.phi - self.trailing.phi),2)  )
        return deltaR

    def getTau1Flavour(self):
        flavour = 0
        if not abs(self.leading.pdgId) == 15:
            return flavour
        if self.__isMuonDecay(self.leading.prunedParticle):     flavour = 1;
        if self.__isElectronDecay(self.leading.prunedParticle): flavour = 2;
        if self.__isHadronicDecay(self.leading.prunedParticle):   flavour = 3;
        return flavour

    def getTau2Flavour(self):
        flavour = 0
        if not abs(self.trailing.pdgId) == 15:
            return flavour
        if self.__isMuonDecay(self.trailing.prunedParticle):     flavour = 1;
        if self.__isElectronDecay(self.trailing.prunedParticle): flavour = 2;
        if self.__isHadronicDecay(self.trailing.prunedParticle):   flavour = 3;
        return flavour


    def __isMuonDecay(self, tau):
        isMuon = False
        for i in range(0, tau.numberOfDaughters()):
            daughter = tau.daughter(i)
            if abs(daughter.pdgId()) == 13:
                isMuon = True
        return isMuon

    def __isElectronDecay(self, tau):
        isElectron = False
        for i in range(0, tau.numberOfDaughters()):
            daughter = tau.daughter(i)
            if abs(daughter.pdgId()) == 11:
                isElectron = True
        return isElectron

    def __isHadronicDecay(self, tau):
        return (not self.__isMuonDecay(tau) and not self.__isElectronDecay(tau))

# class tauDecayMode(Enum):
#     kMuDecay = 1
#     kEleDecay = 2
#     kHadDecay = 3

h = ROOT.TH1F( 'h1', 'test', 100, -10., 10. )

f = ROOT.TFile( 'test.root', 'recreate' )
t = ROOT.TTree( 't1', 'tree with histos' )
#finalState = -1
treeEventId = array( 'i', [ 0 ] )
treeLumiId = array( 'i', [ 0 ] )
treeRunId = array( 'i', [ 0 ] )
ZDecay = array( 'i', [ 0 ] )
HDecayMode = array( 'i', [ 0 ] )
pt_Z = array( 'f', [ 0 ] )
eta_Z = array( 'f', [ 0 ] )
mass_Z = array( 'f', [ 0 ] )
pt_H = array( 'f', [ 0 ] )
eta_H = array( 'f', [ 0 ] )
mass_H = array( 'f', [ 0 ] )
genmass_H = array( 'f', [ 0 ] )
deltaR_Z = array( 'f', [ 0 ] )
deltaR_H = array( 'f', [ 0 ] )
pt_lead_Z = array( 'f', [ 0 ] )
pt_trail_Z = array( 'f', [ 0 ] )
eta_lead_Z = array( 'f', [ 0 ] )
eta_trail_Z = array( 'f', [ 0 ] )
pt_lead_H = array( 'f', [ 0 ] )
pt_trail_H = array( 'f', [ 0 ] )
eta_lead_H = array( 'f', [ 0 ] )
eta_trail_H = array( 'f', [ 0 ] )

t.Branch('eventId', treeEventId , "eventId/I")
t.Branch('lumiId' , treeLumiId  , "lumiId/I" )
t.Branch('runId'  , treeRunId   , "runId/I"  )

t.Branch('pt_Z'  , pt_Z   , "pt_Z/F"  )
t.Branch('eta_Z'  , eta_Z   , "eta_Z/F"  )
t.Branch('mass_Z'  , mass_Z   , "mass_Z/F"  )
t.Branch('pt_H'  , pt_H   , "pt_H/F"  )
t.Branch('eta_H'  , eta_H   , "eta_H/F"  )
t.Branch('mass_H'  , mass_H   , "mass_H/F"  )
t.Branch('GENmass_H'  , genmass_H   , "GENmass_H/F"  )

t.Branch('deltaR_Z'  , deltaR_Z   , "deltaR_Z/F"  )
t.Branch('deltaR_H'  , deltaR_H   , "deltaR_H/F"  )
t.Branch('pt_lead_Z'  , pt_lead_Z   , "pt_lead_Z/F"  )
t.Branch('eta_lead_Z'  , eta_lead_Z   , "eta_lead_Z/F"  )
t.Branch('pt_trail_Z'  , pt_trail_Z   , "pt_trail_Z/F"  )
t.Branch('eta_trail_Z'  , eta_trail_Z   , "eta_trail_Z/F"  )
t.Branch('pt_lead_H'  , pt_lead_H   , "pt_lead_H/F"  )
t.Branch('eta_lead_H'  , eta_lead_H   , "eta_lead_H/F"  )
t.Branch('pt_trail_H'  , pt_trail_H   , "pt_trail_H/F"  )
t.Branch('eta_trail_H'  , eta_trail_H   , "eta_trail_H/F"  )

t.Branch( 'ZfinalState', ZDecay, 'ZfinalState/I' )
t.Branch( 'HfinalState', HDecayMode, 'HfinalState/I' )

events = Events (["/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/20245004-44C7-E611-A5AF-A0369F3016EC.root",
                  "/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/220CB8B6-53C7-E611-BC88-0CC47AD99052.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2411BE9B-74C8-E611-84FA-02163E019D73.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2A9B92A1-89C6-E611-A50D-002590E2D9FE.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4C20BEC8-D9C7-E611-B74B-FA163E2D1FE1.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/505A0021-F7C7-E611-9687-D067E5F91E51.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/600D3FD8-48C7-E611-B341-90B11C0BD210.root",
#"/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/88867835-A7C7-E611-B26A-008CFA11113C.root",
])

handlePruned, labelPruned  = Handle ("std::vector<reco::GenParticle>"), "prunedGenParticles"
handlePacked, labelPacked  = Handle ("std::vector<pat::PackedGenParticle> "), "packedGenParticles"
triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"
# loop over events
count= 0
for iev,event in enumerate(events):
    event.getByLabel (labelPacked, handlePacked)
    event.getByLabel (labelPruned, handlePruned)
    event.getByLabel(triggerBitLabel, triggerBits)
    event.getByLabel(triggerPrescaleLabel, triggerPrescales)
    # get the product
    packed = handlePacked.product()
    pruned = handlePruned.product()

    # print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    # print "\n === TRIGGER PATHS ==="
    # names = event.object().triggerNames(triggerBits.product())
    # for i in xrange(triggerBits.product().size()):
    #     print "Trigger ", names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)")

    treeEventId[0] = event.eventAuxiliary().event();
    treeLumiId[0]  = event.eventAuxiliary().luminosityBlock();
    treeRunId[0]   = event.eventAuxiliary().run();

    for p in pruned:
        # Z at GEN level
        if abs(p.pdgId()) == 23 :
            ZDecay[0] = -1
            if p.status() == 62:
                Zcandidate = genParticle(p)
                ZDecay[0] = Zcandidate.decayMode
                pt_Z [0] = Zcandidate.pt
                eta_Z [0] = Zcandidate.phi
                mass_Z[0] = Zcandidate.mass
                deltaR_Z[0] = Zcandidate.daughtersDeltaR()
                pt_lead_Z[0]   = Zcandidate.leading.pt
                pt_trail_Z[0]  = Zcandidate.trailing.pt
                eta_lead_Z[0]  = Zcandidate.leading.eta
                eta_trail_Z[0] = Zcandidate.trailing.eta
                # print Zcandidate #"  PdgId : %s   pt : %s  eta : %s   phi : %s   numberOfDaughters: %s" %(Zcandidate.pdgId,Zcandidate.pt,Zcandidate.eta,Zcandidate.phi,Zcandidate.numberOfDaughters)
                # print "     daughters"
                #for fermion in Zcandidate.decayProducts:
                    #print fermion
        # Higgs (25)
        #finalState = 1
        # isMuDecay = False
        # isEleDecay = False
        # isHadDecay = False
        if abs(p.pdgId()) == 25 :
            if p.status() == 62:
                HiggsCandidate = genParticle(p)
                #print HiggsCandidate.motherMass
                # print "="*30
                # print HiggsCandidate
                # #print "PdgId : %s   pt : %s  eta : %s   phi : %s   numberOfDaughters: %s" %(p.pdgId(),p.pt(),p.eta(),p.phi(),p.numberOfDaughters())
                # print " ---- daughters"
                # for tau in HiggsCandidate.decayProducts:
                #     print"---------"
                #     print tau
                #     print"---------"
                #     for i in range(0, tau.prunedParticle.numberOfDaughters()):
                #         dau = baseParticle( tau.prunedParticle.daughter(i))
                #         print "\t {0}".format(dau)
                # print HiggsCandidate.tausFinalState
                HDecayMode[0] = HiggsCandidate.tausFinalState
                pt_H [0] = HiggsCandidate.pt
                eta_H [0] = HiggsCandidate.phi
                mass_H[0] = HiggsCandidate.motherMass
                genmass_H[0] = HiggsCandidate.mass
                deltaR_H[0] = HiggsCandidate.daughtersDeltaR()
                pt_lead_H[0]   = HiggsCandidate.leading.pt
                pt_trail_H[0]  = HiggsCandidate.trailing.pt
                eta_lead_H[0]  = HiggsCandidate.leading.eta
                eta_trail_H[0] = HiggsCandidate.trailing.eta
                # print "="*30

    t.Fill()
f.Write()
f.Close()
            #    for pa in packed:
            #            mother = pa.mother(0)
            #            if mother and isAncestor(p,mother) :
            #                  print "     PdgId : %s   pt : %s  eta : %s   phi : %s" %(pa.pdgId(),pa.pt(),pa.eta(),pa.phi())
