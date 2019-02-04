import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('') )
from RecoJets.JetProducers.PileupJetIDParams_cfi import cutbased as pu_jetid


###### Electron VID
from RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff import *

if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_loose,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_loose.isPOGApproved
if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_medium,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_medium.isPOGApproved
if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_tight,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_tight.isPOGApproved

myVidElectronId = cms.PSet(
    loose = cutBasedElectronID_Spring15_25ns_V1_standalone_loose,
    medium = cutBasedElectronID_Spring15_25ns_V1_standalone_medium,
    tight = cutBasedElectronID_Spring15_25ns_V1_standalone_tight
)
#######

#from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jFullNoQG as mySignalMVA
from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jFull as mySignalMVA
#from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jBase as mySignalMVA

datapileup_latest = cms.vdouble(0, 238797, 837543, 2.30843e+06, 3.12475e+06, 4.47619e+06, 5.99591e+06, 7.0009e+06, 1.28917e+07, 3.52617e+07, 7.87012e+07, 1.76946e+08, 3.6009e+08, 6.02766e+08, 8.76519e+08, 1.17447e+09, 1.48906e+09, 1.75935e+09, 1.94393e+09, 2.04917e+09, 2.10158e+09, 2.13279e+09, 2.1491e+09, 2.12899e+09, 2.06265e+09, 1.96288e+09, 1.84187e+09, 1.70414e+09, 1.55452e+09, 1.39949e+09, 1.24353e+09, 1.08882e+09, 9.37305e+08, 7.92044e+08, 6.56718e+08, 5.34467e+08, 4.27127e+08, 3.35106e+08, 2.57725e+08, 1.93751e+08, 1.41831e+08, 1.00671e+08, 6.90139e+07, 4.55401e+07, 2.88475e+07, 1.75063e+07, 1.01626e+07, 5.63778e+06, 2.98728e+06, 1.512e+06, 731845, 339822, 152545, 67404.8, 30489.7, 15152.1, 8975.91, 6496.15, 5434.81, 4889.96, 4521.72, 4208.46, 3909.76, 3614.27, 3320.72, 3031.1, 2748.24, 2474.98, 2213.82, 1966.82, 1735.55, 1521.11, 1324.15, 1144.9, 983.22, 838.668, 710.534, 597.91, 499.739, 414.866, 342.082, 280.162, 227.901, 184.137, 147.773, 117.789, 93.2547, 73.3322, 57.2764, 44.4338, 34.2379, 26.2034, 19.9188, 15.0392, 11.2782, 8.40063, 6.21495, 4.56686, 3.33312, 2.41623)

from os import path as path

theLumiMask = path.expandvars("")


process.summuryNtupler = cms.EDAnalyzer('BaseNTuplizer',
    isMC = cms.bool(True),
    genEventInfoProduct   = cms.InputTag("generator"),
    PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
    lheEventProducts      = cms.InputTag("externalLHEProducer"),
)

process.syncNtupler = cms.EDAnalyzer('ZHTauTauAnalyzer',
#runProcess = cms.PSet(
    dtag  = cms.string("ZHToTauTau_M125_13TeV_powheg_pythia8"),
    outfile = cms.string("ZHToTauTau_M125_13TeV_powheg_pythia8.root"),
    isMC = cms.bool(True),
    triggerstudy = cms.bool(False),
    xsec = cms.double(0.05785),
    suffix = cms.string(""),
    cprime = cms.double(-1),
    brnew = cms.double(-1),
    mctruthmode = cms.int32(0),
    jacknife = cms.vint32(0,0),
    saveSummaryTree = cms.bool(True),
    runSystematics = cms.bool(False),
    runSVfit = cms.bool(False),
    resonance = cms.double(1),
    weightsFile = cms.vstring(""),
    puWeightsFile = cms.vstring(""),
    dirName = cms.string("dataAnalyzer"),
    useMVA = cms.bool(True),
    tmvaInput = mySignalMVA,
    muscleDir =  cms.string('${CMSSW_BASE}/src/UserCode/llvv_fwk/data/jec/'),
    jecDir = cms.string('${CMSSW_BASE}/src/UserCode/llvv_fwk/data/jec/25ns/'),
    datapileup = datapileup_latest,
    datapileupSingleLep = datapileup_latest,
    debug = cms.bool(False),
    lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange(),
    pujetidparas = cms.PSet(pu_jetid),
    electronidparas = cms.PSet(myVidElectronId),
    maxevents = cms.int32(-1), # set to -1 when running on grid.
    vtxSrc                = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho                   = cms.InputTag("fixedGridRhoFastjetAll"),
    muonSrc               = cms.InputTag("slimmedMuons"),
    electronSrc           = cms.InputTag("slimmedElectrons"),
    recHitCollectionEBSrc = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEESrc = cms.InputTag("reducedEgamma","reducedEERecHits"),
    tauSrc                = cms.InputTag("slimmedTaus"),
    packCandSrc           = cms.InputTag("packedPFCandidates"),
    jetSrc                = cms.InputTag("slimmedJets"),
    pfMETSrc              = cms.InputTag("slimmedMETs"),
    triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
    metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
    genParticleSrc        = cms.InputTag("prunedGenParticles"),
    genEventInfoProduct   = cms.InputTag("generator"),
    PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
    lheEventProducts      = cms.InputTag("externalLHEProducer"),
)

#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

runOnData = not process.syncNtupler.isMC

# Latest JEC
if runOnData:
  process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
  #process.source.lumisToProcess = LumiList.LumiList(filename = '../json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()
else:
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

inputfile = cms.untracked.vstring("file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/20245004-44C7-E611-A5AF-A0369F3016EC.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/220CB8B6-53C7-E611-BC88-0CC47AD99052.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2411BE9B-74C8-E611-84FA-02163E019D73.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2A9B92A1-89C6-E611-A50D-002590E2D9FE.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4C20BEC8-D9C7-E611-B74B-FA163E2D1FE1.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/505A0021-F7C7-E611-9687-D067E5F91E51.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/600D3FD8-48C7-E611-B341-90B11C0BD210.root",
#"file:/storage/data/cms/store/mc/RunIISummer16MiniAODv2/ZHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/88867835-A7C7-E611-B26A-008CFA11113C.root",
)

process.source = cms.Source ("PoolSource", fileNames = inputfile)

process.TFileService = cms.Service("TFileService", fileName = process.syncNtupler.outfile )

process.p = cms.Path(process.summuryNtupler +
                     process.syncNtupler
                     )
