import FWCore.ParameterSet.Config as cms

process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

isMC=False

gt='80X_dataRun2_Prompt_v16'          ### '80X_dataRun2_2016SeptRepro_v7' 
outName='Data'
if isMC :
    gt='80X_mcRun2_asymptotic_2016_TrancheIV_v6'
    outName='MC'

process.GlobalTag.globaltag = gt
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PFchs    = cms.EDAnalyzer('JetCorrectorDBReader',
                                         payloadName    = cms.untracked.string('AK4PFchs'),
                                         globalTag      = cms.untracked.string(outName),
                                         printScreen    = cms.untracked.bool(False),
                                         createTextFile = cms.untracked.bool(True)
                                         )

process.p = cms.Path(process.readAK5PFchs)
