from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'test_summuryTTree_analysis'
config.General.workArea = 'ZHToTauTau'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/users/c/c/ccaputo/AZH/Full-FW/CMSSW_8_0_26_patch1/src/UserCode/llvv_fwk/crab/ZHToTauTau_M125_13TeV_powheg_pythia8_0_cfg.py'

config.Data.inputDataset = '/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 300
config.Data.outLFNDirBase = '/store/user/%s/provaCRAB_TTree/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_test_MC_analysis'

config.Site.whitelist   = ['T2_BE_UCL']
config.Site.storageSite = 'T2_BE_UCL'
