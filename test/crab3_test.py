from CRABClient.UserUtilities import config
config = config()

#BParking set number, and Run letter of 2022 for the inputDataset
# 40 different datasets total
#0 to 7
setnum = 0
#C to G
runlet = 'C'

config.General.requestName = 'test1_%d%s'%(setnum, runlet)
config.General.workArea = 'crab_MiniAnalyzer'
config.General.transferOutputs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ntuplizer_cfg.py'
config.JobType.maxMemoryMB = 5000
##allow to run for up to 2 full days (maximum allowed, I think)
#config.JobType.maxJobRuntimeMin = 2880
#config.JobType.maxJobRuntimeMin = 1440
config.JobType.numCores = 1 

config.Data.inputDataset = '/ParkingDoubleMuonLowMass%d/Run2022%s-PromptReco-v1/MINIAOD'%(setnum, runlet)
#config.Data.inputDataset = '/EphemeralZeroBias0/Run2022C-v1/RAW'
#config.Data.inputDataset = '/DoubleElectron_Pt-1To300-gun/Run3Summer21DR-FlatPU0to70FEVT_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/BuToKee_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/Run3Summer21DRPremix-rndm_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.userInputFiles = [ 'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/GSDR_Run3_Signal_Samples/EtaTo2Mu2E_1020174_DIGIRAWHLT.root ',

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 1 #10 #2
#config.Data.totalUnits = 100000 #10000
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.ignoreLocality = True

config.Data.outLFNDirBase = '/store/user/bgreenbe/BParking2022'
config.Site.whitelist = ['T2_US*','T2_CH*']
#config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T3_US_FNALLPC'

##need to use this to allow the use of this weird CMSSW release
#config.JobType.allowUndistributedCMSSW = True

