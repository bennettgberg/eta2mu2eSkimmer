from CRABClient.UserUtilities import config
config = config()

#BParking set number, and Run letter of 2022 for the inputDataset
# 40 different datasets total
#0 to 7
setnum = 1
#C to G
runlet = 'C'

config.General.requestName = 'test3_%d%s'%(setnum, runlet)
config.General.workArea = 'crab_MiniAnalyzer'
config.General.transferOutputs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_TM_ntuplizer_cfg.py'
config.JobType.maxMemoryMB = 5000
##allow to run for up to 2 full days (maximum allowed, I think)
#
#config.JobType.maxJobRuntimeMin = 2880
config.JobType.numCores = 1 

#config.Data.inputDataset = '/ParkingDoubleMuonLowMass%d/Run2022%s-PromptReco-v1/MINIAOD'%(setnum, runlet)
config.Data.userInputFiles = 'root://cmsxrootd.fnal.gov//store/user/bgreenbe/MiniTest/Run3_2022_BParking_MINItest.root'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1 #10 #2
#config.Data.totalUnits = 100000 #10000
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.ignoreLocality = True

config.Data.outLFNDirBase = '/store/user/marlow/TrueMuonium/BParking2022'
config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T3_US_FNALLPC'


