from CRABClient.UserUtilities import config
config = config()

isMC = True
#signal or MC background?
isSig = True

#BParking set number, and Run letter of 2022 for the inputDataset
# 40 different datasets total
#0 to 7
setnum = 0
#C to G
runlet = 'G'

#v1 or v2?
version = 1
if runlet == 'D':
    version = 2

if not isMC:
    config.General.requestName = 'test6_%d%s'%(setnum, runlet)
else:
    if isSig:
        #config.General.requestName = 'test0MC_EtaTo2Mu2E'
        config.General.requestName = 'test11MC_EtaTo2Mu2E'
    else:
        config.General.requestName = 'test0MC_EtaToMuMuGamma'
config.General.workArea = 'crab_MiniAnalyzer'
config.General.transferOutputs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ntuplizer_cfg.py'
#parameters to pass into the cfg file
config.JobType.pyCfgParams=["data=%d"%(0 if isMC else 1)]
config.JobType.maxMemoryMB = 5000
##allow to run for up to 2 full days (maximum allowed, I think)
# Not allowed with automatic splitting.
#config.JobType.maxJobRuntimeMin = 2880
#config.JobType.maxJobRuntimeMin = 1440
config.JobType.numCores = 1 

if not isMC:
    config.Data.inputDataset = '/ParkingDoubleMuonLowMass%d/Run2022%s-PromptReco-v%d/MINIAOD'%(setnum, runlet, version)
#config.Data.inputDataset = '/EphemeralZeroBias0/Run2022C-v1/RAW'
#config.Data.inputDataset = '/DoubleElectron_Pt-1To300-gun/Run3Summer21DR-FlatPU0to70FEVT_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/BuToKee_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/Run3Summer21DRPremix-rndm_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.userInputFiles = [ 'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/GSDR_Run3_Signal_Samples/EtaTo2Mu2E_1020174_DIGIRAWHLT.root ',
else:
    if isSig:
        config.Data.userInputFiles = [
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_0_MINIAOD.root ',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_9_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_1_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_25_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_46_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_7_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_8_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_20_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_55_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_5_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_17_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_15_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_64_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_28_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_24_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_42_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_27_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_12_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_11_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_58_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_13_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_59_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_63_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_50_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_6_MINIAOD.root ',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_21_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_34_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_14_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_23_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_57_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_41_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_4_MINIAOD.root ',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_36_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_35_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_16_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_37_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_26_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_49_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_22_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_40_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_32_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_51_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_52_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_19_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_61_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_29_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_44_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_56_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_33_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_38_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_60_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_48_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_54_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_30_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_45_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_39_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_47_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_2_MINIAOD.root ',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_31_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_62_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_43_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_53_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_18_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_10_MINIAOD.root',
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_3_MINIAOD.root '
        ]
    else:
        config.Data.userInputFiles = [ 
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700) 
        ]

config.Data.inputDBS = 'global'
if isMC:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1 #10 #2
    config.Data.totalUnits = 100000 #10000
else:
    config.Data.splitting = 'Automatic'
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

