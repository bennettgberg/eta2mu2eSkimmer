from CRABClient.UserUtilities import config
config = config()

isMC = True
#signal or MC background?
isSig = False

#EtaToMuMuGamma sample
isMuMu = True

#doing test with central MC production?
central = False

#BParking set number, and Run letter of 2022 for the inputDataset
# 40 different datasets total
#0 to 7
setnum = 0
#C to G
runlet = 'C'

#v1 or v2?
version = 1
if runlet == 'D':
    version = 2

if not isMC:
    config.General.requestName = 'test8_%d%s'%(setnum, runlet)
    #config.General.requestName = 'test1_disptau'
else:
    if isSig:
        #config.General.requestName = 'test0MC_EtaTo2Mu2E'
        config.General.requestName = 'test37MC_EtaTo2Mu2E2'
    elif isMuMu:
        config.General.requestName = 'test37MC_EtaToMuMu'
    else:
        config.General.requestName = 'test37MC_EtaToMuMuGamma3'
if central:
    config.General.requestName = 'test22MC_CentralJPsi'
config.General.workArea = 'crab_MiniAnalyzer'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ntuplizer_cfg.py'
#parameters to pass into the cfg file
config.JobType.pyCfgParams=["data=%d"%(0 if isMC else 1)]
#config.JobType.maxMemoryMB = 5000
##allow to run for up to 2 full days (maximum allowed, I think)
# Not allowed with automatic splitting.
#config.JobType.maxJobRuntimeMin = 2880
#config.JobType.maxJobRuntimeMin = 1440
config.JobType.numCores = 1 

if not isMC:
    config.Data.inputDataset = '/ParkingDoubleMuonLowMass%d/Run2022%s-PromptReco-v%d/MINIAOD'%(setnum, runlet, version)
    #config.Data.userInputFiles = [ 'root://cmsxrootd.fnal.gov//store/user/bgreenbe/MiniTest/Run3_2022_BParking_MINItest.root' ]
#config.Data.inputDataset = '/EphemeralZeroBias0/Run2022C-v1/RAW'
#config.Data.inputDataset = '/DoubleElectron_Pt-1To300-gun/Run3Summer21DR-FlatPU0to70FEVT_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/Run3Summer21DRPremix-120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/BuToKee_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/Run3Summer21DRPremix-rndm_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-DIGI-RAW'
#config.Data.userInputFiles = [ 'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/GSDR_Run3_Signal_Samples/EtaTo2Mu2E_1020174_DIGIRAWHLT.root ',
else:
    if isSig:
        inputFiles = []
        for i in range(65):
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD/EtaTo2Mu2E_2022Test_%d_MINIAOD_2022.root'%i) 
        for i in range(300):
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD_2/EtaTo2Mu2E_2022Test_%d_MINIAOD_2022.root'%i)
        config.Data.userInputFiles = inputFiles #[
        #    #'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD/EtaTo2Mu2E_2022Test_%d_MINIAOD_2022.root'%i for i in range(65)
        #    'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD_2/EtaTo2Mu2E_2022Test_%d_MINIAOD_2022.root'%i for i in range(300)
        #]
    elif isMuMu:
        config.Data.userInputFiles = [
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMu/Run3_2022_MINIAOD/EtaToMuMu_2022Test_MINIAOD_%d.root'%i for i in range(50)
        ]
    elif not central:
        inputFiles = []
        for i in range(700):
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i)
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_2/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i)
        for i in range(5020):
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i)  

        #list of file numbers for the 4th set of EtaToMuMuGamma MC files
        nums4 = []
        for i in range(126):
            if i == 111: continue
            nums4.append(i)
        for i in range(909, 2000):
            if i in [1804, 1805, 1807, 1808, 1809, 1812, 1817, 1818, 1821, 1824, 1826, 1838, 1839]: continue
            nums4.append(i)
        for i in nums4:
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/lpcdisptau/eta2mu2e/EtaToMuMuGamma/Run3_2022_MINIAOD_4/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i) 
        #Dan
        for i in range(305):
            if i == 12 or i == 15: continue
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/lpcdisptau/eta2mu2e/EtaToMuMuGamma/Run3_2022_MINIAOD_Dan/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i) 
        #Dan_2
        for i in range(269)
            inputFiles.append('root://cmsxrootd.fnal.gov//store/user/lpcdisptau/eta2mu2e/EtaToMuMuGamma/Run3_2022_MINIAOD_Dan_2/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i) 
        config.Data.userInputFiles = inputFiles # [ 
#            #'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700) 
#            #'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_2/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700) 
#            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(3500) 
#        ]
    else:
        #is central
        config.Data.userInputFiles = [
            'root://cmsxrootd.fnal.gov//store/mc/Run3Summer22EEMiniAODv3/LambdaBToX3872Lambda_X3872ToJPsiRho_JPsiTo2Mu_RhoTo2Pi_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/004dd08c-4890-4cc9-8acf-6ef929bf443a.root']

config.Data.inputDBS = 'global'
if isMC:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1 #10 #2
    config.Data.totalUnits = 100000 #10000
else:
    #config.Data.splitting = 'Automatic'
    ### for test only
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1 #2
    config.Data.totalUnits = 1000000 #10000
    ###  ###
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.ignoreLocality = True

if isMC:
    config.Data.outLFNDirBase = '/store/user/bgreenbe/BParking2022'
else:
    config.Data.outLFNDirBase = '/store/group/lpcdisptau/eta2mu2e'
config.Site.whitelist = ['T2_US*','T2_CH*']
#config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T3_US_FNALLPC'

##need to use this to allow the use of this weird CMSSW release
#config.JobType.allowUndistributedCMSSW = True

