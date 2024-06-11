import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# setup 'analysis' options
options = VarParsing.VarParsing('analysis')

options.register('data',
        False,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.bool,
        "Run on data (1) or MC (0)"
        )
#use electron triggers instead of muon triggers (for calculating trigger eff in data)
options.register('elTrig',
        False,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.bool,
        "Run with electron triggers (1) or muon triggers (0)"
        )
options.register('test',
        0,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
        "Run for a test (1) or not (0)"
        )
options.parseArguments()

year = 2022
#data = False
data = options.data
elTrig = options.elTrig

#test file name
#BParking data
#testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/MiniTest/Run3_2022_BParking_MINItest.root"
#testfname = "root://xrootd-cms.infn.it//store/data/Run2022C/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/357/271/00000/012c38a8-143b-46a9-8846-3523c8d19862.root"
testfname = "root://xrootd-cms.infn.it//store/data/Run2022C/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/355/870/00000/f96cbcab-4a63-4806-a420-de4dbe586c89.root"
if not data:
    #signal sample
    #testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD/EtaTo2Mu2E_2022Test_0_MINIAOD_2022.root"
    #testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_2022_MINIAOD_3/EtaTo2Mu2E_2022Test_7_MINIAOD_2022.root"
    #testfname = "root://cmseos.fnal.gov//store/user/lpcdisptau/eta2mu2e/EtaToMuMuGamma/Run3_2022_MINIAOD_Exp3/EtaToMuMuGamma_2022Test_MINIAOD_7.root"
    testfname = "root://cmseos.fnal.gov//store/user/lpcdisptau/eta2mu2e/EtaToMuMuGamma/Run3_2022_MINIAOD_Exp3/EtaToMuMuGamma_2022Test_MINIAOD_9.root"
    ##resonant bkg sample
    #testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_1.root"
    #testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_Exp/EtaToMuMuGamma_2022Test_MINIAODTEST_303.root"
    ##EtaToMuMu sample
    #testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaToMuMu/Run3_2022_MINIAOD/EtaToMuMu_2022Test_MINIAOD_0.root"
    ##central mc sample
    #testfname = "/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/0009bf67-3c7a-4e81-b50a-e3914b3d2ffa.root"

outputFile = 'test.root'

if options.test: 
    options.inputFiles = testfname
    options.outputFile = outputFile

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.EventContent.EventContent_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if data:
    globaltag = '124X_dataRun3_Prompt_v4'
else:
    globaltag = '120X_mcRun3_2021_realistic_v6'

process.GlobalTag.globaltag = globaltag

process.MessageLogger = cms.Service("MessageLogger",
    destinations   =  cms.untracked.vstring('messages', 'cerr'),
    statistics     =  cms.untracked.vstring('messages', 'cerr'),
    debugModules   = cms.untracked.vstring('*'),
    categories     = cms.untracked.vstring('FwkReport'),
    messages       = cms.untracked.PSet(
        extension = cms.untracked.string('.txt'),
        threshold =  cms.untracked.string('WARNING')
        ),
    cerr           = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING'),
        WARNING = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(10)
            ),
        INFO = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(10)
            ),
        FwkReport = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(10000)
            )
        )
    )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(1)
    )
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1000)
    #input = cms.untracked.int32(100)
    #input = cms.untracked.int32(-1)
    input = cms.untracked.int32(options.maxEvents)
    )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles), #testfname),
    skipBadFiles = cms.untracked.bool(True)
    )
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile),
    closeFileFast = cms.untracked.bool(True),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
    )
process.Timing = cms.Service("Timing",
    summaryOnly = cms.untracked.bool(True),
    useJobReport = cms.untracked.bool(True)
    )

## Electron and photon VID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

from eta2mu2e.eta2mu2eSkimmer.eta2mu2eAnalyzer_cfi import eta2mu2eAnalyzer
process.ntuples = eta2mu2eAnalyzer.clone(
    isData = cms.bool(data),
    useElTrig = cms.bool(elTrig),
    primary_vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
)

process.commonSequence = cms.Sequence(
    process.ntuples
    )

process.p = cms.Path(
    process.commonSequence
)
