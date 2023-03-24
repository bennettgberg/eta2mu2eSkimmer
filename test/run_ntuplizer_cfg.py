import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# setup 'analysis' options
options = VarParsing.VarParsing('analysis')

options.parseArguments()

year = 2022
data = True #False

#test file name
#BParking data
testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/MiniTest/Run3_2022_BParking_MINItest.root"
if not data:
    #signal sample
    testfname = "root://cmseos.fnal.gov//store/user/bgreenbe/EtaTo2Mu2E/Run3_MiniAOD/EtaTo2Mu2E_10218787_MINIAOD.root"

process = cms.Process("USER")

outputFile = 'test.root'
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
    fileNames = cms.untracked.vstring(testfname),
    skipBadFiles = cms.untracked.bool(True)
    )
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFile),
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
    primary_vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
)

process.commonSequence = cms.Sequence(
    process.ntuples
    )

process.p = cms.Path(
    process.commonSequence
)
