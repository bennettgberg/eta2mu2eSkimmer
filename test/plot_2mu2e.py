import sys
import printEvent
import printEvent_backup

#allow events that pass the trigger only?
trg_only = True
if not trg_only:
    input("Are you *sure* you want trg_only False???")

#instead of looking at trigger bits, do toy MC simulation to find effect of trigger turn-on threshold uncty
# how much to vary the turn-on threshold by? 0: nominal (closure); over 100: off (regular MC)
trg_var = 999 #0.08 #999 #0.07
if abs(trg_var < 100):
    import random

#set this true to REJECT ANY event that have a reconstructed photon!
reject_photon = False
#set this true to reject any event that has a valid eta->mumu (.52 to .58 GeV invar. mass)
reject_etamumu = False
#set this to true to reject any event that has a valid eta->mumugamma (.52 to .58 GeV invar. mass)
reject_mmg = False
#require the conversion veto and nMissingHits <= 3 on electrons?
basic_cuts = False #True
#require electron_ID?
require_elID = False
#require that the electron DOESN'T pass elID? (for a study)
invert_elID = True
#require muon ID to be greater than 0?
require_muID = True
#cut on electron pT? -1 for no cut
elpTcut = 2
#cut on electron |pseudoRapidity|
elEtacut = 2.5
#also cut out the electrons in the transition from barrel to endcap: 1.44 < |eta| < 1.57 ?
cut_trans = True
#cut on muon pT
mupTcut = 3
muleadpTcut = 4
#cut on muon |pseudoRapidity|
muEtacut = 2.4
#cut out .04 < M_ee < .09, or nah?
cut_mee = False
#include pileup corrections?
do_pileup = True
#include trigger efficiency correction from the D0 meson analysis (flat plateau correction)?
do_trigCor = False #True
#do the trigger turn on correction (L1 from BPH-24-005, HLT derived ourselves)?
do_TTOC = True 
#use the new event weights (calculated from DG/Cheb4 2mu fits)?
new_wt = 4 #2
#-1: preEE, +1: postEE, 0: all together
EEperiod = -1
#use WP80 elID instead of WP90?
tightEl = False
#use medium muID instead of loose?
tightMu = False
#fractional shift (scale factor) to apply to all muon pT values, for MC ONLY 
musf = -0.0014
elsf = 0.001
#add in the (vertical) electron sf (from the json file)
do_ElIDsf = False
do_ElRecsf = False
#add in (vertical) muon sf?
do_MuIDsf = True

peakmin = .52
peakmax = .58

#apply Golden JSON?
golden = True

#what test number to label the output files with
testnum = 4787

isMC = False
#use the central MC just to test the triggers (not really useful anymore)
central = False
#is EtaToMuMu MC? (will be set by arguments)
isMuMu = False

nPrint = 5
#true if running a synchronization test (so print out each event, diff nbins, etc) -- now set in args
syncTest = False

#set to true if running sync over just one single file instead of a whole subset of data
singleFile = False

#year only used for data-- for MC the lumi is sum of 2022 only for now
year = 2022
#argument is telling which files to analyze (should be about 100 files each)
if len(sys.argv) < 2:
    #print("Error: please provide an argument as an integer betwixt 0 and 35") 
    print("No argument provided so running in signal MC mode") 
    isMC = True
    isSig = True
    arg = -1
elif sys.argv[1] == "sync" and singleFile:
    syncTest = True
    isMC = False
    isSig = False
    print("Running syncTest!")
    arg = -1
elif sys.argv[1] == "bkg":
    isMC = True
    isSig = False
    print("Running in bkg MC mode")
    if len(sys.argv) > 2:
        arg = int(sys.argv[2])
    else:
        arg = -1
elif sys.argv[1] == "central":
    isMC = True
    isSig = False
    central = True
    arg = -1
    print("Running test of central MC for trigger debugging purposes.")
elif sys.argv[1] == "mumu":
    isMC = True
    isSig = False
    print("Running in EtaToMuMu MC mode")
    arg = -1
    isMuMu = True
elif "2023" in sys.argv[1]:
    year = 2023
    if "sig" in sys.argv[1] or "mmee" in sys.argv[1]:
        isMC = True
        isSig = True
    elif "bkg" in sys.argv[1] or "mmg" in sys.argv[1]:
        isMC = True
        isSig = False
    elif "mumu" in sys.argv[1]:
        isMC = True
        isMuMu = True
        isSig = False
    else:
        isSig = False
        isMC = False
    arg = -1
elif len(sys.argv) < 4:
    print("Error: for data you must specify the run letter and set number.\nUsage: python plot_2mu2e.py [let] [setnum] [subnum]")
    print("let: C to G ; setnum: 0 to 7 ; subnum: 0 to 31") 
    print("Example: python plot_2mu2e.py C 0 12") 
    sys.exit()
else:
    let = sys.argv[1]
    num = int(sys.argv[2]) 
    arg = int(sys.argv[3])

import glob
import ROOT
import time
import os
import array

#if syncTest:
#    import printEvent
#    nPrint = 5

if syncTest and isMC:
    print("Error: syncTest only runnable on data rn.")
    sys.exit()

start = time.time()

#lets = ['C']

#nums = [0]
mu_mass = .105658
el_mass = .000511
etamass = .547862
pi_mass = .13957

#scale factors are for MC ONLY
if not isMC:
    musf = 0
    elsf = 0
    do_ElIDsf = False
    do_ElRecsf = False
    do_MuIDsf = False
    do_TTOC = False

if do_ElIDsf or do_ElRecsf or do_MuIDsf:
    #open the json file
    from correctionlib import _core
if do_ElIDsf or do_ElRecsf:
    evaluator = _core.CorrectionSet.from_file('electron.json.gz')
if do_MuIDsf:
    if EEperiod == -1:
        mvaluator = _core.CorrectionSet.from_file('muon_preEE.json.gz')
    else:
        mvaluator = _core.CorrectionSet.from_file('muon_postEE.json.gz')

#if true, add only the vertex with the highest chi2 prob to the histogram
singleVert = False #True #not syncTest
##maximum reduced chi2 on the vertex that is allowed to be kept (-1 for no cut, 2.70554 for chi2 prob>.1 for 2-lep vertices; 1.84727 for 4-lepton) 
#rChi2Cut = -1 #10.0 #2.6 # -- used only on test36 and before!!
#-1 for no cut
vProbCut = 0.1
#use low pt electrons too?
useLowPt = True #not syncTest

#if True, genmatch the Onia photon tracks to the LowPt electrons and discard any that are a match
useOnia = False
#minimum dR for it to be considered a successful genmatch
dRcut = 0.03 #0.13 #0.08

#parameters for invar mass histogram
#nbins = 350
nbins = 550
#xmin = .45
xmin = .25
xmax = .8
if syncTest:
    nbins = 100
    xmin = .25
    xmax = 1.0

#parse an input file that lists event numbers and invariant masses; store them in a list
def readEvents(eventsFile):
    ef = open(eventsFile, "r")
    events = []
    runs = []
    for line in ef:
        words = line.split()
        evt = int(words[1])
        run = int(words[0])
        events.append(evt)
        runs.append(run)
    return runs,events

#list of vertex types 
#2pat::Electron; 2pat::Muon-2pat::Electron; 2pat::Muon-Photon; 2pat::Muon
#vtypes = ["elel", "mmelel", "mmg", "mumu"]
vtypes = ["mmelel", "mumu", "elel"]
#vtypes = ["mmg", "mumu", "mmelel"] 
if useLowPt:
    #2-low-pT pat::Electron; mu-mu-2-low-pT pat::Electron
    #vtypes.append("lplp")
    vtypes.append("mmlplp")
if isMuMu and "mumu" not in vtypes:
    print("WARNING: adding mumu vertices to vtypes.")
    vtypes.append("mumu")
if isMC and not isSig and not isMuMu and len(vtypes) > 1:
    print("WARNING: changing vtypes to just one type for resonant background MC!")
    if useLowPt:
        vtypes = ["mmlplp"]
    else:
        vtypes = ["mmelel"]
    #vtypes = ["elel"]
    #vtypes = ["mumu"]

#if True, ONLY mmelel vertices made up of 2 good 2-lepton vertices are allowed
# ie if no good mu pair or no good el pair, cannot use the 4-lepton pair
mmelelExclusive = False # True

#list of particle types
#ptypes = ["PCTrack", "patElectron", "lowpTElectron"]
ptypes = ["patElectron"]
if useLowPt:
    ptypes.append("lowpTElectron")

#just keep track of the total number of events processed.
hNevt = ROOT.TH1D("hNevt", "Total number of events processed", 1, 0, 2) 

#number of primary vertices distribution
hnPV = ROOT.TH1D("hnPV", "Number of primary vertices", 200, 0, 200) 
#all events passing the L1 trigger HLT_Mu0_L1DoubleMu, as a function of subleading muon pT
hpTL1Pass = ROOT.TH1F("hpTL1Pass", "", 500, 0, 50)
#events ALSO passing the HLT trigger HLT_DoubleMu4_3_LowMass
hpTHLTPass = ROOT.TH1F("hpTHLTPass", "", 500, 0, 50)
#L1Pass and ALSO is a good mumu event
hpTGoodL1Pass = ROOT.TH1F("hpTGoodL1Pass", "", 500, 0, 50)
#HLTPass and ALSO is a good mumu event
hpTGoodHLTPass = ROOT.TH1F("hpTGoodHLTPass", "", 500, 0, 50)
#pT of all electrons
hpTElAll = ROOT.TH1F("hpTElAll", "", 100, 0.0, 50.0) 
hpTElAll0 = ROOT.TH1F("hpTElAll0", "electron #Delta R < 0.0005", 100, 0.0, 50.0) 
hpTElAll1 = ROOT.TH1F("hpTElAll1", "0.0005 < electron #Delta R < 0.0025", 100, 0.0, 50.0) 
hpTElAll2 = ROOT.TH1F("hpTElAll2", "0.0025 < electron #Delta R < 0.005", 100, 0.0, 50.0) 
hpTElAll3 = ROOT.TH1F("hpTElAll3", "electron #Delta R > 0.005", 100, 0.0, 50.0) 
#pT of all electrons that pass the electron ID (individually)
hpTElIDPass = ROOT.TH1F("hpTElIDPass", "", 100, 0.0, 50.0) 
hpTElIDPass0 = ROOT.TH1F("hpTElIDPass0", "electron #Delta R < 0.0005", 100, 0.0, 50.0) 
hpTElIDPass1 = ROOT.TH1F("hpTElIDPass1", "0.0005 < electron #Delta R < 0.0025", 100, 0.0, 50.0) 
hpTElIDPass2 = ROOT.TH1F("hpTElIDPass2", "0.0025 < electron #Delta R < 0.005", 100, 0.0, 50.0) 
hpTElIDPass3 = ROOT.TH1F("hpTElIDPass3", "electron #Delta R > 0.005", 100, 0.0, 50.0) 
#deltaR between reco'd electrons of mmelel vertices
hdREl = {}  
hdRvPtEl = {}  
#invar mass, pT dists (dictionary for ease of adding new vertex types)
hM = {}
#invar mass with INVERTED electron charge requirement (ie same-sign electrons req'd instead of oppo-sign) 
hMSSe = ROOT.TH1F("hMSSe", "#mu#muee invar. mass with same-sign electrons", nbins, xmin, xmax) 
hMOSe = ROOT.TH1F("hMOSe", "#mu#muee invar. mass with oppo-sign electrons, filled same way as hMSSe", nbins, xmin, xmax)
hMee = ROOT.TH1F("hMee", "Invar. mass of elel vertices", 200, 0.0, 1.0)
hMeePeak = ROOT.TH1F("hMeePeak", "Invar. mass of elel vertices with mmelel in the eta mass window", 200, 0.0, 1.0)
if isMC: # and not isMuMu:
    hMeedR = ROOT.TH1F("hMeedR", "min #DeltaR between reco. electrons in mmelel vertices and gen electrons or photon", 1000, 0.0, 1.0)
    hMmmdR = ROOT.TH1F("hMmmdR", "min #DeltaR between reco. muons in peak and gen muons", 1000, 0.0, 1.0)
    hMeeGM = ROOT.TH1F("hMeeGM", "Invar. mass of mmelel vertices with e^{-},e^{+} gen-matched to gen electrons or photon (#DeltaR<%f)"%dRcut, 200, 0.0, 1.0)
    hEldR = ROOT.TH1F("hEldR", "min #DeltaR b/t gen electron and reco electron of same charge or photon", 1000, 0.0, 1.0) 
    
#reduced chi2 of the vertices filling the histogram
hRchi2 = {}
#mass histogram with all weights 1
hMNoWt = {}
#mass hist of events that failed the trigger (for some reason some data events are doing this??) 
hMFailedTrig = {}
#PV histogram for good events (passing all selections!)-- just for closure test
hnPVgood = {}
#histogram of events combining event N's dimuon with event N-1's dielectron, to test for combinatorial bkg
hMComb = ROOT.TH1F("hMComb", "Invar. mass of event N's dimuon, event N-1's dielectron", 2000, 0.0, 10.0) 
#control of hMComb, filled exactly the same way except with all this event.
hMReal = ROOT.TH1F("hMReal", "Control plot for hMComb: event N's dimuon+dielectron", 2000, 0.0, 10.0) 
#subleading pt for events passing L1
hL1TC = ROOT.TH1F("hL1TC", "|#eta| < 1.2", 200, 0.0, 20.0)
hL1TF = ROOT.TH1F("hL1TF", "|#eta| > 1.2", 200, 0.0, 20.0)
#subleading pt for events passing HLT (and therefore ofc also L1)
hHLTC = ROOT.TH1F("hHLTC", "|#eta| < 1.2", 200, 0.0, 20.0)
hHLTF = ROOT.TH1F("hHLTF", "|#eta| > 1.2", 200, 0.0, 20.0)
#pT histogram 
hpT = {}
#pT of events IN eta mass window (.52-.58 GeV) ONLY
hpTSig = {}
#pT of events JUST OUTSIDE eta mass window (.45-.49,.62-.67 GeV) ONLY
hpTSide = {}
#pT of electrons in the vertex
hpTEl = {}
hpTElNoWt = {}
#pT of muons in the vertex
hpTMu = {}
hpTMuF = {}
hpTMuNoWt = {}
#pseudorapidity of muons in the vertex
hEtaMu = {}
#pseudorapidity histogram
hEta = {}
##dR histogram for positive muons (MC only)
#hdRP = {}
##"      "        negative   "
#hdRN = {}
#separate hists for high vs. low Vxy values (cutoff 1.2 cm)
#hMhiVxy = {}
#hMloVxy = {}
#hDxyHiVxy = {}
#hDxyLoVxy = {}
#hSigmaVxyHiVxy = {}
#hSigmaVxyLoVxy = {}
hVxy = {}
hNcand = {}
#2-d hist of invar. mass vs. pT
hMvsPt = {}
#2-d hist of electron scale factors as a function of pT, eta
hElRecsf = {}
hElRecsfUp = {}
hElRecsfDn = {}
hElIDsf = {}
hElIDsfUp = {}
hElIDsfDn = {}
hEletapt = {}
#sublead pT vs dR of muon pair
hsubVdR = {}
hL1all = {}
#L1 seeds passed
hL1pass = {}

#hsubpTL1 = ROOT.TH2F("hsubpTL1", "sublead p_{T} for events passing L1", 200, 0.0, 50.0) 
#hsubpTHLT = ROOT.TH2F("hsubpTL1", "sublead p_{T} for events passing L1", 200, 0.0, 50.0) 
#2-d hist of ee invar mass vs mmee
hMeeVsMmmee = ROOT.TH2F("hMeeVsMmmee", "", 200, 0.0, 1.0, 200, 0.0, 1.0) 
if isMC and not isSig:
    hMeeNoG = {}
    hMeeG = {}
if useOnia:
    hVertMatched = {}
    hVertNoMatch = {}
#plot of the mu-mu mass only of mmelel vertices (just to see what it looks like)
hMNoEl = ROOT.TH1F("hMNoEl", "Invar. mass of muons ONLY in mmelel signal window", 200, 0.0, 1.0) 
hMNoMu = ROOT.TH1F("hMNoMu", "Invar. mass of electrons ONLY in mmelel signal window", 200, 0.0, 1.0) 
hMNoElLSide = ROOT.TH1F("hMNoElLSide", "Invar. mass of muons ONLY to LEFT of mmelel signal window (lower sideband)", 200, 0.0, 1.0) 
hMNoMuLSide = ROOT.TH1F("hMNoMuLSide", "Invar. mass of electrons ONLY to LEFT of mmelel signal window (lower sideband)", 200, 0.0, 1.0) 
hMNoElSide = ROOT.TH1F("hMNoElSide", "Invar. mass of muons ONLY OUTSIDE of mmelel signal window (lower or upper sideband)", 200, 0.0, 1.0) 
hMNoMuSide = ROOT.TH1F("hMNoMuSide", "Invar. mass of electrons ONLY OUTSIDE of mmelel signal window (lower or upper sideband)", 200, 0.0, 1.0) 
hMNoElRSide = ROOT.TH1F("hMNoElRSide", "Invar. mass of muons ONLY to RIGHT of mmelel signal window (upper sideband)", 200, 0.0, 1.0) 
hMNoMuRSide = ROOT.TH1F("hMNoMuRSide", "Invar. mass of electrons ONLY to RIGHT of mmelel signal window (upper sideband)", 200, 0.0, 1.0) 
#invariant mass distribution of just electrons in the mu-mu-e-e, but assuming pion mass instead of electron mass
#hMNoMuPiM = ROOT.TH1F("hMNoMuPiM", "Invar. mass of electrons ONLY in mmelel signal window, assuming pion mass instead", 200, 0.0, 1.0) 
for vtype in vtypes:
    if isMC and vtype not in ["mumu", "mmg"]:
        #hElRecsf[vtype] = ROOT.TH2F("hElRecsf"+vtype, "Electron SF as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        #hElRecsfUp[vtype] = ROOT.TH2F("hElRecsfUp"+vtype, "Electron SF up variation as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        #hElRecsfDn[vtype] = ROOT.TH2F("hElRecsfDn"+vtype, "Electron SF down variation as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        #hElRecsf[vtype] = ROOT.TH2F("hElRecsf"+vtype, "Electron SF as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        #hElRecsfUp[vtype] = ROOT.TH2F("hElRecsfUp"+vtype, "Electron SF up variation as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        #hElRecsfDn[vtype] = ROOT.TH2F("hElRecsfDn"+vtype, "Electron SF down variation as a function of #eta, p_{T}", 48, -2.4, 2.4, 50, 0, 50.0) 
        hElRecsf[vtype] = ROOT.TH1F("hElRecsf"+vtype, "Electron SF as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hElRecsfUp[vtype] = ROOT.TH1F("hElRecsfUp"+vtype, "Electron SF up variation as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hElRecsfDn[vtype] = ROOT.TH1F("hElRecsfDn"+vtype, "Electron SF down variation as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hElIDsf[vtype] = ROOT.TH1F("hElIDsf"+vtype, "Electron SF as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hElIDsfUp[vtype] = ROOT.TH1F("hElIDsfUp"+vtype, "Electron SF up variation as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hElIDsfDn[vtype] = ROOT.TH1F("hElIDsfDn"+vtype, "Electron SF down variation as a function of #eta, p_{T}", 48, -2.4, 2.4) 
        hEletapt[vtype] = ROOT.TH2F("hEletapt"+vtype, "Electron #eta, pt distributions", 480, -2.4, 2.4, 50, 0, 50.0) 
    if vtype in ["mumu", "elel", "pcpc", "lplp"]:
        hM[vtype] = ROOT.TH1F("hM"+vtype, "Invar. mass with "+vtype+" vertices", 1000, 0.0, 1.0) 
    else:
        hM[vtype] = ROOT.TH1F("hM"+vtype, "Invar. mass with "+vtype+" vertices", nbins, xmin, xmax) 
    hM[vtype].Sumw2()
    hdREl[vtype] = ROOT.TH1F("hdREl"+vtype, "", 1000, 0.0, 0.1)
    hdRvPtEl[vtype] = ROOT.TH2F("hdRvPtEl"+vtype, "#DeltaR vs. p_{T} for electrons in "+vtype+" vertices", 1000, 0.0, 0.1, 1000, 0., 100.)
    hL1all[vtype] = ROOT.TH1D("hL1all"+vtype, "Events surviving all cuts for "+vtype+" vertices that have L1 seed info", 1, 0, 1)
    hL1pass[vtype] = ROOT.TH1D("hL1pass"+vtype, "L1 pass rates for each seed, amongst events surviving all cuts for "+vtype+" vertices", 8, 0, 8)
    hRchi2[vtype] = ROOT.TH1F("hRchi2"+vtype, "Reduced #chi^{2} with "+vtype+" vertices", 500, 0.0, 5.0)
    hRchi2[vtype].Sumw2()
    hMNoWt[vtype] = ROOT.TH1F("hMNoWt"+vtype, "UNWEIGHTED Invar. mass with "+vtype+" vertices", nbins, xmin, xmax) 
    hMNoWt[vtype].Sumw2()
    hMFailedTrig[vtype] = ROOT.TH1F("hMFailedTrig"+vtype, "Invar. mass of events FAILING trigger with "+vtype+" vertices", nbins, xmin, xmax) 
    hMFailedTrig[vtype].Sumw2()
    hnPVgood[vtype] = ROOT.TH1D("hnPVgood"+vtype, "Number of primary vertices in good "+vtype+" vertices", 200, 0, 200) 
    hpT[vtype] = ROOT.TH1F("hpT"+vtype, "pT with "+vtype+" vertices", 500, 0., 100.) 
    hpTSig[vtype] = ROOT.TH1F("hpTSig"+vtype, "p_{T}^{2#mu2e} w/ .52<m_{2#mu2e}<.58 GeV w/ "+vtype+" vertices", 500, 0., 100.) 
    hpTSide[vtype] = ROOT.TH1F("hpTSide"+vtype, "p_{T}^{2#mu2e} w/ 0.61<m_{2#mu2e}<0.8 GeV w/ "+vtype+" vertices", 500, 0., 100.) 
    hpT[vtype].Sumw2()
    hpTEl[vtype] = ROOT.TH1F("hpTEl"+vtype, "Weighted pT of electrons in "+vtype+" vertices", 500, 0., 100.) 
    hpTElNoWt[vtype] = ROOT.TH1F("hpTElNoWt"+vtype, "Unweighted pT of electrons in "+vtype+" vertices", 500, 0., 100.) 
    hpTMu[vtype] = ROOT.TH1F("hpTMu"+vtype, "Weighted p_{T} of muons in "+vtype+" vertices", 500, 0., 100.) 
    hpTMuF[vtype] = ROOT.TH1F("hpTMu"+vtype+"F", "WeightedpT of muons in "+vtype+" vertices with |#eta|>1.5", 500, 0., 100.) 
    hpTMuNoWt[vtype] = ROOT.TH1F("hpTMuNoWt"+vtype, "Unweighted pT of muons in "+vtype+" vertices", 500, 0., 100.) 
    hEtaMu[vtype] = ROOT.TH1F("hEtaMu"+vtype, "Pseudorapidity of muons in "+vtype+" vertices", 2000, -10., 10.) 
    #pseudorapidity distribution of the reconstructed eta mesons
    hEta[vtype] = ROOT.TH1F("hEta"+vtype, "pseudorapidity with "+vtype+" vertices", 2000, -10., 10.) 
    hEta[vtype].Sumw2()
    hNcand[vtype] = ROOT.TH1D("hNcand"+vtype, "Number of candidates in [.52, .58] for "+vtype+" vertices", 20, 0, 20) 
    #hdRP[vtype] = ROOT.TH1F("hdRP"+vtype, "dR b/t Onia and lowPtelectrons Pos", 1000, 0, 10.)
    #hdRN[vtype] = ROOT.TH1F("hdRN"+vtype, "dR b/t Onia and lowPtelectrons Neg", 1000, 0, 10.)
    #hMhiVxy[vtype] = ROOT.TH1F("hMhiVxy"+vtype, "Invar. mass with "+vtype+" vertices, Vxy>1.2", nbins, xmin, xmax) 
    #hMloVxy[vtype] = ROOT.TH1F("hMloVxy"+vtype, "Invar. mass with "+vtype+" vertices, Vxy<1.2", nbins, xmin, xmax) 
    #hDxyHiVxy[vtype] = ROOT.TH1F("hDxyHiVxy"+vtype, "Dxy with "+vtype+" vertices, Vxy>1.2", 10000, -.5, .5) 
    #hDxyLoVxy[vtype] = ROOT.TH1F("hDxyLoVxy"+vtype, "Dxy with "+vtype+" vertices, Vxy<1.2", 10000, -.5, .5) 
    #hSigmaVxyHiVxy[vtype] = ROOT.TH1F("hSigmaVxyHiVxy"+vtype, "#sigmaVxy with "+vtype+" vertices, Vxy>1.2", 10000, 0.0, 10.0) 
    #hSigmaVxyLoVxy[vtype] = ROOT.TH1F("hSigmaVxyLoVxy"+vtype, "#sigmaVxy with "+vtype+" vertices, Vxy<1.2", 10000, 0.0, 10.0) 
    #hMvsPt[vtype] = ROOT.TH2F("hMvsPt"+vtype, "Invar. mass as a function of pT", 1000, 0.0, 100.0, 350, .45, .8) 
    hMvsPt[vtype] = ROOT.TH2F("hMvsPt"+vtype, "Invar. mass as a function of pT", 100, 0.0, 100.0, 800, .2, 1.0) 
    hsubVdR[vtype] = ROOT.TH2F("hsubVdR"+vtype, "sublead p_{T} vs. #DeltaR for vertices of type "+vtype, 500, 0.0, 100.0, 1000, 0.0, 1.0) 
    #hMvsPt[vtype].Sumw2()
    hVxy[vtype] = ROOT.TH1F("hVxy"+vtype, "Vxy for all "+vtype+" vertices", 1000, 0.0, 10.0) 
    hVxy[vtype].Sumw2()
    #dielectron mass withOUT photon saved
    if isMC and not isSig:
        hMeeNoG[vtype] = ROOT.TH1F("hMeeNoG"+vtype, "Invar. mass with "+vtype+" vertices for evts with NO gen photon", 1000, 0, 1.0)
        hMeeNoG[vtype].Sumw2()
        hMeeG[vtype] = ROOT.TH1F("hMeeG"+vtype, "Invar. mass with "+vtype+" for evts with a gen photon", 1000, 0, 1.0)
        hMeeG[vtype].Sumw2()
    if useOnia:
        hVertMatched[vtype] = ROOT.TH2F("hVertMatched"+vtype, "vx vs. vy for matched "+vtype+" tracks", 12000, -60, 60, 12000, -60, 60)
        hVertMatched[vtype].Sumw2()
        hVertNoMatch[vtype] = ROOT.TH2F("hVertNoMatch"+vtype, "vx vs. vy for unmatched "+vtype+" tracks", 12000, -60, 60, 12000, -60, 60)
        hVertNoMatch[vtype].Sumw2()
        if isMC:
            hpTVertMatched[vtype] = ROOT.TH1F("hpTVertMatched"+vtype, "Gen pT of etas in which "+vtype+" elecs were matched to converted photons", 500, 0., 100.) 
            hpTVertMatched[vtype].Sumw2()

if mmelelExclusive:
    #histogram of invariant masses that have valid mmelel vertices but either invalid mu-mu or invalid el-el vertex
    hMdeleted = ROOT.TH1F("hMdeleted", "DELETED invar. mass with mmelel vertices", nbins, xmin, xmax) 
    hMdeleted.Sumw2()
    hRchi2deleted = ROOT.TH1F("hRchi2deleted", "DELETED reduced chi2 with mmelel vertices", 500, 0.0, 5.0) 
    hRchi2deleted.Sumw2()
    
if isMC:
    hMudRvPt = ROOT.TH2F("hMudRvPt", "min #DeltaR b/t gen,reco muon of same charge Vs. gen muon p_{T}", 1000, 0.0, 1.0, 1000, 0.0, 100.0) 
    hMudRvPtEta = ROOT.TH2F("hMudRvPtEta", "min #DeltaR b/t gen,reco muon of same charge Vs. gen eta p_{T}", 1000, 0.0, 1.0, 1000, 0.0, 100.0) 
    #histogram of pT for all gen eta mesons
    hpTGenAll = ROOT.TH1F("hpTGenAll", "pT of all gen Eta Mesons", 500, 0., 100.)
    hpTGenAll.Sumw2()
    #should all be exactly .548 GeV (just a cross-check)
    hMGenAll = ROOT.TH1F("hMGenAll", "mass of all gen Eta Mesons", nbins, xmin, xmax)
    hMGenAll.Sumw2()
    #ONLY with muon pt > 5 GeV (cuz otherwise it's a very biased comparison for eta) 
    hEtaGenAll = ROOT.TH1F("hEtaGenAll", "pseudorapidity of all gen Eta Mesons", 2000, -10., 10.) 
    hEtaGenAll.Sumw2()
    #gen eta mesons passing at least one trigger
    hpTGenTrig = ROOT.TH1F("hpTGenTrig", "pT of gen Eta Mesons that pass the trigger", 500, 0., 100.)
    hpTGenTrig.Sumw2()
    hEtaGenTrig = ROOT.TH1F("hEtaGenTrig", "pseudorapidity gen Eta Mesons that pass the trigger", 2000, -10., 10.) 
    hEtaGenTrig.Sumw2()
    hEtaGenTrigL0 = ROOT.TH1F("hEtaGenTrigL0", "pseudorapidity of gen Eta Mesons that pass HLT_Mu3_PFJet40", 2000, -10., 10.)
    hEtaGenTrigPFL1 = ROOT.TH1F("hEtaGenTrigPFL1", "pseudorapidity of gen Eta Mesons that pass HLT_Mu3_PFJet40 and HLT_Mu4_L1DoubleMu", 2000, -10., 10.)
    hEtaGenTrigL1 = ROOT.TH1F("hEtaGenTrigL1", "pseudorapidity of gen Eta Mesons that pass HLT_Mu4_L1DoubleMu", 2000, -10., 10.)
    hEtaGenTrigHLTL1 = ROOT.TH1F("hEtaGenTrigHLTL1", "pseudorapidity of gen Eta Mesons that pass HLT_Mu4_L1DoubleMu AND HLT_DoubleMu4_3_LowMass", 2000, -10., 10.)
    hpTGenReco = {}
    hpTGenNotReco = {}
    hpTGenAcc = {}
    #matched in dR instead of mass range
    hpTGenAccdR = {}
    hEtaGenReco = {}
    hEtaGenAcc = {}
    hvxy_gm = {}
    hEventWeight = {}
    hEvtWtVsPt = {}
    hsubVdRGenAllC = ROOT.TH2F("hsubVdRGenAllC", "sublead p_{T} vs. #DeltaR for Gen muons w/|#eta|<1.5", 500, 0.0, 100.0, 1000, 0.0, 1.0) 
    hsubVdRGenTrigC = ROOT.TH2F("hsubVdRGenTrigC", "sublead p_{T} vs. #DeltaR for muons passing trigger w/|#eta|<1.5", 500, 0.0, 100.0, 1000, 0.0, 1.0) 
    hsubVdRGenAllF = ROOT.TH2F("hsubVdRGenAllF", "sublead p_{T} vs. #DeltaR for all Gen muons w/|#eta|>1.5", 500, 0.0, 100.0, 1000, 0.0, 1.0) 
    hsubVdRGenTrigF = ROOT.TH2F("hsubVdRGenTrigF", "sublead p_{T} vs. #DeltaR for muons passing trigger w/|#eta|>1.5", 500, 0.0, 100.0, 1000, 0.0, 1.0) 
    #2d hist to compare 2 different xsection measurements
    #hxs0 = []
    #hxs1 = []
    #upper and lower values of hM considering evt weight uncertainties
    hMRecUp = {}
    hMRecDn = {}
    hMIDUp = {}
    hMIDDn = {}
    hMTTOCUp = {}
    hMTTOCDn = {}
    hMMuUp = {}
    hMMuDn = {}
    #number of modified electron efficiency points
    nMod = 12
    hMMod = {}
    #this modification is ON TOP of musf (if it is non-zero)
    hMModmu = {}
    hpTModmu = {}
    sf = [0.0 for s in range(nMod)]
    sfmu = [0.0 for s in range(nMod)]
    for s in range(nMod):
        ss = s if s < 6 else s+1
        sf[s] = -1.2 + 0.2*ss
        sfmu[s] = -0.3 + 0.05*ss
        hMMod[s] = ROOT.TH1F("hMMod"+str(s), "mmelel Invar. mass with electron p_{T} modified by "+str(sf[s])+"\%", 1000, 0.0, 1.0) 
        hMModmu[s] = ROOT.TH1F("hMModmu"+str(s), "mmelel Invar. mass with muon p_{T} modified by "+str(sfmu[s])+"\%", 1000, 0.0, 1.0) 
        hpTModmu[s] = ROOT.TH1F("hpTModmu"+str(s), "p_{T} mod'd by "+str(sfmu[s])+"\%, evts in sig peak .52-.58 GeV for"+str(vtype)+" verts", 500, 0.0, 100.0) 
    for vtype in vtypes:
        if vtype in ["mumu", "elel", "pcpc", "lplp"]:
            hMRecUp[vtype] = ROOT.TH1F("hMRecUp"+vtype, "Invar. mass with "+vtype+" vertices, electron reco SF varied up", 1000, 0.0, 1.0) 
            hMRecDn[vtype] = ROOT.TH1F("hMRecDn"+vtype, "Invar. mass with "+vtype+" vertices, electron reco SF varied down", 1000, 0.0, 1.0) 
            hMIDUp[vtype] = ROOT.TH1F("hMIDUp"+vtype, "Invar. mass with "+vtype+" vertices, electron ID SF varied up", 1000, 0.0, 1.0) 
            hMIDDn[vtype] = ROOT.TH1F("hMIDDn"+vtype, "Invar. mass with "+vtype+" vertices, electron ID SF varied down", 1000, 0.0, 1.0) 
            hMTTOCUp[vtype] = ROOT.TH1F("hMTTOCUp"+vtype, "Invar. mass with "+vtype+" vertices, trigger turn-on SF varied up", 1000, 0.0, 1.0) 
            hMTTOCDn[vtype] = ROOT.TH1F("hMTTOCDn"+vtype, "Invar. mass with "+vtype+" vertices, trigger turn-on SF varied down", 1000, 0.0, 1.0) 
            hMMuUp[vtype] = ROOT.TH1F("hMMuUp"+vtype, "Invar. mass with "+vtype+" vertices, muon ID SF varied up", 1000, 0.0, 1.0) 
            hMMuDn[vtype] = ROOT.TH1F("hMMuDn"+vtype, "Invar. mass with "+vtype+" vertices, muon ID SF varied down", 1000, 0.0, 1.0) 
        else:
            hMRecUp[vtype] = ROOT.TH1F("hMRecUp"+vtype, "Invar. mass with "+vtype+" vertices, electron reco SF varied up", nbins, xmin, xmax) 
            hMRecDn[vtype] = ROOT.TH1F("hMRecDn"+vtype, "Invar. mass with "+vtype+" vertices, electron reco SF varied down", nbins, xmin, xmax) 
            hMIDUp[vtype] = ROOT.TH1F("hMIDUp"+vtype, "Invar. mass with "+vtype+" vertices, electron ID SF varied up", nbins, xmin, xmax) 
            hMIDDn[vtype] = ROOT.TH1F("hMIDDn"+vtype, "Invar. mass with "+vtype+" vertices, electron ID SF varied down", nbins, xmin, xmax) 
            hMTTOCUp[vtype] = ROOT.TH1F("hMTTOCUp"+vtype, "Invar. mass with "+vtype+" vertices, trigger turn-on SF varied up", nbins, xmin, xmax) 
            hMTTOCDn[vtype] = ROOT.TH1F("hMTTOCDn"+vtype, "Invar. mass with "+vtype+" vertices, trigger turn-on SF varied down", nbins, xmin, xmax) 
            hMMuUp[vtype] = ROOT.TH1F("hMMuUp"+vtype, "Invar. mass with "+vtype+" vertices, muon ID SF varied up", nbins, xmin, xmax) 
            hMMuDn[vtype] = ROOT.TH1F("hMMuDn"+vtype, "Invar. mass with "+vtype+" vertices, muon ID SF varied down", nbins, xmin, xmax) 
        #histogram of surviving event weights
        hEventWeight[vtype] = ROOT.TH1F("hEventWeight"+vtype, "Event weights surviving all cuts", 1000, 0.0, 100.0)
        hEventWeight[vtype].Sumw2()
        hEvtWtVsPt[vtype] = ROOT.TH2F("hEvtWtVsPt"+vtype, "Event weights surviving all cuts as a function of pT", 100, 0.0, 100.0, 100, 0.0, 100.0) 
        hEvtWtVsPt[vtype].Sumw2()
        #gen eta mesons reco'd using packed candidates
        hpTGenReco[vtype] = ROOT.TH1F("hpTGenReco"+vtype, "pT of gen Eta Mesons that are reconstructed with "+vtype+" vertices", 500, 0., 100.)
        hpTGenNotReco[vtype] = ROOT.TH1F("hpTGenNotReco"+vtype, "pT of gen Eta Mesons that are NOT reconstructed with "+vtype+" vertices", 500, 0., 100.)
        hpTGenReco[vtype].Sumw2()
        #ONLY with muon pt > 5 GeV (cuz otherwise it's a very biased comparison for eta) 
        hEtaGenReco[vtype] = ROOT.TH1F("hEtaGenReco"+vtype, "pseudorapidity of gen Eta Mesons that are reco'd w/"+vtype+" vertices", 2000, -10., 10.) 
        hEtaGenReco[vtype].Sumw2()
        #gen eta mesons passing at least one trigger AND reco'd with packed cand's
        hpTGenAcc[vtype] = ROOT.TH1F("hpTGenAcc"+vtype, "pT of gen Eta Mesons that are accepted (trig+reco) with "+vtype+" vertices", 500, 0., 100.)
        hpTGenAcc[vtype].Sumw2()
        hpTGenAccdR[vtype] = ROOT.TH1F("hpTGenAccdR"+vtype, "p_{T} of accepted gen #eta mesons (#DeltaR<"+str(dRcut)+") with "+vtype+" vertices", 500, 0., 100.)
        hpTGenAccdR[vtype].Sumw2()
        #ONLY with muon pt > 5 GeV (cuz otherwise it's a very biased comparison for eta) 
        hEtaGenAcc[vtype] = ROOT.TH1F("hEtaGenAcc"+vtype, "pseudorapidity of gen Eta Mesons that are accepted (trig+reco) w/"+vtype+" vertices", 2000, -10., 10.) 
        hEtaGenAcc[vtype].Sumw2()
        #vxy values for ONLY genmatched (LowPt)electron vertices
        hvxy_gm[vtype] = ROOT.TH1F("hvxygm"+vtype, "Vxy for genmatched "+vtype+" vertices", 1000, 0.0, 10.0) 
        hvxy_gm[vtype].Sumw2()
    hGenMudR = ROOT.TH1F("hGenMudR", "dR b/t gen muons", 1000, 0.0, 1.0) 
    hGenMudR.Sumw2()
    #leading muon pT
    hGenMupT0 = ROOT.TH1F("hGenMupT0", "pT of lead gen muon", 10000, 0.0, 100.0)
    #hGenMupT0.Sumw2()
    #subleading muon pT
    hGenMupT1 = ROOT.TH1F("hGenMupT1", "pT of subleading gen muon", 10000, 0.0, 100.0)
    hGenMupTAllC = ROOT.TH1F("hGenMupTAllC", "pT of all gen muons with -1.5<#eta^{#mu}<1.5", 10000, 0.0, 100.0)
    hGenMupTAllF = ROOT.TH1F("hGenMupTAllF", "pT of all gen muons with |#eta^{#mu}|>1.5", 10000, 0.0, 100.0)
    #hGenMupT1.Sumw2()
    hGenMupTRec = ROOT.TH1F("hGenMupTRec", "p_{T} of gen-matched gen muons", 10000, 0.0, 100.0) 
    #electron pT
    hGenElpTAll = ROOT.TH1F("hGenElpTAll", "p_{T} of gen electrons", 10000, 0.0, 100.0) 
    hGenElpTRec = ROOT.TH1F("hGenElpTRec", "p_{T} of gen-matched gen electrons", 10000, 0.0, 100.0) 
    #electron dR
    hGenEldRvPtAll = ROOT.TH2F("hGenEldRvPtAll", "gen dR vs. p_{T} for all gen electrons", 8000, 0.0, 0.2, 100, 0.0, 50.0) 
    hGenEldRvPtRec = ROOT.TH2F("hGenEldRvPtRec", "gen dR vs. p_{T} gen-matched gen electrons", 8000, 0.0, 0.2, 100, 0.0, 50.0) 
    #hGenElpT.Sumw2()

    if isSig:
        hEldRvPtEta = ROOT.TH2F("hEldRvPtEta", "min #DeltaR b/t gen,reco electron of same charge Vs. gen eta p_{T}", 1000, 0.0, 1.0, 1000, 0.0, 100.0) 
        #dR between reconstructed tracks and gen electrons
        hdR = {}
        for ptype in ptypes:
            hdR[ptype] = ROOT.TH1F("hdR"+ptype, "dR b/t reco "+ptype+" and gen electrons", 10000, 0., 1.) 
            hdR[ptype].Sumw2()

    #open the xsec file to get the pT-dependent xsec's (so can get weighted overall efficiencies)
    #f_xsec = ROOT.TFile.Open("xsecs.root")
    #f_xsec = ROOT.TFile.Open("xsec2022_nominal.root")
    ##h_xsec = f_xsec.Get("corr_xsec")
    #h_xsec = f_xsec.Get("hXsecCor")
    if do_pileup:
        PUFile = ROOT.TFile.Open("pileup_corrections%s_%d.root"%("_newWt" if new_wt else "", year))
        if isSig:
            PUCor = PUFile.Get("hSigCorr%d"%year)
        elif isMuMu:
            PUCor = PUFile.Get("hRefCorr%d"%year)
        else:
            if not new_wt:
                PUCor = PUFile.Get("hBkgCorr%d"%year)
            else:
                PUCor = PUFile.Get("hSigCorr%d"%year) 
    #if do_trigCor:
    #    trigfile = ROOT.TFile.Open("trigger_corrections_%sMC2022.root"%("mumu" if isMuMu else "sig")) 
    #    #2d hist
    #    trigCor = trigfile.Get("%sMCcorrection"%("mumu" if isMuMu else "sig"))
    if do_TTOC:
        trigfileF = ROOT.TFile.Open("trigger_correctionsF.root")
        trigfileC = ROOT.TFile.Open("trigger_correctionsC.root")
        triggerCorrectionF = trigfileF.Get("trigger_correctionsF")     
        triggerCorrectionC = trigfileC.Get("trigger_correctionsC")     
        print("Opened triggerCorrectionF: %s"%triggerCorrectionF.Print())
        print("Opened triggerCorrectionC: %s"%triggerCorrectionC.Print())
    
    #total weight of all MC events
    all_weight = 0.
    all_weightUp = 0.
    all_weightDn = 0.
    #total weight of MC events passing trigger
    trg_weight = 0.
    trg_weightUp = 0.
    trg_weightDn = 0.
    rec_weight = {}
    rec_weightUp = {}
    rec_weightDn = {}
    #how many different selections to calculate the acceptance weights for
    nselections = 9
    acc_weight = {}
    acc_weightUp = {}
    acc_weightDn = {}
    accdR_weight = {}
    accdR_weightUp = {}
    accdR_weightDn = {}
    for vtype in vtypes:
        #total weight of MC events passing reco
        rec_weight[vtype] = 0.
        rec_weightUp[vtype] = 0.
        rec_weightDn[vtype] = 0.
        #total weight of events accepted (trg and reco)
        acc_weight[vtype] = [0. for s in range(nselections+1)]
        acc_weightUp[vtype] = 0.
        acc_weightDn[vtype] = 0.
        accdR_weight[vtype] = 0.
        accdR_weightUp[vtype] = 0.
        accdR_weightDn[vtype] = 0.

inc_vertM = False
if inc_vertM:
    #mass hist as saved per vertex (should be same as hM)
    hMV = ROOT.TH1F("hMV", "Invar. mass with tracks in vertex", 500, .4, .9) 
    hpTV = ROOT.TH1F("hpTV", "pT with tracks in vertex", 500, 0, 100.) 

if arg < 27:
    subdir = "0000"
else:
    subdir = "0001"

#open file to write the event nums and masses
if syncTest:
    #syncFname = "syncTest_%d%s_%d_test%d_mmelel.txt"%(num, let, arg, testnum)
    if singleFile:
        syncFname = "syncTest_test%d_mumu.txt"%(testnum)
    else:
        syncFname = "syncTest_%d%s_%d_test%d_mumu.txt"%(num, let, arg, testnum)
    syncFile = open(syncFname, "w") 
    #dan_runs,dan_events = readEvents("danEvents.txt")
#if not isMC and not singleFile:
#    trigFname = "triggerFailed_test%d_%d%s_%d.txt"%(testnum, num, let, arg)
#    trigF = open(trigFname, "w") 

#count the number of error events
#nerr = 0

#how frequently to print out the event number
printevery = 10000
#printevery = 1
#for let in lets:
#    for num in nums:
        #inname = "test_{0:s}*.root".format(dig)
        #path = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ParkingDoubleMuonLowMass%d/test2_%d%s/*/*/test_*.root"%(num, num, let)
        #path = "/eos/uscms/store/user/bgreenbe/BParking2022/ParkingDoubleMuonLowMass%d/test2_%d%s/*/*/test_*.root"%(num, num, let)
        #path = "/eos/uscms/store/user/bgreenbe/BParking2022/ParkingDoubleMuonLowMass{0:d}/test2_{0:d}{1:s}/*/{2:s}/test_{3:s}*.root".format(num, let, subdir, dig)
        #path = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ParkingDoubleMuonLowMass{0:d}/test2_{0:d}{1:s}/*/{2:s}/{3:s}".format(num, let, subdir, inname)

#to be used for hMComb
old_diel = None
diel = None

#fit for the electron efficiency
def elEff(elpt):
    #if elpt > 27:
    #    return elEff(27)
    #p0                        =      5.43551  
    #p1                        =      1.47295  
    #p2                        =      26.8049  
    #p3                        =       9.7099  
    #p4                        =    0.0245243  
    #p5                        = -2.19139e-06  
    #p6                        =     -161.112
    #return p0 / (p1*elpt + ROOT.TMath.Exp(-(elpt - p2)/p3)) - p4/elpt + p5*(elpt-p6)**2
    p0                        =       203.27 
    p1                        =      2.38773 
    p2                        =     0.331089 
    p3                        =      4.79575 
    p4                        =       184.13 
    p5                        =      1081.15 
    return p0 / (ROOT.TMath.Exp(ROOT.TMath.Exp(-(elpt - p1)/p2) + p3)+p4+p5/elpt)

#fit for the muon efficiency
def muEff(mupt):
    ##dR cut 0.13
    #p0     =      2.01098e+01
    #p1     =     -3.42298e-01
    #p2     =      9.82785e-01
    #p3     =      5.64684e+00
    #p4     =     -2.62466e+02
    #p5     =      1.18172e+01
    #dR cut 0.03
    p0                        =      93.0934  
    p1                        =      1.74031  
    p2                        =      1.05146  
    p3                        =      5.07463  
    p4                        =     -64.5288  
    p5                        =      49.4846  
    return p0 / (ROOT.TMath.Exp(ROOT.TMath.Exp(-(mupt - p1)/p2) + p3)+p4+p5/mupt)

#randomly assign a trigger pass/fail value according to the fitted model as a function of subleading muon pT (threshold varied by trg_var)
def toyTrig(subpt): #, subeta=9999):
    #if subeta == 9999:
    #    plat = 0.7336 #0.8754
    #    w = 0.1714 #0.4117 #* (1+trg_var)
    #    #p50 = 3.508 * (1+trg_var)
    #    p50 = 3.517 * (1+trg_var)
    #elif abs(subeta) < 1.5:
    #    plat = 0.8060
    #    w = 0.1581
    #    p50 = 3.6065 * (1+trg_var)
    #else:
    #    plat = 0.6062
    #    w = 0.1350
    #    p50 = 3.2678 * (1+trg_var)
    #for eta->2mu
    if isMuMu:
        plat = 0.2867
        p50 = 3.345 * (1+trg_var)
        w = 0.4625
    else:
        #signal
        plat = 2.77351e-01 #0.7539
        p50 = 3.38611  * (1+trg_var) #0.7190
        w = 3.89334e-01 #0.3887
    #alpha = 953.9
    #efficiency at this pT value
    #toyEff = plat / (1 + ROOT.TMath.Exp(-(subpt - p50) / w)) 
    #toyEff = plat / (1 + ROOT.TMath.Exp(-(subpt - p50) / w))**alpha 
    toyEff = plat / (ROOT.TMath.Exp(ROOT.TMath.Exp(-(subpt - p50) / w) - 1))
    #old fit doesn't work well anymore because now using only a single trigger path.
    #new fit: [0]/((x-[1])/[2]+[3]/(x-[4])+TMath::Exp(-(x-[5])/[6]))
    #p0 = 4.256
    #p1 = -1.696
    #p2 = 2.119
    #p3 = 6.342
    #p4 = 2.598
    #p5 = 3.641
    #p6 = 0.171
    #if subpt > 8:
    #    subpt = 8
    #toyEff = p0 / ( (subpt - p1) / p2 + p3 / (subpt - p4) + ROOT.TMath.Exp(-(subpt - p5) / p6) )

    #randomly generate a number 0-1, if it's lower than the efficiency, we pass.
    ran = random.random()
    return (ran < toyEff)

#process all the vertices of this type.
# vtype: string for vertex type: elel, lplp, mmelel, mmlplp -- for elel and lplp also use mumu vertices
# singleVert: true for only ONE vertex, false for all
# useOnia: true to remove any tracks matched to Onia conversion candidates, false to not
# g: gen event (if MC only, obviously)
#def process_vertices(e, vtype, singleVert, useOnia, xsec, evt_weight, g=None, genEtaPt=0, passedTrig=True):
def process_vertices(e, vtype, singleVert, useOnia, xsec, evt_weight, g=None, gen_eta=None, passedTrig=True, recod=False):

    #print("vtx %s Run %d ls %d evt %d"%(vtype, e.run, e.lumi_sec, e.evt)) 

    #is good eta or nah
    good = False
    bestj = -1
    if "elel" in vtype:
        #fill the comb bkg hist with electron pair from last event, EVEN if there's no electrons in this event!! 
        for mV in range(len(e.Vertex_mumu_muP)):
            #form the test for combinatorial bkg by combining this event's dimuon pair with last event's dielectron
            vec_muP = ROOT.TLorentzVector()
            vec_muN = ROOT.TLorentzVector()
            mp = ord(e.Vertex_mumu_muP[mV])
            mn = ord(e.Vertex_mumu_muN[mV])
            vec_muP.SetPtEtaPhiM(e.Muon_pt[mp]*(1+musf), e.Muon_eta[mp], e.Muon_phi[mp], mu_mass)
            vec_muN.SetPtEtaPhiM(e.Muon_pt[mn]*(1+musf), e.Muon_eta[mn], e.Muon_phi[mn], mu_mass)
            dimu = vec_muP + vec_muN
            global old_diel
            if old_diel != None:
                comb = dimu + old_diel
                hMComb.Fill(comb.M(), evt_weight)
                #now set old_diel to None so that the same one isn't used again in the future.
                old_diel = None

        if ord(e.nGoodElectron) < 2:
            return False
        lptstr = ""
    elif "lplp" in vtype:
        if ord(e.nGoodLowPtElectron) < 2:
            return False
        lptstr = "LowPt"
    elif vtype in ["mumu", "mmg"]:
        if ord(e.nGoodMuon) < 2:
            return False
        lptstr = ""
    else:
        sys.exit("Error: unrecognized vertex type %s"%(vtype)) 
    nvert = None
    #the string that this type of vertex was called in the ntuples
    vstr = vtype
    if vtype == "mmg":
        vstr = "mumu"
    nvert = eval("len(e.Vertex_%s_vz)"%(vstr))

    try:
        vtx_vchi2 = eval("e.Vertex_%s_chi2"%(vstr))
        vtx_vndof = eval("e.Vertex_%s_ndof"%(vstr))
        vtx_vrechi2 = [vtx_vchi2[a] / vtx_vndof[a] for a in range(nvert)] 
        vtx_vprob = [ROOT.TMath.Prob(vtx_vchi2[a], vtx_vndof[a]) for a in range(nvert)] 
    except:
        vtx_vrechi2 = eval("e.Vertex_%s_reduced_chi2"%(vstr))
    if "elel" in vstr or "lplp" in vstr:
        vtx_veleP = eval("e.Vertex_%s_eleP"%(vstr))
        vtx_veleN = eval("e.Vertex_%s_eleN"%(vstr))
    vtx_svxy = eval("e.Vertex_%s_sigmaVxy"%(vstr))
    try:
        vtx_vvx = eval("e.Vertex_%s_vx"%(vstr))
        vtx_vvy = eval("e.Vertex_%s_vy"%(vstr))
    except:
        vtx_vxy = eval("e.Vertex_%s_vxy"%(vstr)) 
    if vstr in ["elel", "lplp"]:
        vtx_vmuP = e.Vertex_mumu_muP
        vtx_vmuN = e.Vertex_mumu_muN
    else:
        vtx_vmuP = eval("e.Vertex_%s_muP"%(vstr))
        vtx_vmuN = eval("e.Vertex_%s_muN"%(vstr))

    #delet invalid mmelel vertices
    if mmelelExclusive and vtype == "mmelel":
        #loop thru all vertices
        j = 0
        while j < nvert:
            #check if this mu pair is in the list of mu pairs
            found_eVtx = False
            for k in range(len(e.Vertex_elel_eleP)):
                if ord(e.Vertex_elel_eleP[k]) == ord(vtx_veleP[j]) and ord(e.Vertex_elel_eleN[k]) == ord(vtx_veleN[j]) :
                    found_eVtx = True
                    break
            found_muVtx = False
            #don't bother wasting time on the mu vtx loop if electron vtx wasn't even found
            if found_eVtx:
                for k in range(len(e.Vertex_mumu_muP)):
                    if ord(e.Vertex_mumu_muP[k]) == ord(vtx_vmuP[j]) and ord(e.Vertex_mumu_muN[k]) == ord(vtx_vmuN[j]) :
                        found_muVtx = True
                        break

            #if the electron vertex or muon vertex was not found, delet this mmelel vertex.
            if not (found_eVtx and found_muVtx):
                #before deletting, fill the 'deletted' invariant mass histogram
                vec_elP = ROOT.TLorentzVector()
                vec_elN = ROOT.TLorentzVector()
                vec_muP = ROOT.TLorentzVector()
                vec_muN = ROOT.TLorentzVector()
                vec_elP.SetPtEtaPhiM(e.Electron_pt[ord(vtx_veleP[j])]*(1+elsf), e.Electron_eta[ord(vtx_veleP[j])], e.Electron_phi[ord(vtx_veleP[j])], el_mass)
                vec_elN.SetPtEtaPhiM(e.Electron_pt[ord(vtx_veleN[j])]*(1+elsf), e.Electron_eta[ord(vtx_veleN[j])], e.Electron_phi[ord(vtx_veleN[j])], el_mass)
                vec_muP.SetPtEtaPhiM(e.Muon_pt[ord(vtx_vmuP[j])]*(1+musf), e.Muon_eta[ord(vtx_vmuP[j])], e.Muon_phi[ord(vtx_vmuP[j])], mu_mass)
                vec_muN.SetPtEtaPhiM(e.Muon_pt[ord(vtx_vmuN[j])]*(1+musf), e.Muon_eta[ord(vtx_vmuN[j])], e.Muon_phi[ord(vtx_vmuN[j])], mu_mass)
                mdeleted = (vec_elP + vec_elN + vec_muP + vec_muN).M()
                hMdeleted.Fill(mdeleted)
                hRchi2deleted.Fill(vtx_vrechi2[j]) 
                #print("EVT " + str(i) + ": removing " + str(j) + " before: " + str(vtx_vrechi2)) 
                try:
                    vtx_vchi2.erase(vtx_vchi2.begin() + j)
                    vtx_vndof.erase(vtx_vndof.begin() + j)
                    vtx_vrechi2.remove(j)
                    vtx_vprob.remove(j)
                except:
                    vtx_vrechi2.erase(vtx_vrechi2.begin() + j)
                vtx_veleP.erase(vtx_veleP.begin() + j)
                vtx_veleN.erase(vtx_veleN.begin() + j)
                vtx_vmuP.erase(vtx_vmuP.begin() + j)
                vtx_vmuN.erase(vtx_vmuN.begin() + j)
                vtx_vvx.erase(vtx_vvx.begin() + j)
                vtx_vvy.erase(vtx_vvy.begin() + j)
                vtx_svxy.erase(vtx_svxy.begin()+j)
                j -= 1
                nvert -= 1
            j += 1

    #if requested, genmatch the converted photon tracks to the lowPt electrons, discard any that are a match
    if useOnia and vtype != "mmg":
        #find the best genmatches
        for oni in range(ord(e.nOnia)):
            #make the LorentzVectors for the pos and neg onia tracks
            vec_oni0 = ROOT.TLorentzVector()
            vec_oni1 = ROOT.TLorentzVector()
            vec_oni0.SetPtEtaPhiM(e.Onia_pt0[oni], e.Onia_eta0[oni], e.Onia_phi0[oni], el_mass)
            vec_oni1.SetPtEtaPhiM(e.Onia_pt1[oni], e.Onia_eta1[oni], e.Onia_phi1[oni], el_mass) 
            mindRP = 9999.0
            minelP = -1
            mindRN = 9999.0
            minelN = -1
            for j in range(nvert):
                vec_elP = ROOT.TLorentzVector()
                vec_elN = ROOT.TLorentzVector()
                elP = ord(vtx_veleP[j])
                elN = ord(vtx_veleN[j])
                elptP = eval("e.%sElectron_pt[elP]"%(lptstr))*(1+elsf)
                elptN = eval("e.%sElectron_pt[elN]"%(lptstr))*(1+elsf)
                eletaP = eval("e.%sElectron_eta[elP]"%(lptstr))
                eletaN = eval("e.%sElectron_eta[elN]"%(lptstr))
                elphiP = eval("e.%sElectron_phi[elP]"%(lptstr))
                elphiN = eval("e.%sElectron_phi[elN]"%(lptstr))
                vec_elP.SetPtEtaPhiM(elptP, eletaP, elphiP, el_mass)
                vec_elN.SetPtEtaPhiM(elptN, eletaN, elphiN, el_mass)
                if e.Onia_charge0[oni] == 1:
                    dRP = vec_elP.DeltaR(vec_oni0)
                    dRN = vec_elN.DeltaR(vec_oni1)
                else:
                    dRP = vec_elP.DeltaR(vec_oni1)
                    dRN = vec_elN.DeltaR(vec_oni0)
                if dRP < mindRP:
                    mindRP = dRP
                    minelP = elP
                    #print("New mindRP: " + str(mindRP))
                if dRN < mindRN:
                    mindRN = dRN
                    minelN = elN
                    #print("New mindRN: " + str(mindRN))

            vx = e.Onia_vx[oni]
            vy = e.Onia_vy[oni]
            #now if the dR is less than the cutoff, remove this vertex from the list of vertices!
            if mindRP < dRcut or mindRN < dRcut:
                #mindR = min(mindRP, mindRN)
                #print("mindR: " + str(mindR) + ": less than dRcut!! elP: " + str(minelP))
                #fill the Matched hist
                hVertMatched.Fill(vx, vy)
                if isMC:
                    #hpTVertMatched.Fill(genEtaPt, evt_weight)
                    hpTVertMatched.Fill(gen_eta.Pt(), evt_weight)
                #remove ALL vertices containing the mindR, both pos and neg
                j = 0
                while j < nvert:
                    if (ord(vtx_veleP[j]) == minelP and mindRP < dRcut) or (ord(vtx_veleN[j]) == minelN and mindRN < dRcut):
                        #print("EVT " + str(i) + ": removing " + str(j) + " before: " + str(vtx_vrechi2)) 
                        try:
                            vtx_vchi2.erase(vtx_vchi2.begin() + j)
                            vtx_vndof.erase(vtx_vndof.begin() + j)
                            vtx_vrechi2.remove(j)
                            vtx_vprob.remove(j)
                        except:
                            vtx_vrechi2.erase(vtx_vrechi2.begin() + j)
                        vtx_veleP.erase(vtx_veleP.begin() + j)
                        vtx_veleN.erase(vtx_veleN.begin() + j)
                        if vstr in ["mmelel", "mmlplp", "mmee"]:
                            vtx_vmuP.erase(vtx_vmuP.begin() + j)
                            vtx_vmuN.erase(vtx_vmuN.begin() + j)
                        vtx_vvx.erase(vtx_vvx.begin() + j)
                        vtx_vvy.erase(vtx_vvy.begin() + j)
                        vtx_svxy.erase(vtx_svxy.begin()+j)
                        j -= 1
                        nvert -= 1
                    j += 1
            #elif mindRP < 9999:
                #print("mindRP: " + str(mindRP)) 
            else:
                hVertNoMatch.Fill(vx, vy) 

            #hdRP[vtype].Fill(mindRP)
            #hdRN[vtype].Fill(mindRN)

    danEvent = False
    #if syncTest and vtype == "mumu" and e.evt in dan_events:
    #    idx = dan_events.index(e.evt)
    #    if dan_runs[idx] == e.run:
    #        print("**danEvent**")
    #        danEvent = True
    #        try:
    #            printEvent.printEvent(e)
    #        except AttributeError:
    #            printEvent_backup.printEvent(e)

    #fill the pT, M histograms after looping thru all vertex candidates? or in the middle
    # (for the 4-lepton vertices and mumu, only executing the loop once anyway so minus well just fill in the middle)
    # JK now executing multiple times so fill at end for 4-leptons!!!
    #fillAtEnd = (vstr not in ["mmelel", "mmlplp", "mmee", "mumu"] and singleVert)
    fillAtEnd =  singleVert

    if vProbCut > 0:
        #These reduced chi2 cuts correspond to vProb cutoff of 0.1 ONLY -- MUST BE MODIFIED IF YOU CHANGE vProbCut!!!!
        rChi2Cut = 1.84727
        if vtype in ["mumu", "elel", "lplp", "pcpc"]:
            rChi2Cut = 2.70554
    else:
        rChi2Cut = -1

    Vxy = -1
    if singleVert and nvert > 0:
        #bestj, bestChi2 = min(enumerate(e.Vertex_lplp_reduced_chi2), key=lambda x: x[1])
        #bestj, bestChi2 = min(enumerate(e.Vertex_mmlplp_reduced_chi2), key=lambda x: x[1])
        #bestj, bestChi2 = min(enumerate(e.Vertex_mmlplp_reduced_chi2), key=lambda x: x[1])
        #no longer really used for anything!!
        #bestj, bestChi2 = min(enumerate(vtx_vrechi2), key=lambda x: x[1])
        try:
            bestj, bestprob = max(enumerate(vtx_vprob), key=lambda x: x[1])
            bestChi2 = vtx_vrechi2[bestj]
        except:
            bestj, bestChi2 = min(enumerate(vtx_vrechi2), key=lambda x: x[1])
    bestm = 99999
    bestm2mu = 99999
    bestm2el = 99999
    bestpt = -1
    if isMC:
        acc_filled = [False for s in range(nselections+1)] 
    accdR_filled = False
    elptP = -99
    elptN = -99
    #how many candidates are filled into the histogram in this event?
    ncand = 0
    for j in range(nvert):
        if singleVert and j != bestj : continue
        #try a cut on the chi2 value? or on vxy??
        #print("evt = %d, j = %d, nvert = %d"%(i, j, nvert)) 
        #chi2 of 2.6 corresponds to prob of .1
        #if rChi2Cut > 0 and vtx_vrechi2[j] > rChi2Cut: continue
        try:
            failed = vProbCut > 0 and vtx_vprob[j] < vProbCut
        except:
            failed = rChi2Cut > 0 and vtx_vrechi2[j] > rChi2Cut
        if failed: continue

        if danEvent:
            print(vtype + " vtx " + str(j)) 

        if vstr != "mumu":
            vec_elP = ROOT.TLorentzVector()
            vec_elN = ROOT.TLorentzVector()
            vec_piP = ROOT.TLorentzVector()
            vec_piN = ROOT.TLorentzVector()
            elP = ord(vtx_veleP[j])
            pt = eval("e.%sElectron_pt[elP]"%(lptstr))*(1+elsf) 
            eta = eval("e.%sElectron_eta[elP]"%(lptstr)) 
            phi = eval("e.%sElectron_phi[elP]"%(lptstr))
            vec_elP.SetPtEtaPhiM(pt, eta, phi, el_mass)
            vec_piP.SetPtEtaPhiM(pt, eta, phi, pi_mass)
            elN = ord(vtx_veleN[j])
            pt = eval("e.%sElectron_pt[elN]"%(lptstr))*(1+elsf) 
            eta = eval("e.%sElectron_eta[elN]"%(lptstr)) 
            phi = eval("e.%sElectron_phi[elN]"%(lptstr))
            vec_elN.SetPtEtaPhiM(pt, eta, phi, el_mass)
            vec_piN.SetPtEtaPhiM(pt, eta, phi, pi_mass)
            elIDP = eval("ord(e.%sElectron_id[elP])"%(lptstr))
            elIDN = eval("ord(e.%sElectron_id[elN])"%(lptstr))
            #loose ID
            looseID_p = elIDP & 0b00110000 
            looseID_n = elIDN & 0b00110000 
            #WP90 (nominal selection)
            WP90ID_p = elIDP & 0b00000100 
            WP90ID_n = elIDN & 0b00000100 
            #if WP90ID_p == 0 or WP90ID_n == 0:
            #if WP90ID_p == 0 and WP90ID_n == 0:
            #WP80
            WP80ID_p = elIDP & 0b00000001
            WP80ID_n = elIDN & 0b00000001
            ##require looseID on BOTH electrons
            #if looseID_p == 0 or looseID_n == 0:
            #if WP80ID_p == 0 or WP80ID_n == 0:
            #require elID on BOTH electrons
            if basic_cuts and "elel" in vtype:
                vetoP = ord(e.Electron_convVeto[elP]) != 0
                vetoN = ord(e.Electron_convVeto[elN]) != 0
                if vetoP or vetoN:
                    continue
                nMissP = ord(e.Electron_nMissingHits[elP])
                nMissN = ord(e.Electron_nMissingHits[elN])
                if nMissP > 3 or nMissN > 3:
                #if nMissP > 2 or nMissN > 2:
                #if nMissP > 1 or nMissN > 1:
                #if nMissP > 0 or nMissN > 0:
                    continue
            #for purposes of electron efficiency uncertainty calculation, need to move this further down.
            #if elpTcut > -1 and (e.Electron_pt[elP] < elpTcut or e.Electron_pt[elN] < elpTcut):
            #    continue
            if "elel" in vtype and elEtacut > -1 and (abs(e.Electron_eta[elP]) > elEtacut or abs(e.Electron_eta[elN]) > elEtacut):
                continue
            if "elel" in vtype and cut_trans and ( (abs(e.Electron_eta[elP]) > 1.44 and abs(e.Electron_eta[elP]) < 1.57) or (abs(e.Electron_eta[elN]) > 1.44 and abs(e.Electron_eta[elN]) < 1.57) ):
                continue
            if vtype == "mmelel":
                hpTElAll.Fill(e.Electron_pt[elP]*(1+elsf), evt_weight)
                hpTElAll.Fill(e.Electron_pt[elN]*(1+elsf), evt_weight) 
                if WP90ID_p > 0:
                    hpTElIDPass.Fill(e.Electron_pt[elP]*(1+elsf), evt_weight)
                if WP90ID_n > 0:
                    hpTElIDPass.Fill(e.Electron_pt[elN]*(1+elsf), evt_weight)
            if (require_elID and "lplp" not in vtype) and (WP90ID_p == 0 or WP90ID_n == 0):
                continue
            if (invert_elID and "elel" in vtype):
                vetoP = ord(e.Electron_convVeto[elP]) != 0
                vetoN = ord(e.Electron_convVeto[elN]) != 0
                if not vetoP and not vetoN: 
                    continue
            if (tightEl or "lplp" in vtype) and require_elID and (WP80ID_p == 0 or WP80ID_n == 0):
                continue
            if invert_elID and "lplp" in vtype and (WP80ID_p != 0 or WP80ID_n != 0):
                continue
            try:
                vxy = (vtx_vvx[j]**2 + vtx_vvy[j]**2)**0.5
            except:
                vxy = vtx_vxy[j] 
            Vxy = vxy
            #to be used for hMComb
            global diel
            diel = vec_elP + vec_elN
            if isMC:
                if isSig:
                    gm = False
                    #find the dR betwixt the track and the genElectron of the same charge
                    for k in range(g.nGenPart):
                        if not abs(g.GenPart_pdgId[k]) == 11 : continue
                        if ord(g.GenPart_charge[k]) == 1 :
                            reco_el = vec_elP
                        else:
                            reco_el = vec_elN
                        gen_ele = ROOT.TLorentzVector()
                        gen_ele.SetPtEtaPhiM(g.GenPart_pt[k], g.GenPart_eta[k], g.GenPart_phi[k], el_mass)
                        dr = reco_el.DeltaR(gen_ele) 
                        if "lplp" in vstr:
                            hdR["lowpTElectron"].Fill(dr)
                        elif "elel" in vstr:
                            hdR["patElectron"].Fill(dr) 
                        if dr < dRcut:
                            gm = True 
                            #fill half here and half for the other ele
                            hvxy_gm[vtype].Fill(vxy, 0.5)
                else:
                    #not signal -- background MC
                    #has photon or nah
                    hasG = False
                    for k in range(g.nGenPart):
                        if g.GenPart_pdgId[k] == 22:
                            hasG = True
                            #break
                        ##if it's the original eta meson, see if it's gen-matched?
                        #if g.GenPart_pdgId[k] == 221:
                        #    gen_eta = ROOT.TLorentzVector()
                        #    gen_eta.SetPtEtaPhiM(g.GenPart_pt[k], g.GenPart_eta[k], g.GenPart_phi[k], etamass)
                    if vstr == "lplp":
                        mee = (vec_elP + vec_elN).M()
                        if hasG:
                            hMeeG[vtype].Fill(mee, evt_weight)
                        else:
                            hMeeNoG[vtype].Fill(mee, evt_weight)
                #if isSig and gm:
                #    #hpTGenReco[vtype].Fill(genEtaPt)
                #    hpTGenReco[vtype].Fill(gen_eta.Pt())
                #    rec_weight[vtype] += xsec
                #    if passedTrig:
                #        #hpTGenAcc[vtype].Fill(genEtaPt)
                #        hpTGenAcc[vtype].Fill(gen_eta.Pt())
                #        acc_weight[vtype] += xsec
        nvertMu = len(vtx_vmuP)
        bestjj = -1
        trig_filled = False
        for jj in range(nvertMu):
            #for 4-lepton vertices MUST use the muons corresponding to these electrons!
            if vstr in ["mmelel", "mmlplp", "mmee", "mumu"] and jj != j: continue
            if vstr == "mumu":
                try:
                    vxy = vtx_vxy[jj]
                except:
                    vxy = (vtx_vvx[jj]**2 + vtx_vvy[jj]**2)**0.5 
            else:
                #for all other vertex types, use the vxy including the electrons
                vxy = Vxy
            #charge stored as: 1 for positive; 255 for negative
            mp = ord(vtx_vmuP[jj]) 
            vec_muP = ROOT.TLorentzVector()
            #print("mp: " + str(mp) + "; nGoodMuon: " + str(ord(e.nGoodMuon)) + "; len(mp): " + str(len(vtx_vmuP)))
            #vec_muP.SetPtEtaPhiM(e.Muon_pt[mp], e.Muon_eta[mp], e.Muon_phi[mp], mu_mass)
            vec_muP.SetPtEtaPhiM(e.Muon_pt[mp]*(1+musf), e.Muon_eta[mp], e.Muon_phi[mp], mu_mass)
            mm = ord(vtx_vmuN[jj]) 
            # < 4 for requiring tight muID
            try:
                if tightMu:
                    failed_muID = require_muID and (ord(e.Muon_id[mp]) < 4 or ord(e.Muon_id[mm]) < 4)
                else:
                    failed_muID = require_muID and (ord(e.Muon_id[mp]) == 0 or ord(e.Muon_id[mm]) == 0)
            except TypeError:
                if tightMu:
                    failed_muID = require_muID and (e.Muon_id[mp] < 4 or e.Muon_id[mm] < 4)
                else:
                    failed_muID = require_muID and (e.Muon_id[mp] == 0 or e.Muon_id[mm] == 0)
            if failed_muID:
                continue
            #move these down further so can measure efficiency uncty!!
            #if mupTcut > -1 and (e.Muon_pt[mp] < mupTcut or e.Muon_pt[mm] < mupTcut):
            #    continue
            #if muleadpTcut > -1 and (e.Muon_pt[mp] < muleadpTcut and e.Muon_pt[mm] < muleadpTcut):
            #    continue
            if muEtacut > -1 and (abs(e.Muon_eta[mp]) > muEtacut or abs(e.Muon_eta[mm]) > muEtacut):
                continue
            #print("mm: " + str(mm) + "; nGoodMuon: " + str(ord(e.nGoodMuon)) + "; len(mm): " + str(len(vtx_vmuN))) 
            #make the 4-vector, fill the hists
            vec_muN = ROOT.TLorentzVector()
            vec_muN.SetPtEtaPhiM(e.Muon_pt[mm]*(1+musf), e.Muon_eta[mm], e.Muon_phi[mm], mu_mass)
            if isMC:
                #adjust the event weight (slightly) according to the scaled muon efficiency!
                Rmusf = (muEff(e.Muon_pt[mm]*(1+musf)) / muEff(e.Muon_pt[mm]) * muEff(e.Muon_pt[mp]*(1+musf))/muEff(e.Muon_pt[mp]))
                my_evt_weight = Rmusf*evt_weight 
                if vtype not in ["mumu", "mmg", "mmlplp"]:
                    try:
                        Relsf = (elEff(e.Electron_pt[elP]*(1+elsf)) / elEff(e.Electron_pt[elN]) * elEff(e.Electron_pt[elP]*(1+elsf))/elEff(e.Electron_pt[elP]))
                    except ZeroDivisionError:
                        Relsf = 0
                    my_evt_weight *= Relsf
                    if do_ElRecsf:
                        recoNameP = "RecoBelow20"
                        recoNameN = "RecoBelow20"
                        ptp = e.Electron_pt[elP]*(1+elsf)
                        if ptp > 20 and ptp < 75:
                            recoNameP = "Reco20to75"
                        elif ptp > 75:
                            recoNameP = "RecoAbove75"
                        ptn = e.Electron_pt[elN]*(1+elsf)
                        if ptn > 20 and ptn < 75:
                            recoNameN = "Reco20to75"
                        elif ptn > 75:
                            recoNameN = "RecoAbove75"
                        valsfPRec = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sf", recoNameP, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNRec = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sf", recoNameN, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        valsfPRecUp = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfup", recoNameP, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNRecUp = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfup", recoNameN, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        valsfPRecDn = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfdown", recoNameP, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNRecDn = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfdown", recoNameN, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        #valsfPRec = 0.98
                        #valsfPRecUp = 1.01
                        #valsfPRecDn = 0.95
                        #if (abs(e.Electron_eta[elP]) > 1.5 and abs(e.Electron_eta[elP]) < 2) or abs(e.Electron_eta[elP]) < 1:
                        #    valsfPRec = 0.99
                        #    valsfPRecUp = 1.02
                        #    valsfPRecDn = 0.96
                        #valsfNRec = 0.98
                        #valsfNRecUp = 1.01
                        #valsfNRecDn = 0.95
                        #if (abs(e.Electron_eta[elN]) > 1.5 and abs(e.Electron_eta[elN]) < 2) or abs(e.Electron_eta[elN]) < 1:
                        #    valsfNRec = 0.99
                        #    valsfNRecUp = 1.02
                        #    valsfNRecDn = 0.96
                        my_evt_weight *= (valsfPRec*valsfNRec)
                        evt_weightRecUp = my_evt_weight * (valsfPRecUp/valsfPRec) * (valsfNRecUp/valsfNRec)
                        evt_weightRecDn = my_evt_weight * (valsfPRecDn/valsfPRec) * (valsfNRecDn/valsfNRec)
                        if recoNameP == "RecoBelow20":
                            #print("valsfPRec: " + str(valsfPRec)) 
                            hElRecsf[vtype].SetBinContent(hElRecsf[vtype].FindBin(e.Electron_eta[elP]), valsfPRec)
                            hElRecsfUp[vtype].SetBinContent(hElRecsfUp[vtype].FindBin(e.Electron_eta[elP]), valsfPRecUp)
                            hElRecsfDn[vtype].SetBinContent(hElRecsfDn[vtype].FindBin(e.Electron_eta[elP]), valsfPRecDn)
                        if recoNameN == "RecoBelow20":
                            hElRecsf[vtype].SetBinContent(hElRecsf[vtype].FindBin(e.Electron_eta[elN]), valsfNRec)
                            hElRecsfUp[vtype].SetBinContent(hElRecsfUp[vtype].FindBin(e.Electron_eta[elN]), valsfNRecUp)
                            hElRecsfDn[vtype].SetBinContent(hElRecsfDn[vtype].FindBin(e.Electron_eta[elN]), valsfNRecDn)
                    if do_ElIDsf:
                        wptxt = "wp90noiso"
                        if tightEl:
                            wptxt = "wp80noiso"
                        valsfPID = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sf", wptxt, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNID = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sf", wptxt, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        valsfPIDUp = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfup", wptxt, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNIDUp = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfup", wptxt, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        valsfPIDDn = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfdown", wptxt, e.Electron_eta[elP], max(e.Electron_pt[elP]*(1+elsf), 10.0))
                        valsfNIDDn = evaluator["Electron-ID-SF"].evaluate("2022Re-recoE+PromptFG", "sfdown", wptxt, e.Electron_eta[elN], max(e.Electron_pt[elN]*(1+elsf), 10.0))
                        ##I double the uncertainties shown in the Run2 paper to be extra conservative!!
                        #if abs(e.Electron_eta[elP]) < 0.80:
                        #    valsfPID = 0.98
                        #    valsfPIDUp = 1.02
                        #    valsfPIDDn = 0.94
                        #elif abs(e.Electron_eta[elP]) < 1.44:
                        #    valsfPID = 0.99
                        #    valsfPIDUp = 1.01
                        #    valsfPIDDn = 0.97
                        #elif abs(e.Electron_eta[elP]) > 1.57 and abs(e.Electron_eta[elP]) < 2.00:
                        #    valsfPID = 0.96
                        #    valsfPIDUp = 1.02
                        #    valsfPIDDn = 0.90
                        #elif abs(e.Electron_eta[elP]) > 2.00 and abs(e.Electron_eta[elP]) < 2.50:
                        #    valsfPID = 1.00
                        #    valsfPIDUp = 1.04
                        #    valsfPIDDn = 0.96
                        #else:
                        #    print("Warning: Run %d evt %d electron eta value %f, not correcting."%(e.run, e.evt, e.Electron_eta[elP])) 
                        #    valsfPID = 1.00
                        #    valsfPIDUp = 1.00
                        #    valsfPIDDn = 1.00
                        #if abs(e.Electron_eta[elN]) < 0.80:
                        #    valsfNID = 0.98
                        #    valsfNIDUp = 1.02
                        #    valsfNIDDn = 0.94
                        #elif abs(e.Electron_eta[elN]) < 1.44:
                        #    valsfNID = 0.99
                        #    valsfNIDUp = 1.01
                        #    valsfNIDDn = 0.97
                        #elif abs(e.Electron_eta[elN]) > 1.57 and abs(e.Electron_eta[elN]) < 2.00:
                        #    valsfNID = 0.96
                        #    valsfNIDUp = 1.02
                        #    valsfNIDDn = 0.90
                        #elif abs(e.Electron_eta[elN]) > 2.00 and abs(e.Electron_eta[elN]) < 2.50:
                        #    valsfNID = 1.00
                        #    valsfNIDUp = 1.04
                        #    valsfNIDDn = 0.96
                        #else:
                        #    print("Warning: Run %d evt %d electron eta value %f, not correcting."%(e.run, e.evt, e.Electron_eta[elN])) 
                        #    valsfNID = 1.00
                        #    valsfNIDUp = 1.00
                        #    valsfNIDDn = 1.00
                        my_evt_weight *= (valsfPID*valsfNID)
                        evt_weightRecUp *= (valsfPID*valsfNID)
                        evt_weightRecDn *= (valsfPID*valsfNID)
                        evt_weightIDUp = my_evt_weight * (valsfPIDUp/valsfPID) * (valsfNIDUp/valsfNID)
                        evt_weightIDDn = my_evt_weight * (valsfPIDDn/valsfPID) * (valsfNIDDn/valsfNID)
                        if e.Electron_pt[elP]*(1+elsf) < 10:
                            hElIDsf[vtype].SetBinContent(hElIDsf[vtype].FindBin(e.Electron_eta[elP]), valsfPID)
                            hElIDsfUp[vtype].SetBinContent(hElIDsfUp[vtype].FindBin(e.Electron_eta[elP]), valsfPIDUp)
                            hElIDsfDn[vtype].SetBinContent(hElIDsfDn[vtype].FindBin(e.Electron_eta[elP]), valsfPIDDn)
                        if e.Electron_pt[elN]*(1+elsf) < 10:
                            hElIDsf[vtype].SetBinContent(hElIDsf[vtype].FindBin(e.Electron_eta[elN]), valsfNID)
                            hElIDsfUp[vtype].SetBinContent(hElIDsfUp[vtype].FindBin(e.Electron_eta[elN]), valsfNIDUp)
                            hElIDsfDn[vtype].SetBinContent(hElIDsfDn[vtype].FindBin(e.Electron_eta[elN]), valsfNIDDn)
                        hEletapt[vtype].Fill(e.Electron_eta[elP], e.Electron_pt[elP]*(1+elsf), my_evt_weight)
                        hEletapt[vtype].Fill(e.Electron_eta[elN], e.Electron_pt[elN]*(1+elsf), my_evt_weight)
                    
                if do_MuIDsf:
                    valsfPMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mp], max(e.Muon_pt[mp]*(1+musf), 3.0), "nominal")
                    valsfNMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mm], max(e.Muon_pt[mm]*(1+musf), 3.0), "nominal")
                    valsfPUpMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mp], max(e.Muon_pt[mp]*(1+musf), 3.0), "systup")
                    valsfNUpMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mm], max(e.Muon_pt[mm]*(1+musf), 3.0), "systup")
                    valsfPDnMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mp], max(e.Muon_pt[mp]*(1+musf), 3.0), "systdown")
                    valsfNDnMu = mvaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(e.Muon_eta[mm], max(e.Muon_pt[mm]*(1+musf), 3.0), "systdown")
                    my_evt_weight *= (valsfPMu*valsfNMu) 
                    evt_weightMuUp = my_evt_weight * (valsfPUpMu/valsfPMu) * (valsfNUpMu/valsfNMu)
                    evt_weightMuDn = my_evt_weight * (valsfPDnMu/valsfPMu) * (valsfNDnMu/valsfNMu)
                    #print("valsfPMu: %f, valsfNMu: %f"%(valsfPMu, valsfNMu)) 
                if do_pileup:
                    #pileup correction
                    puCor = PUCor.GetBinContent(PUCor.FindBin(e.nPV))
                    my_evt_weight *= puCor
            else:
                my_evt_weight = 1
            if vtype == "elel":
                hMReal.Fill((diel+dimu).M(), my_evt_weight)
            if vtype == "mmg":
                #find the best photon
                bestg = -1
                bestM = 99999
                vec_gams = [ ROOT.TLorentzVector() for gg in range(ord(e.nGoodPhoton)) ]
                for gg in range(ord(e.nGoodPhoton)): 
                    vec_gams[gg].SetPtEtaPhiM(e.Photon_pt[gg], e.Photon_eta[gg], e.Photon_phi[gg], 0)
                    meta = (vec_muP + vec_muN + vec_gams[gg]).M()
                    dif = abs(etamass - meta)
                    if dif < abs(bestM - etamass):
                        bestM = meta
                        bestg = gg
                if bestg == -1: continue
                vec_eta = vec_muP + vec_muN + vec_gams[bestg]
            elif vtype == "mumu":
                vec_eta = vec_muP + vec_muN
                #for mumu vertices only, also fill the hMSSe and hMOSe histograms (for comb. bkg test)
                for aa in range(ord(e.nGoodElectron)):
                    for bb in range(aa+1, ord(e.nGoodElectron)):
                        #form the 4-lepton invariant mass
                        vec_aa = ROOT.TLorentzVector()
                        vec_bb = ROOT.TLorentzVector()
                        vec_aa.SetPtEtaPhiM(e.Electron_pt[aa]*(1+elsf), e.Electron_eta[aa], e.Electron_phi[aa], el_mass)
                        vec_bb.SetPtEtaPhiM(e.Electron_pt[bb]*(1+elsf), e.Electron_eta[bb], e.Electron_phi[bb], el_mass)
                        vec_rand = vec_muP + vec_muN + vec_aa + vec_bb
                        if e.Electron_charge[aa] == e.Electron_charge[bb]:
                            hMSSe.Fill(vec_rand.M(), my_evt_weight)
                        else:
                            hMOSe.Fill(vec_rand.M(), my_evt_weight)
            else:
                vec_eta = vec_muP + vec_muN + vec_elP + vec_elN
                mee = (vec_elP + vec_elN).M()
                ##Cut on Mee to get rid of most resonant bkg!
                if cut_mee:
                    if mee > 0.04 and mee < 0.09:
                    #    #print("continuing!!!! evt=%d, mee=%f"%(e.evt, mee))
                        continue
            pt = vec_eta.Pt() 
            if vec_muP.Pt() < vec_muN.Pt():
                subpt = e.Muon_pt[mp] #vec_muP.Pt()
                subeta = vec_muP.PseudoRapidity()
            else:
                subpt = e.Muon_pt[mm] #vec_muN.Pt()
                subeta = vec_muN.PseudoRapidity()
            m = vec_eta.M()
            mudR = vec_muP.DeltaR(vec_muN)

            #make sure you're only filling if this is the true 2nd highest pT muon!
            himupt = max(e.Muon_pt)
            #find the second highest muon pt
            hi2mupt = 0.0
            for mpt in e.Muon_pt:
                if mpt < himupt and mpt > hi2mupt:
                    hi2mupt = mpt
            #be sure to avoid rounding errors...
            if vtype == "mumu" and mudR < 0.1 and not trig_filled and abs(hi2mupt-subpt) < 0.01:
                trig_filled = True
                if int(str(testnum)[:2]) > 38: 
                    trig = e.Triggers_fired1
                else:
                    trig = e.Triggers_fired0
                l1pass = (trig & (1<<25)) != 0
                hlpass = (trig & (1<<11) != 0)
                if l1pass:
                    if abs(subeta) < 1.2:
                        hL1TC.Fill(subpt, my_evt_weight)
                    else:
                        hL1TF.Fill(subpt, my_evt_weight)
                if hlpass:
                    if subpt < 2.5:
                        print("aaaa werid event!! subpt: %f, himupt: %f"%(subpt, himupt)) 
                        printEvent.printEvent(e, isMC, g if isMC else None)
                    if abs(subeta) < 1.2:
                        hHLTC.Fill(subpt, my_evt_weight)
                    else:
                        hHLTF.Fill(subpt, my_evt_weight)
            
            #L1 PART of these SFs from https://indico.cern.ch/event/1398079/contributions/5876242/attachments/2826267/4937141/DoubleMu4_3_LowMass_Trigger_SFs_v1.pdf
            #HLT part derived ourselves!
            if do_TTOC:   
                if vec_muP.Pt() > vec_muN.Pt() :
                    etaProbe = vec_muN.PseudoRapidity()
                else:
                    etaProbe = vec_muP.PseudoRapidity()
                #L1 part, then HLT part
                if abs(etaProbe) < 0.9:
                    if subpt < 6.5:
                        my_evt_weight *= 0.940464
                        #my_evt_weight *= 0.636192
                        ttocL1 = 0.00513163
                        #ttocHLT = 0.0824157 
                    elif subpt < 10:
                        my_evt_weight *= 0.995195
                        #my_evt_weight *= 0.595192
                        ttocL1 = 0.0175254
                        #ttocHLT = 0.0420128 
                    elif subpt < 15:
                        my_evt_weight *= 0.936098
                        #my_evt_weight *= 0.934145
                        ttocL1 = 0.0147551 
                        #ttocHLT = 0.042563  
                    elif subpt < 20:
                        my_evt_weight *= 0.869899
                        #my_evt_weight *= 0.963782
                        ttocL1 = 0.00579038
                        #ttocHLT = 0.00625994
                    else:
                        my_evt_weight *= 0.844248
                        #my_evt_weight *= 0.953331
                        ttocL1 = 0.00607634
                        #ttocHLT = 0.00583198
                elif abs(etaProbe) < 1.2:
                    #0.9 < probe eta < 1.2
                    if subpt < 6.5:
                        my_evt_weight *= 0.974163
                        #my_evt_weight *= 0.932408
                        ttocL1 = 0.112774  
                        #ttocHLT = 0.28918   
                    elif subpt < 10:
                        my_evt_weight *= 0.983871
                        #my_evt_weight *= 0.732941
                        ttocL1 = 0.0106754 
                        #ttocHLT = 0.0907367 
                    elif subpt < 15:
                        my_evt_weight *= 0.946724
                        #my_evt_weight *= 0.929410
                        ttocL1 = 0.0814474 
                        #ttocHLT = 0.0475792 
                    elif subpt < 20:
                        my_evt_weight *= 0.854219
                        #my_evt_weight *= 0.936802
                        ttocL1 = 0.0135183 
                        #ttocHLT = 0.0109483 
                    else:
                        my_evt_weight *= 0.725096
                        #my_evt_weight *= 0.957847
                        ttocL1 = 0.0606656 
                        #ttocHLT = 0.0500664 
                else:
                    #probe eta > 1.2
                    if subpt < 6.5:
                        my_evt_weight *= 0.890303
                        #my_evt_weight *= 0.730769
                        ttocL1 = 0.0562114 
                        #ttocHLT = 0.205952  
                    elif subpt < 10:
                        my_evt_weight *= 0.968846
                        #my_evt_weight *= 0.70067
                        ttocL1 = 0.0553319 
                        #ttocHLT = 0.0825737 
                    elif subpt < 15:
                        my_evt_weight *= 0.829491
                        #my_evt_weight *= 0.905507
                        ttocL1 = 0.058834  
                        #ttocHLT = 0.0222346 
                    elif subpt < 20:
                        my_evt_weight *= 0.785405
                        #my_evt_weight *= 0.934546
                        ttocL1 = 0.0220265 
                        #ttocHLT = 0.031868  
                    else:
                        my_evt_weight *= 0.760657
                        #my_evt_weight *= 0.974615
                        ttocL1 = 0.0133128 
                        #ttocHLT = 0.051794  
                if abs(etaProbe) < 1.2:
                    tcbin = triggerCorrectionC.FindBin(min(subpt, triggerCorrectionC.GetXaxis().GetXmax()))
                    my_evt_weight *= triggerCorrectionC.GetBinContent(tcbin)
                    ttocHLT = triggerCorrectionC.GetBinError(tcbin)
                else:
                    tcbin = triggerCorrectionF.FindBin(min(subpt, triggerCorrectionF.GetXaxis().GetXmax()))
                    my_evt_weight *= triggerCorrectionF.GetBinContent(tcbin)
                    ttocHLT = triggerCorrectionF.GetBinError(tcbin)
                evt_weightTTOCUp = my_evt_weight * ( 1 + ttocHLT)
                evt_weightTTOCUp *= ( 1 + ttocL1)
                evt_weightTTOCDn = my_evt_weight * ( 1 - ttocHLT)
                evt_weightTTOCDn *= ( 1 - ttocL1)

            if vtype == "mmelel":
                hMeeVsMmmee.Fill(mee, m, my_evt_weight) 

            if danEvent:
                print(vtype + " vtx " + str(jj) + " mass: %f"%m) 

            if isMC:
                gm = False
                #if isSig:
                dr = vec_eta.DeltaR(gen_eta)
                #else:
                #    #for bkg MC, need to see if the eta meson is reco'd with mu-mu-gamma
                #    dr = 9999
                #    for kk in range(ord(e.nGoodPhoton)):
                #        vec_phot = ROOT.TLorentzVector()
                #        vec_phot.SetPtEtaPhiM(e.Photon_pt[kk], e.Photon_eta[kk], e.Photon_phi[kk], 0)
                #        vec_bkgEta = vec_muP + vec_muN + vec_phot
                #        dr = min(dr, vec_bkgEta.DeltaR(gen_eta)) 
                #TESTING TO SEE WHAT DIFFERENCE THIS MAKES (SHOULD BE NEGLIGIBLE)
                #if dr < dRcut:
                #    gm = True
                if m > peakmin and m < .58:
                    gm = True
                if gm: 
                    #hpTGenReco[vtype].Fill(genEtaPt)
                    #hpTGenReco[vtype].Fill(gen_eta.Pt())
                    #hEtaGenReco[vtype].Fill(gen_eta.PseudoRapidity())
                    #rec_weight[vtype] += evt_weight
                    if passedTrig and not acc_filled[1]:
                        #hpTGenAcc[vtype].Fill(genEtaPt)
                        hpTGenAcc[vtype].Fill(gen_eta.Pt(), evt_weight)
                        if gen_eta.Pt() > 25:
                            hEtaGenAcc[vtype].Fill(gen_eta.PseudoRapidity(), my_evt_weight)
                        acc_weight[vtype][1] += my_evt_weight
                        acc_filled[1] = True
                    #now see if it also passes the stricter cuts
                    if passedTrig and vtype not in ["mmg", "mumu", "mmlplp"]:
                        if e.Electron_pt[elP]*(1+elsf) > elpTcut and e.Electron_pt[elN]*(1+elsf) > elpTcut and \
                            ((elEtacut < 0) or (abs(e.Electron_eta[elP]) < elEtacut and abs(e.Electron_eta[elN]) < elEtacut)) and \
                            e.Muon_pt[mp]*(1+musf) > mupTcut and e.Muon_pt[mm]*(1+musf) > mupTcut and \
                            ((muEtacut < -1) or (abs(e.Muon_eta[mp]) < muEtacut and abs(e.Muon_eta[mm]) < muEtacut)) and \
                            (e.Muon_pt[mp]*(1+musf) > muleadpTcut or e.Muon_pt[mm]*(1+musf) > muleadpTcut) :
                            if (not acc_filled[2]):
                                acc_weight[vtype][2] += my_evt_weight
                                acc_filled[2] = True
                            if (not acc_filled[3]) and (WP90ID_p > 0 and WP90ID_n > 0):
                                acc_weight[vtype][3] += my_evt_weight
                                acc_filled[3] = True
                                if (not acc_filled[7]) and (mee < 0.04 or mee > 0.09):
                                    acc_weight[vtype][7] += my_evt_weight
                                    acc_filled[7] = True
                            if (not acc_filled[4]) and (WP90ID_p > 0 or WP90ID_n > 0):
                                acc_weight[vtype][4] += my_evt_weight
                                acc_filled[4] = True
                                if (not acc_filled[8]) and (mee < 0.04 or mee > 0.09):
                                    acc_weight[vtype][8] += my_evt_weight
                                    acc_filled[8] = True
                            if (not acc_filled[5]) and (WP80ID_p > 0 and WP80ID_n > 0):
                                acc_weight[vtype][5] += my_evt_weight
                                acc_filled[5] = True
                                if (not acc_filled[9]) and (mee < 0.04 or mee > 0.09):
                                    acc_weight[vtype][9] += my_evt_weight
                                    acc_filled[9] = True
                            if (not acc_filled[6]) and (mee < 0.04 or mee > 0.09):
                                acc_weight[vtype][6] += my_evt_weight
                                acc_filled[6] = True
                #if vtype == "mmlplp":
                #    print("vtype=%s recod=%s passedTrig=%s accdR_filled=%s"%(vtype, recod, passedTrig, accdR_filled))
                if recod and passedTrig and not accdR_filled:
                    #make sure it passed the rest of the event selections
                    if e.Muon_pt[mp]*(1+musf) > mupTcut and e.Muon_pt[mm]*(1+musf) > mupTcut and \
                        (muEtacut < 0 or (abs(e.Muon_eta[mp]) < muEtacut and abs(e.Muon_eta[mm]) < muEtacut)) and \
                        (e.Muon_pt[mp]*(1+musf) > muleadpTcut or e.Muon_pt[mm]*(1+musf) > muleadpTcut) and \
                        (not failed_muID):

                        if vtype in ["mumu", "mmg"] or ("elel" in vtype and e.Electron_pt[elP]*(1+elsf) > elpTcut and e.Electron_pt[elN]*(1+elsf) > elpTcut and \
                            (elEtacut < 0 or (abs(e.Electron_eta[elP]) < elEtacut and abs(e.Electron_eta[elN]) < elEtacut)) and \
                            ((not require_elID) or (WP90ID_p > 0 and WP90ID_n > 0))) or "lplp" in vtype :

                            #fill hpTGenAccdR with my_evt_weight to properly correct the acceptance for all the corrections!
                            hpTGenAccdR[vtype].Fill(gen_eta.Pt(), my_evt_weight) 
                            accdR_weight[vtype] += my_evt_weight
                            accdR_filled = True
            #danEvent = False

            #if syncTest and vtype == "mumu" and m < .58 and m > .52:
            ##if syncTest and vtype == "mumu" and e.evt in dan_events:
            #    danEvent = True
            global nPrint
            if isMC and vtype == "mmelel" and m > peakmin and m < .58 and nPrint > 0:
                #global nPrint
                try:
                    printEvent.printEvent(e, isMC, g if isMC else None)
                except AttributeError:
                    printEvent_backup.printEvent(e, isMC, g if isMC else None)
                print(vtype + " vtx " + str(jj) + " mass: %f"%m + " evtWeight: %f"%my_evt_weight) 
                nPrint -= 1

            if not fillAtEnd:
                #do the nominal muons pass the pt cuts?
                mu_fail = (mupTcut > -1 and (vec_muP.Pt() < mupTcut or vec_muN.Pt() < mupTcut)) or (muleadpTcut > -1 and vec_muP.Pt() < muleadpTcut and vec_muN.Pt() < muleadpTcut)
                #do the nominal electrons pass the pt cuts?
                el_fail = ("elel" in vtype and elpTcut > -1 and (e.Electron_pt[elP]*(1+elsf) < elpTcut or e.Electron_pt[elN]*(1+elsf) < elpTcut)) 
                if vtype == "mumu" and m > 0.45 and m < 0.65:
                    #see if can fill hpTGoodL1Pass, hpTGoodHLTPass
                    #if vec_muP.Pt() > 5 and vec_muN.Pt() > 5 and (e.Triggers_fired0 & (1 << 25)):
                    thistrig = e.Triggers_fired0 if int(str(testnum)[:2]) <= 38 else e.Triggers_fired1
                    if thistrig & (1 << 25):
                        musubpt = min(vec_muP.Pt(), vec_muN.Pt())
                        hpTGoodL1Pass.Fill(musubpt, my_evt_weight ) 
                        #if (e.Triggers_fired0 & (1 << 11)) :
                        if thistrig & (1 << 11) :
                            hpTGoodHLTPass.Fill(musubpt, my_evt_weight)
                if vstr != "mumu":
                    elptP = vec_elP.Pt()
                    elptN = vec_elN.Pt()
                muptP = vec_muP.Pt()
                muptN = vec_muN.Pt()
                #for electron efficiency uncty calculations
                if isMC:
                    if vtype == "mmelel" and not mu_fail:
                        for s in range(nMod):
                            R = elEff(elptP*(1.0 + sf[s]/100.0)) / elEff(elptP) * elEff(elptN*(1.0 + sf[s]/100.0)) / elEff(elptN)
                            evt_weightMod = my_evt_weight * R
                            vec_elPmod = ROOT.TLorentzVector()
                            vec_elPmod.SetPtEtaPhiM(elptP*(1.0 + sf[s]/100.0), vec_elP.PseudoRapidity(), vec_elP.Phi(), el_mass)
                            vec_elNmod = ROOT.TLorentzVector()
                            vec_elNmod.SetPtEtaPhiM(elptN*(1.0 + sf[s]/100.0), vec_elN.PseudoRapidity(), vec_elN.Phi(), el_mass)
                            if not (elpTcut > -1 and (vec_elPmod.Pt() < elpTcut or vec_elNmod.Pt() < elpTcut)):
                                mmod = (vec_elPmod + vec_elNmod + vec_muP + vec_muN).M()
                                hMMod[s].Fill(mmod, evt_weightMod)
                    if ((isMuMu and vtype == "mumu") or (isMC and vtype == "mmelel" and not isMuMu)) and not el_fail:
                        for s in range(nMod):
                            Rmu = muEff(muptP*(1.0 + sfmu[s]/100.0)) / muEff(muptP) * muEff(muptN*(1.0 + sfmu[s]/100.0)) / muEff(muptN)
                            evt_weightModmu = my_evt_weight * Rmu
                            vec_muPmod = ROOT.TLorentzVector()
                            vec_muPmod.SetPtEtaPhiM(muptP*(1.0 + sfmu[s]/100.0), vec_muP.PseudoRapidity(), vec_muP.Phi(), mu_mass)
                            vec_muNmod = ROOT.TLorentzVector()
                            vec_muNmod.SetPtEtaPhiM(muptN*(1.0 + sfmu[s]/100.0), vec_muN.PseudoRapidity(), vec_muN.Phi(), mu_mass)
                            if not ((mupTcut > -1 and (vec_muPmod.Pt() < mupTcut or vec_muNmod.Pt() < mupTcut)) or (muleadpTcut > -1 and vec_muPmod.Pt() < muleadpTcut and vec_muNmod.Pt() < muleadpTcut)):
                                if "elel" in vtype:
                                    etamodmu = vec_muPmod + vec_muNmod + vec_elP + vec_elN
                                    #print("muptP,muptN=%f,%f, Rmu=%f"%(muptP, muptN, Rmu)) 
                                    #print("mmodmu: %f; regular m: %f"%(mmodmu, m)) 
                                else:
                                    etamodmu = vec_muPmod + vec_muNmod
                                mmodmu = etamodmu.M()
                                ptmodmu = etamodmu.Pt()
                                hMModmu[s].Fill(mmodmu, evt_weightModmu)
                                if mmodmu > peakmin and mmodmu < peakmax:
                                    hpTModmu[s].Fill(ptmodmu, evt_weightModmu)
                if not (el_fail or mu_fail):
                    if my_evt_weight != 1 and not isMC:
                        print("*****ERROR: SOMEONE IS F***ING WITH EVENT WEIGHTS IN DATA*****")
                        sys.exit()
                    hpT[vtype].Fill(pt, my_evt_weight)
                    if m > peakmin and m < peakmax:
                        hpTSig[vtype].Fill(pt, my_evt_weight)
                    #elif (m > (peakmin-2*.03) and m < (peakmin-.03)) or (m > (peakmax+.03) and m < (peakmax+2*.03)):
                    #sidebands are strictly greater than peak now!
                    elif m > (peakmax+.03) and m < 0.80:
                        hpTSide[vtype].Fill(pt, my_evt_weight)
                    if vstr != "mumu":
                        hpTEl[vtype].Fill(elptP*(1+elsf), my_evt_weight)
                        hpTEl[vtype].Fill(elptN*(1+elsf), my_evt_weight)
                        hpTElNoWt[vtype].Fill(elptP*(1+elsf))
                        hpTElNoWt[vtype].Fill(elptN*(1+elsf))
                        eldR = vec_elP.DeltaR(vec_elN)
                        hdREl[vtype].Fill(eldR, my_evt_weight)
                        hdRvPtEl[vtype].Fill(eldR, elptP*(1+elsf), my_evt_weight)
                        hdRvPtEl[vtype].Fill(eldR, elptN*(1+elsf), my_evt_weight)
                        if eldR > 0.005:
                            hpTElAll3.Fill(elptP*(1+elsf), my_evt_weight)
                            hpTElAll3.Fill(elptN*(1+elsf), my_evt_weight)
                            if WP90ID_p != 0:
                                hpTElIDPass3.Fill(elptP*(1+elsf), my_evt_weight)
                            if WP90ID_n != 0:
                                hpTElIDPass3.Fill(elptN*(1+elsf), my_evt_weight)
                        elif eldR > 0.0025:
                            hpTElAll2.Fill(elptP*(1+elsf), my_evt_weight)
                            hpTElAll2.Fill(elptN*(1+elsf), my_evt_weight)
                            if WP90ID_p != 0:
                                hpTElIDPass2.Fill(elptP*(1+elsf), my_evt_weight)
                            if WP90ID_n != 0:
                                hpTElIDPass2.Fill(elptN*(1+elsf), my_evt_weight)
                        elif eldR > 0.0005:
                            hpTElAll1.Fill(elptP*(1+elsf), my_evt_weight)
                            hpTElAll1.Fill(elptN*(1+elsf), my_evt_weight)
                            if WP90ID_p != 0:
                                hpTElIDPass1.Fill(elptP*(1+elsf), my_evt_weight)
                            if WP90ID_n != 0:
                                hpTElIDPass1.Fill(elptN*(1+elsf), my_evt_weight)
                        else:
                            hpTElAll0.Fill(elptP*(1+elsf), my_evt_weight)
                            hpTElAll0.Fill(elptN*(1+elsf), my_evt_weight)
                            if WP90ID_p != 0:
                                hpTElIDPass0.Fill(elptP*(1+elsf), my_evt_weight)
                            if WP90ID_n != 0:
                                hpTElIDPass0.Fill(elptN*(1+elsf), my_evt_weight)
                    hpTMu[vtype].Fill(muptP*(1+musf), my_evt_weight)
                    hpTMu[vtype].Fill(muptN*(1+musf), my_evt_weight)
                    if abs(vec_muP.PseudoRapidity()) > 1.5:
                        hpTMuF[vtype].Fill(muptP*(1+musf), my_evt_weight)
                    if abs(vec_muN.PseudoRapidity()) > 1.5:
                        hpTMuF[vtype].Fill(muptN*(1+musf), my_evt_weight)
                    hpTMuNoWt[vtype].Fill(muptP*(1+musf))
                    hpTMuNoWt[vtype].Fill(muptN*(1+musf))
                    hEtaMu[vtype].Fill(vec_muP.PseudoRapidity())
                    hEtaMu[vtype].Fill(vec_muN.PseudoRapidity())
                    hEta[vtype].Fill(vec_eta.PseudoRapidity(), my_evt_weight)
                    if isMC:
                        #fill event weight histograms only for ACCEPTED events (correct mass reco'd)
                        if m > peakmin and m < peakmax:
                            hEventWeight[vtype].Fill(my_evt_weight)
                            hEvtWtVsPt[vtype].Fill(gen_eta.Pt(), my_evt_weight)
                    ##test: see what happens if fill hist ONLY for high-vxy bois
                    #if vxy > 1.2: continue
                    ###### TEST #####
                    hM[vtype].Fill(m, my_evt_weight)
                    if m > peakmin and m < peakmax:
                        ncand += 1
                        try:
                            for l1b in range(8):
                                if ord(e.L1_fired) & (1<<l1b) != 0:
                                    hL1pass[vtype].Fill(l1b+0.5)
                            hL1all[vtype].Fill(0.5)
                        except AttributeError:
                            pass
                    if isMC and vtype not in ["mumu", "mmg", "mmlplp"]: 
                        if do_ElRecsf:
                            hMRecUp[vtype].Fill(m, evt_weightRecUp)
                            hMRecDn[vtype].Fill(m, evt_weightRecDn)
                        if do_ElIDsf:
                            hMIDUp[vtype].Fill(m, evt_weightIDUp)
                            hMIDDn[vtype].Fill(m, evt_weightIDDn)
                    if isMC and do_TTOC:
                        hMTTOCUp[vtype].Fill(m, evt_weightTTOCUp)
                        hMTTOCDn[vtype].Fill(m, evt_weightTTOCDn)
                    if isMC and do_MuIDsf:
                        hMMuUp[vtype].Fill(m, evt_weightMuUp)
                        hMMuDn[vtype].Fill(m, evt_weightMuDn)
                    if not passedTrig:
                        hMFailedTrig[vtype].Fill(m, my_evt_weight) 
                    hMNoWt[vtype].Fill(m)
                    hnPVgood[vtype].Fill(e.nPV, my_evt_weight)
                    if not (year == 2023 and vtype == "mumu"):
                        hMvsPt[vtype].Fill(pt, m, my_evt_weight) 
                    hsubVdR[vtype].Fill(subpt, mudR)
                    hVxy[vtype].Fill(vxy)
                    hRchi2[vtype].Fill(vtx_vrechi2[j])
                    if vtype == "mmelel" and m > peakmin and m < peakmax:
                        hMNoEl.Fill( (vec_muP+vec_muN).M(), my_evt_weight)
                        mee = (vec_elP+vec_elN).M()
                        hMNoMu.Fill(mee, my_evt_weight)
                        #hMNoMuPiM.Fill( (vec_piP+vec_piN).M(), my_evt_weight)
                    if vtype == "mmelel" and m > (peakmin-2*.03) and m < (peakmin-.03) :
                        hMNoElLSide.Fill( (vec_muP+vec_muN).M(), my_evt_weight)
                        hMNoMuLSide.Fill( (vec_elP+vec_elN).M(), my_evt_weight)
                    if vtype == "mmelel" and m > (peakmax+.03) and m < (peakmax+2*.03):
                        hMNoElRSide.Fill( (vec_muP+vec_muN).M(), my_evt_weight)
                        hMNoMuRSide.Fill( (vec_elP+vec_elN).M(), my_evt_weight)
                    if vtype == "mmelel" and ((m > (peakmax+.03) and m < (peakmax+2*.03)) or (m > (peakmin-2*.03) and m < (peakmin-.03))):
                        hMNoElSide.Fill( (vec_muP+vec_muN).M(), my_evt_weight)
                        hMNoMuSide.Fill( (vec_elP+vec_elN).M(), my_evt_weight)
                        
                    #if vtype == "elel":
                    if isMC and (((isSig or not isMuMu) and vtype == "mmelel") or (isMuMu and vtype == "mumu")):
                        if not isMuMu:
                            mee = (vec_elP+vec_elN).M()
                            hMee.Fill(mee, my_evt_weight)
                            if m > peakmin and m < peakmax:
                                hMeePeak.Fill(mee, my_evt_weight)
                        if isMC: # and not isMuMu:
                            #min dR for pos, neg electrons
                            mindrpE = 9999.9
                            mindrnE = 9999.9
                            #min dR for pos, neg muons
                            mindrpM = 9999.9
                            mindrnM = 9999.9
                            #see if can genmatch-- loop thru all gen electrons/photons AND muons
                            for ge in range(g.nGenPart):
                                gv = ROOT.TLorentzVector()
                                gv.SetPtEtaPhiM(g.GenPart_pt[ge], g.GenPart_eta[ge], g.GenPart_phi[ge], g.GenPart_mass[ge]) 
                                #if (isSig and abs(g.GenPart_pdgId[ge]) == 11) or ((not isSig) and g.GenPart_pdgId[ge] == 22):
                                if abs(g.GenPart_pdgId[ge]) in [11, 22]:
                                    drpE = vec_elP.DeltaR(gv)
                                    drnE = vec_elN.DeltaR(gv)
                                    #valid dR if it's the same charge electron, or if it's a photon (in case of resBkg)
                                    if drpE < mindrpE and ((not isSig) or g.GenPart_charge == 1):
                                        mindrpE = drpE
                                    if drnE < mindrnE and ((not isSig) or g.GenPart_charge != 1):
                                        mindrnE = drnE
                                elif abs(g.GenPart_pdgId[ge]) == 13:
                                    drpM = vec_muP.DeltaR(gv)
                                    drnM = vec_muN.DeltaR(gv)
                                    if g.GenPart_charge == 1 and drpM < mindrpM:
                                        mindrpM = drpM
                                    if g.GenPart_charge != 1 and drnM < mindrnM:
                                        mindrnM = drnM
                            if vtype == "mmelel":
                                hMeedR.Fill(mindrpE, my_evt_weight)
                                hMeedR.Fill(mindrnE, my_evt_weight)
                                if mindrpE < dRcut and mindrnE < dRcut:
                                    hMeeGM.Fill(mee, my_evt_weight) 
                            hMmmdR.Fill(mindrpM, my_evt_weight)
                            hMmmdR.Fill(mindrnM, my_evt_weight)
                #if syncTest and vtype == "mmelel" and m < 1.0:
                #if danEvent and m > .52 and m < .58:
                if syncTest and vtype == "mumu" and m > peakmin and m < peakmax:
                #if syncTest and vtype == "mumu" and e.evt in dan_events:
                    syncFile.write("%d %f\n"%(e.evt, m)) 
                    #global nPrint
                    if nPrint > 0:
                        try:
                            printEvent.printEvent(e, isMC, g if isMC else None)
                        except AttributeError:
                            printEvent_backup.printEvent(e, isMC, g if isMC else None)
                        print("mass: %f"%m) 
                        nPrint -= 1
            #good eta if in the right mass range
            if m > peakmin and m < peakmax:
                good = True
            if not singleVert:
                Vxy = vxy
                #hMNoWt[vtype].Fill(m)
                #hMvsPt[vtype].Fill(pt, m, evt_weight) 
                #if vxy > 1.2:
                #    hMhiVxy[vtype].Fill(m, evt_weight)
                #else:
                #    hMloVxy[vtype].Fill(m, evt_weight)
            else:
                if abs(m-etamass) < abs(bestm-etamass):
                    bestm = m
                    if vtype == "mmelel":
                        bestm2mu = (vec_muP+vec_muN).M()
                        bestm2el = (vec_elP+vec_elN).M()
                    bestjj = jj
                    bestpt = pt
                    besteta = vec_eta.PseudoRapidity()
                    Vxy = vxy
                    if vstr != "mumu":
                        elptP = vec_elP.Pt()
                        elptN = vec_elN.Pt()
                    muptP = vec_muP.Pt()
                    muptN = vec_muN.Pt()
                    muEtaP = vec_muP.PseudoRapidity()
                    muEtaN = vec_muN.PseudoRapidity()
        if fillAtEnd and bestjj > -1:
            hpT[vtype].Fill(pt, evt_weight)
            if vstr != "mumu":
                hpTEl[vtype].Fill(elptP, evt_weight)
                hpTEl[vtype].Fill(elptN, evt_weight)
                hpTElNoWt[vtype].Fill(elptP)
                hpTElNoWt[vtype].Fill(elptN)
            hpTMu[vtype].Fill(muptP, evt_weight)
            hpTMu[vtype].Fill(muptN, evt_weight)
            hpTMuNoWt[vtype].Fill(muptP)
            hpTMuNoWt[vtype].Fill(muptN)
            hEtaMu[vtype].Fill(muEtaP, evt_weight)
            hEtaMu[vtype].Fill(muEtaN, evt_weight)
            
            hEta[vtype].Fill(besteta, evt_weight)
            if isMC:
                #fill event weight histograms ONLY for accepted events! (correct mass range)
                if bestm > peakmin and bestm < peakmax:
                    hEventWeight[vtype].Fill(evt_weight) 
                    hEvtWtVsPt[vtype].Fill(gen_eta.Pt(), evt_weight)
            #if danEvent and bestm > .52 and bestm < .58:
            if syncTest and vtype == "mumu" and bestm > peakmin and bestm < peakmax:
                syncFile.write("%d %f\n"%(e.evt, bestm)) 
                #try:
                #    printEvent.printEvent(e)
                #except AttributeError:
                #    printEvent_backup.printEvent(e)
                print("mass: %f"%bestm) 
            #why m instead of bestm???
            #hM[vtype].Fill(m, evt_weight)
            #hMNoWt[vtype].Fill(m)
            #hMvsPt[vtype].Fill(pt, m, evt_weight)
            hM[vtype].Fill(bestm, evt_weight)
            if bestm > peakmin and bestm < peakmax:
                ncand += 1
            hMNoWt[vtype].Fill(bestm)
            hnPVgood[vtype].Fill(e.nPV, evt_weight)
            if isMC:
                if do_ElRecsf:
                    hMRecUp[vtype].Fill(bestm, evt_weightRecUp)
                    hMRecDn[vtype].Fill(bestm, evt_weightRecDn)
                if do_ElIDsf:
                    hMIDUp[vtype].Fill(bestm, evt_weightIDUp)
                    hMIDDn[vtype].Fill(bestm, evt_weightIDDn)
            if isMC and do_MuIDsf:
                hMMuUp[vtype].Fill(bestm, evt_weightMuUp)
                hMMuDn[vtype].Fill(bestm, evt_weightMuDn)
            if not passedTrig:
                print("filling failedTrig!!!! mass %f, weight %f"%(bestm, evt_weight))
                hMFailedTrig[vtype].Fill(bestm, evt_weight) 
            if not (year == 2023 and vtype == "mumu"):
                hMvsPt[vtype].Fill(bestpt, bestm, evt_weight)
            hVxy[vtype].Fill(Vxy)
            hRchi2[vtype].Fill(bestChi2)
            if vtype == "mmelel" and m > peakmin and m < peakmax:
                hMNoEl.Fill(bestm2mu, evt_weight)
                hMNoMu.Fill(bestm2el, evt_weight)
            #if vtype == "mmelel" and m < .52 and m > .45 :
            #    hMNoElLSide.Fill(bestm2mu, evt_weight)
            #    hMNoMuLSide.Fill(bestm2el, evt_weight)
            #if vtype == "mmelel" and m > .58 and m < .8 :
            #    hMNoElRSide.Fill(bestm2mu, evt_weight)
            #    hMNoMuRSide.Fill(bestm2el, evt_weight)
            if vtype == "mmelel" and ((m > (peakmax+.03) and m < (peakmax+2*.03)) or (m > (peakmin-2*.03) and m < (peakmin-2*.03))) :
                hMNoElSide.Fill(bestm2mu, evt_weight)
                hMNoMuSide.Fill(bestm2el, evt_weight)
            #if vxy > 1.2:
            #    #hMhiVxy[vtype].Fill(m, evt_weight)
            #    hMhiVxy[vtype].Fill(bestm, evt_weight)
            #else:
            #    hMloVxy[vtype].Fill(bestm, evt_weight)
        hNcand[vtype].Fill(ncand)
    return good

#return true if at least two oppositely charged lepnames are found in the event, False otherwise
def found_oppo(e, lepname):
    found_pos = False
    found_neg = False
    nlep = eval("ord(e.nGood%s)"%lepname) 
    for l in range(nlep):
        charge = eval("ord(e.%s_charge[l])"%lepname)
        if charge == 1:
            found_pos = True
        else:
            found_neg = True
    return (found_pos and found_neg)

#return the sub-leading muon pT (lower amongst the highest-pt muons of each charge)
def find_subpt(e):
    if ord(e.nGoodMuon) < 2:
        return -1
    c0 = ord(e.Muon_charge[0])
    imu = 1
    while imu < ord(e.nGoodMuon) and ord(e.Muon_charge[imu]) == c0:
        imu += 1
    if imu == ord(e.nGoodMuon):
        #no oppo-charge muon
        return -1
    return e.Muon_pt[imu]

if golden and not isMC:
    import json
    goldname="Cert_Collisions2022_355100_362760_Golden.json"
    goldf = open(goldname, "r")
    goldl = json.load(goldf)
    goldList = dict(goldl) 

#return True if it's a GOOD event, in the golden JSON.
def checkGolden(run, lusec):
    if str(run) not in goldList:
        return False
    for lrange in goldList[str(run)]:
        if lusec > lrange[0] and lusec < lrange[1]:
            return True
    return False

#now open each file individually again, and process the events
def process_file(fname, singleVert, useOnia, hWeights):
    global all_weight, trg_weight
    print("Opening file %s"%fname) 
    f = ROOT.TFile.Open(fname)
    t = f.Get("ntuples/recoT") 
    if "output_" in fname:
        isCentral = True
    else:
        isCentral = False
    if not t:
        t = f.Get("recoT")  
    nTot = t.GetEntries()
    if isMC:
        gt = f.Get("ntuples/genT") 
        if not gt:
            gt = f.Get("genT") 
        nEntries = t.GetEntries()
        if year == 2022:
            if EEperiod == -1:
                #preEE
                lumi = 9.76
            elif EEperiod == 1:
                #postEE
                lumi = 28.25
            else:
                lumi = 38.48 #fb^-1 (this is the lumi for all CMS in 2022)
        elif year == 2023:
            lumi = 28.89
        #lumi = 38.48 + 28.89 #fb^-1; for 2022+2023
        if isSig:
            sys.path.insert(1, "tm_analysis/analysis/python/utils")
            import blinding
            bratio = blinding.get_blinded_BR_EtaTo2Mu2E()
        elif isMuMu:
            bratio = 5.8e-6
        else:
            bratio = 3.1e-4

    for i,e in enumerate(t):
        ##debugging
        #if i > 100000: break
        if i%printevery == (printevery-1):
            print("Event %d/%d"%(i+1, nTot)) 

        #fill the number of events processed tracker
        hNevt.Fill(1)

        #check the golden JSON (if necessary)
        if golden and not isMC:
            isGold = checkGolden(e.run, e.lumi_sec)
            if not isGold:
                continue

        #if syncTest and e.evt in dan_events:
        #    idx = dan_events.index(e.evt)
        #    if dan_runs[idx] == e.run:
        #        print("**danEvent**")
        #        try:
        #            printEvent.printEvent(e)
        #        except AttributeError:
        #            printEvent_backup.printEvent(e)
        #if i < 100:
        #    print("Run %d ls %d evt %d"%(e.run, e.lumi_sec, e.evt)) 
        passedTrig = False
        #true if rejected for having photon
        failedPhot = (reject_photon and ord(e.nGoodPhoton) > 0)
        if failedPhot: 
            continue
        trig0 = e.Triggers_fired0
        trig1 = e.Triggers_fired1
        if year == 2022 and int(str(testnum)[:2]) > 38:
            trig2 = e.Triggers_fired2

        #if trig0 > 0 or ord(trig1) > 0:
        #trying to accept only a few special triggers
        #if trig0 & ((1<<27) + (1<<26)) > 0:
        #if year == 2022 and (trig0 & ((1<<11) + (1<<12) + (1<<17) + (1<<25) + (1<<26) + (1<<28)) > 0):
        #if year == 2022 and (trig0 & (1<<25) > 0):
        if int(str(testnum)[:2]) > 38 or ((not isMC) and (let == "Dv1" or (let == "F" and num == 1))) or (isMC and (not isSig) and (not isMuMu) and isCentral) :
            #new trigger words
            if year == 2022 and (trig1 & (1<<11) > 0):
                passedTrig = True
            elif year == 2023: 
                try:
                    passedTrig = (trig1 & (1<<5) > 0)
                except TypeError:
                    passedTrig = (ord(trig1) & (1<<5) > 0)
            elif trg_only:
                if not isMC:
                    pass
                    #print("Error!!! Data event failed trigger??? Run %d LS %d Event %d"%(e.run, e.lumi_sec, e.evt))
                    #try: 
                    #    printEvent.printEvent(e)
                    #except AttributeError:
                    #    printEvent_backup.printEvent(e)
                    #if not singleFile:
                    #    trigF.write("%d\n"%e.evt) 
        else:
            #old words
            if year == 2022 and (trig0 & (1<<11) > 0):
                passedTrig = True
            #elif year == 2023 and (ord(trig0) & (1<<5) > 0):
            elif year == 2023 :
                if 'v' in let and let != 'Cv4':
                    passedTrig = (trig1 & (1<<11) > 0)
                else:
                    passedTrig = (trig0 & (1<<11) > 0)
        if testnum == 38104:
            passedTrig = (trig0 != 0 or trig1 != 0) 

        #weight is 1 for data, different for MC
        evt_weight = 1.0
        if isMC:
            #get the gen tree event.
            gt.GetEntry(i)
            g = gt
            #construct the eta meson
            genEtaPt = 0
            genEtaEta = -9999
            pTmu = -1
            pTel = -1
            old_matched = False
            gen_eta = ROOT.TLorentzVector()
            #first muon encountered
            old_vec = ROOT.TLorentzVector()
            old_elV = ROOT.TLorentzVector()
            if not isSig and 221 not in g.GenPart_pdgId:
                #print("No gen eta meson??? Evt %d; pdgIds: %s"%(i, str(g.GenPart_pdgId)))
                continue
            #else:
            #    print("eta meson found!! Evt %d; pdgIds: %s"%(i, str(g.GenPart_pdgId)))

            #recod will be false if any gen particle doesn't get genmatched to a reco'd particle
            recodmumu = True
            recodelel = True if isSig or isMuMu else False
            mindrmup = 9999
            mindrmun = 9999
            mindrelp = 9999
            mindreln = 9999
            genmudR = -1
            gensubpt = -1
            gensubeta = 9999
            for j in range(g.nGenPart):
                gid = g.GenPart_pdgId[j]
                gpt = g.GenPart_pt[j]
                if abs(gid) == 221:
                    #eta meson
                    genEtaPt = gpt
                    genEtaEta = g.GenPart_eta[j]
                    #print("genEtaPt=%f"%genEtaPt) 
                    gen_eta.SetPtEtaPhiM(gpt, g.GenPart_eta[j], g.GenPart_phi[j], g.GenPart_mass[j]) 
                    #break
                else:
                    foundmatch = False
                    geta = g.GenPart_eta[j] 
                    gphi = g.GenPart_phi[j]
                    gmass = g.GenPart_mass[j]
                    new_vec = ROOT.TLorentzVector()
                    new_vec.SetPtEtaPhiM(gpt, geta, gphi, gmass) 
                    if genEtaPt == 0 and gid in [-11, 11, -13, 13, 22]:
                        gen_eta = gen_eta + new_vec
                    #elif gid != 990:
                    #    print("Unrecognized pdgId: %d"%(gid)) 
                    if abs(gid) == 13:
                        #muon
                        if abs(g.GenPart_eta[j]) < 1.5:
                            hGenMupTAllC.Fill(g.GenPart_pt[j]) 
                        else:
                            hGenMupTAllF.Fill(g.GenPart_pt[j]) 
                        #if this is the 2nd muon, compute the dR and fill the histograms
                        if pTmu == -1:
                            pTmu = g.GenPart_pt[j]
                            old_vec = new_vec
                        else:
                            genmudR = old_vec.DeltaR(new_vec)
                            hGenMudR.Fill(genmudR)
                            if pTmu > g.GenPart_pt[j]:
                                #pTmu is lead pT, this one is sublead
                                hGenMupT0.Fill(pTmu)
                                hGenMupT1.Fill(g.GenPart_pt[j]) 
                                gensubpt = g.GenPart_pt[j]
                                gensubeta = g.GenPart_eta[j]
                            else:
                                hGenMupT0.Fill(g.GenPart_pt[j])
                                hGenMupT1.Fill(pTmu)
                                gensubpt = pTmu
                                gensubeta = old_vec.PseudoRapidity()

                            if abs(gensubeta) < 1.5:
                                hsubVdRGenAllC.Fill(gensubpt, genmudR) 
                            else:
                                hsubVdRGenAllF.Fill(gensubpt, genmudR) 
                            #here can add in the trigger efficiency correction
                            # below 5 GeV not enough statistics.
                            if do_trigCor and gensubpt > 5:  
                                #xybin = trigCor.FindBin(gensubpt, genmudR)
                                #evt_weight *= trigCor.GetBinContent( trigCor.FindBin(gensubpt, genmudR) )
                                #evt_weight *= ( trigCor.GetBinContent( xybin ) ) # + trigCor.GetBinError( xybin ) ) 
                                #flat trigger correction found from AN2023-026 (D0->mumu analysis)
                                evt_weight *= 0.943
                            
                        #see if this gen particle has a reco gen match -- if not, recod=False
                        for k in range(len(e.Muon_eta)):
                            ch = ord(e.Muon_charge[k])
                            #NEGATIVE muons have POSITIVE pdgId!!!
                            if ch == 1 and gid > 0: continue
                            elif ch != 1 and gid < 0: continue
                            rec_vec = ROOT.TLorentzVector()
                            rec_vec.SetPtEtaPhiM(e.Muon_pt[k]*(1+musf), e.Muon_eta[k], e.Muon_phi[k], mu_mass) 
                            rec_dr = rec_vec.DeltaR(new_vec) 
                            if rec_dr < dRcut and not foundmatch:
                                hGenMupTRec.Fill(g.GenPart_pt[j]) 
                                foundmatch = True
                            if ch == 1 and rec_dr < mindrmup:
                                mindrmup = rec_dr
                            elif ch != 1 and rec_dr < mindrmun:
                                mindrmun = rec_dr
                                #break
                            hMudRvPt.Fill(rec_dr, gpt) 
                    elif abs(gid) == 11:
                        hGenElpTAll.Fill(g.GenPart_pt[j]) 
                        if not useLowPt:
                            if pTel == -1:
                                pTel = g.GenPart_pt[j]
                                old_elV = new_vec
                            else:
                                geneldR = old_elV.DeltaR(new_vec)
                                #fill the dR hist twice for the two electrons
                                hGenEldRvPtAll.Fill(geneldR, pTel, evt_weight)
                                hGenEldRvPtAll.Fill(geneldR, g.GenPart_pt[j], evt_weight)
                                if old_matched:
                                    hGenEldRvPtRec.Fill(geneldR, pTel, evt_weight) 
                                pTel = -1
                            #see if this gen particle has a reco gen match -- if not, recod=False
                            for k in range(len(e.Electron_eta)):
                                ch = ord(e.Electron_charge[k])
                                #electrons have a POSITIVE ID, positrons have NEGATIVE
                                if ch == 1 and gid > 0: continue
                                elif ch != 1 and gid < 0: continue
                                rec_vec = ROOT.TLorentzVector()
                                rec_vec.SetPtEtaPhiM(e.Electron_pt[k]*(1+elsf), e.Electron_eta[k], e.Electron_phi[k], el_mass) 
                                rec_dr = rec_vec.DeltaR(new_vec) 
                                #for electrons need higher dRcut!
                                if rec_dr < dRcut and not foundmatch:
                                    foundmatch = True
                                    hGenElpTRec.Fill(g.GenPart_pt[j]) 
                                    if pTel == -1:
                                        #pTel would only be -1 here if both gen electrons have been looped thru so we have a geneldR
                                        hGenEldRvPtRec.Fill(geneldR, g.GenPart_pt[j], evt_weight)
                                    else:
                                        #we'll fill this when we find the partner electron (need it to know what the dR is)
                                        old_matched = True
                                    #break
                                if ch == 1 and rec_dr < mindrelp:
                                    mindrelp = rec_dr
                                elif ch != 1 and rec_dr < mindreln:
                                    mindreln = rec_dr
                                #hEldRvPt.Fill(rec_dr, gpt) 
                                hEldR.Fill(rec_dr)
                        else:
                            #use lowpt electrons instead
                            for k in range(len(e.LowPtElectron_eta)):
                                ch = ord(e.LowPtElectron_charge[k])
                                if ch == 1 and gid > 0: continue
                                elif ch != 1 and gid < 0: continue
                                rec_vec = ROOT.TLorentzVector()
                                rec_vec.SetPtEtaPhiM(e.LowPtElectron_pt[k]*(1+elsf), e.LowPtElectron_eta[k], e.LowPtElectron_phi[k], el_mass) 
                                rec_dr = rec_vec.DeltaR(new_vec) 
                                #for electrons need higher dRcut!
                                if rec_dr < dRcut and not foundmatch:
                                    foundmatch = True
                                    hGenElpTRec.Fill(g.GenPart_pt[j]) 
                                    #break
                                if ch == 1 and rec_dr < mindrelp:
                                    mindrelp = rec_dr
                                elif ch != 1 and rec_dr < mindreln:
                                    mindreln = rec_dr
                                #hEldRvPt.Fill(rec_dr, gpt) 
                                hEldR.Fill(rec_dr)
                    #photon
                    elif gid == 22:
                        if not isSig and not isMuMu:
                            #set it to true for now (will be set to false later if not actually recod)
                            recodelel = True
                        #for photons, need to find the pos and neg electron match!
                        foundpos = False
                        foundneg = False
                        for k in range(len(e.Electron_eta)):
                            ch = ord(e.Electron_charge[k])
                            if ch == 1 and foundpos: continue
                            elif ch != 1 and foundneg: continue
                            rec_vec = ROOT.TLorentzVector()
                            rec_vec.SetPtEtaPhiM(e.Electron_pt[k]*(1+elsf), e.Electron_eta[k], e.Electron_phi[k], el_mass) 
                            rec_dr = rec_vec.DeltaR(new_vec) 
                            if rec_dr < dRcut:
                                if ch == 1:
                                    foundpos = True
                                else:
                                    foundneg = True
                                if foundpos and foundneg:
                                    foundmatch = True
                                    #break
                            hEldR.Fill(rec_dr)
                    if abs(gid) == 13 and not foundmatch:
                        recodmumu = False
                    elif (abs(gid) == 11 or gid == 22) and not foundmatch:
                        recodelel = False
                        #if g.GenPart_pt[j] > 20:
                        #    print("NO e MATCH!!!")
                        #    if int(str(testnum)[:2]) > 38:
                        #        printEvent.printEvent(e, True, g)
                        #    else:
                        #        printEvent_backup.printEvent(e, True, g)
                #only 4 gen particles max probably
                if j == (g.nGenPart-1) and genEtaPt == 0:
                    genEtaPt = gen_eta.Pt()
                    genEtaEta = gen_eta.PseudoRapidity()
                    #if genEtaPt > 65:
                    #    print("**genEtaPt: %f; gid: %d nGenPart:%d"%(genEtaPt, gid, g.nGenPart)) 
            hMudRvPtEta.Fill(mindrmup, genEtaPt) 
            hMudRvPtEta.Fill(mindrmun, genEtaPt) 
            if isSig:
                hEldRvPtEta.Fill(mindrelp, genEtaPt) 
                hEldRvPtEta.Fill(mindreln, genEtaPt) 
            if abs(trg_var) < 100:
                passedTrig = toyTrig(gensubpt)
            #if genEtaPt > 65:
            #    print("genEtaPt: %f; nGenPart:%d"%(genEtaPt, g.nGenPart)) 
            #xbin = h_xsec.FindBin( genEtaPt )
            #xsec0 = h_xsec.GetBinContent( xbin ) * h_xsec.GetBinWidth( xbin )
            #testing this new xsec measurement???
            #xsec1 = 6.2615311e+15 / genEtaPt**5.8244956
            #fitted value of parameter 0
            #xsec_p0 = 6.6657161e+14
            #xsec_p0 = 7.1656026e+14
            xsec_p0 = 7.60902e+14
            #xsec_p0 = 4.35099e+14
            #fitted uncertainty on parameter 0
            #xsec_unc0 = 7.3097364e+13
            #xsec_unc0 = 5.3862360e+14
            xsec_unc0 = 5.73166e+14
            #xsec_unc0 = 3.51337e+14
            #xsec_p1 = 4.9552
            #xsec_p1 = 5.1121328
            xsec_p1 = 5.00655
            #xsec_p1 = 5.17427
            #xsec_unc1 = 0.037738932
            #xsec_unc1 = 0.22696893
            xsec_unc1 = 0.226186
            #xsec_unc1 = 0.250471 
            #now using new weights!! from DG/Cheb4 2mu fits
            #if require_muID:
            if new_wt == 1:
                #medium ID
                #xsec_p0 = 1.6093980e+15
                #xsec_p1 = 5.2268649
                #loose ID
                xsec_p0 = 7.22777e+14
                xsec_p1 = 5.01568
            #elif new_wt == 1:
            #    xsec_p0 = 1.34311e+15
            #    xsec_unc0 = 8.30077e+14
            #    xsec_p1 = 5.22827
            #    xsec_unc1 = 1.80061e-01 
            elif new_wt == 2:
                #DG, Cheb3!
                #xsec_p0 = 9.4012018e+14
                #xsec_p1 = 5.1077085
                xsec_p0 = 9.31211e+14
                xsec_p1 = 5.1051
            elif new_wt == 3:
                #xsec_p0 = 8.4903893e+14
                #xsec_p1 = 5.0395503
                #corrected for 2mu!!
                xsec_p0 = 9.82423e+14
                xsec_p1 = 5.03434
            xsec1 = xsec_p0 / genEtaPt**xsec_p1
            xsecUp = (xsec_p0 + xsec_unc0) / genEtaPt**(xsec_p1 - xsec_unc1) 
            xsecDn = (xsec_p0 - xsec_unc0) / genEtaPt**(xsec_p1 + xsec_unc1) 
            #print("xsec0: %f, xsec1: %f"%(xsec0, xsec1))
            #hxs0.append(xsec0) 
            #hxs1.append(xsec1)
            xsec = xsec1
            #xsec = xsec0
            wbin = hWeights.FindBin(genEtaPt)
            ptWeight = 1.0 * hWeights.GetBinContent( wbin )
            if ptWeight == 0:
                print("ptWeight = 0 !!!!! genEtaPt = %f, wbin = %d"%(genEtaPt, wbin)) 
            if bratio == 0:
                print("bratio = 0!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            
            evt_weight *= (xsec * bratio * lumi / ptWeight)
            ##Add corrections!
            #if do_pileup:
            #    #pileup correction
            #    puCor = PUCor.GetBinContent(PUCor.FindBin(e.nPV))
            #    evt_weight *= puCor

            if year != 2023:
                hpTGenAll.Fill(genEtaPt, evt_weight)
                hMGenAll.Fill(gen_eta.M(), evt_weight) 

            if genEtaPt > 0: #25:
                hEtaGenAll.Fill(genEtaEta, evt_weight)

            #uncertainties
            #evt_weightUp = xsecUp/xsec * evt_weight
            #evt_weightDn = xsecDn/xsec * evt_weight
            if ord(e.nGoodElectron) > 1 and ord(e.nGoodMuon) > 1 and evt_weight > 50 and genEtaPt > 20:
                print("***high event weight: %f***"%(evt_weight))
                print("event: %d, genEtaPt: %f, xsec: %f, ptWeight: %f"%(i, genEtaPt, xsec, ptWeight)) 
            #this can happen a lot because of pileup corrections
            #if ord(e.nGoodElectron) > 1 and ord(e.nGoodMuon) > 1 and genEtaPt < 12 and evt_weight < 0.01:
            #    print("******low event weight: %f***"%(evt_weight))
            #    print("event: %d, genEtaPt: %f, xsec: %f, ptWeight: %f, lumi: %f"%(i, genEtaPt, xsec, ptWeight, lumi)) 
            #    print("dsigma/dpT: %f, binwidth: %f"%(h_xsec.GetBinContent( xbin ), h_xsec.GetBinWidth( xbin ))) 
            all_weight += evt_weight
            if passedTrig and year != 2023:
                hpTGenTrig.Fill(genEtaPt, evt_weight)
                trg_weight += evt_weight
            if passedTrig:
                if genEtaPt > 0: #25:
                    hEtaGenTrig.Fill(genEtaEta, evt_weight) 
                if abs(gensubeta) < 1.5:
                    hsubVdRGenTrigC.Fill(gensubpt, genmudR) 
                else:
                    hsubVdRGenTrigF.Fill(gensubpt, genmudR) 
            if int(str(testnum)[:2]) > 38:
                if (year == 2022 and (ord(trig2) & (1<<2))) or (year == 2023 and (trig1 & (1<<27))):
                    hEtaGenTrigL0.Fill(genEtaEta, evt_weight) 
                    if (year == 2022 and (trig1 & (1<<25))) or (year == 2023 and (trig1 & (1<<26))):
                        hEtaGenTrigPFL1.Fill(genEtaEta, evt_weight) 
                if (year == 2022 and (trig1 & (1<<25))) or (year == 2023 and (trig1 & (1<<26))):
                    hEtaGenTrigL1.Fill(genEtaEta, evt_weight) 
                    if (year == 2022 and (trig1 & (1<<11))) or (year == 2023 and (trig1 & (1<<5))):
                        hEtaGenTrigHLTL1.Fill(genEtaEta, evt_weight) 
            else:
                #if (year == 2022 and (trig0 & (1<<25))) or (year == 2023 and (trig0 & (1<<26))):
                if (year == 2022 and (trig0 & (1<<25))) or (year == 2023 and (trig0 & (1<<25))):
                    hEtaGenTrigL1.Fill(genEtaEta, evt_weight) 
                    #if (year == 2022 and (trig0 & (1<<11))) or (year == 2023 and (trig0 & (1<<5))):
                    if (year == 2022 and (trig0 & (1<<11))) or (year == 2023 and (trig0 & (1<<11))):
                        hEtaGenTrigHLTL1.Fill(genEtaEta, evt_weight) 
            
        #fill in the Primary vertices histogram with the appropriate event weight
        hnPV.Fill(e.nPV, evt_weight)

        #See if can fill the trigger eff histograms (only for events with subleading pt > 5 GeV
        #find the subleading pt
        mu_subpt = find_subpt(e)
        if int(str(testnum)[:2]) > 38: 
            #see if passed the L1 one
            l1pass = False
            hlpass = False
            if year == 2022 :
                l1pass = (trig1 & (1<<25)) != 0
            elif year == 2023 :
                try:
                    l1pass = (ord(trig1) & (1<<26))
                except TypeError:
                    l1pass = False
            if year == 2022 :
                hlpass = (trig1 & (1<<11) > 0)
            elif year == 2023:
                try:
                    hlpass = (trig1 & (1<<5) != 0)
                except TypeError:
                    hlpass = False
            if l1pass:
                hpTL1Pass.Fill(mu_subpt, evt_weight)
            if hlpass:
                hpTHLTPass.Fill(mu_subpt, evt_weight)
                
        else:
            #see if passed the L1 one
            #if (year == 2022 and (trig0 & (1<<25))) or (year == 2023 and (trig0 & (1<<26))):
            if (year == 2022 and (trig0 & (1<<25))) or (year == 2023 and (trig0 & (1<<25))):
                hpTL1Pass.Fill(mu_subpt, evt_weight)
                #if (year == 2022 and (trig0 & (1<<11))) or (year == 2023 and (trig0 & (1<<5))):
                if (year == 2022 and (trig0 & (1<<11))) or (year == 2023 and (trig0 & (1<<11))):
                    hpTHLTPass.Fill(mu_subpt, evt_weight)

        goodmumu = False
        goodmmg = False
        for vtype in vtypes:
            if isMC:
                ##see if reco criteria is fulfilled: just the two (or four) opposite-charged leptons were reconstructed
                recod = False
                if vtype == "mumu" and recodmumu:
                    recod = True
                elif vtype == "elel" and recodelel:
                    recod = True
                elif vtype in ["mmelel", "mmlplp"] and recodmumu and recodelel:
                    recod = True
                #new reco eff measure: req gen matching!! -- done above
                ##require opposite charge pairs even for background sample (where memory consumption can get extreme)?
                #req_bkgOppo = True
                #if ord(e.nGoodMuon) > 1:
                #    if vtype == "mumu":
                #        #for some reason bkg takes way too much memory when I check the charges....whack
                #        if isSig or isMuMu:
                #            recod = found_oppo(e, "Muon")
                #        else:
                #            recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #    elif vtype == "mmelel" and ord(e.nGoodElectron) > 1:
                #        if isSig or isMuMu:
                #            recod = found_oppo(e, "Muon") and found_oppo(e, "Electron")
                #        else:
                #            recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #    elif vtype == "mmlplp" and ord(e.nGoodLowPtElectron) > 1:
                #        if isSig or isMuMu:
                #            recod = found_oppo(e, "Muon") and found_oppo(e, "LowPtElectron")
                #        else:
                #            recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #    elif vtype == "mmg" and ord(e.nGoodPhoton) > 0:
                #        if isSig or isMuMu:
                #            recod = found_oppo(e, "Muon") 
                #        else:
                #            recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #if vtype == "elel" and ord(e.nGoodElectron) > 1:
                #    if isSig or isMuMu:
                #        recod = found_oppo(e, "Electron")
                #    else:
                #        recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #elif vtype == "lplp" and ord(e.nGoodLowPtElectron) > 1:
                #    if isSig or isMuMu:
                #        recod = found_oppo(e, "LowPtElectron")
                #    else:
                #        recod = True if not req_bkgOppo else found_oppo(e, "Muon")
                #    #for 2023 it's calculated in the ntuplizer, this would be wrong
                if recod:
                    if gen_eta.Pt() > 0: #25:
                        hEtaGenReco[vtype].Fill(gen_eta.PseudoRapidity(), evt_weight)
                    #if year != 2023:
                    #reco def is just enough particles now; doesn't need to be in right mass range
                    hpTGenReco[vtype].Fill(gen_eta.Pt(), evt_weight)
                    rec_weight[vtype] += evt_weight
                else:
                    hpTGenNotReco[vtype].Fill(gen_eta.Pt(), evt_weight)
                #else:
                #    continue

                #must continue AFTER filling the reco eff, to make sure it has no regard to trigger eff!!!
                if trg_only and not passedTrig:
                    continue
                if vtype in ["mmlplp", "mmelel"] and ((reject_etamumu and goodmumu) or (reject_mmg and goodmmg)):
                    #print("event %d rejected for etaToMuMu!"%i) 
                    continue
                if vtype == "mumu" and reject_mmg and goodmmg:
                    continue

                isGood = process_vertices(e, vtype, singleVert, useOnia, xsec, evt_weight, g, gen_eta, passedTrig, recod)
                #isGood = False
            else:
                #must continue AFTER filling the reco eff, to make sure it has no regard to trigger eff!!!
                if trg_only and not passedTrig:
                    continue

                isGood = process_vertices(e, vtype, singleVert, useOnia, 0, 1.0, passedTrig=passedTrig)
            if vtype == "mumu":
                goodmumu = isGood
            elif vtype == "mmg":
                goodmmg = isGood

        if diel != None:
            global old_diel
            old_diel = diel


    #rm the file now that you're done with it.
    #if not (isMC and not isSig):
    if (not isMC) and (not syncTest) and year != 2023:
        os.system("rm %s"%fname)

#finish the processing and write to an output file
def finish_processing(foutname):
    succ = True
    try:
        fout = ROOT.TFile.Open(foutname, "recreate") 
    except OSError:
        #print("Error! Could not open %s in the current directory so trying in /afs/cern.ch/work/b/bgreenbe/public"%foutname) 
        print("Error! Could not open %s :("%foutname) 
        fout = ROOT.TFile.Open("/afs/cern.ch/work/b/bgreenbe/public/%s"%foutname, "recreate") 
        succ = False
        return succ
        
    hNevt.Write()
    hnPV.Write()
    hpTL1Pass.Write()
    hpTHLTPass.Write()
    hpTGoodL1Pass.Write()
    hpTGoodHLTPass.Write()
    hMSSe.Write()
    hMOSe.Write()
    hMee.Write()
    hMeePeak.Write()
    hMComb.Write()
    hMReal.Write()
    hL1TC.Write()
    hL1TF.Write()
    hHLTC.Write()
    hHLTF.Write()
    if isMC: #and not isMuMu:
        if not isMuMu:
            hMeedR.Write()
            hMeeGM.Write()
            hEldR.Write()
        hMmmdR.Write()
    hMeeVsMmmee.Write()
    hpTElAll.Write()
    hpTElAll0.Write()
    hpTElAll1.Write()
    hpTElAll2.Write()
    hpTElAll3.Write()
    hpTElIDPass.Write()
    hpTElIDPass0.Write()
    hpTElIDPass1.Write()
    hpTElIDPass2.Write()
    hpTElIDPass3.Write()
    for vtype in vtypes:
        if isMC and vtype not in ["mumu", "mmg"]:
            hElRecsf[vtype].Write()
            hElRecsfUp[vtype].Write()
            hElRecsfDn[vtype].Write()
            hElIDsf[vtype].Write()
            hElIDsfUp[vtype].Write()
            hElIDsfDn[vtype].Write()
            hEletapt[vtype].Write()
        hM[vtype].Write()
        hMNoWt[vtype].Write()
        hMFailedTrig[vtype].Write()
        hdREl[vtype].Write()
        hdRvPtEl[vtype].Write()
        hL1all[vtype].Write()
        hL1pass[vtype].Write()
        hnPVgood[vtype].Write()
        #if vtype != "mumu" and vtype != "mmg":
        hpT[vtype].Write()
        hpTSig[vtype].Write()
        hpTSide[vtype].Write()
        hpTEl[vtype].Write()
        hpTElNoWt[vtype].Write()
        hpTMu[vtype].Write()
        hpTMuF[vtype].Write()
        hpTMuNoWt[vtype].Write()
        hEtaMu[vtype].Write()
        hEta[vtype].Write()
        hNcand[vtype].Write()
        #hMhiVxy[vtype].Write()
        #hMloVxy[vtype].Write()
        #hDxyHiVxy[vtype].Write()
        #hDxyLoVxy[vtype].Write()
        #hSigmaVxyHiVxy[vtype].Write()
        #hSigmaVxyLoVxy[vtype].Write()
        #hdRP[vtype].Write()
        #hdRN[vtype].Write()
        hMvsPt[vtype].Write()
        hsubVdR[vtype].Write()
        hVxy[vtype].Write()
        hRchi2[vtype].Write()
        if isMC and not isSig:
            hMeeG[vtype].Write()
            hMeeNoG[vtype].Write()
        if useOnia:
            hVertMatched[vtype].Write()
            hVertNoMatch[vtype].Write()
            if isMC:
                hpTVertMatched[vtype].Write()
    if mmelelExclusive:
        hMdeleted.Write()
        hRchi2deleted.Write()
    hMNoEl.Write()
    hMNoMu.Write()
    hMNoElLSide.Write()
    hMNoMuLSide.Write()
    hMNoElRSide.Write()
    hMNoMuRSide.Write()
    hMNoElSide.Write()
    hMNoMuSide.Write()
    #hMNoMuPiM.Write()
    if inc_vertM:
        hMV.Write()
        hpTV.Write()
    if isMC:
        hMudRvPt.Write()
        hMudRvPtEta.Write()
        hpTGenAll.Write()
        hMGenAll.Write()
        hEtaGenAll.Write()
        for s in range(nMod):
            hMMod[s].Write()
            hMModmu[s].Write()
            hpTModmu[s].Write()
        sel = {}
        for vtype in vtypes:
            hMRecUp[vtype].Write()
            hMRecDn[vtype].Write()
            hMIDUp[vtype].Write()
            hMIDDn[vtype].Write()
            hMTTOCUp[vtype].Write()
            hMTTOCDn[vtype].Write()
            hMMuUp[vtype].Write()
            hMMuDn[vtype].Write()
            hEventWeight[vtype].Write()
            hEvtWtVsPt[vtype].Write()
            hpTGenReco[vtype].Write()
            hpTGenNotReco[vtype].Write()
            hpTGenAcc[vtype].Write()
            hpTGenAccdR[vtype].Write()
            hEtaGenReco[vtype].Write()
            hEtaGenAcc[vtype].Write()
            hvxy_gm[vtype].Write()
            #write out the acceptances as a TGraph
            acc_weight[vtype][0] = all_weight
            sel[vtype] = ROOT.TGraph( nselections+1, array.array('f', [s for s in range(nselections+1)]), array.array('f', acc_weight[vtype]) ) 
            sel[vtype].SetName("selAcc"+vtype)
            sel[vtype].Write()
            print("wrote the tg for " + vtype)
        hpTGenTrig.Write()
        hEtaGenTrig.Write()
        hEtaGenTrigL0.Write()
        hEtaGenTrigPFL1.Write()
        hEtaGenTrigL1.Write()
        hEtaGenTrigHLTL1.Write()
        hGenMudR.Write()
        hGenMupT0.Write()
        hGenMupT1.Write()
        hGenMupTAllC.Write()
        hGenMupTAllF.Write()
        hGenMupTRec.Write()
        hGenElpTAll.Write()
        hGenElpTRec.Write()
        hGenEldRvPtAll.Write()
        hGenEldRvPtRec.Write()
        hsubVdRGenAllC.Write()
        hsubVdRGenTrigC.Write()
        hsubVdRGenAllF.Write()
        hsubVdRGenTrigF.Write()
        if isSig:
            hEldRvPtEta.Write()
            for ptype in ptypes:
                hdR[ptype].Write()
        #xsComp = ROOT.TGraph(len(hxs0), array.array('d', hxs0), array.array('d', hxs1))
        #xsComp.Write()
    fout.Close()
    print("time: %d seconds"%(time.time()-start)) 
    #hM.Draw()
    #input("hhhhhhhh")
    if isMC:
        print("**Weighted Efficiencies**")
        print("Trigger eff: %f%%"%(trg_weight/all_weight*100))
        for vtype in vtypes:
            print("Reco eff (%s): %f%%"%(vtype, rec_weight[vtype]/all_weight*100))
        print("**Overall weighted acceptances (trig + reco)**") 
        for vtype in vtypes:
            for s in range(nselections):
                print("Acc (%s) selection %d: %f%%"%(vtype, s, acc_weight[vtype][s]/all_weight*100)) 
            print("AccdR (%s): %f%%"%(vtype, accdR_weight[vtype]/all_weight*100)) 

    #return True if saved file in the current dir, False otherwise
    return succ

if syncTest and singleFile:
    foutname = "bparking_syncTest_test%d.root"%(testnum)
elif year == 2023 and not isMC:
    foutname = "bparking_2023datatest%d_ALL.root"%testnum
elif year == 2023 and isMC and isSig:
    foutname = "bparking_2023sigMCtest%d.root"%testnum
elif year == 2023 and isMC and isMuMu:
    foutname = "bparking_2023mumuMCtest%d.root"%testnum
elif year == 2023 and isMC:
    foutname = "bparking_2023bkgMCtest%d.root"%testnum 
elif not isMC:
    foutname = "bparking_test%d_%s%d_%d.root"%(testnum, let, num, arg)
else:
    if isSig:
        foutname = "bparking_sigMCtest%d.root"%testnum
    elif isMuMu:
        foutname = "bparking_mumuMCtest%d.root"%testnum
    elif central:
        foutname = "bparking_centralMCtest%d.root"%testnum
    else:
        if arg == -1:
            foutname = "bparking_bkgMCtest%d.root"%testnum
        else:
            foutname = "bparking_bkgMCtest%d_%d.root"%(testnum, arg)
#if year == 2022:
foutname = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking%d/ultraskimmed/"%year + foutname

if isMC:
    if isSig and year == 2022:
        if EEperiod == -1:
            flistname = "sigMCpreEEList.txt"
        else:
            flistname = "sigMCList.txt"
    elif year == 2023 and isSig:
        flistname = "sigMC2023List.txt"
    elif isMuMu and year == 2022:
        if EEperiod == -1:
            flistname = "mumuMCpreEEList.txt"
        else:
            flistname = "mumuMCList.txt"
    elif isMuMu and year == 2023:
        flistname = "mumuMC2023List.txt"
    elif central:
        flistname = "centralMCList.txt"
    else:
        if EEperiod == -1:
            flistname = "bkgMCpreEEList.txt"
            if arg > -1:
                flist2name = flistname
        else:
            flistname = "%sbkgMCList.txt"%("bkgMCReducedLPTLists/" if useLowPt else "")
            if arg > -1:
                flist2name = "%sbkgMCList%d.txt"%("bkgMCReducedLPTLists/" if useLowPt else "", arg)
elif syncTest and singleFile:
    flistname = "syncList.txt"
elif year == 2023 and not isMC:
    flistname = "data2023List.txt"
else:
    #flistname = "flist_%s%d_%d.txt"%(let, num, arg)
    flistname = "flist_%d%s_%d.txt"%(num, let, arg)
    #flistname = "flist_whack.txt"
    #print("WARNING: USING WHACK FLIST!!!!!")
    if not os.path.exists(flistname):
        flistname = "filelists/" + flistname

#initialize hWeights histogram
if isMC:
    newbins = [i*1.0 for i in range(6, 31)]
    for i in range(32, 41, 2):
        newbins.append(1.0*i)
    for i in range(45, 56, 5):
        newbins.append(1.0*i)
    newbins.append(70.0)
    newbins.append(100.0)
    newnptbins = len(newbins)-1
    hWeights = ROOT.TH1F("hWeights", "hWeights", newnptbins, array.array('d', newbins)) 
else:
    hWeights = None

all_fnames = []
fl = open(flistname, "r")
for lnum,line in enumerate(fl):
    ##debugging
    #if lnum > 3: break
    #fname = line.strip('/eos/uscms')
    #get rid of the /eos/uscms
    path = line.strip() #[10:]
    #fullpath = "root://cmsxrootd.fnal.gov/" + path
    if singleFile:
        fullpath = path
    else:
        fullpath = "root://cmseos.fnal.gov/" + path
    print("Copying file %s"%(fullpath)) 
    #print("WARNING: NOT DOING xrdcp!!!")
    #if not (isMC and not isSig):
    if (not isMC) and (not syncTest) and year != 2023:
        try:
            os.system("xrdcp %s ."%fullpath)
            fname = path.split('/')[-1]
        except:
            try:
                print("2nd try at copying file")
                os.system("xrdcp %s ."%fullpath)
                fname = path.split('/')[-1]
            except:
                fname = "root://cmseos.fnal.gov/" + fullpath
    elif syncTest:
        fname = fullpath #path
    else:
        fname = fullpath
    #fname = fullpath
    if isMC:
        f = ROOT.TFile.Open(fname)
        hW = f.Get("ntuples/allGenPtEta")
        if not hW:
            hW = f.Get("allGenPtEta")
        print("hWeight entries: %d"%(hW.GetEntries())) 
        #hW.Rebin(100)
        #now need to merge the bins exactly like the xsec histogram!!
        for pb in range(hW.GetNbinsX()+1):
            hWeights.Fill( hW.GetBinCenter(pb), hW.GetBinContent(pb) ) 
            if year == 2023:
                hpTGenAll.Fill( hW.GetBinCenter(pb), hW.GetBinContent(pb) ) 
        #for 2023 now also fill the trigger efficiency plot
        #if year == 2023:
        #    #get the trig eff plot
        #    hTrig = f.Get("ntuples/trigGenPtEta") 
        #    for pb in range(hTrig.GetNbinsX()+1):
        #        hpTGenTrig.Fill(hTrig.GetBinCenter(pb), hTrig.GetBinContent(pb)) 
        #    #now get the reco eff plot
        #    hReco = f.Get("ntuples/recoGenPtEta")
        #    for pb in range(hReco.GetNbinsX()+1):
        #        hpTGenReco["mmelel"].Fill(hReco.GetBinCenter(pb), hReco.GetBinContent(pb)) 
        #    #now get the 2-d mumu M vs pT histogram
        #    if "mumu" in vtypes:
        #        h2dmumu = f.Get("ntuples/allMvsPt")
        #        for pb in range(h2dmumu.GetNbinsX()+1):
        #            for mb in range(h2dmumu.GetNbinsY()+1):
        #                hMvsPt["mumu"].Fill(h2dmumu.GetXaxis().GetBinCenter(pb), h2dmumu.GetYaxis().GetBinCenter(mb), h2dmumu.GetBinContent(pb, mb)) 
        f.Close()

    if not (isMC and not isSig and arg > -1):
        all_fnames.append(fname)

#for bkgMC with a nonzero arg, get the list of files to run over from a separate filelist
if isMC and not isSig and arg > -1:
    fl2 = open(flist2name, "r")
    for line in fl2:
        path = line.strip()
        fullpath = "root://cmseos.fnal.gov/" + path
        all_fnames.append(fullpath)

nfiles = len(all_fnames)
for fnum,fname in enumerate(all_fnames):
    print("Starting file %d / %d"%(fnum, nfiles)) 
    process_file(fname, singleVert, useOnia, hWeights)
succ = finish_processing(foutname)
#if succ:
#    os.system( "xrdcp -f %s root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/"%foutname )
#else:
#    os.system( "xrdcp -f /afs/cern.ch/work/b/bgreenbe/public/%s root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/"%foutname )
print("Done copying :)") 
if syncTest:
    syncFile.close()
    os.system( "xrdcp -f %s root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/"%syncFname )
#if not isMC and not singleFile:
#    trigF.close()
#    os.system( "xrdcp -f %s root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/"%trigFname )
    
