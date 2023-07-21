import sys

#allow events that pass the trigger only?
trg_only = True

#set this true to REJECT ANY event that have a reconstructed photon!
reject_photon = False
#set this true to reject any event that has a valid eta->mumu (.53 to .57 GeV invar. mass)
reject_etamumu = False
#require electron_ID to be greater than 0?
require_elID = True

#what test number to label the output files with
testnum = 324

isMC = False
#use the central MC just to test the triggers (not really useful anymore)
central = False
#argument is telling which files to analyze (should be about 100 files each)
if len(sys.argv) < 2:
    #print("Error: please provide an argument as an integer betwixt 0 and 35") 
    print("No argument provided so running in signal MC mode") 
    isMC = True
    isSig = True
    arg = -1
elif sys.argv[1] == "bkg":
    isMC = True
    isSig = False
    print("Running in bkg MC mode")
    arg = -1
elif sys.argv[1] == "central":
    isMC = True
    isSig = False
    central = True
    arg = -1
    print("Running test of central MC for trigger debugging purposes.")
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

start = time.time()

#lets = ['C']

#nums = [0]
mu_mass = .105658
el_mass = .000511
etamass = .547862

#if true, add only the vertex with the highest chi2 prob to the histogram
singleVert = True #False
#use low pt electrons too?
useLowPt = True

#if True, genmatch the Onia photon tracks to the LowPt electrons and discard any that are a match
useOnia = False
#minimum dR for it to be considered a successful genmatch
dRcut = .002 #.05

#parameters for invar mass histogram
nbins = 350
xmin = .45
xmax = .8

#list of vertex types 
#2pat::Electron; 2pat::Muon-2pat::Electron; 2pat::Muon-Photon; 2pat::Muon
vtypes = ["elel", "mmelel", "mmg", "mumu"]
if useLowPt:
    #2-low-pT pat::Electron; mu-mu-2-low-pT pat::Electron
    vtypes.append("lplp")
    vtypes.append("mmlplp")

#list of particle types
ptypes = ["PCTrack", "patElectron", "lowpTElectron"]

#invar mass, pT dists (dictionary for ease of adding new vertex types)
hM = {}
#mass histogram with all weights 1
hMNoWt = {}
#pT histogram
hpT = {}
#dR histogram for positive muons (MC only)
hdRP = {}
#"      "        negative   "
hdRN = {}
#separate hists for high vs. low Vxy values (cutoff 1.2 cm)
hMhiVxy = {}
hMloVxy = {}
hDxyHiVxy = {}
hDxyLoVxy = {}
hSigmaVxyHiVxy = {}
hSigmaVxyLoVxy = {}
#2-d hist of invar. mass vs. pT
hMvsPt = {}
if isMC and not isSig:
    hMeeNoG = {}
    hMeeG = {}
if useOnia:
    hVertMatched = {}
    hVertNoMatch = {}
for vtype in vtypes:
    hM[vtype] = ROOT.TH1F("hM"+vtype, "Invar. mass with "+vtype+" vertices", nbins, xmin, xmax) 
    hMNoWt[vtype] = ROOT.TH1F("hMNoWt"+vtype, "UNWEIGHTED Invar. mass with "+vtype+" vertices", nbins, xmin, xmax) 
    hpT[vtype] = ROOT.TH1F("hpT"+vtype, "pT with "+vtype+" vertices", 500, 0., 100.) 
    hdRP[vtype] = ROOT.TH1F("hdRP"+vtype, "dR b/t Onia and lowPtelectrons Pos", 1000, 0, 10.)
    hdRN[vtype] = ROOT.TH1F("hdRN"+vtype, "dR b/t Onia and lowPtelectrons Neg", 1000, 0, 10.)
    hMhiVxy[vtype] = ROOT.TH1F("hMhiVxy"+vtype, "Invar. mass with "+vtype+" vertices, Vxy>1.2", nbins, xmin, xmax) 
    hMloVxy[vtype] = ROOT.TH1F("hMloVxy"+vtype, "Invar. mass with "+vtype+" vertices, Vxy<1.2", nbins, xmin, xmax) 
    hDxyHiVxy[vtype] = ROOT.TH1F("hDxyHiVxy"+vtype, "Dxy with "+vtype+" vertices, Vxy>1.2", 10000, -.5, .5) 
    hDxyLoVxy[vtype] = ROOT.TH1F("hDxyLoVxy"+vtype, "Dxy with "+vtype+" vertices, Vxy<1.2", 10000, -.5, .5) 
    hSigmaVxyHiVxy[vtype] = ROOT.TH1F("hSigmaVxyHiVxy"+vtype, "#sigmaVxy with "+vtype+" vertices, Vxy>1.2", 10000, 0.0, 10.0) 
    hSigmaVxyLoVxy[vtype] = ROOT.TH1F("hSigmaVxyLoVxy"+vtype, "#sigmaVxy with "+vtype+" vertices, Vxy<1.2", 10000, 0.0, 10.0) 
    #hMvsPt[vtype] = ROOT.TH2F("hMvsPt"+vtype, "Invar. mass as a function of pT", 1000, 0.0, 100.0, 350, .45, .8) 
    hMvsPt[vtype] = ROOT.TH2F("hMvsPt"+vtype, "Invar. mass as a function of pT", 1000, 0.0, 100.0, 10000, .2, 10) 
    #dielectron mass withOUT photon saved
    if isMC and not isSig:
        hMeeNoG[vtype] = ROOT.TH1F("hMeeNoG"+vtype, "Invar. mass with "+vtype+" vertices for evts with NO gen photon", 1000, 0, 1.0)
        hMeeG[vtype] = ROOT.TH1F("hMeeG"+vtype, "Invar. mass with "+vtype+" for evts with a gen photon", 1000, 0, 1.0)
    if useOnia:
        hVertMatched[vtype] = ROOT.TH2F("hVertMatched"+vtype, "vx vs. vy for matched "+vtype+" tracks", 12000, -60, 60, 12000, -60, 60)
        hVertNoMatch[vtype] = ROOT.TH2F("hVertNoMatch"+vtype, "vx vs. vy for unmatched "+vtype+" tracks", 12000, -60, 60, 12000, -60, 60)
        if isMC:
            hpTVertMatched[vtype] = ROOT.TH1F("hpTVertMatched"+vtype, "Gen pT of etas in which "+vtype+" elecs were matched to converted photons", 500, 0., 100.) 
if isMC:
    #histogram of pT for all gen eta mesons
    hpTGenAll = ROOT.TH1F("hpTGenAll", "pT of all gen Eta Mesons", 500, 0., 100.)
    #gen eta mesons passing at least one trigger
    hpTGenTrig = ROOT.TH1F("hpTGenTrig", "pT of gen Eta Mesons that pass the trigger", 500, 0., 100.)
    hpTGenReco = {}
    hpTGenAcc = {}
    hvxy_gm = {}
    for vtype in vtypes:
        #gen eta mesons reco'd using packed candidates
        hpTGenReco[vtype] = ROOT.TH1F("hpTGenReco"+vtype, "pT of gen Eta Mesons that are reconstructed with "+vtype+" vertices", 500, 0., 100.)
        #gen eta mesons passing at least one trigger AND reco'd with packed cand's
        hpTGenAcc[vtype] = ROOT.TH1F("hpTGenAcc"+vtype, "pT of gen Eta Mesons that are accepted (trig+packed cand reco) with "+vtype+" vertices", 500, 0., 100.)
        #vxy values for ONLY genmatched (LowPt)electron vertices
        hvxy_gm[vtype] = ROOT.TH1F("hvxygm"+vtype, "Vxy for genmatched "+vtype+" vertices", 1000, 0.0, 10.0) 
    hGenMudR = ROOT.TH1F("hGenMudR", "dR b/t gen muons", 1000, 0.0, 1.0) 
    #leading muon pT
    hGenMupT0 = ROOT.TH1F("hGenMupT0", "pT of lead gen muon", 10000, 0.0, 100.0)
    #subleading muon pT
    hGenMupT1 = ROOT.TH1F("hGenMupT1", "pT of subleading gen muon", 10000, 0.0, 100.0)

    if isSig:
        #dR between reconstructed tracks and gen electrons
        hdR = {}
        for ptype in ptypes:
            hdR[ptype] = ROOT.TH1F("hdR"+ptype, "dR b/t reco "+ptype+" and gen electrons", 10000, 0., 1.) 

    #open the xsec file to get the pT-dependent xsec's (so can get weighted overall efficiencies)
    f_xsec = ROOT.TFile.Open("xsecs.root")
    h_xsec = f_xsec.Get("corr_xsec")
    
    #total weight of all MC events
    all_weight = 0.
    #total weight of MC events passing trigger
    trg_weight = 0.
    rec_weight = {}
    acc_weight = {}
    for vtype in vtypes:
        #total weight of MC events passing reco
        rec_weight[vtype] = 0.
        #total weight of events accepted (trg and reco)
        acc_weight[vtype] = 0.

inc_vertM = False
if inc_vertM:
    #mass hist as saved per vertex (should be same as hM)
    hMV = ROOT.TH1F("hMV", "Invar. mass with tracks in vertex", 500, .4, .9) 
    hpTV = ROOT.TH1F("hpTV", "pT with tracks in vertex", 500, 0, 100.) 

if arg < 27:
    subdir = "0000"
else:
    subdir = "0001"

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

#process all the vertices of this type.
# vtype: string for vertex type: elel, lplp, mmelel, mmlplp -- for elel and lplp also use mumu vertices
# singleVert: true for only ONE vertex, false for all
# useOnia: true to remove any tracks matched to Onia conversion candidates, false to not
# g: gen event (if MC only, obviously)
def process_vertices(e, vtype, singleVert, useOnia, xsec, evt_weight, g=None, genEtaPt=0, passedTrig=True):
    #is good eta or nah
    good = False
    bestj = -1
    if "elel" in vtype:
        lptstr = ""
    elif "lplp" in vtype:
        lptstr = "LowPt"
    elif vtype in ["mumu", "mmg"]:
        lptstr = ""
    else:
        sys.exit("Error: unrecognized vertex type %s"%(vtype)) 
    nvert = None
    #the string that this type of vertex was called in the ntuples
    vstr = vtype
    if vtype == "mmg":
        vstr = "mumu"
    nvert = eval("len(e.Vertex_%s_reduced_chi2)"%(vstr))

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
                elptP = eval("e.%sElectron_pt[elP]"%(lptstr))
                elptN = eval("e.%sElectron_pt[elN]"%(lptstr))
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
                    hpTVertMatched.Fill(genEtaPt, evt_weight)
                #remove ALL vertices containing the mindR, both pos and neg
                j = 0
                while j < nvert:
                    if (ord(vtx_veleP[j]) == minelP and mindRP < dRcut) or (ord(vtx_veleN[j]) == minelN and mindRN < dRcut):
                        #print("EVT " + str(i) + ": removing " + str(j) + " before: " + str(vtx_vrechi2)) 
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

            hdRP[vtype].Fill(mindRP)
            hdRN[vtype].Fill(mindRN)

    #fill the pT, M histograms after looping thru all vertex candidates? or in the middle
    fillAtEnd = (vstr not in ["mmelel", "mmlplp", "mmee", "mumu"] and singleVert)

    if singleVert and nvert > 0:
        #bestj, bestChi2 = min(enumerate(e.Vertex_lplp_reduced_chi2), key=lambda x: x[1])
        #bestj, bestChi2 = min(enumerate(e.Vertex_mmlplp_reduced_chi2), key=lambda x: x[1])
        #bestj, bestChi2 = min(enumerate(e.Vertex_mmlplp_reduced_chi2), key=lambda x: x[1])
        bestj, bestChi2 = min(enumerate(vtx_vrechi2), key=lambda x: x[1])
    for j in range(nvert):
        if singleVert and j != bestj : continue
        #try a cut on the chi2 value? or on vxy??
        #print("evt = %d, j = %d, nvert = %d"%(i, j, nvert)) 
        if vstr != "mumu":
            vec_elP = ROOT.TLorentzVector()
            vec_elN = ROOT.TLorentzVector()
            elP = ord(vtx_veleP[j])
            pt = eval("e.%sElectron_pt[elP]"%(lptstr)) 
            eta = eval("e.%sElectron_eta[elP]"%(lptstr)) 
            phi = eval("e.%sElectron_phi[elP]"%(lptstr))
            vec_elP.SetPtEtaPhiM(pt, eta, phi, el_mass)
            elN = ord(vtx_veleN[j])
            if require_elID:
                elIDP = eval("ord(e.%sElectron_id[elP])"%(lptstr))
                elIDN = eval("ord(e.%sElectron_id[elN])"%(lptstr))
                if elIDP == 0 or elIDN == 0:
                    continue
            pt = eval("e.%sElectron_pt[elN]"%(lptstr)) 
            eta = eval("e.%sElectron_eta[elN]"%(lptstr)) 
            phi = eval("e.%sElectron_phi[elN]"%(lptstr))
            vec_elN.SetPtEtaPhiM(pt, eta, phi, el_mass)
            try:
                vxy = (vtx_vvx[j]**2 + vtx_vvy[j]**2)**0.5
            except:
                vxy = vtx_vxy[j] 
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
                            break
                    if vstr == "lplp":
                        mee = (vec_elP + vec_elN).M()
                        if hasG:
                            hMeeG[vtype].Fill(mee, evt_weight)
                        else:
                            hMeeNoG[vtype].Fill(mee, evt_weight)
                if isSig and gm:
                    hpTGenReco[vtype].Fill(genEtaPt)
                    rec_weight[vtype] += xsec
                    if passedTrig:
                        hpTGenAcc[vtype].Fill(genEtaPt)
                        acc_weight[vtype] += xsec
        nvertMu = len(vtx_vmuP)
        bestm = 99999
        bestpt = -1
        bestjj = -1
        for jj in range(nvertMu):
            #for 4-lepton vertices MUST use the muons corresponding to these electrons!
            if vstr in ["mmelel", "mmlplp", "mmee", "mumu"] and jj != j: continue
            #charge stored as: 1 for positive; 255 for negative
            mp = ord(vtx_vmuP[jj]) 
            vec_muP = ROOT.TLorentzVector()
            #print("mp: " + str(mp) + "; nGoodMuon: " + str(ord(e.nGoodMuon)) + "; len(mp): " + str(len(vtx_vmuP)))
            vec_muP.SetPtEtaPhiM(e.Muon_pt[mp], e.Muon_eta[mp], e.Muon_phi[mp], mu_mass)
            mm = ord(vtx_vmuN[jj]) 
            #print("mm: " + str(mm) + "; nGoodMuon: " + str(ord(e.nGoodMuon)) + "; len(mm): " + str(len(vtx_vmuN))) 
            #make the 4-vector, fill the hists
            vec_muN = ROOT.TLorentzVector()
            vec_muN.SetPtEtaPhiM(e.Muon_pt[mm], e.Muon_eta[mm], e.Muon_phi[mm], mu_mass)
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
            else:
                vec_eta = vec_muP + vec_muN + vec_elP + vec_elN
            pt = vec_eta.Pt() 
            m = vec_eta.M()
            if not fillAtEnd:
                hpT[vtype].Fill(pt, evt_weight)
                ##test: see what happens if fill hist ONLY for high-vxy bois
                #if vxy > 1.2: continue
                ###### TEST #####
                hM[vtype].Fill(m, evt_weight)
            #good eta if in the right mass range
            if m > .52 and m < .58:
                good = True
            if not singleVert:
                hMNoWt[vtype].Fill(m)
                hMvsPt[vtype].Fill(pt, m, evt_weight) 
                if vxy > 1.2:
                    hMhiVxy[vtype].Fill(m, evt_weight)
                else:
                    hMloVxy[vtype].Fill(m, evt_weight)
            else:
                if abs(m-etamass) < abs(bestm-etamass):
                    bestm = m
                    bestvert = jj
                    bestpt = pt
        if fillAtEnd and bestjj > -1:
            hpT[vtype].Fill(pt, evt_weight)
            hM[vtype].Fill(m, evt_weight)
            hMNoWt[vtype].Fill(m)
            hMvsPt[vtype].Fill(pt, m, evt_weight)
            if vxy > 1.2:
                hMhiVxy[vtype].Fill(m, evt_weight)
            else:
                hMloVxy[vtype].Fill(m, evt_weight)

    return good

def process_file(fname, singleVert, useOnia):
    global all_weight, trg_weight
    print("Opening file %s"%fname) 
    f = ROOT.TFile.Open(fname)
    t = f.Get("ntuples/recoT") 
    nTot = t.GetEntries()
    if isMC:
        gt = f.Get("ntuples/genT") 
        #but needs to be one single file for this to work correctly!!!???
        sow = t.GetEntries()
        lumi = 38.48 #fb^-1 (this is the lumi for all CMS in 2022)
        if isSig:
            sys.path.insert(1, "tm_analysis/analysis/python/utils")
            import blinding
            bratio = blinding.get_blinded_BR_EtaTo2Mu2E()
        else:
            bratio = 3.1e-4

    for i,e in enumerate(t):
        #debugging
        #if i > 202: break
        if i%printevery == (printevery-1):
            print("Event %d/%d"%(i+1, nTot)) 

        passedTrig = False
        #true if rejected for having photon
        failedPhot = (reject_photon and ord(e.nGoodPhoton) > 0)
        if failedPhot: 
            continue
        trig0 = e.Triggers_fired0
        trig1 = e.Triggers_fired1
        if trig0 > 0 or ord(trig1) > 0:
        #trying to accept only 2 special triggers
        #if trig0 & ((1<<27) + (1<<26)) > 0:
            passedTrig = True
        elif trg_only:
            continue

        #weight is 1 for data, different for MC
        evt_weight = 1.0
        if isMC:
            #get the gen tree event.
            gt.GetEntry(i)
            g = gt
            #construct the eta meson
            genEtaPt = 0
            pTmu = -1
            gen_eta = ROOT.TLorentzVector()
            #first muon encountered
            old_vec = ROOT.TLorentzVector()
            for j in range(g.nGenPart):
                gid = g.GenPart_pdgId[j]
                gpt = g.GenPart_pt[j]
                if abs(gid) == 221:
                    #eta meson
                    genEtaPt = gpt
                    #break
                else:
                    geta = g.GenPart_eta[j] 
                    gphi = g.GenPart_phi[j]
                    gmass = g.GenPart_mass[j]
                    new_vec = ROOT.TLorentzVector()
                    new_vec.SetPtEtaPhiM(gpt, geta, gphi, gmass) 
                    if genEtaPt == 0:
                        gen_eta = gen_eta + new_vec
                    if abs(gid) == 13:
                        #muon
                        #if this is the 2nd muon, compute the dR and fill the histograms
                        if pTmu == -1:
                            pTmu = g.GenPart_pt[j]
                            old_vec = new_vec
                        else:
                            mudR = old_vec.DeltaR(new_vec)
                            hGenMudR.Fill(mudR)
                            if pTmu > g.GenPart_pt[j]:
                                #pTmu is lead pT, this one is sublead
                                hGenMupT0.Fill(pTmu)
                                hGenMupT1.Fill(g.GenPart_pt[j]) 
                            else:
                                hGenMupT0.Fill(g.GenPart_pt[j])
                                hGenMupT1.Fill(pTmu)
                #only 4 gen particles max probably
                if j == (g.nGenPart-1) and genEtaPt == 0:
                    genEtaPt = gen_eta.Pt()

            hpTGenAll.Fill(genEtaPt)
            xbin = h_xsec.FindBin( genEtaPt )
            xsec = h_xsec.GetBinContent( xbin ) * h_xsec.GetBinWidth( xbin )
            all_weight += xsec
            
            evt_weight = xsec * bratio * lumi / sow
            if passedTrig:
                hpTGenTrig.Fill(genEtaPt)
                trg_weight += xsec
            
        goodmumu = False
        for vtype in vtypes:
            if vtype in ["mmlplp", "mmelel"] and reject_etamumu and goodmumu:
                #print("event %d rejected for etaToMuMu!"%i) 
                continue
            if isMC:
                isGood = process_vertices(e, vtype, singleVert, useOnia, xsec, evt_weight, g, genEtaPt, passedTrig)
            else:
                #print("processing vertices: evt %d"%i) 
                isGood = process_vertices(e, vtype, singleVert, useOnia, 0, 1.0)
            if vtype == "mumu":
                goodmumu = isGood

    #rm the file now that you're done with it.
    if not (isMC and not isSig):
        os.system("rm %s"%fname)

#finish the processing and write to an output file
def finish_processing(foutname):
    fout = ROOT.TFile.Open(foutname, "recreate") 
    for vtype in vtypes:
        hM[vtype].Write()
        hMNoWt[vtype].Write()
        hpT[vtype].Write()
        hMhiVxy[vtype].Write()
        hMloVxy[vtype].Write()
        hDxyHiVxy[vtype].Write()
        hDxyLoVxy[vtype].Write()
        hSigmaVxyHiVxy[vtype].Write()
        hSigmaVxyLoVxy[vtype].Write()
        hdRP[vtype].Write()
        hdRN[vtype].Write()
        hMvsPt[vtype].Write()
        if isMC and not isSig:
            hMeeG[vtype].Write()
            hMeeNoG[vtype].Write()
        if useOnia:
            hVertMatched[vtype].Write()
            hVertNoMatch[vtype].Write()
            if isMC:
                hpTVertMatched[vtype].Write()
    if inc_vertM:
        hMV.Write()
        hpTV.Write()
    if isMC:
        hpTGenAll.Write()
        for vtype in vtypes:
            hpTGenReco[vtype].Write()
            hpTGenAcc[vtype].Write()
            hvxy_gm[vtype].Write()
        hpTGenTrig.Write()
        hGenMudR.Write()
        hGenMupT0.Write()
        hGenMupT1.Write()
        if isSig:
            for ptype in ptypes:
                hdR[ptype].Write()
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
            print("Acc (%s): %f%%"%(vtype, acc_weight[vtype]/all_weight*100)) 

if not isMC:
    foutname = "bparking_test%d_%s%d_%d.root"%(testnum, let, num, arg)
else:
    if isSig:
        foutname = "bparking_sigMCtest%d.root"%testnum
    elif central:
        foutname = "bparking_centralMCtest%d.root"%testnum
    else:
        foutname = "bparking_bkgMCtest%d.root"%testnum

if isMC:
    if isSig:
        flistname = "sigMCList.txt"
    elif central:
        flistname = "centralMCList.txt"
    else:
        flistname = "bkgMCList.txt"
else:
    flistname = "flist_%s%d_%d.txt"%(let, num, arg)
    #flistname = "flist_whack.txt"
    #print("WARNING: USING WHACK FLIST!!!!!")
    if not os.path.exists(flistname):
        flistname = "filelists/" + flistname

fl = open(flistname, "r")
for lnum,line in enumerate(fl):
    ##debugging
    #if lnum > 3: break
    #fname = line.strip('/eos/uscms')
    #get rid of the /eos/uscms
    path = line.strip()#[10:]
    #fullpath = "root://cmsxrootd.fnal.gov/" + path
    fullpath = "root://cmseos.fnal.gov/" + path
    print("Copying file %s"%(fullpath)) 
    #print("WARNING: NOT DOING xrdcp!!!")
    if not (isMC and not isSig):
        os.system("xrdcp %s ."%fullpath)
        fname = path.split('/')[-1]
    else:
        fname = fullpath
    #fname = fullpath
    process_file(fname, singleVert, useOnia)
    finish_processing(foutname)
os.system( "xrdcp -f %s root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/"%foutname )
